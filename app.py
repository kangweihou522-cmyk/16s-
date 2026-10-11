#!/usr/bin/env python3
"""Local-only 16S workstation. Python standard library; QIIME2/R run in WSL."""
import csv, gzip, hashlib, io, json, os, re, secrets, shutil, signal, subprocess, threading, time, uuid, webbrowser, zipfile, platform, mimetypes
from pathlib import Path
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs, quote

ROOT = Path(__file__).resolve().parent
TOKEN = secrets.token_urlsafe(32)
JOBS = {}
LOCK = threading.Lock()
IUPAC = dict(zip('ACGTRYSWKMBDHVN', ['A','C','G','T','AG','CT','CG','AT','GT','AC','CGT','AGT','ACT','ACG','ACGT']))
PRIMERS = {'341F': 'CCTACGGGNGGCWGCAG', '338F': 'ACTCCTACGGGAGGCAGCAG', '806R': 'GGACTACHVGGGTWTCTAAT'}

def path(s):
    s = str(s).strip()
    if not s: raise ValueError('路径不能为空')
    if re.match(r'^[A-Za-z]:[\\/]', s):
        s = '/mnt/' + s[0].lower() + '/' + s[3:].replace('\\','/')
    return Path(s).expanduser().resolve()

def integer(c, key, low, high):
    v = int(c[key])
    if not low <= v <= high: raise ValueError(f'{key} 必须在 {low}–{high} 之间')
    return v

def read_sample_sheet(p):
    with p.open(encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f)
        header=reader.fieldnames or []
        allowed=[{'SampleID','Group'},{'SampleID','Group','forward','reverse'}]
        if len(header)!=len(set(header)) or set(header) not in allowed:
            raise ValueError('分组 CSV 列应为 SampleID,Group 或 SampleID,Group,forward,reverse；列名不能重复')
        rows=[]
        for line,r in enumerate(reader,start=2):
            if None in r or any(r.get(k) is None for k in header):
                width='两列' if len(header)==2 else '四列'
                raise ValueError(f'CSV 第{line}行必须恰好{width}')
            sid,group=r['SampleID'],r['Group']
            if not sid or not group or any(x in sid+group for x in '\t\r\n'):
                raise ValueError(f'第{line}行样本名和分组不能为空或含制表符/换行')
            if sid!=sid.strip() or group!=group.strip(): raise ValueError('样本名/分组两端有空格，请清理 CSV')
            if len(header)==4:
                if not r['forward']: raise ValueError(f'{sid} 的 forward 不能为空')
                if any(x in r['forward']+r['reverse'] for x in '\t\r\n'): raise ValueError('测序文件路径不能包含制表符/换行')
                if r['forward']!=r['forward'].strip() or r['reverse']!=r['reverse'].strip(): raise ValueError('测序文件路径两端有空格，请清理 CSV')
            rows.append(r)
    if not rows or len({r['SampleID'] for r in rows})!=len(rows): raise ValueError('CSV 为空或 SampleID 重复')
    return rows

def metadata(p):
    return [(r['SampleID'],r['Group']) for r in read_sample_sheet(p)]

def mapped_fastq(folder,value,sid,direction):
    if re.match(r'^[A-Za-z]:[\\/]',value) or value.startswith('/'):
        candidate=path(value)
    else:
        candidate=(folder/value.replace('\\','/')).resolve()
    if not candidate.is_relative_to(folder): raise ValueError(f'{sid} 的 {direction} 文件必须位于选择的原始目录内')
    if not re.search(r'\.(fastq|fq)(\.gz)?$',candidate.name,re.I): raise ValueError(f'{sid} 的 {direction} 不是支持的 FASTQ 文件')
    if not candidate.is_file(): raise ValueError(f'{sid} 的 {direction} 文件不存在或不是文件：{value}')
    return candidate

def reads(p, limit=1000):
    op = gzip.open if p.name.lower().endswith('.gz') else open
    with op(p, 'rt', encoding='ascii') as f:
        for _ in range(limit):
            h = f.readline()
            if not h: return
            s, plus, q = f.readline().rstrip(), f.readline(), f.readline().rstrip()
            if not h.startswith('@') or not plus.startswith('+') or not s or len(q)!=len(s): raise ValueError(f'FASTQ 记录损坏：{p.name}')
            if any(ch not in 'ACGTN' for ch in s.upper()): raise ValueError(f'FASTQ 含不支持的碱基：{p.name}')
            yield h.split()[0].removesuffix('/1').removesuffix('/2'), s.upper()

def match_primer(seq, primer):
    n = len(primer)
    return any(sum(seq[k+i] not in IUPAC[b] for i,b in enumerate(primer)) <= 1 for k in range(min(9,max(0,len(seq)-n+1))))

def discover(c):
    folder=path(c['input']); sheet=read_sample_sheet(path(c['metadata']))
    rows=[(r['SampleID'],r['Group']) for r in sheet]
    if not folder.is_dir(): raise ValueError('输入目录不存在')
    mode=c['mode']
    if mode not in ('paired','single'): raise ValueError('未知数据格式')
    explicit='forward' in sheet[0]
    found={}
    if explicit:
        used={}
        for row in sheet:
            sid=row['SampleID']
            if mode=='paired' and not row['reverse']: raise ValueError(f'{sid} 的 reverse 不能为空，双端需要两份文件')
            if mode=='single' and row['reverse']: raise ValueError(f'{sid} 的 reverse 应留空，当前选择的是单端')
            files={}
            for direction,column in [('1','forward')]+([('2','reverse')] if mode=='paired' else []):
                file=mapped_fastq(folder,row[column],sid,column)
                st=file.stat(); identity=(st.st_dev,st.st_ino) if st.st_ino else str(file)
                if identity in used: raise ValueError(f'文件重复引用：{sid}.{column} 与 {used[identity]} 使用同一 FASTQ')
                used[identity]=f'{sid}.{column}'
                files[direction]=file
            found[sid]=files
    else:
        for p in sorted(folder.iterdir()):
            if not re.search(r'\.(fastq|fq)(\.gz)?$',p.name,re.I) or not p.is_file(): continue
            stem = re.sub(r'\.(fastq|fq)(\.gz)?$','',p.name,flags=re.I)
            if mode=='paired':
                candidates=[]
                for known,_ in rows:
                    m0=re.fullmatch(r'(?:raw[._-])?'+re.escape(known)+r'(?:[._-]raw)?(?:_S\d+)?(?:_L\d{3})?[_\.-]R?([12])(?:_\d{3})?',stem,re.I)
                    if m0: candidates.append((known,m0.group(1)))
                if len(candidates)>1: raise ValueError(f'文件名对应多个 SampleID: {p.name}')
                m = re.fullmatch(r'(.+?)(?:_S\d+)?(?:_L\d{3})?[_\.-]R?([12])(?:_\d{3})?',stem,re.I)
                if not m: raise ValueError(f'不能识别配对文件 {p.name}；请使用 SampleID_R1/R2.fastq.gz 命名')
                sid, direction = candidates[0] if candidates else (m.group(1),m.group(2))
            else:
                candidates=[known for known,_ in rows if re.fullmatch(r'(?:raw[._-])?'+re.escape(known)+r'(?:[._-]raw)?',stem,re.I)]
                if len(candidates)>1: raise ValueError(f'文件名对应多个 SampleID: {p.name}')
                sid,direction=(candidates[0] if candidates else stem),'1'
            if direction in found.setdefault(sid,{}): raise ValueError(f'{sid} 有多个同方向文件；多 lane 请先经确认合并')
            found[sid][direction]=p
        expected = {r[0] for r in rows}
        if set(found)!=expected: raise ValueError(f'样本不匹配：缺文件 {sorted(expected-set(found))}；缺分组 {sorted(set(found)-expected)}')
    result=[]
    for sid,group in rows:
        files = found[sid]
        if mode=='paired' and set(files)!={'1','2'}: raise ValueError(f'{sid} 缺少 R1 或 R2')
        seqs = {k:list(reads(v)) for k,v in files.items()}
        if any(not x for x in seqs.values()): raise ValueError(f'{sid} 文件为空')
        if mode=='paired' and ([x[0] for x in seqs['1']]!=[x[0] for x in seqs['2']]): raise ValueError(f'{sid} 前1000条 R1/R2 read ID 不对应')
        evidence={}
        for k,rr in seqs.items():
            evidence['R'+k]={'sampled_reads':len(rr),'min_length':min(len(s) for _,s in rr),'max_length':max(len(s) for _,s in rr),
                'primer_match_fraction':{name:round(sum(match_primer(s,pr) for _,s in rr)/len(rr),4) for name,pr in PRIMERS.items()}}
        result.append({'SampleID':sid,'Group':group,'files':{k:str(v) for k,v in files.items()},'evidence':evidence,'mapping':'explicit' if explicit else 'automatic'})
    return result

def find_conda(c):
    executable = c.get('conda','').strip() or shutil.which('conda')
    if not executable:
        executable = next((str(p) for p in [Path.home()/'miniconda3/bin/conda',Path.home()/'anaconda3/bin/conda',Path('/opt/conda/bin/conda')] if p.is_file()),None)
    if not executable: raise ValueError('找不到 conda；请填写 WSL 中 conda 可执行文件的完整路径')
    return executable

def conda_cmd(c):
    executable=find_conda(c)
    env = c.get('env','').strip()
    if not env: raise ValueError('请填写已安装的 QIIME 2 conda 环境名')
    return [executable,'run','--no-capture-output','-n',env]

def r_command(c,*args):
    renv=c.get('r_env') or os.environ.get('WORKSTATION_R_ENV','16s-workstation-r')
    if not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*',renv): raise ValueError('R环境名只能包含字母、数字、点、下划线和短横线')
    r_options=['--version'] if args==('--version',) else ['--vanilla',*map(str,args)]
    return [find_conda(c),'run','--no-capture-output','-n',renv,'/usr/bin/env','-u','R_HOME','-u','LD_LIBRARY_PATH','-u','LD_PRELOAD',
        'R_LIBS=/dev/null','R_LIBS_USER=/dev/null','R_LIBS_SITE=/dev/null','R_PROFILE=/dev/null','R_PROFILE_USER=/dev/null','R_ENVIRON=/dev/null','R_ENVIRON_USER=/dev/null',
        'Rscript',*r_options]

def r_execution_env(c=None):
    process_env=os.environ.copy()
    for key in ['R_HOME','LD_LIBRARY_PATH','LD_PRELOAD','DYLD_LIBRARY_PATH','PKG_CONFIG_PATH','LIBRARY_PATH','CPATH','C_INCLUDE_PATH','CPLUS_INCLUDE_PATH','FONTCONFIG_PATH','FONTCONFIG_FILE']:
        process_env.pop(key,None)
    for key in ['R_LIBS','R_LIBS_USER','R_LIBS_SITE']: process_env[key]='/dev/null'
    for key in ['R_PROFILE','R_PROFILE_USER','R_ENVIRON','R_ENVIRON_USER']: process_env[key]='/dev/null'
    return process_env

def command(job,args,cwd=None,process_env=None):
    if job.get('cancel'): raise RuntimeError('用户已停止任务')
    with open(job['log'],'a',encoding='utf8') as f:
        f.write('\n$ '+ __import__('shlex').join(map(str,args))+'\n'); f.flush()
        proc=subprocess.Popen(list(map(str,args)),cwd=cwd,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,env=process_env)
        job['pid']=proc.pid
        code=proc.wait(); job.pop('pid',None)
    with Path(job['log']).with_name('commands.jsonl').open('a',encoding='utf8') as f:
        f.write(json.dumps({'args':list(map(str,args)),'cwd':str(cwd),'finished_at':time.time(),'exit_code':code},ensure_ascii=False)+'\n')
    if code: raise RuntimeError(f'命令退出码 {code}，请查看日志；已生成的结果会保留')
    if job.get('cancel'): raise RuntimeError('用户已停止任务')

def q(job,c,*args): command(job,conda_cmd(c)+['qiime',*args],job['run'])

def marker(run):
    (run/'workstation_run.json').write_text(json.dumps({'application':'16S-workstation','format':2}),encoding='utf8')

def approved_run(run):
    run=path(run)
    m=run/'workstation_run.json'
    if not m.is_file() or json.loads(m.read_text()) != {'application':'16S-workstation','format':2}: raise ValueError('只能查看本工作站创建的结果目录')
    return run

def approved_file(filename):
    filename=path(filename)
    if not filename.is_file(): raise ValueError('结果文件不存在')
    for parent in filename.parents:
        if (parent/'workstation_run.json').is_file():
            approved_run(parent)
            return filename
    raise ValueError('只能下载本工作站生成的结果文件')

def result_files(run):
    run=approved_run(run); files=[]
    for p in sorted(run.rglob('*')):
        try:
            if not p.is_file() or not p.resolve().is_relative_to(run): continue
            files.append({'name':str(p.relative_to(run)),'path':str(p),'size':p.stat().st_size,'extension':p.suffix.lower()})
        except OSError: continue
        if len(files)>=10000: break
    return {'run':str(run),'files':files}

def inventory(run):
    result=result_files(run)
    with (run/'outputs.tsv').open('w',encoding='utf8',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['path','bytes','type'])
        w.writerows((x['name'],x['size'],x['extension']) for x in result['files'])

def export_artifact(job,c,input_path,output_path):
    q(job,c,'tools','export','--input-path',input_path,'--output-path',output_path)

def preflight(job,c):
    job['step']='检查 QIIME 2、插件、R 与分类器'
    result={'python':platform.python_version(),'cpu_threads':os.cpu_count(),'checks':[]}
    for name in ['conda','QIIME 2']:
        if job.get('cancel'): raise RuntimeError('用户已停止环境检查')
        try:
            args=[find_conda(c),'env','list','--json'] if name=='conda' else conda_cmd(c)+['qiime','info']
            proc=subprocess.run(args,text=True,capture_output=True,timeout=180)
            with Path(job['log']).open('a',encoding='utf8') as f:f.write('\n'+name+'\n'+proc.stdout+proc.stderr)
            item={'component':name,'ok':proc.returncode==0,'details':(proc.stdout+proc.stderr)[-8000:]}
            result['checks'].append(item)
            if name=='conda' and proc.returncode==0: result['environments']=[Path(x).name for x in json.loads(proc.stdout).get('envs',[])]
        except Exception as e:result['checks'].append({'component':name,'ok':False,'details':str(e)})
    for plugin in ['dada2','cutadapt','feature-classifier','phylogeny','diversity','taxa','feature-table']:
        if job.get('cancel'): raise RuntimeError('用户已停止环境检查')
        try:
            proc=subprocess.run(conda_cmd(c)+['qiime',plugin,'--help'],text=True,capture_output=True,timeout=90)
            result['checks'].append({'component':'plugin '+plugin,'ok':proc.returncode==0,'details':'' if proc.returncode==0 else (proc.stdout+proc.stderr)[-2000:]})
        except Exception as e:result['checks'].append({'component':'plugin '+plugin,'ok':False,'details':str(e)})
    try:
        proc=subprocess.run(r_command(c,str(ROOT/'check_r.R'),c.get('font','Times New Roman')),text=True,capture_output=True,timeout=90,env=r_execution_env(c))
        result['checks'].append({'component':'专用 R 环境与绘图依赖','ok':proc.returncode==0,'details':(proc.stdout+proc.stderr)[-6000:]})
    except Exception as e:result['checks'].append({'component':'专用 R 环境与绘图依赖','ok':False,'details':str(e)})
    if c.get('classifier','').strip():
        try:
            classifier=path(c['classifier'])
            proc=subprocess.run(conda_cmd(c)+['qiime','tools','peek',str(classifier)],text=True,capture_output=True,timeout=90)
            detail=proc.stdout+proc.stderr
            result['checks'].append({'component':'用户选择的分类器','ok':proc.returncode==0 and 'TaxonomicClassifier' in detail,'details':detail[-3000:]})
        except Exception as e:result['checks'].append({'component':'用户选择的分类器','ok':False,'details':str(e)})
    else:result['checks'].append({'component':'用户选择的分类器','ok':False,'details':'尚未选择本地 .qza；质量检查不需要分类器，分类阶段必需。'})
    if job.get('cancel'): raise RuntimeError('用户已停止环境检查')
    job['result']=result
    job['step']='环境检查完成；请按报告处理未通过项（版本兼容性需实际分类时核验）'

def write_metadata(run, rows, filename='metadata.tsv'):
    with (run/filename).open('w',encoding='utf8',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['sample-id','Group']); w.writerow(['#q2:types','categorical']); w.writerows(rows)

def check_overlap(c):
    f=integer(c,'trunc_f',1,1000); r=integer(c,'trunc_r',1,1000)
    trimf=integer(c,'trim_f',0,500); trimr=integer(c,'trim_r',0,500)
    length=integer(c,'amplicon',1,2000)
    # length must describe sequences entering DADA2 (after cutadapt if enabled).
    if f<=trimf or r<=trimr: raise ValueError('截断位置必须大于左端 trim 长度')
    overlap=min(f,length-trimr)-max(trimf,length-r)
    if overlap<12: raise ValueError('按填写的扩增子长度计算，两端重叠不足12bp；请调整参数')
    return f,r,trimf,trimr

def preprocess(job,c,quality_only=False):
    job['step']='检查输入及运行依赖'
    samples=discover(c); out=path(c['output'])
    raw=path(c['input'])
    if out==raw or out in raw.parents: raise ValueError('输出根目录可放在原始目录的子目录中，但不能是原始目录本身或其父目录')
    if out.exists() and not out.is_dir(): raise ValueError('输出位置必须是目录，不能是已有文件')
    integer(c,'threads',1,max(1,os.cpu_count() or 1))
    if not quality_only:
        integer(c,'min_frequency',1,1000000)
        if c.get('max_depth','').strip(): integer(c,'max_depth',2,100000000)
    if not quality_only and c['mode']=='paired': check_overlap(c)
    elif not quality_only:
        if integer(c,'trunc_f',1,1000)<=integer(c,'trim_f',0,500): raise ValueError('截断位置须大于 trim')
    primer_mode=c.get('primer_mode')
    if primer_mode not in ('removed','341F','338F','custom'): raise ValueError('请根据检测证据确认引物处理方式后再运行')
    forward=c.get('forward','').strip().upper() if primer_mode=='custom' else PRIMERS.get(primer_mode,'')
    reverse=c.get('reverse','').strip().upper() if primer_mode=='custom' else PRIMERS['806R']
    if primer_mode!='removed' and (not forward or not reverse or any(b not in IUPAC for b in forward+reverse)): raise ValueError('引物序列为空或不是有效 IUPAC 序列')
    classifier=path(c['classifier']) if c.get('classifier','').strip() else Path('/nonexistent-classifier')
    if not quality_only and (not classifier.is_file() or classifier.suffix!='.qza'): raise ValueError('请选择与 QIIME 环境兼容的 V3–V4 分类器 .qza')
    reuse=not quality_only and bool(c.get('quality_run','').strip())
    if reuse:
        run=path(c['quality_run'])
        if not (run/'quality_complete.json').is_file() or (run/'preprocess_started.json').exists(): raise ValueError('质量目录未完成或已经开始过预处理；请创建新的质量检查任务')
        old=json.loads((run/'parameters.json').read_text(encoding='utf8'))
        keys=['input','metadata','mode','primer_mode','forward','reverse','env','conda']
        if any(c.get(k)!=old.get(k) for k in keys): raise ValueError('数据、引物或环境改变，请重新执行质量检查')
        input_inventory=json.loads((run/'input_inventory.json').read_text())
        for x in input_inventory:
            stat=Path(x['path']).stat()
            if stat.st_size!=x['size'] or stat.st_mtime_ns!=x['mtime_ns']: raise ValueError('原始文件已改变，请重新质量检查')
        oldmeta=[(x['SampleID'],x['Group']) for x in old['samples']]
        if oldmeta!=[(x['SampleID'],x['Group']) for x in samples]: raise ValueError('分组已改变，请重新质量检查')
        if [(x['SampleID'],x['files']) for x in old['samples']]!=[(x['SampleID'],x['files']) for x in samples]: raise ValueError('样本与测序文件的映射已改变，请重新质量检查')
        (run/'preprocess_started.json').write_text('{}')
    else:
        out.mkdir(parents=True,exist_ok=True)
        run=out/('run_'+time.strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:8]); run.mkdir()
        if not quality_only: (run/'preprocess_started.json').write_text('{}')
    marker(run)
    job['run']=str(run); job['log']=str(run/'workflow.log')
    saved=dict(c); saved['samples']=samples
    (run/'parameters.json').write_text(json.dumps(saved,ensure_ascii=False,indent=2),encoding='utf8')
    (run/'input_inventory.json').write_text(json.dumps([{'path':v,'size':Path(v).stat().st_size,'mtime_ns':Path(v).stat().st_mtime_ns} for x in samples for v in x['files'].values()],indent=2),encoding='utf8')
    write_metadata(run,[(x['SampleID'],x['Group']) for x in samples])
    with (run/'manifest.tsv').open('w',encoding='utf8',newline='') as f:
        w=csv.writer(f,delimiter='\t')
        w.writerow(['sample-id','forward-absolute-filepath','reverse-absolute-filepath'] if c['mode']=='paired' else ['sample-id','absolute-filepath'])
        for x in samples: w.writerow([x['SampleID'],x['files']['1']]+([x['files']['2']] if c['mode']=='paired' else []))
    q(job,c,'info')
    if not quality_only: command(job,r_command(c,'--version'),run,process_env=r_execution_env(c))
    if not quality_only:
        command(job,r_command(c,str(ROOT/'check_r.R'),c['font']),run,process_env=r_execution_env(c))
        q(job,c,'tools','peek',classifier)
    job['step']='导入 FASTQ 与质量概览'
    paired=c['mode']=='paired'
    source='demux.qza' if primer_mode=='removed' else 'trimmed.qza'
    if not reuse:
        q(job,c,'tools','import','--type','SampleData[PairedEndSequencesWithQuality]' if paired else 'SampleData[SequencesWithQuality]', '--input-path','manifest.tsv','--output-path','demux.qza','--input-format','PairedEndFastqManifestPhred33V2' if paired else 'SingleEndFastqManifestPhred33V2')
        q(job,c,'demux','summarize','--i-data','demux.qza','--o-visualization','demux.qzv')
        if primer_mode!='removed':
            job['step']='去除引物（未匹配引物的 reads 丢弃，数量见日志）'
            args=['cutadapt','trim-paired' if paired else 'trim-single','--i-demultiplexed-sequences','demux.qza','--p-front-f' if paired else '--p-front',forward,'--p-discard-untrimmed','--p-cores',c['threads'],'--o-trimmed-sequences','trimmed.qza']
            if paired: args += ['--p-front-r',reverse]
            q(job,c,*args); source='trimmed.qza'
            q(job,c,'demux','summarize','--i-data',source,'--o-visualization','trimmed.qzv')
    if quality_only:
        inventory(run)
        (run/'quality_complete.json').write_text('{}')
        job['step']='质量检查完成；查看质量图后填写截断参数，再继续预处理'
        return
    job['step']='DADA2 降噪'
    args=['dada2','denoise-paired' if paired else 'denoise-single','--i-demultiplexed-seqs',source,'--p-n-threads',c['threads'],'--o-table','dada2-table.qza','--o-representative-sequences','dada2-rep-seqs.qza','--o-denoising-stats','dada2-stats.qza']
    args += ['--p-trunc-len-f' if paired else '--p-trunc-len',c['trunc_f'],'--p-trim-left-f' if paired else '--p-trim-left',c['trim_f']]
    if paired: args+=['--p-trunc-len-r',c['trunc_r'],'--p-trim-left-r',c['trim_r']]
    q(job,c,*args)
    q(job,c,'metadata','tabulate','--m-input-file','dada2-stats.qza','--o-visualization','dada2-stats.qzv')
    job['step']='分类与过滤'
    q(job,c,'feature-classifier','classify-sklearn','--i-classifier',classifier,'--i-reads','dada2-rep-seqs.qza','--p-n-jobs','1','--o-classification','taxonomy.qza')
    q(job,c,'metadata','tabulate','--m-input-file','taxonomy.qza','--o-visualization','taxonomy.qzv')
    q(job,c,'feature-table','summarize','--i-table','dada2-table.qza','--m-sample-metadata-file','metadata.tsv','--o-visualization','dada2-table.qzv')
    q(job,c,'taxa','filter-table','--i-table','dada2-table.qza','--i-taxonomy','taxonomy.qza','--p-exclude','Archaea,mitochondria,chloroplast','--o-filtered-table','filtered-step1.qza')
    q(job,c,'feature-table','filter-features','--i-table','filtered-step1.qza','--p-min-frequency',c['min_frequency'],'--o-filtered-table','filtered-table.qza')
    q(job,c,'feature-table','filter-seqs','--i-data','dada2-rep-seqs.qza','--i-table','filtered-table.qza','--o-filtered-data','filtered-rep-seqs.qza')
    q(job,c,'feature-table','summarize','--i-table','filtered-table.qza','--m-sample-metadata-file','metadata.tsv','--o-visualization','filtered-table.qzv')
    q(job,c,'taxa','barplot','--i-table','filtered-table.qza','--i-taxonomy','taxonomy.qza','--m-metadata-file','metadata.tsv','--o-visualization','taxa-barplot.qzv')
    job['step']='系统发育树与抽平曲线'
    q(job,c,'phylogeny','align-to-tree-mafft-fasttree','--i-sequences','filtered-rep-seqs.qza','--p-n-threads',c['threads'],'--o-alignment','aligned.qza','--o-masked-alignment','masked.qza','--o-tree','unrooted-tree.qza','--o-rooted-tree','rooted-tree.qza')
    q(job,c,'tools','export','--input-path','filtered-table.qza','--output-path','exported-table')
    q(job,c,'tools','export','--input-path','taxonomy.qza','--output-path','exported-taxonomy')
    q(job,c,'tools','export','--input-path','filtered-rep-seqs.qza','--output-path','exported-sequences')
    command(job,conda_cmd(c)+['biom','convert','-i','exported-table/feature-table.biom','-o','exported-table/feature-table.tsv','--to-tsv'],run)
    depth_counts=totals(run)
    curve_depth=int(c['max_depth']) if c.get('max_depth','').strip() else int(max(depth_counts.values()))
    if curve_depth<2: raise ValueError('过滤后序列数不足，无法生成抽平曲线')
    saved['effective_curve_max_depth']=curve_depth
    (run/'parameters.json').write_text(json.dumps(saved,ensure_ascii=False,indent=2),encoding='utf8')
    q(job,c,'diversity','alpha-rarefaction','--i-table','filtered-table.qza','--i-phylogeny','rooted-tree.qza','--m-metadata-file','metadata.tsv','--p-max-depth',curve_depth,'--p-metrics','shannon','observed_features','faith_pd','--o-visualization','alpha-rarefaction.qzv')
    export_artifact(job,c,'dada2-stats.qza','exports/denoising-stats')
    export_artifact(job,c,'rooted-tree.qza','exports/phylogeny')
    inventory(run)
    (run/'preprocess_complete.json').write_text(json.dumps({'completed_at':time.time()}))
    job['step']='预处理完成；查看曲线和样本深度后选择抽平深度'

def totals(run):
    table=run/'exported-table/feature-table.tsv'
    with table.open(encoding='utf8') as f:
        lines=[l for l in f if not l.startswith('# Constructed')]
    rows=list(csv.reader(lines,delimiter='\t'))
    if not rows or len(rows[0])<2: raise ValueError('特征表为空')
    ids=rows[0][1:]; counts={sid:0 for sid in ids}
    for row in rows[1:]:
        if len(row)!=len(ids)+1: raise ValueError('特征表列数不一致')
        for sid,v in zip(ids,row[1:]): counts[sid]+=float(v)
    return counts

def retention(run,depth):
    c=json.loads((run/'parameters.json').read_text(encoding='utf8')); counts=totals(run)
    return [{'SampleID':x['SampleID'],'Group':x['Group'],'count':counts.get(x['SampleID'],0),'retained':counts.get(x['SampleID'],0)>=depth} for x in c['samples']]

def alpha_skip_reason(filename,rows):
    with Path(filename).open(encoding='utf8') as f: data=list(csv.reader(f,delimiter='\t'))
    if len(data)<2: raise ValueError('导出的 alpha 向量为空')
    values={r[0]:float(r[1]) for r in data[1:] if r and not r[0].startswith('#q2:')}
    if not all(sid in values for sid,_ in rows): raise ValueError('导出的 alpha 向量与保留样本不匹配')
    groups={g for _,g in rows}
    if len(set(values.values()))==1: return '指标在所有样本中完全相同，KW 无法检验'
    gs=sorted(groups)
    for i,g in enumerate(gs):
        for h in gs[i+1:]:
            if len({values[sid] for sid,group in rows if group in (g,h)})==1: return f'{g} 与 {h} 的指标完全相同，QIIME 两两 KW 无法检验'
    return ''

def distance_skip_reason(filename):
    with Path(filename).open(encoding='utf8') as f: data=list(csv.reader(f,delimiter='\t'))
    if len(data)<3: raise ValueError('导出的距离矩阵为空或样本不足')
    values=[float(v) for r in data[1:] for v in r[1:]]
    if not values: raise ValueError('距离矩阵为空')
    return '所有距离为0，群落之间没有可用于置换检验的距离差异' if all(v==0 for v in values) else ''

def diversity(job,c):
    run=path(c['run'])
    if not (run/'preprocess_complete.json').is_file(): raise ValueError('请选择本工作站已完成预处理的 run 目录')
    saved=json.loads((run/'parameters.json').read_text(encoding='utf8'))
    depth=integer(c,'depth',1,100000000); perms=integer(c,'permutations',99,99999)
    kept=retention(run,depth); rows=[(x['SampleID'],x['Group']) for x in kept if x['retained']]
    groups={g for _,g in rows}; orig={x['Group'] for x in kept}
    if len(groups)<2 or groups!=orig: raise ValueError('此深度使某个组全部丢失或不足两组，不能进行分组检验')
    if len(rows)<=len(groups) or any(sum(g==gg for _,g in rows)<2 for gg in groups): raise ValueError('每组至少需要2个保留样本；样本不足不能运行完整分组检验')
    sub=run/('diversity_'+str(depth)+'_'+uuid.uuid4().hex[:8]); sub.mkdir()
    marker(sub)
    job['run']=str(sub); job['log']=str(sub/'workflow.log'); job['step']='单次抽平与核心多样性'
    write_metadata(sub,rows); (sub/'retention.json').write_text(json.dumps(kept,ensure_ascii=False,indent=2),encoding='utf8')
    (sub/'parameters.json').write_text(json.dumps({'depth':depth,'permutations':perms,'preprocess':str(run),'config':saved},ensure_ascii=False,indent=2),encoding='utf8')
    # QIIME's random rarefaction may differ between runs. Preserve this table; do not re-rarefy for plots.
    q(job,saved,'diversity','core-metrics-phylogenetic','--i-phylogeny',run/'rooted-tree.qza','--i-table',run/'filtered-table.qza','--p-sampling-depth',depth,'--m-metadata-file','metadata.tsv','--output-dir','core-metrics-results')
    (sub/'alpha').mkdir(); (sub/'beta').mkdir()
    q(job,saved,'feature-table','summarize','--i-table','core-metrics-results/rarefied_table.qza','--m-sample-metadata-file','metadata.tsv','--o-visualization','core-metrics-results/rarefied_table.qzv')
    export_artifact(job,saved,'core-metrics-results/rarefied_table.qza','exported-rarefied-table')
    command(job,conda_cmd(saved)+['biom','convert','-i','exported-rarefied-table/feature-table.biom','-o','exported-rarefied-table/feature-table.tsv','--to-tsv'],sub)
    alpha_metrics=['observed_features','shannon','evenness','faith_pd','chao1','simpson','goods_coverage']
    for metric in ['chao1','simpson','goods_coverage']:
        q(job,saved,'diversity','alpha','--i-table','core-metrics-results/rarefied_table.qza','--p-metric',metric,'--o-alpha-diversity',f'core-metrics-results/{metric}_vector.qza')
    for metric in alpha_metrics:
        export_artifact(job,saved,f'core-metrics-results/{metric}_vector.qza',f'exports/alpha/{metric}')
    for metric in ['weighted_unifrac','unweighted_unifrac','bray_curtis','jaccard']:
        export_artifact(job,saved,f'core-metrics-results/{metric}_pcoa_results.qza',f'exports/pcoa/{metric}')
        export_artifact(job,saved,f'core-metrics-results/{metric}_distance_matrix.qza',f'exports/distance/{metric}')
    (sub/'qiime_outputs_complete.json').write_text(json.dumps({'completed_at':time.time()}))

    job['step']='α Kruskal–Wallis 与 β PERMANOVA / PERMDISP'
    stats=[]
    for metric in alpha_metrics:
        reason=alpha_skip_reason(sub/f'exports/alpha/{metric}/alpha-diversity.tsv',rows)
        stats.append({'metric':metric,'method':'KW and pairwise KW/BH','status':'skipped' if reason else 'pending','reason':reason})
        (sub/'statistics_status.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf8')
        if reason: continue
        q(job,saved,'diversity','alpha-group-significance','--i-alpha-diversity',f'core-metrics-results/{metric}_vector.qza','--m-metadata-file','metadata.tsv','--o-visualization',f'alpha/{metric}-group-significance.qzv')
        export_artifact(job,saved,f'alpha/{metric}-group-significance.qzv',f'exports/statistics/alpha_{metric}')
        stats[-1]['status']='executed'
        (sub/'statistics_status.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf8')
    for metric in ['weighted_unifrac','unweighted_unifrac','bray_curtis','jaccard']:
        reason=distance_skip_reason(sub/f'exports/distance/{metric}/distance-matrix.tsv')
        for method in ['permanova','permdisp']:
            stats.append({'metric':metric,'method':method,'status':'skipped' if reason else 'pending','reason':reason})
            (sub/'statistics_status.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf8')
            if reason: continue
            q(job,saved,'diversity','beta-group-significance','--i-distance-matrix',f'core-metrics-results/{metric}_distance_matrix.qza','--m-metadata-file','metadata.tsv','--m-metadata-column','Group','--p-method',method,'--p-pairwise','--p-permutations',perms,'--o-visualization',f'beta/{metric}-{method}.qzv')
            export_artifact(job,saved,f'beta/{metric}-{method}.qzv',f'exports/statistics/beta_{metric}_{method}')
            stats[-1]['status']='executed'
            (sub/'statistics_status.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf8')
    inventory(sub)
    (sub/'qiime_statistics_complete.json').write_text(json.dumps({'completed_at':time.time()}))
    job['step']='R 绘图：读取保存的 QIIME 2 指标、PCoA 与分类表'
    command(job,r_command(saved,str(ROOT/'render.R'),str(sub),saved['font'],str(run),saved.get('plot_scripts','')),sub,process_env=r_execution_env(saved))
    inventory(sub)
    (sub/'complete.json').write_text(json.dumps({'completed_at':time.time(),'note':'保存本次随机抽平表；重跑可能不同。'}),encoding='utf8')
    job['step']='分析完成：SVG、统计表、QZA/QZV 与抽平特征表已保存'

def worker(job,action,c):
    try:
        if action=='preflight': preflight(job,c)
        elif action=='inspect':
            job['result']=discover(c); job['step']='输入检查完成（检测前1000条，不代表全文件完整性）'
        elif action=='quality': preprocess(job,c,quality_only=True)
        elif action=='preprocess': preprocess(job,c)
        elif action=='diversity': diversity(job,c)
        job['status']='done'
    except Exception as e:
        job['status']='cancelled' if job.get('cancel') else 'failed'; job['error']=str(e)
    finally:
        if job.get('run'):
            inventory(Path(job['run']))
            Path(job['run'],'job_status.json').write_text(json.dumps({k:v for k,v in job.items() if k not in ('pid',)},ensure_ascii=False,indent=2),encoding='utf8')

def launch(action,c):
    with LOCK:
        if any(j['status']=='running' for j in JOBS.values()): raise ValueError('已有任务运行中，请等待或停止后再开始')
        jid=uuid.uuid4().hex
        job={'id':jid,'action':action,'status':'running','step':'准备任务','log':str(ROOT/('job_'+jid+'.log'))}; JOBS[jid]=job
    threading.Thread(target=worker,args=(job,action,c),daemon=True).start(); return jid

class Handler(BaseHTTPRequestHandler):
    def log_message(self,*args): pass
    def send(self,obj,status=200):
        body=json.dumps(obj,ensure_ascii=False).encode(); self.send_response(status); self.send_header('Content-Type','application/json; charset=utf-8'); self.send_header('Cache-Control','no-store'); self.send_header('X-Content-Type-Options','nosniff'); self.end_headers(); self.wfile.write(body)
    def authorized(self):
        if self.headers.get('Host') not in (f'127.0.0.1:{self.server.server_port}',f'localhost:{self.server.server_port}'): return False
        return self.headers.get('X-Workstation-Token')==TOKEN
    def do_GET(self):
        if self.path=='/':
            if self.headers.get('Host') not in (f'127.0.0.1:{self.server.server_port}',f'localhost:{self.server.server_port}'):
                self.send({'error':'Host rejected'},403); return
            body=(ROOT/'index.html').read_text(encoding='utf8').replace('__TOKEN__',TOKEN).encode()
            self.send_response(200); self.send_header('Content-Type','text/html; charset=utf-8'); self.send_header('Cache-Control','no-store'); self.send_header('Content-Security-Policy',"frame-ancestors 'none'"); self.end_headers(); self.wfile.write(body); return
        if not self.authorized(): self.send({'error':'访问被拒绝'},403); return
        try:
            u=urlparse(self.path); p=parse_qs(u.query)
            if u.path=='/api/status':
                j=JOBS[p['id'][0]]; result=dict(j)
                logfile=Path(j['log']); result['log_tail']=''
                if logfile.exists():
                    with logfile.open('rb') as f:
                        f.seek(max(0,logfile.stat().st_size-50000)); result['log_tail']=f.read().decode('utf8',errors='replace')[-20000:]
                self.send(result)
            elif u.path=='/api/browse':
                directory=path(p.get('path',['/mnt'])[0])
                if not directory.is_dir(): raise ValueError('目录不存在')
                entries=[]
                for x in sorted(directory.iterdir(),key=lambda x:x.name.lower()):
                    try: entries.append({'name':x.name,'path':str(x),'directory':x.is_dir()})
                    except OSError: continue
                    if len(entries)>=2000: break
                entries.sort(key=lambda x:(not x['directory'],x['name'].lower()))
                self.send({'path':str(directory),'parent':str(directory.parent),'entries':entries,'limited':len(entries)==2000})
            elif u.path=='/api/results': self.send(result_files(path(p['run'][0])))
            elif u.path=='/api/file':
                filename=approved_file(p['path'][0]); size=filename.stat().st_size
                self.send_response(200)
                self.send_header('Content-Type',mimetypes.guess_type(filename)[0] or 'application/octet-stream')
                self.send_header('Content-Length',str(size))
                self.send_header('Content-Disposition',"attachment; filename*=UTF-8''"+quote(filename.name))
                self.send_header('Content-Security-Policy',"default-src 'none'; sandbox")
                self.send_header('X-Content-Type-Options','nosniff'); self.end_headers()
                with filename.open('rb') as f: shutil.copyfileobj(f,self.wfile,1024*1024)
            elif u.path=='/api/retention': self.send(retention(path(p['run'][0]),int(p['depth'][0])))
            else: self.send({'error':'Not found'},404)
        except Exception as e: self.send({'error':str(e)},400)
    def do_POST(self):
        if not self.authorized(): self.send({'error':'访问被拒绝'},403); return
        try:
            length=int(self.headers.get('Content-Length','0'))
            if not 0<length<=100000: raise ValueError('请求大小异常')
            c=json.loads(self.rfile.read(length))
            if self.path=='/api/stop':
                j=JOBS[c['id']]; j['cancel']=True
                if j.get('pid'):
                    try: os.killpg(j['pid'],signal.SIGTERM)
                    except ProcessLookupError: pass
                self.send({'ok':True}); return
            action=self.path.removeprefix('/api/')
            if action not in ('preflight','inspect','quality','preprocess','diversity'): raise ValueError('未知操作')
            self.send({'id':launch(action,c)})
        except Exception as e: self.send({'error':str(e)},400)

def main():
    import argparse
    parser=argparse.ArgumentParser(); parser.add_argument('--port',type=int,default=8765); parser.add_argument('--no-browser',action='store_true'); args=parser.parse_args()
    server=ThreadingHTTPServer(('127.0.0.1',args.port),Handler)
    print(f'16S 工作站已启动：http://localhost:{args.port}\n关闭终端会停止工作站。',flush=True)
    if not args.no_browser:
        try: subprocess.Popen(['cmd.exe','/c','start','',f'http://localhost:{args.port}'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        except OSError: webbrowser.open(f'http://localhost:{args.port}')
    try: server.serve_forever()
    except KeyboardInterrupt: pass
    finally:
        for j in JOBS.values():
            if j.get('pid'):
                try: os.killpg(j['pid'],signal.SIGTERM)
                except ProcessLookupError: pass
        server.server_close()
if __name__=='__main__': main()
