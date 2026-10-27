#!/usr/bin/env python3
"""Local-only 16S workstation. Python standard library; QIIME2/R run in WSL."""
import csv, gzip, hashlib, io, json, os, re, secrets, shutil, signal, subprocess, threading, time, uuid, webbrowser, zipfile, platform, mimetypes
from pathlib import Path
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs, quote
from datetime import datetime, timezone, timedelta
from run_layout import initialize_layout, artifact_path, mapped_argument, quality_counts, denoising_counts, table_summary, stage_report
from quality_advice import recommend_parameters
from depth_advice import read_curves, recommend_depth
from execution_record import prepare_execution, capture_r_sources, command_explanation, finalize_execution, now
from group_style import resolve_group_style, save_group_style
from plot_settings import resolve_plot_settings, save_plot_settings, ui_schema
from sequence_structure import inspect_samples, save_report as save_sequence_report, export_clean_fastq
from deployment_context import SERVICE_ID, load_context, frontend_defaults, public_context
from recovery import persist_task, seal_checkpoint, ensure_task_run, recovery_plan, prepare_retry

ROOT = Path(__file__).resolve().parent
DEPLOYMENT = load_context(ROOT)

def module_api():
    from types import SimpleNamespace
    names=('ROOT','artifact_path','q','command','r_command','r_execution_env','find_conda','conda_cmd',
           'path','begin_stage','finish_stage','stage_report','export_artifact','probe',
           'approved_run','approved_file','marker','render_saved_diversity','render_saved_differential',
           'discover','JOBS')
    return SimpleNamespace(**{name:globals()[name] for name in names})

def service_identity():
    revision='local'
    try:
        result=subprocess.run(['git','-C',str(ROOT),'rev-parse','--short','HEAD'],text=True,capture_output=True,timeout=3)
        if result.returncode==0: revision=result.stdout.strip()
    except (OSError,subprocess.TimeoutExpired): pass
    if DEPLOYMENT:
        revision = DEPLOYMENT.get('version', revision)
    return {'revision':revision,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'directory':str(ROOT),'pid':os.getpid(),'deployment':public_context(DEPLOYMENT),
            'started_at':datetime.now(timezone(timedelta(hours=8))).isoformat(timespec='seconds')}

RUNTIME = service_identity()
TOKEN = secrets.token_urlsafe(32)
INSTANCE_ID = secrets.token_hex(16)
JOBS = {}
LOCK = threading.Lock()
STOPPING = False

def jobs_busy():
    # Final command/source archives and inventories are part of a running task.
    return any(j.get('status')=='running' or j.get('_execution') is not None or j.get('_task') is not None for j in JOBS.values())
IUPAC = dict(zip('ACGTRYSWKMBDHVN', ['A','C','G','T','AG','CT','CG','AT','GT','AC','CGT','AGT','ACT','ACG','ACGT']))
PRIMERS = {'341F': 'CCTACGGGNGGCWGCAG', '338F': 'ACTCCTACGGGAGGCAGCAG', '806R': 'GGACTACHVGGGTWTCTAAT'}

def path(s):
    s = str(s).strip()
    if not s: raise ValueError('路径不能为空')
    if re.match(r'^[A-Za-z]:[\\/]', s):
        s = '/mnt/' + s[0].lower() + '/' + s[3:].replace('\\','/')
    return Path(s).expanduser().resolve()

def integer(c, key, low, high):
    names={'trunc_f':'R1 截断位置','trunc_r':'R2 截断位置','trim_f':'R1 左端额外 trim','trim_r':'R2 左端额外 trim','threads':'线程数','min_frequency':'ASV 全表最低总丰度','max_depth':'抽平曲线最大深度','amplicon':'预计扩增子长度','depth':'本次抽平深度','permutations':'β 置换次数'}
    label=names.get(key,key); raw=c.get(key)
    if raw is None or not str(raw).strip(): raise ValueError(f'请填写{label}（整数）；如需自动值，请先获取本批结果的推荐')
    if isinstance(raw,bool) or not re.fullmatch(r'[+-]?\d+',str(raw).strip()): raise ValueError(f'{label}必须为整数，当前输入为 {raw!r}')
    v = int(str(raw).strip())
    if not low <= v <= high: raise ValueError(f'{label}必须在 {low}–{high} 之间')
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
        used={}
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
            real=p.resolve()
            if not real.is_relative_to(folder): raise ValueError(f'{sid} 的 R{direction} 文件必须位于选择的原始目录内')
            st=real.stat(); identity=(st.st_dev,st.st_ino) if st.st_ino else str(real)
            if identity in used: raise ValueError(f'文件重复引用：{sid}.R{direction} 与 {used[identity]} 使用同一 FASTQ')
            used[identity]=f'{sid}.R{direction}'
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
    env = c.get('env','qiime2-amplicon-2025.7').strip()
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
    args=list(map(str,args)); capture_r_sources(job,args,cwd)
    entry={'args':args,'cwd':str(Path(cwd or Path.cwd()).resolve()),'stage':job.get('stage_title',job.get('step','未指定阶段')),
           'stage_parameters':dict(job.get('stage_parameters',{})),'started_at':now(),'exit_code':None,'status':'running'}
    with open(job['log'],'a',encoding='utf8') as f:
        f.write('\n# '+command_explanation(args)+'\n$ '+ __import__('shlex').join(args)+'\n'); f.flush()
        try:
            proc=subprocess.Popen(args,cwd=cwd,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,env=process_env)
            job['pid']=proc.pid
            code=proc.wait(); entry['exit_code']=code
            entry['status']='done' if code==0 else ('cancelled' if job.get('cancel') else 'failed')
        except OSError as error:
            entry.update(status='launch_failed',launch_error=str(error))
            raise
        finally:
            job.pop('pid',None); entry['finished_at']=now()
            with Path(job['log']).with_name('commands.jsonl').open('a',encoding='utf8') as commands:
                commands.write(json.dumps(entry,ensure_ascii=False)+'\n')
            if '_execution' in job: job['_execution']['commands'].append(entry)
    if code: raise RuntimeError(f'命令退出码 {code}，请查看日志；已生成的结果会保留')
    if job.get('cancel'): raise RuntimeError('用户已停止任务')

def q(job,c,*args):
    run=Path(job['run']); values=list(map(str,args))
    for index in range(1,len(values)):
        option=values[index-1]
        is_path=option.startswith(('--i-','--o-')) or (option.startswith('--m-') and option.endswith('-file')) or option in ('--input-path','--output-path','--output-dir')
        if not is_path: continue
        values[index]=mapped_argument(run,values[index])
        if option.startswith('--o-') or option in ('--output-path','--output-dir'):
            # QIIME output directories must themselves remain absent until the action creates them.
            output=Path(values[index]); output=output if output.is_absolute() else run/output
            output.parent.mkdir(parents=True,exist_ok=True)
    command(job,conda_cmd(c)+['qiime',*values],run)

def begin_stage(job,folder,title,parameters=None):
    job['stage_folder']=folder; job['stage_title']=title; job['stage_parameters']=parameters or {}
    job['step']=title
    stage_report(job['run'],folder,title,status='running',parameters=parameters)
    persist_task(job)

def finish_stage(job,results=None,notes=()):
    stage_report(job['run'],job['stage_folder'],job['stage_title'],parameters=job.get('stage_parameters'),results=results,notes=notes)
    if job['stage_folder']=='06_rarefaction': seal_checkpoint(job['run'],'abundance')
    job.pop('stage_folder',None); job.pop('stage_title',None); job.pop('stage_parameters',None)
    persist_task(job)

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

def probe(job,args,timeout,process_env=None):
    """Record dependency probes too, before a check-result directory exists."""
    args=list(map(str,args)); capture_r_sources(job,args)
    entry={'args':args,'cwd':str(Path.cwd()),'stage':job.get('step','运行环境检查'),
           'stage_parameters':{'timeout_seconds':timeout},'started_at':now(),'exit_code':None,'status':'running'}
    with Path(job['log']).open('a',encoding='utf8') as f:
        f.write('\n# '+command_explanation(args)+'\n$ '+__import__('shlex').join(args)+'\n')
    try:
        proc=subprocess.run(args,text=True,capture_output=True,timeout=timeout,env=process_env)
        entry.update(exit_code=proc.returncode,status='done' if proc.returncode==0 else 'failed')
        return proc
    except Exception as error:
        entry.update(status='timeout' if isinstance(error,subprocess.TimeoutExpired) else 'launch_failed',launch_error=str(error))
        raise
    finally:
        entry['finished_at']=now()
        if '_execution' in job: job['_execution']['commands'].append(entry)

def preflight(job,c):
    requested_font=resolve_plot_settings(c)['general']['font_family']
    if requested_font!='original': c=dict(c,font=requested_font)
    job['step']='检查 QIIME 2、插件、R 与分类器'
    result={'python':platform.python_version(),'cpu_threads':os.cpu_count(),'checks':[]}
    for name in ['conda','QIIME 2']:
        if job.get('cancel'): raise RuntimeError('用户已停止环境检查')
        try:
            args=[find_conda(c),'env','list','--json'] if name=='conda' else conda_cmd(c)+['qiime','info']
            proc=probe(job,args,180)
            with Path(job['log']).open('a',encoding='utf8') as f:f.write('\n'+name+'\n'+proc.stdout+proc.stderr)
            item={'component':name,'ok':proc.returncode==0,'details':(proc.stdout+proc.stderr)[-8000:]}
            result['checks'].append(item)
            if name=='conda' and proc.returncode==0: result['environments']=[Path(x).name for x in json.loads(proc.stdout).get('envs',[])]
        except Exception as e:result['checks'].append({'component':name,'ok':False,'details':str(e)})
    for plugin in ['dada2','cutadapt','feature-classifier','phylogeny','diversity','taxa','feature-table']:
        if job.get('cancel'): raise RuntimeError('用户已停止环境检查')
        try:
            proc=probe(job,conda_cmd(c)+['qiime',plugin,'--help'],90)
            result['checks'].append({'component':'plugin '+plugin,'ok':proc.returncode==0,'details':'' if proc.returncode==0 else (proc.stdout+proc.stderr)[-2000:]})
        except Exception as e:result['checks'].append({'component':'plugin '+plugin,'ok':False,'details':str(e)})
    try:
        proc=probe(job,r_command(c,str(ROOT/'check_r.R'),c.get('font','Times New Roman')),90,process_env=r_execution_env(c))
        result['checks'].append({'component':'专用 R 环境与绘图依赖','ok':proc.returncode==0,'details':(proc.stdout+proc.stderr)[-6000:]})
    except Exception as e:result['checks'].append({'component':'专用 R 环境与绘图依赖','ok':False,'details':str(e)})
    if c.get('classifier','').strip():
        try:
            classifier=path(c['classifier'])
            proc=probe(job,conda_cmd(c)+['qiime','tools','peek',str(classifier)],90)
            detail=proc.stdout+proc.stderr
            result['checks'].append({'component':'用户选择的分类器','ok':proc.returncode==0 and 'TaxonomicClassifier' in detail,'details':detail[-3000:]})
        except Exception as e:result['checks'].append({'component':'用户选择的分类器','ok':False,'details':str(e)})
    else:result['checks'].append({'component':'用户选择的分类器','ok':False,'details':'尚未选择本地 .qza；质量检查不需要分类器，分类阶段必需。'})
    if job.get('cancel'): raise RuntimeError('用户已停止环境检查')
    job['result']=result
    job['step']='环境检查完成；请按报告处理未通过项（版本兼容性需实际分类时核验）'

def write_metadata(run, rows, filename='metadata.tsv'):
    with artifact_path(run,filename).open('w',encoding='utf8',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['sample-id','Group']); w.writerow(['#q2:types','categorical']); w.writerows(rows)

def check_overlap(c):
    f=integer(c,'trunc_f',1,1000); r=integer(c,'trunc_r',1,1000)
    trimf=integer(c,'trim_f',0,500); trimr=integer(c,'trim_r',0,500)
    if f<=trimf or r<=trimr: raise ValueError('截断位置必须大于左端 trim 长度')
    if not str(c.get('amplicon','')).strip(): return f,r,trimf,trimr
    length=integer(c,'amplicon',1,2000)
    # Unknown biological length stays unknown; actual DADA2 merging enforces overlap.
    overlap=min(f,length-trimr)-max(trimf,length-r)
    if overlap<12: raise ValueError('按填写的扩增子长度计算，两端重叠不足12bp；请调整参数')
    return f,r,trimf,trimr

def export_quality(job,c,visualization,destination):
    folder=artifact_path(job['run'],destination)
    if not (folder/'per-sample-fastq-counts.tsv').is_file():
        export_artifact(job,c,visualization,destination)
    return quality_counts(folder)

def quality_advice_for_run(run,amplicon=None):
    run=approved_run(run)
    saved=json.loads((run/'parameters.json').read_text(encoding='utf8'))
    logical='quality-raw' if saved['primer_mode']=='removed' else 'quality-trimmed'
    folder=artifact_path(run,logical)
    advice=recommend_parameters(folder/'forward-seven-number-summaries.tsv',
        folder/'reverse-seven-number-summaries.tsv' if saved['mode']=='paired' else None,
        amplicon=amplicon or None,cpu_count=os.cpu_count())
    advice['quality_source']=str(folder)
    advice['primer_mode']=saved['primer_mode']
    return advice

def save_quality_advice(job,c):
    try:
        advice=quality_advice_for_run(job['run'],c.get('amplicon'))
    except (ValueError,FileNotFoundError) as e:
        advice={'parameters':{},'evidence':{},'warnings':[str(e)],'requires_review':True,'method':'质量表无法生成可靠推荐，请检查质量输出并手动设置'}
    (Path(job['run'])/'parameter_recommendation.json').write_text(json.dumps(advice,ensure_ascii=False,indent=2),encoding='utf8')
    job['result']={'quality_advice':advice}
    return advice

def refresh_quality_advice(job,c):
    run=approved_run(c['quality_run'])
    if not (run/'quality_complete.json').is_file(): raise ValueError('请选择已完成的质量检查目录')
    saved=json.loads((run/'parameters.json').read_text(encoding='utf8'))
    job['run']=str(run); job['log']=str(run/'workflow.log')
    export_quality(job,saved,'demux.qzv','quality-raw')
    if saved['primer_mode']!='removed': export_quality(job,saved,'trimmed.qzv','quality-trimmed')
    advice=save_quality_advice(job,c)
    stage_report(run,'08_parameter_advice','自动参数推荐',results=advice,notes=['此操作只导出已有质量文件并生成推荐，不运行 DADA2。','扩增子长度不能由逐位置质量统计推断；未知时留空。'])
    job['step']='参数推荐完成；可手动修改，实际保留率需 DADA2 输出核验'

def save_check_results(job,c,action):
    if not c.get('output','').strip(): return
    out=path(c['output'])
    if c.get('input','').strip():
        raw=path(c['input'])
        if out==raw or out in raw.parents: raise ValueError('检查结果目录不能是原始数据目录本身或其父目录')
    out.mkdir(parents=True,exist_ok=True)
    run=out/('check_'+time.strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:8]); run.mkdir(); marker(run)
    previous=Path(job.get('log',run/'workflow.log'))
    if previous.is_file(): shutil.copyfile(previous,run/'workflow.log')
    if '_execution' in job:
        (run/'commands.jsonl').write_text(''.join(json.dumps(entry,ensure_ascii=False)+'\n' for entry in job['_execution']['commands']),encoding='utf8')
    job['run']=str(run); job['log']=str(run/'workflow.log')
    if action=='inspect' and isinstance(job.get('result'),dict) and 'samples' in job['result']:
        save_sequence_report(run/'00_checks',job['result'])
    (run/'parameters.json').write_text(json.dumps({'config':c},ensure_ascii=False,indent=2),encoding='utf8')
    stage_report(run,'00_checks','输入与配对检查' if action=='inspect' else '运行环境检查',results={'inspection':job.get('result')},notes=['输入检查只抽查前 1000 条 reads，不等于全文件完整性检测。' if action=='inspect' else '请查看每项 ok 状态；检查结束不表示所有依赖均通过，也不表示已分析测序数据。'])

def preprocess(job,c,quality_only=False):
    plotting=resolve_plot_settings(c)
    if plotting['general']['font_family']!='original': c=dict(c,font=plotting['general']['font_family'])
    c=dict(c)
    c.setdefault('env','qiime2-amplicon-2025.7'); c.setdefault('r_env','16s-workstation-r')
    if c.get('input_layout')=='merged' and not quality_only:
        raise ValueError('已提前拼接的FASTQ不能在本工作站当作真正单端执行DADA2；可检查质量，建议索取未拼接clean R1/R2')
    data_state=c.get('data_state') or ('clean' if c.get('primer_mode')=='removed' else 'raw')
    if data_state not in ('raw','clean','unknown'): raise ValueError('请选择原始、已clean或状态未知')
    c['data_state']=data_state
    c['primer_mode']='removed' if data_state=='clean' else ('auto' if c.get('primer_mode') in ('unknown','removed',None) else c['primer_mode'])
    auto=c.get('parameter_mode')=='auto'
    if auto or not str(c.get('threads','')).strip(): c['threads']=str(min(4,os.cpu_count() or 1))
    if not str(c.get('min_frequency','')).strip(): c['min_frequency']='1'
    job['step']='检查输入及运行依赖'
    samples=discover(c); out=path(c['output'])
    raw=path(c['input'])
    if out==raw or out in raw.parents: raise ValueError('输出根目录可放在原始目录的子目录中，但不能是原始目录本身或其父目录')
    if out.exists() and not out.is_dir(): raise ValueError('输出位置必须是目录，不能是已有文件')
    integer(c,'threads',1,max(1,os.cpu_count() or 1))
    if not quality_only:
        integer(c,'min_frequency',1,1000000)
        if c.get('max_depth','').strip(): integer(c,'max_depth',3,100000000)
    if not quality_only and not auto and c['mode']=='paired': check_overlap(c)
    elif not quality_only and not auto:
        if integer(c,'trunc_f',1,1000)<=integer(c,'trim_f',0,500): raise ValueError('截断位置须大于 trim')
    primer_mode=c.get('primer_mode')
    if primer_mode not in ('removed','341F','338F','custom','auto'): raise ValueError('未知引物处理方式')
    forward=c.get('forward','').strip().upper() if primer_mode=='custom' else PRIMERS.get(primer_mode,'')
    reverse=c.get('reverse','').strip().upper() if primer_mode=='custom' else PRIMERS['806R']
    sequence_report=None; cleanup_plan=None
    if data_state!='clean':
        sequence_report=inspect_samples(samples,c,reads)
        if not sequence_report['reliable']:
            job['result']={'sequence_structure':sequence_report}
            save_check_results(job,c,'inspect')
            save_sequence_report(Path(job['run'])/'00_checks',sequence_report)
            raise ValueError('序列结构未确认，已保存检测报告；请补充真实引物或依据来源确认已clean。未检出不等于clean')
        cleanup_plan=sequence_report['recommended_primers']
    classifier=path(c['classifier']) if c.get('classifier','').strip() else Path('/nonexistent-classifier')
    if not quality_only and (not classifier.is_file() or classifier.suffix!='.qza'): raise ValueError('请选择与 QIIME 环境兼容的 V3–V4 分类器 .qza')
    reuse=not quality_only and bool(c.get('quality_run','').strip())
    if reuse:
        run=path(c['quality_run'])
        if not (run/'quality_complete.json').is_file() or (run/'preprocess_started.json').exists(): raise ValueError('质量目录未完成或已经开始过预处理；请创建新的质量检查任务')
        old=json.loads((run/'parameters.json').read_text(encoding='utf8'))
        keys=['input','metadata','mode','primer_mode','forward','reverse','env','conda','input_layout']
        defaults={'input_layout':'unmerged','forward':'','reverse':''}
        if any((c.get(k) or defaults.get(k))!=(old.get(k) or defaults.get(k)) for k in keys): raise ValueError('数据、引物或环境改变，请重新执行质量检查')
        if old.get('data_state',data_state)!=data_state: raise ValueError('输入clean状态改变，请重新质量检查')
        if old.get('batch_id') and c.get('batch_id')!=old['batch_id']: raise ValueError('质量目录属于其他批次，请打开对应历史任务')
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
        initialize_layout(run,'preprocess')
        if not quality_only: (run/'preprocess_started.json').write_text('{}')
    marker(run)
    job['run']=str(run); job['log']=str(run/'workflow.log')
    saved=dict(c, **resolve_group_style(c,samples),plot_settings=plotting); saved['samples']=samples
    saved['cleanup_plan']=cleanup_plan
    if sequence_report: save_sequence_report(run/'00_inputs/sequence_detection',sequence_report)
    save_group_style(run,resolve_group_style(saved,samples))
    save_plot_settings(run,plotting)
    (run/'parameters.json').write_text(json.dumps(saved,ensure_ascii=False,indent=2),encoding='utf8')
    (run/'input_inventory.json').write_text(json.dumps([{'path':v,'size':Path(v).stat().st_size,'mtime_ns':Path(v).stat().st_mtime_ns} for x in samples for v in x['files'].values()],indent=2),encoding='utf8')
    write_metadata(run,[(x['SampleID'],x['Group']) for x in samples])
    with artifact_path(run,'manifest.tsv').open('w',encoding='utf8',newline='') as f:
        w=csv.writer(f,delimiter='\t')
        w.writerow(['sample-id','forward-absolute-filepath','reverse-absolute-filepath'] if c['mode']=='paired' else ['sample-id','absolute-filepath'])
        for x in samples: w.writerow([x['SampleID'],x['files']['1']]+([x['files']['2']] if c['mode']=='paired' else []))
    begin_stage(job,'00_inputs','输入、分组及软件记录',{'mode':c['mode'],'qiime_env':c.get('env',''),'r_env':c.get('r_env','16s-workstation-r')})
    q(job,c,'info')
    if not quality_only: command(job,r_command(c,'--version'),run,process_env=r_execution_env(c))
    if not quality_only:
        command(job,r_command(c,str(ROOT/'check_r.R'),c['font']),run,process_env=r_execution_env(c))
        q(job,c,'tools','peek',classifier)
    groups={g:sum(x['Group']==g for x in samples) for g in sorted({x['Group'] for x in samples})}
    finish_stage(job,{'sample_count':len(samples),'groups':groups,'mapping':[{'SampleID':x['SampleID'],'Group':x['Group'],'files':x['files']} for x in samples]},['每个样本的输入文件及分组已记录；原始 FASTQ 保持不变。','软件版本与实际命令见 workflow.log 和 commands.jsonl。','本工作站接收已拆分到样本的 FASTQ，不进行 barcode 拆样。'])
    paired=c['mode']=='paired'
    source='demux.qza' if primer_mode=='removed' else 'trimmed.qza'
    if not reuse:
        begin_stage(job,'01_import','FASTQ 导入与原始质量统计',{'encoding':'Phred33','quality_sample_size':10000})
        q(job,c,'tools','import','--type','SampleData[PairedEndSequencesWithQuality]' if paired else 'SampleData[SequencesWithQuality]', '--input-path','manifest.tsv','--output-path','demux.qza','--input-format','PairedEndFastqManifestPhred33V2' if paired else 'SingleEndFastqManifestPhred33V2')
        q(job,c,'demux','summarize','--i-data','demux.qza','--o-visualization','demux.qzv')
        raw_counts=export_quality(job,c,'demux.qzv','quality-raw')
        finish_stage(job,raw_counts,['样本序列数量为完整文件计数；逐位置质量分布来自随机抽样，不代表每条 read 的完整质量。','已保存 per-sample-fastq-counts.tsv、forward-seven-number-summaries.tsv；双端还保存 reverse 表。'])
        if primer_mode!='removed':
            begin_stage(job,'02_primers','匹配引物及其前方附加序列去除、clean质量检查',{'primers':cleanup_plan,'error_rate':0.1,'full_primer_overlap':True,'discard_untrimmed':True,'downstream_unknown_sequence_removed':False})
            # Each adapter has its own full-length overlap, never a shorter global
            # reverse overlap inherited from the forward primer's length.
            fronts=lambda d:[p+';min_overlap='+str(len(p)) for p in cleanup_plan[d]]
            args=['cutadapt','trim-paired' if paired else 'trim-single','--i-demultiplexed-sequences','demux.qza','--p-front-f' if paired else '--p-front',*fronts('1'),'--p-discard-untrimmed','--p-cores',c['threads'],'--o-trimmed-sequences','trimmed.qza','--p-error-rate','0.1','--p-match-adapter-wildcards']
            if paired: args += ['--p-front-r',*fronts('2')]
            q(job,c,*args); source='trimmed.qza'
            q(job,c,'demux','summarize','--i-data',source,'--o-visualization','trimmed.qzv')
            trimmed_counts=export_quality(job,c,'trimmed.qzv','quality-trimmed')
            q(job,c,'tools','export','--input-path','trimmed.qza','--output-path','02_primers/qiime_fastq_export')
            clean=export_clean_fastq(run/'02_primers/qiime_fastq_export',run/'02_primers/clean_fastq',samples,c,reads)
            saved['clean_fastq']=clean
            (run/'parameters.json').write_text(json.dumps(saved,ensure_ascii=False,indent=2),encoding='utf8')
            finish_stage(job,{'before':raw_counts,'after':trimmed_counts,'clean_fastq':clean},['按匹配序列切除完整引物及其前方附加序列；barcode身份与引物后的未知附加序列不能据此确认。','cutadapt 允许错误率0.1，支持IUPAC；与前1000条保守检测阈值不同。','未匹配的reads会被丢弃，数量见处理前后表。','clean FASTQ、样本文件对应CSV和残留检查已保存，原始数据只读。'])
            if clean['residual_warning']: raise ValueError('clean FASTQ残留引物匹配超过5%抽样阈值，请检查残留报告；暂不进入DADA2')
        else:
            stage_report(run,'02_primers','引物处理',status='skipped',parameters={'primer_mode':'removed'},notes=['用户确认 FASTQ 已去除引物，本步不重复切除；额外 trim 默认 0。'])
    else:
        export_quality(job,c,'demux.qzv','quality-raw')
        if primer_mode!='removed': export_quality(job,c,'trimmed.qzv','quality-trimmed')
    advice=save_quality_advice(job,c)
    if sequence_report: job['result']['sequence_structure']=sequence_report
    if saved.get('clean_fastq'): job['result']['clean_fastq']=saved['clean_fastq']
    stage_report(run,'08_parameter_advice','自动参数推荐',results=advice,notes=['推荐采用质量分布和抽样长度的保守启发式，不是 QIIME 2 官方固定阈值，也不是保证最优参数。','质量表无法计算每条 read 的预期错误数或推断真实扩增子长度。'])
    seal_checkpoint(run,'quality')
    if quality_only:
        inventory(run)
        (run/'quality_complete.json').write_text('{}')
        job['step']='质量检查完成；质量表和参数推荐已保存，可继续预处理'
        return
    if auto:
        if not advice['parameters']: raise ValueError('质量统计不足，不能自动设置截断；请检查推荐说明或切换手动参数')
        c.update(advice['parameters'])
    if paired: check_overlap(c)
    elif integer(c,'trunc_f',1,1000)<=integer(c,'trim_f',0,500): raise ValueError('截断位置须大于 trim')
    saved.update(c); saved['parameter_source']='automatic_quality_heuristic' if auto else 'user'
    saved['dada2_qc']={'max_ee':2,'trunc_q':2,'min_overlap':12 if paired else None,'max_merge_mismatch':0 if paired else None,'chimera_method':'consensus'}
    (run/'parameters.json').write_text(json.dumps(saved,ensure_ascii=False,indent=2),encoding='utf8')
    begin_stage(job,'03_dada2','DADA2 质控、降噪、拼接及去嵌合体',{k:c.get(k,'') for k in ('trunc_f','trunc_r','trim_f','trim_r','amplicon','threads')} | saved['dada2_qc'])
    args=['dada2','denoise-paired' if paired else 'denoise-single','--i-demultiplexed-seqs',source,'--p-n-threads',c['threads'],'--o-table','dada2-table.qza','--o-representative-sequences','dada2-rep-seqs.qza','--o-denoising-stats','dada2-stats.qza']
    args += ['--p-trunc-len-f' if paired else '--p-trunc-len',c['trunc_f'],'--p-trim-left-f' if paired else '--p-trim-left',c['trim_f']]
    args+=['--p-chimera-method','consensus','--p-trunc-q','2']
    if paired: args+=['--p-trunc-len-r',c['trunc_r'],'--p-trim-left-r',c['trim_r'],'--p-max-ee-f','2','--p-max-ee-r','2','--p-min-overlap','12','--p-max-merge-mismatch','0']
    else: args+=['--p-max-ee','2']
    q(job,c,*args)
    q(job,c,'metadata','tabulate','--m-input-file','dada2-stats.qza','--o-visualization','dada2-stats.qzv')
    export_artifact(job,c,'dada2-stats.qza','exports/denoising-stats')
    export_artifact(job,c,'dada2-table.qza','exported-dada2-table')
    command(job,conda_cmd(c)+['biom','convert','-i',str(artifact_path(run,'exported-dada2-table/feature-table.biom')),'-o',str(artifact_path(run,'exported-dada2-table/feature-table.tsv')),'--to-tsv'],run)
    dada_summary=denoising_counts(artifact_path(run,'exports/denoising-stats'))
    from qc_report import write_dada2_retention
    retention_report=write_dada2_retention(run/'03_dada2',dada_summary,paired=paired)
    job.setdefault('result',{})['dada2_retention']=retention_report
    finish_stage(job,{'denoising':dada_summary,'retention_report':retention_report,'asv_table':table_summary(artifact_path(run,'exported-dada2-table/feature-table.tsv'))},['input 是输入数量，filtered 是质量过滤后数量，denoised 是降噪后数量，merged 是拼接后数量，non-chimeric 是去嵌合体后数量；双端按读对/合并序列计数。','逐样本保留率见 dada2_retention.csv/json；低保留提示仅供排查，不自动删样，阈值不是官方合格线。','按预期错误数过滤并学习错误模型；不采用逐 read 的 Q20/10bp 滑窗质控。','双端在降噪后拼接，至少 12bp 重叠且允许 0 个错配。','consensus 在各样本检测嵌合体并汇总去除；不额外使用 USEARCH/gold。','扩增子长度未知时没有预先验证重叠，请核对每个样本 merged 和 non-chimeric 保留率。','官方依据：https://github.com/qiime2/q2-dada2/blob/2025.7.0/q2_dada2/_denoise.py'])
    begin_stage(job,'04_taxonomy_filter','物种分类与特征过滤',{'classifier':str(classifier),'min_frequency':c['min_frequency'],'exclude':['Archaea','mitochondria','chloroplast']})
    q(job,c,'feature-classifier','classify-sklearn','--i-classifier',classifier,'--i-reads','dada2-rep-seqs.qza','--p-n-jobs','1','--o-classification','taxonomy.qza')
    q(job,c,'metadata','tabulate','--m-input-file','taxonomy.qza','--o-visualization','taxonomy.qzv')
    q(job,c,'feature-table','summarize','--i-table','dada2-table.qza','--m-sample-metadata-file','metadata.tsv','--o-visualization','dada2-table.qzv')
    q(job,c,'taxa','filter-table','--i-table','dada2-table.qza','--i-taxonomy','taxonomy.qza','--p-exclude','Archaea,mitochondria,chloroplast','--o-filtered-table','filtered-step1.qza')
    q(job,c,'feature-table','filter-features','--i-table','filtered-step1.qza','--p-min-frequency',c['min_frequency'],'--o-filtered-table','filtered-table.qza')
    q(job,c,'feature-table','filter-seqs','--i-data','dada2-rep-seqs.qza','--i-table','filtered-table.qza','--o-filtered-data','filtered-rep-seqs.qza')
    q(job,c,'feature-table','summarize','--i-table','filtered-table.qza','--m-sample-metadata-file','metadata.tsv','--o-visualization','filtered-table.qzv')
    q(job,c,'taxa','barplot','--i-table','filtered-table.qza','--i-taxonomy','taxonomy.qza','--m-metadata-file','metadata.tsv','--o-visualization','taxa-barplot.qzv')
    for src,dest in [('filtered-table.qza','exported-table'),('taxonomy.qza','exported-taxonomy'),('filtered-rep-seqs.qza','exported-sequences')]: export_artifact(job,c,src,dest)
    command(job,conda_cmd(c)+['biom','convert','-i',str(artifact_path(run,'exported-table/feature-table.biom')),'-o',str(artifact_path(run,'exported-table/feature-table.tsv')),'--to-tsv'],run)
    finish_stage(job,{'before':table_summary(artifact_path(run,'exported-dada2-table/feature-table.tsv')),'after':table_summary(artifact_path(run,'exported-table/feature-table.tsv'))},['默认最低总丰度为 1，不额外删除低丰度 ASV；手动提高会影响丰富度与 Chao1。','当前流程以细菌为目标，排除古菌、叶绿体、线粒体。','本阶段输出为未抽平表；分类组成图使用它。'])
    begin_stage(job,'05_phylogeny','系统发育树构建',{'method':'MAFFT + masked alignment + FastTree','threads':c['threads']})
    q(job,c,'phylogeny','align-to-tree-mafft-fasttree','--i-sequences','filtered-rep-seqs.qza','--p-n-threads',c['threads'],'--o-alignment','aligned.qza','--o-masked-alignment','masked.qza','--o-tree','unrooted-tree.qza','--o-rooted-tree','rooted-tree.qza')
    export_artifact(job,c,'rooted-tree.qza','exports/phylogeny')
    finish_stage(job,{'input_asv_count':table_summary(artifact_path(run,'exported-table/feature-table.tsv'))['asv_count']},['使用过滤后的代表序列构建系统发育树，供 Faith PD 和 UniFrac 使用。'])
    depth_counts=totals(run)
    maximum_count=int(max(depth_counts.values(),default=0))
    requested_curve_depth=integer(c,'max_depth',3,100000000) if c.get('max_depth','').strip() else None
    curve_depth=min(requested_curve_depth,maximum_count) if requested_curve_depth is not None else maximum_count
    if curve_depth<3: raise ValueError(f'过滤后最大样本计数为 {maximum_count}；抽平曲线最大深度至少为 3，才能生成两个或更多深度点')
    curve_warnings=[]
    if requested_curve_depth is not None and requested_curve_depth>maximum_count:
        curve_warnings.append(f'请求的曲线最大深度 {requested_curve_depth} 超过过滤后最大样本计数 {maximum_count}，本次使用 {curve_depth}；原请求值已保留。')
    saved['requested_curve_max_depth']=requested_curve_depth
    saved['effective_curve_max_depth']=curve_depth
    saved['curve_depth_warnings']=curve_warnings
    (run/'parameters.json').write_text(json.dumps(saved,ensure_ascii=False,indent=2),encoding='utf8')
    curve_metrics=['shannon','simpson','chao1','faith_pd','observed_features','goods_coverage','pielou_e']
    curve_steps=min(10,curve_depth-1)
    begin_stage(job,'06_rarefaction','多深度抽平曲线',{'min_depth':1,'max_depth':curve_depth,'requested_max_depth':requested_curve_depth,'effective_max_depth':curve_depth,'maximum_sample_count':maximum_count,'warnings':curve_warnings,'steps':curve_steps,'iterations':10,'source':'filtered unrarefied table','metrics':curve_metrics})
    q(job,c,'diversity','alpha-rarefaction','--i-table','filtered-table.qza','--i-phylogeny','rooted-tree.qza','--m-metadata-file','metadata.tsv','--p-min-depth','1','--p-max-depth',curve_depth,'--p-steps',curve_steps,'--p-iterations','10','--p-metrics',*curve_metrics,'--o-visualization','alpha-rarefaction.qzv')
    export_artifact(job,c,'alpha-rarefaction.qzv','exported-rarefaction')
    depth_advice=depth_advice_for_run(run)
    (run/'depth_recommendation.json').write_text(json.dumps(depth_advice,ensure_ascii=False,indent=2),encoding='utf8')
    job.setdefault('result',{})['depth_advice']=depth_advice
    if curve_warnings: job['result']['curve_depth_warnings']=curve_warnings
    finish_stage(job,{'sample_depths':depth_counts,'max_depth':curve_depth,'requested_max_depth':requested_curve_depth,'effective_max_depth':curve_depth,'warnings':curve_warnings,'metrics':curve_metrics,'depth_advice':depth_advice},['此处只生成多个深度的诊断曲线，原始过滤后特征表不被覆盖。','正式抽平深度推荐为严格低于本批最低过滤后 reads 的最大整百数；可人工修改，不读取 P 值。','曲线的三点平台判断独立显示，阈值为工作站启发式，不是官方唯一最佳深度。','Pielou 均匀度的 QIIME 曲线参数名为 pielou_e，核心分析保存的向量名为 evenness。',*curve_warnings])
    from abundance_tables import run_abundance
    abundance=run_abundance(job,saved,module_api())
    job.setdefault('result',{})['abundance']=abundance
    inventory(run)
    (run/'preprocess_complete.json').write_text(json.dumps({'completed_at':time.time()}))
    job['step']='预处理完成；查看曲线和样本深度后选择抽平深度'

def totals(run):
    return table_summary(artifact_path(run,'exported-table/feature-table.tsv'))['per_sample']

def require_analysis_source(run,module):
    """A verified external source is distinct from completed FASTQ preprocessing."""
    if (run/'preprocess_complete.json').is_file(): return
    marker=run/'external_inputs_complete.json'
    if marker.is_file():
        data=json.loads(marker.read_text(encoding='utf8'))
        if module in data.get('supported_modules',[]): return
        raise ValueError('此导入任务尚未准备 '+module+' 所需文件；请在该模块重新选择完整输入')
    raise ValueError('请选择已完成预处理或该模块已校验的外部输入；失败任务请先检查恢复边界')

def depth_advice_for_run(run):
    run=approved_run(run)
    config=json.loads((run/'parameters.json').read_text(encoding='utf8'))
    curves=read_curves(artifact_path(run,'exported-rarefaction'),artifact_path(run,'alpha-rarefaction.qzv'))
    return recommend_depth(config['samples'],totals(run),curves)

def retention(run,depth):
    c=json.loads((run/'parameters.json').read_text(encoding='utf8')); counts=totals(run)
    return [{'SampleID':x['SampleID'],'Group':x['Group'],'count':counts.get(x['SampleID'],0),'retained':counts.get(x['SampleID'],0)>=depth} for x in c['samples']]

def alpha_skip_reason(filename,rows):
    with Path(filename).open(encoding='utf8') as f: data=list(csv.reader(f,delimiter='\t'))
    if len(data)<2: raise ValueError('导出的 alpha 向量为空')
    values={r[0]:float(r[1]) for r in data[1:] if r and not r[0].startswith('#q2:')}
    if not all(sid in values for sid,_ in rows): raise ValueError('导出的 alpha 向量与保留样本不匹配')
    if len(set(values.values()))==1: return '指标在所有样本中完全相同，KW 无法检验'
    # QIIME omits an identical pair itself; other pairs and the overall KW
    # remain estimable and must not be discarded with that one comparison.
    return ''

def distance_skip_reason(filename):
    with Path(filename).open(encoding='utf8') as f: data=list(csv.reader(f,delimiter='\t'))
    if len(data)<3: raise ValueError('导出的距离矩阵为空或样本不足')
    values=[float(v) for r in data[1:] for v in r[1:]]
    if not values: raise ValueError('距离矩阵为空')
    return '所有距离为0，群落之间没有可用于置换检验的距离差异' if all(v==0 for v in values) else ''

def diversity(job,c):
    run=path(c['run'])
    require_analysis_source(run,'diversity')
    saved=json.loads((run/'parameters.json').read_text(encoding='utf8'))
    if c.get('batch_id') and saved.get('batch_id') and c['batch_id']!=saved['batch_id']: raise ValueError('预处理结果属于其他批次，请先打开对应历史批次')
    saved=dict(saved,**resolve_group_style(c,saved['samples'],fallback=saved),plot_settings=resolve_plot_settings(c,fallback=saved),
               **{key:c[key] for key in ('conda','env','r_env','font') if c.get(key)})
    depth=integer(c,'depth',1,100000000); perms=integer(c,'permutations',99,99999)
    kept=retention(run,depth); rows=[(x['SampleID'],x['Group']) for x in kept if x['retained']]
    groups={g for _,g in rows}; orig={x['Group'] for x in kept}
    if len(groups)<2 or groups!=orig: raise ValueError('此深度使某个组全部丢失或不足两组，不能进行分组检验')
    if len(rows)<=len(groups) or any(sum(g==gg for _,g in rows)<2 for gg in groups): raise ValueError('每组至少需要2个保留样本；样本不足不能运行完整分组检验')
    analysis_root=run/'07_diversity'; analysis_root.mkdir(exist_ok=True)
    sub=analysis_root/('diversity_'+str(depth)+'_'+uuid.uuid4().hex[:8]); sub.mkdir()
    initialize_layout(sub,'diversity')
    marker(sub)
    job['run']=str(sub); job['log']=str(sub/'workflow.log'); job['step']='单次抽平与核心多样性'
    style_file=save_group_style(sub,resolve_group_style(saved,saved['samples']))
    plot_file=save_plot_settings(sub,saved['plot_settings'])
    write_metadata(sub,rows); artifact_path(sub,'retention.json').write_text(json.dumps(kept,ensure_ascii=False,indent=2),encoding='utf8')
    (sub/'parameters.json').write_text(json.dumps({'depth':depth,'permutations':perms,'preprocess':str(run),'config':saved},ensure_ascii=False,indent=2),encoding='utf8')
    # QIIME's random rarefaction may differ between runs. Preserve this table; do not re-rarefy for plots.
    begin_stage(job,'01_core_metrics','正式单次抽平与核心多样性',{'sampling_depth':depth,'retained_samples':len(rows)})
    q(job,saved,'diversity','core-metrics-phylogenetic','--i-phylogeny',artifact_path(run,'rooted-tree.qza'),'--i-table',artifact_path(run,'filtered-table.qza'),'--p-sampling-depth',depth,'--m-metadata-file','metadata.tsv','--output-dir','core-metrics-results')
    artifact_path(sub,'alpha').mkdir(parents=True,exist_ok=True); artifact_path(sub,'beta').mkdir(parents=True,exist_ok=True)
    q(job,saved,'feature-table','summarize','--i-table','core-metrics-results/rarefied_table.qza','--m-sample-metadata-file','metadata.tsv','--o-visualization','core-metrics-results/rarefied_table.qzv')
    finish_stage(job,{'retained_samples':len(rows),'dropped_samples':len(kept)-len(rows),'retention':kept},['每次核心分析仅抽平一次；保存本次抽平表，后续多样性指标和 R 绘图不重复抽平。','不足本次深度的样本被排除；重跑的随机抽平结果可能不同。'])
    begin_stage(job,'02_exports','多样性指标、距离与坐标导出')
    export_artifact(job,saved,'core-metrics-results/rarefied_table.qza','exported-rarefied-table')
    command(job,conda_cmd(saved)+['biom','convert','-i',str(artifact_path(sub,'exported-rarefied-table/feature-table.biom')),'-o',str(artifact_path(sub,'exported-rarefied-table/feature-table.tsv')),'--to-tsv'],sub)
    alpha_metrics=['observed_features','shannon','evenness','faith_pd','chao1','simpson','goods_coverage']
    for metric in ['chao1','simpson','goods_coverage']:
        q(job,saved,'diversity','alpha','--i-table','core-metrics-results/rarefied_table.qza','--p-metric',metric,'--o-alpha-diversity',f'core-metrics-results/{metric}_vector.qza')
    for metric in alpha_metrics:
        export_artifact(job,saved,f'core-metrics-results/{metric}_vector.qza',f'exports/alpha/{metric}')
    for metric in ['weighted_unifrac','unweighted_unifrac','bray_curtis','jaccard']:
        export_artifact(job,saved,f'core-metrics-results/{metric}_pcoa_results.qza',f'exports/pcoa/{metric}')
        export_artifact(job,saved,f'core-metrics-results/{metric}_distance_matrix.qza',f'exports/distance/{metric}')
    (sub/'qiime_outputs_complete.json').write_text(json.dumps({'completed_at':time.time()}))
    finish_stage(job,{'alpha_metrics':alpha_metrics,'beta_metrics':['weighted_unifrac','unweighted_unifrac','bray_curtis','jaccard'],'rarefied_table':table_summary(artifact_path(sub,'exported-rarefied-table/feature-table.tsv'))},['七项 alpha 指标与四套距离/PCoA 均来自本次保存的抽平表。'])

    group_sizes={g:sum(group==g for _,group in rows) for g in sorted(groups)}
    begin_stage(job,'03_statistics','α Kruskal–Wallis 与 β PERMANOVA / PERMDISP',{'permutations':perms,'group_column':'Group','group_sizes':group_sizes,'design':'independent samples; Group only; unrestricted permutations','pairwise_fdr_scope':'BH within each alpha metric or beta metric/method; overall P values unadjusted'})
    stats=[]
    for metric in alpha_metrics:
        reason=alpha_skip_reason(artifact_path(sub,f'exports/alpha/{metric}/alpha-diversity.tsv'),rows)
        stats.append({'metric':metric,'method':'KW and pairwise KW/BH','status':'skipped' if reason else 'pending','reason':reason})
        artifact_path(sub,'statistics_status.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf8')
        if reason: continue
        q(job,saved,'diversity','alpha-group-significance','--i-alpha-diversity',f'core-metrics-results/{metric}_vector.qza','--m-metadata-file','metadata.tsv','--o-visualization',f'alpha/{metric}-group-significance.qzv')
        export_artifact(job,saved,f'alpha/{metric}-group-significance.qzv',f'exports/statistics/alpha_{metric}')
        stats[-1]['status']='executed'
        artifact_path(sub,'statistics_status.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf8')
    for metric in ['weighted_unifrac','unweighted_unifrac','bray_curtis','jaccard']:
        distance_reason=distance_skip_reason(artifact_path(sub,f'exports/distance/{metric}/distance-matrix.tsv'))
        for method in ['permanova','permdisp']:
            reason=distance_reason
            if not reason and method=='permdisp' and min(group_sizes.values())<3:
                reason='至少一组仅2个样本；到组中位点的距离退化，跳过 PERMDISP，不生成 P 值'
            stats.append({'metric':metric,'method':method,'status':'skipped' if reason else 'pending','reason':reason})
            artifact_path(sub,'statistics_status.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf8')
            if reason: continue
            q(job,saved,'diversity','beta-group-significance','--i-distance-matrix',f'core-metrics-results/{metric}_distance_matrix.qza','--m-metadata-file','metadata.tsv','--m-metadata-column','Group','--p-method',method,'--p-pairwise','--p-permutations',perms,'--o-visualization',f'beta/{metric}-{method}.qzv')
            export_artifact(job,saved,f'beta/{metric}-{method}.qzv',f'exports/statistics/beta_{metric}_{method}')
            stats[-1]['status']='executed'
            artifact_path(sub,'statistics_status.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf8')
    inventory(sub)
    (sub/'qiime_statistics_complete.json').write_text(json.dumps({'completed_at':time.time()}))
    finish_stage(job,{'executed':sum(x['status']=='executed' for x in stats),'skipped':sum(x['status']=='skipped' for x in stats),'tests':stats},['检验采用 QIIME 2 保存的结果；常数指标、全零距离或 PERMDISP 组内样本不足会注明原因并跳过，不生成虚构 P 值。','仅支持独立样本的 Group 单因素比较；未调整笼位、重复测量或协变量。','BH 只校正每项 alpha 或每项 beta 距离/方法内的两两比较；总体 P 值未跨指标校正。']+(['至少一组少于5个样本；KW 的卡方近似可能不可靠，结果应视为探索性。'] if min(group_sizes.values())<5 else []))
    begin_stage(job,'04_R_figures','R 绘图：读取保存的 QIIME 2 结果',{'font':saved['font'],'r_env':saved.get('r_env','16s-workstation-r')})
    command(job,r_command(saved,str(ROOT/'render.R'),str(sub),saved['font'],str(run),saved.get('plot_scripts',''),str(style_file),'',str(plot_file)),sub,process_env=r_execution_env(saved))
    figure_root=artifact_path(sub,'R_original_style')
    summary_file=figure_root/'render_summary.json'
    if not summary_file.is_file(): raise ValueError('R 命令结束但没有生成绘图摘要，请查看运行日志')
    summary=json.loads(summary_file.read_text(encoding='utf8'))
    finish_stage(job,{'render':summary,'figures':{ext:len(list(figure_root.rglob('*.'+ext))) for ext in ('svg','pdf','png')}},['多样性图读取保存的抽平结果；分类组成图读取过滤后未抽平表。','R 只绘图和整理已保存的统计结果，不重新抽平、计算 PCoA 或进行检验。'])
    inventory(sub)
    (sub/'complete.json').write_text(json.dumps({'completed_at':time.time(),'note':'保存本次随机抽平表；重跑可能不同。'}),encoding='utf8')
    job['step']='分析完成：SVG、统计表、QZA/QZV 与抽平特征表已保存'
    job['result']={'diversity_run':str(sub),'group_style':resolve_group_style(saved,saved['samples'])}

def group_style_input(metadata=None,run=None):
    if run:
        source=approved_run(run)
        data=json.loads((source/'parameters.json').read_text(encoding='utf8'))
        config=data.get('config',data)
        samples=config['samples']
        style_file=source/'group_style.json'
        fallback=json.loads(style_file.read_text(encoding='utf8')) if style_file.is_file() else config
    else:
        source=path(metadata)
        samples=read_sample_sheet(source)
        fallback={}
    style=resolve_group_style({},samples,fallback=fallback)
    return {**style,'source':str(source),'sample_counts':{group:sum(row['Group']==group for row in samples) for group in style['group_order']}}

def plot_settings_input(run):
    source=approved_run(run)
    filename=source/'plot_settings.json'
    if filename.is_file():
        settings=resolve_plot_settings({'plot_settings':json.loads(filename.read_text(encoding='utf8'))})
        recorded=True
    else:
        parameters=json.loads((source/'parameters.json').read_text(encoding='utf8'))
        config=parameters.get('config',parameters)
        settings=resolve_plot_settings({},fallback=config); recorded='plot_settings' in config
    return {'plot_settings':settings,'recorded':recorded,'source':str(source)}

def render_saved_diversity(job,c):
    analysis=approved_run(c.get('diversity_run') or c.get('result_run',''))
    if not (analysis/'qiime_statistics_complete.json').is_file():
        raise ValueError('请选择已保存完整多样性指标和检验结果的目录；不能用预处理或LEfSe目录重绘')
    previous=json.loads((analysis/'parameters.json').read_text(encoding='utf8'))
    parent=approved_run(previous['preprocess'])
    saved=previous['config']
    style=resolve_group_style(c,saved['samples'],fallback=saved)
    config=dict(saved,**style,font=c.get('font') or saved['font'],plot_settings=resolve_plot_settings(c,fallback=saved),
                **{key:c[key] for key in ('conda','env','r_env') if c.get(key)})
    output=analysis/'05_replots'; output.mkdir(exist_ok=True)
    run=output/('figures_'+time.strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:8]); run.mkdir(); marker(run)
    job.update(run=str(run),log=str(run/'workflow.log'))
    style_file=save_group_style(run,style)
    plot_file=save_plot_settings(run,config['plot_settings'])
    target_args=[]
    if c.get('plot_target'):
        target_file=run/'plot_target.json'
        target_file.write_text(json.dumps(c['plot_target'],ensure_ascii=False,indent=2),encoding='utf8')
        target_args=[str(target_file)]
    (run/'parameters.json').write_text(json.dumps({'extension':'diversity-render','preprocess':str(parent),'diversity_run':str(analysis),'config':config,'extension_parameters':dict(c,diversity_run=str(analysis)),'statistics_refit':False},ensure_ascii=False,indent=2),encoding='utf8')
    figure_root=run/'04_R_figures'/'R_original_style'
    begin_stage(job,'04_R_figures','仅重绘：保留原抽平、指标及统计结果',{'group_style':style,'statistics_refit':False})
    command(job,r_command(config,str(ROOT/'render.R'),str(analysis),config['font'],str(parent),config.get('plot_scripts',''),str(style_file),str(figure_root),str(plot_file),*target_args),run,process_env=r_execution_env(config))
    summary_file=figure_root/'render_summary.json'
    if not summary_file.is_file(): raise ValueError('R 绘图结束但没有生成绘图摘要')
    summary=json.loads(summary_file.read_text(encoding='utf8'))
    finish_stage(job,{'render':summary,'group_style':style,'statistics_refit':False},['图片保存到独立子目录；复用已有抽平表、多样性指标、坐标及P/q值。'])
    (run/'complete.json').write_text(json.dumps({'statistics_refit':False,'diversity_run':str(analysis)}),encoding='utf8')
    job.update(step='按当前分组顺序与配色重绘完成',result={'diversity_run':str(analysis),'group_style':style,'statistics_refit':False})

def differential_gallery(run):
    figures=[]
    for filename in sorted(Path(run).glob('02_differential/L*/abundance/abundance_heatmap.png')):
        folder=filename.parent
        figures.append({'level':folder.parent.name,**{kind+'_'+ext:str(folder/('abundance_'+kind+'.'+ext))
                       for kind in ('heatmap','bubble') for ext in ('png','svg','pdf') if (folder/('abundance_'+kind+'.'+ext)).is_file()}})
    return figures

def render_saved_differential(job,c):
    from differential_analysis import plot_options
    analysis=approved_run(c.get('differential_run',''))
    source_file=analysis/'differential_options.json'
    if not source_file.is_file(): raise ValueError('请选择差异分析结果目录，不能使用预处理目录')
    source=json.loads(source_file.read_text(encoding='utf8'))
    if source.get('method') not in ('ancombc2','qiime_ancom'): raise ValueError('所选结果不是支持的差异分析方法')
    if c.get('plot_target'):
        stem=c['plot_target'].get('stem','')
        match=re.match(r'^02_differential/(L[2-6]_(?:phylum|class|order|family|genus))/',stem)
        if not match: raise ValueError('单图预览的分类层级无效')
        selected=[entry for entry in source['tables'] if f"L{entry['level']}_{entry['name']}"==match[1]]
        if len(selected)!=1: raise ValueError('所选图片与原始分类层级记录不一致')
        source=dict(source,tables=selected)
    if source['method']=='ancombc2' and not list((analysis/'02_differential').glob('L*/ancombc2_fit.rds')):
        raise ValueError('所选目录没有成功保存的 ANCOM-BC2 模型；请先完成一次差异分析')
    if source['method']=='qiime_ancom' and not list((analysis/'02_differential').glob('L*/ancom_results.tsv')):
        raise ValueError('所选目录没有保存的 QIIME ANCOM 检验结果')
    for entry in source['tables']:
        filename=approved_file(entry['counts_tsv'])
        if entry.get('counts_sha256') and hashlib.sha256(filename.read_bytes()).hexdigest()!=entry['counts_sha256']:
            raise ValueError('保存的差异分析计数已改变，不能与旧模型混用')
    metadata=approved_file(source['metadata_tsv'])
    if source.get('metadata_sha256') and hashlib.sha256(metadata.read_bytes()).hexdigest()!=source['metadata_sha256']:
        raise ValueError('保存的差异分析分组已改变，不能与旧模型混用')
    previous=json.loads((analysis/'parameters.json').read_text(encoding='utf8'))
    parent=approved_run(previous['preprocess']); saved=previous['config']
    style=resolve_group_style(c,saved['samples'],fallback=source.get('group_style',saved))
    presentation=plot_options(c,fallback=source)
    config=dict(saved,**style,font=c.get('font') or source.get('font') or saved.get('font','sans'),plot_settings=resolve_plot_settings(c,fallback=source),
                **{key:c[key] for key in ('conda','env','r_env') if c.get(key)})
    options=dict(source,**presentation,group_style=style,plot_settings=config['plot_settings'],font=config['font'],analysis_run=str(analysis),statistics_refit=False)
    if c.get('plot_target'): options['plot_target']=c['plot_target']
    output=analysis/'03_replots'; output.mkdir(exist_ok=True)
    run=output/('figures_'+time.strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:8]);run.mkdir();marker(run)
    job.update(run=str(run),log=str(run/'workflow.log'))
    save_group_style(run,style)
    save_plot_settings(run,config['plot_settings'])
    (run/'parameters.json').write_text(json.dumps({'extension':'differential-render','preprocess':str(parent),'differential_run':str(analysis),'config':config,'extension_parameters':dict(c,differential_run=str(analysis)),'statistics_refit':False},ensure_ascii=False,indent=2),encoding='utf8')
    options_file=run/'differential_options.json'; options_file.write_text(json.dumps(options,ensure_ascii=False,indent=2),encoding='utf8')
    begin_stage(job,'02_differential','差异菌热图、气泡图及已保存两两比较重绘',{'statistics_refit':False,'group_style':style,**presentation})
    summary_file=run/'02_differential'/'summary.json'
    try:
        command(job,r_command(config,str(ROOT/'differential_ancombc.R'),'--plot-ancom',str(options_file),str(run)),run,process_env=r_execution_env(config))
    except RuntimeError:
        if summary_file.is_file():
            result=json.loads(summary_file.read_text(encoding='utf8'))
            job['result']=dict(result,differential_run=str(analysis),differential_figures=differential_gallery(run),statistics_refit=False)
        raise
    if not summary_file.is_file(): raise ValueError('R 命令结束但没有生成差异菌绘图摘要')
    result=json.loads(summary_file.read_text(encoding='utf8'))
    finish_stage(job,result,['图片来自原始模型保存的检验结果与未抽平计数；改变颜色、顺序或图示范围不重新拟合。'])
    (run/'complete.json').write_text(json.dumps({'statistics_refit':False,'differential_run':str(analysis)}),encoding='utf8')
    job.update(step='差异菌两两比较、丰度热图与气泡图重绘完成',result=dict(result,differential_run=str(analysis),differential_figures=differential_gallery(run),statistics_refit=False))

def extension(job,c,kind):
    """Add analyses to a completed run without reprocessing FASTQ or overwriting it."""
    c={key:value for key,value in c.items() if value is not None and str(value).strip()!=''}
    api=module_api()
    if kind=='lefse-render':
        analysis=approved_run(c.get('lefse_run',''))
        if not (analysis/'lefse_complete.json').is_file(): raise ValueError('请选择已完成的 LEfSe 结果目录，不能用其他任务目录重绘')
        previous=json.loads((analysis/'parameters.json').read_text(encoding='utf8'))
        parent=approved_run(previous['preprocess']); saved=previous['config']
        output=analysis/'04_previews'
    else:
        parent=approved_run(c.get('run',''))
        require_analysis_source(parent,kind)
        saved=json.loads((parent/'parameters.json').read_text(encoding='utf8'))
        output=parent/'11_extensions'
    if kind!='lefse-render' and c.get('batch_id') and saved.get('batch_id') and c['batch_id']!=saved['batch_id']: raise ValueError('扩展分析输入属于其他批次，请先打开对应历史批次')
    saved=dict(saved,**{key:c[key] for key in ('conda','env','r_env','font') if c.get(key)})
    parent_samples=json.loads((parent/'parameters.json').read_text(encoding='utf8'))['samples']
    observed=totals(parent)
    if set(observed)-{sample['SampleID'] for sample in parent_samples}: raise ValueError('过滤后表中存在未列入本批分组的样本')
    selected=[sample for sample in parent_samples if sample['SampleID'] in observed]
    if kind in ('differential','lefse','advanced') and {sample['Group'] for sample in selected}!={sample['Group'] for sample in parent_samples}:
        raise ValueError('过滤后某个组已没有有效样本，不能忽略该组进行差异分析')
    style=resolve_group_style(c,selected,fallback=saved)
    plotting=resolve_plot_settings(c,fallback=saved)
    c=dict(c,**style,plot_settings=plotting)
    saved=dict(saved,samples=selected,**style,plot_settings=plotting)
    output.mkdir(exist_ok=True)
    run=output/(kind+'_'+time.strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:8]); run.mkdir(); marker(run)
    job.update(run=str(run),log=str(run/'workflow.log'))
    save_group_style(run,style)
    save_plot_settings(run,plotting)
    parameters={'extension':kind,'preprocess':str(parent),'config':saved,'extension_parameters':c,
                'samples_absent_from_filtered_table':[sample for sample in parent_samples if sample['SampleID'] not in observed]}
    (run/'parameters.json').write_text(json.dumps(parameters,ensure_ascii=False,indent=2),encoding='utf8')
    write_metadata(run,[(sample['SampleID'],sample['Group']) for sample in saved['samples']])
    if kind=='abundance':
        from abundance_tables import run_abundance
        result=run_abundance(job,saved,api,source_run=parent)
        job['result']={'abundance_run':str(run),'abundance':result}
    elif kind=='advanced':
        from advanced_analysis import run_advanced
        job['result']=run_advanced(job,c,saved,parent,api)
        job['step']='高级分析完成；查看各模块状态、结果和独立图片'
    elif kind=='differential':
        from differential_analysis import run_differential
        result=run_differential(job,c,saved,parent,api)
        if result is not None: job['result']=dict(result,differential_run=str(run),differential_figures=differential_gallery(run))
    elif kind=='lefse':
        from lefse_analysis import run_lefse
        result=run_lefse(job,c,saved,parent,api)
        if result is not None: job['result']=result
    elif kind=='lefse-render':
        from lefse_analysis import render_lefse
        result=render_lefse(job,c,analysis,api)
        if result is not None: job['result']=result
    elif kind=='functional':
        from functional_prediction import run_functional
        result=run_functional(job,c,saved,parent,api)
        if result is not None: job['result']=result
    else: raise ValueError('未知扩展分析')
    completion=run/(kind.replace('-','_')+'_complete.json')
    if not completion.exists(): completion.write_text(json.dumps({'completed_at':time.time(),'preprocess':str(parent)}),encoding='utf8')
    job['step']={'abundance':'ASV 与 L2–L7 丰度及相对丰度导出完成','advanced':'高级分析完成；查看各模块状态、结果和独立图片','differential':'门到属差异菌分析完成','lefse':'LEfSe 分析及图片完成','lefse-render':'按所选 LDA 阈值重绘完成，未重算检验','functional':'PICRUSt2 预测及自选 KO 注释汇总完成'}[kind]

def worker(job,action,c):
    job['_task']={'action':action,'config':dict(c),'instance_id':INSTANCE_ID,'started_at':time.time()}
    try:
        if action=='recover': action,c=prepare_retry(job,c,module_api())
        prepare_execution(job,ROOT,action,c,RUNTIME)
        if action=='preflight': preflight(job,c)
        elif action=='inspect':
            samples=discover(c)
            job['result']=inspect_samples(samples,c,reads)
            job['step']='配对与序列结构检测完成（未检出不等于clean）'
        elif action=='quality': preprocess(job,c,quality_only=True)
        elif action=='recommend': refresh_quality_advice(job,c)
        elif action=='preprocess': preprocess(job,c)
        elif action=='external-inputs':
            from external_inputs import prepare_external
            job['result']=prepare_external(job,c,module_api())
            job['step']='独立模块输入已校验并另存；可直接运行所选分析，无需 FASTQ/DADA2'
        elif action=='resume-abundance':
            from abundance_tables import resume_abundance
            resume_abundance(job,c,module_api())
        elif action=='diversity': diversity(job,c)
        elif action=='diversity-render': render_saved_diversity(job,c)
        elif action=='differential-render': render_saved_differential(job,c)
        elif action=='advanced-render':
            from advanced_analysis import render_saved_advanced
            render_saved_advanced(job,c,module_api())
        elif action=='figure-render':
            from figure_editor import render_figure
            render_figure(job,c,module_api())
        elif action in ('abundance','differential','lefse','lefse-render','functional','advanced'): extension(job,c,action)
        if action in ('inspect','preflight'): save_check_results(job,c,action)
        job['status']='done'
    except Exception as e:
        job['status']='cancelled' if job.get('cancel') else 'failed'; job['error']=str(e)
        if job.get('run') and job.get('stage_folder'):
            stage_report(job['run'],job['stage_folder'],job['stage_title'],status=job['status'],parameters=job.get('stage_parameters'),results={'error':str(e)},notes=['此步骤没有成功完成；已生成的文件和此前成功步骤保留，不能据此宣称完整分析完成。'])
    finally:
        try: ensure_task_run(job,c,module_api())
        except Exception as error:
            job['error']=str(job.get('error',''))+'；保存早期失败任务目录失败：'+str(error)
        if job.get('run'):
            try: finalize_execution(job)
            except Exception as error:
                job['status']='failed'; job['error']=str(job.get('error',''))+'；分析结果保留，但运行代码说明保存失败：'+str(error)
            try:
                inventory(Path(job['run']))
                Path(job['run'],'job_status.json').write_text(json.dumps({k:v for k,v in job.items() if k!='pid' and not k.startswith('_')},ensure_ascii=False,indent=2),encoding='utf8')
            except Exception as error:
                job['status']='failed'; job['error']=str(job.get('error',''))+'；保存任务清单失败：'+str(error)
            try: persist_task(job)
            except Exception as error:
                job['status']='failed'; job['error']=str(job.get('error',''))+'；保存任务恢复记录失败：'+str(error)
        job.pop('_execution',None)
        job.pop('_task',None)

def launch(action,c):
    with LOCK:
        if STOPPING: raise ValueError('工作站正在停止或更新，请等待重新启动')
        if jobs_busy(): raise ValueError('已有任务运行中，请等待或停止后再开始')
        jid=uuid.uuid4().hex
        log_root=ROOT
        if DEPLOYMENT:
            log_root=Path(DEPLOYMENT['root_wsl'])/'日志'/'tasks'/jid
            log_root.mkdir(parents=True,exist_ok=True)
        job={'id':jid,'action':action,'status':'running','step':'准备任务','log':str(log_root/('job_'+jid+'.log'))}; JOBS[jid]=job
    threading.Thread(target=worker,args=(job,action,c),daemon=True).start(); return jid

class Handler(BaseHTTPRequestHandler):
    def log_message(self,*args): pass
    def send(self,obj,status=200):
        body=json.dumps(obj,ensure_ascii=False).encode(); self.send_response(status); self.send_header('Content-Type','application/json; charset=utf-8'); self.send_header('Cache-Control','no-store'); self.send_header('X-Content-Type-Options','nosniff'); self.end_headers(); self.wfile.write(body)
    def authorized(self):
        if self.headers.get('Host') not in (f'127.0.0.1:{self.server.server_port}',f'localhost:{self.server.server_port}'): return False
        return self.headers.get('X-Workstation-Token')==TOKEN
    def do_GET(self):
        guide_files={'/guide/index.html':'static/guide/index.html',
                     '/guide/USER_GUIDE.md':'docs/USER_GUIDE.md',
                     '/guide/USER_GUIDE.docx':'docs/USER_GUIDE.docx',
                     '/guide/WINDOWS_DEPLOYMENT.md':'docs/WINDOWS_DEPLOYMENT.md',
                     '/guide/METHODS_AUDIT.md':'docs/METHODS_AUDIT.md'}
        guide=guide_files.get(urlparse(self.path).path)
        if guide:
            if self.headers.get('Host') not in (f'127.0.0.1:{self.server.server_port}',f'localhost:{self.server.server_port}'):
                self.send({'error':'Host rejected'},403); return
            filename=ROOT/guide
            if not filename.is_file(): self.send({'error':'说明文件尚未构建，请更新完整源码'},404); return
            body=filename.read_bytes()
            self.send_response(200)
            content_type=mimetypes.guess_type(filename.name)[0] or 'application/octet-stream'
            self.send_header('Content-Type',content_type+('; charset=utf-8' if filename.suffix in ('.html','.md') else ''))
            self.send_header('Content-Length',str(len(body)))
            self.send_header('Cache-Control','no-store')
            self.send_header('X-Content-Type-Options','nosniff')
            self.send_header('Content-Security-Policy',"default-src 'none'; style-src 'unsafe-inline'; frame-ancestors 'self'")
            if filename.suffix!='.html':
                self.send_header('Content-Disposition',"attachment; filename*=UTF-8''"+quote(filename.name))
            self.end_headers(); self.wfile.write(body); return
        if self.path=='/api/health':
            if self.headers.get('Host') not in (f'127.0.0.1:{self.server.server_port}',f'localhost:{self.server.server_port}'):
                self.send({'error':'Host rejected'},403); return
            with LOCK:
                busy = jobs_busy()
            self.send({'service_id':SERVICE_ID,'install_root':DEPLOYMENT.get('root_wsl',''),
                       'program_dir':str(ROOT),'pid':os.getpid(),'instance_id':INSTANCE_ID,
                       'port':self.server.server_port,'busy':busy,'accepting_jobs':not STOPPING,
                       'revision':RUNTIME['revision']}); return
        if self.path=='/':
            if self.headers.get('Host') not in (f'127.0.0.1:{self.server.server_port}',f'localhost:{self.server.server_port}'):
                self.send({'error':'Host rejected'},403); return
            body=(ROOT/'index.html').read_text(encoding='utf8').replace('__TOKEN__',TOKEN).replace('__PLOT_SCHEMA__',json.dumps(ui_schema(),ensure_ascii=False).replace('<','\\u003c')).replace('__DEPLOYMENT_DEFAULTS__',json.dumps(frontend_defaults(DEPLOYMENT),ensure_ascii=False).replace('<','\\u003c')).encode()
            self.send_response(200); self.send_header('Content-Type','text/html; charset=utf-8'); self.send_header('Cache-Control','no-store'); self.send_header('Content-Security-Policy',"frame-ancestors 'none'"); self.end_headers(); self.wfile.write(body); return
        if not self.authorized(): self.send({'error':'访问被拒绝'},403); return
        try:
            u=urlparse(self.path); p=parse_qs(u.query)
            if u.path=='/api/runtime': self.send(RUNTIME)
            elif u.path=='/api/deployment': self.send(public_context(DEPLOYMENT))
            elif u.path=='/api/batches':
                from batch_context import list_batches
                self.send(list_batches(p['root'][0],module_api()))
            elif u.path=='/api/batch':
                from batch_context import restore_batch
                self.send(restore_batch(p['run'][0],module_api()))
            elif u.path=='/api/recovery': self.send(recovery_plan(p['run'][0],module_api(),tuple(JOBS.values())))
            elif u.path=='/api/tasks':
                recent=list(JOBS.values())[::-1]
                selected=[j for j in recent if j.get('status')=='running']+[j for j in recent if j.get('status')!='running'][:20]
                fields=('id','action','recovery_action','status','step','run','error')
                self.send({'jobs':[{k:j[k] for k in fields if k in j} for j in selected]})
            elif u.path=='/api/group-style': self.send(group_style_input(metadata=p.get('metadata',[None])[0],run=p.get('run',[None])[0]))
            elif u.path=='/api/plot-settings': self.send(plot_settings_input(p['run'][0]))
            elif u.path=='/api/figures':
                from figure_editor import figure_catalog
                self.send(figure_catalog(p['run'][0],module_api()))
            elif u.path=='/api/figure':
                from figure_editor import describe_figure
                self.send(describe_figure(p['image'][0],module_api()))
            elif u.path=='/api/status':
                j=JOBS[p['id'][0]]; result={key:value for key,value in j.items() if not key.startswith('_')}
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
            elif u.path=='/api/quality-advice': self.send(quality_advice_for_run(p['run'][0],p.get('amplicon',[None])[0]))
            elif u.path=='/api/depth-advice': self.send(depth_advice_for_run(p['run'][0]))
            elif u.path=='/api/file':
                filename=approved_file(p['path'][0]); size=filename.stat().st_size
                self.send_response(200)
                self.send_header('Content-Type',mimetypes.guess_type(filename)[0] or 'application/octet-stream')
                self.send_header('Content-Length',str(size))
                self.send_header('Content-Disposition',"attachment; filename*=UTF-8''"+quote(filename.name))
                self.send_header('Content-Security-Policy',"default-src 'none'; sandbox")
                self.send_header('X-Content-Type-Options','nosniff'); self.end_headers()
                with filename.open('rb') as f: shutil.copyfileobj(f,self.wfile,1024*1024)
            elif u.path=='/api/retention': self.send(retention(path(p['run'][0]),integer({'depth':p.get('depth',[''])[0]},'depth',1,100000000)))
            else: self.send({'error':'Not found'},404)
        except Exception as e: self.send({'error':str(e)},400)
    def do_POST(self):
        global STOPPING
        if not self.authorized(): self.send({'error':'访问被拒绝'},403); return
        try:
            length=int(self.headers.get('Content-Length','0'))
            if not 0<length<=100000: raise ValueError('请求大小异常')
            c=json.loads(self.rfile.read(length))
            if self.path=='/api/station-stop':
                if not DEPLOYMENT:
                    self.send({'error':'请在原启动终端停止未部署的工作站'},409); return
                with LOCK:
                    if jobs_busy():
                        self.send({'error':'分析任务仍在运行，请完成或取消后再停止、更新'},409); return
                    STOPPING=True
                self.send({'stopping':True})
                threading.Thread(target=self.server.shutdown,daemon=True).start(); return
            if self.path=='/api/stop':
                j=JOBS[c['id']]; j['cancel']=True
                if j.get('pid'):
                    try: os.killpg(j['pid'],signal.SIGTERM)
                    except ProcessLookupError: pass
                self.send({'ok':True}); return
            action=self.path.removeprefix('/api/')
            if action=='recover':
                plan=recovery_plan(c.get('recovery_run',''),module_api(),tuple(JOBS.values()))
                if not plan['recoverable']: raise ValueError(plan['boundary'])
            if action not in ('external-inputs','preflight','inspect','quality','recommend','preprocess','resume-abundance','diversity','diversity-render','differential-render','advanced','advanced-render','figure-render','abundance','differential','lefse','lefse-render','functional','recover'): raise ValueError('未知操作')
            self.send({'id':launch(action,c)})
        except Exception as e: self.send({'error':str(e)},400)

def main():
    import argparse
    parser=argparse.ArgumentParser(); parser.add_argument('--port',type=int,default=8765); parser.add_argument('--no-browser',action='store_true'); args=parser.parse_args()
    server=ThreadingHTTPServer(('127.0.0.1',args.port),Handler)
    control=None
    if DEPLOYMENT:
        control=Path(DEPLOYMENT['root_wsl'])/'配置/service.json'
        record={'token':TOKEN,'pid':os.getpid(),'instance_id':INSTANCE_ID,'port':args.port,
                'program_dir':str(ROOT),'install_root':DEPLOYMENT['root_wsl']}
        temporary=control.with_suffix('.tmp')
        temporary.write_text(json.dumps(record),encoding='utf8')
        temporary.chmod(0o600)
        temporary.replace(control)
    print(f'16S 工作站已启动：http://localhost:{args.port}\n服务版本：{RUNTIME["revision"]}\n代码目录：{ROOT}\n关闭终端会停止工作站。',flush=True)
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
        if control and control.is_file():
            try:
                if json.loads(control.read_text())['instance_id']==INSTANCE_ID: control.unlink()
            except (OSError,ValueError,KeyError): pass
if __name__=='__main__': main()
