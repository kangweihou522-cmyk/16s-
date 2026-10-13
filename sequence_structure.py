"""Bounded FASTQ inspection; findings are not proof of barcode identity or clean status."""
from collections import Counter
from functools import lru_cache
import csv
import gzip
import json
from pathlib import Path
import re
import shutil

IUPAC = dict(zip('ACGTRYSWKMBDHVN', ['A','C','G','T','AG','CT','CG','AT','GT','AC','CGT','AGT','ACT','ACG','ACGT']))
PRESETS = {'341F': ['CCTACGGGNGGCWGCAG', 'CCTAYGGGRBGCASCAG'],
           '338F': ['ACTCCTACGGGAGGCAGCAG'], '806R': ['GGACTACHVGGGTWTCTAAT']}
COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')


def primer_list(value):
    if not isinstance(value, str): raise ValueError('引物请填写 IUPAC 序列，多条用换行、逗号或分号分隔')
    values = list(dict.fromkeys(re.split(r'[\s,;]+', value.strip().upper())))
    if not values or any(not v or not 12 <= len(v) <= 80 or any(b not in IUPAC for b in v) or
                         sum(b!='N' for b in v)<12 for v in values):
        raise ValueError('引物须为12–80bp有效 IUPAC 序列，并至少有12个非N位点')
    if len(values)>8: raise ValueError('每个方向最多8条真实引物候选')
    return values


def candidates(config):
    mode=config.get('primer_mode','auto')
    if mode=='custom': return primer_list(config.get('forward','')),primer_list(config.get('reverse',''))
    if mode in ('341F','338F'): return PRESETS[mode],PRESETS['806R']
    if mode in ('auto','unknown','removed'): return PRESETS['341F']+PRESETS['338F'],PRESETS['806R']
    raise ValueError('未知引物处理方式')


@lru_cache(maxsize=64)
def search_pattern(primer):
    # One substitution, complete primer; no partial-match inference or indels.
    bases=['['+IUPAC[b]+']' for b in primer]
    patterns=[''.join(bases)]
    for i in range(len(bases)):
        patterns.append(''.join(bases[:i]+['[ACGT]']+bases[i+1:]))
    return re.compile('(?=('+'|'.join(patterns)+'))')


def inspect_direction(sequences, primers):
    rows=[]
    for primer in primers:
        for orientation,pattern_seq in [('expected',primer),('reverse_complement',primer.translate(COMPLEMENT)[::-1])]:
            starts=Counter(); ends=Counter(); exact=0
            pattern=search_pattern(pattern_seq)
            for sequence in sequences:
                hits=[m.start() for m in pattern.finditer(sequence[:60+len(primer)]) if m.start()<=60]
                if not hits: continue
                start=min(hits); starts[start]+=1; ends[start+len(primer)]+=1
                exact+=all(b in IUPAC[p] for b,p in zip(sequence[start:start+len(primer)],pattern_seq))
            matched=sum(starts.values())
            rows.append({'primer':primer,'orientation':orientation,'matched_reads':matched,
                         'fraction':matched/len(sequences) if sequences else 0,'exact_matches':exact,
                         'start_positions':dict(sorted(starts.items())),'end_positions':dict(sorted(ends.items())),
                         'modal_start':starts.most_common(1)[0][0] if starts else None,
                         'prefix_reads':sum(n for pos,n in starts.items() if pos>0)})
    return sorted(rows,key=lambda x:(-x['fraction'],-len(x['primer']),x['primer']))


def inspect_samples(samples, config, read_fastq):
    forward,reverse=candidates(config); result=[]; messages=[]; recommended={}; reliable=True
    for sample in samples:
        evidence={}
        for direction,filename in sample['files'].items():
            sequences=[seq for _,seq in read_fastq(Path(filename),limit=1000)]
            rows=inspect_direction(sequences,forward if direction=='1' else reverse)
            expected=[r for r in rows if r['orientation']=='expected']
            selected=[r for r in expected if r['fraction']>=.9]
            # Require a stable cutting boundary. Similar primer names may overlap;
            # ambiguity is material when the implied end positions disagree.
            boundaries={max(r['end_positions'],key=r['end_positions'].get) for r in selected}
            ambiguous=len(boundaries)>1
            good=len(sequences)>=20 and bool(selected) and not ambiguous and not any(r['orientation']!='expected' and r['fraction']>=.9 for r in rows)
            reliable &= good
            evidence['R'+direction]={'sampled_reads':len(sequences),'candidates':rows,
                'reliable':good,'status':'candidate_detected' if good else 'ambiguous' if ambiguous else 'not_confirmed',
                'prefix_interpretation':'引物前附加序列，不能仅凭序列断言为barcode',
                'recommended_primers':[r['primer'] for r in selected] if good else []}
            if good:
                recommended.setdefault(direction,set()).update(r['primer'] for r in selected)
        result.append({'SampleID':sample['SampleID'],'Group':sample['Group'],'files':sample['files'],'structure':evidence})
    if not reliable: messages.append('存在未检出、方向异常或切除边界不明确的样本；不能自动判定已clean。可补充真实引物，或依据公司说明确认已clean后再运行。')
    messages.append('仅检测前1000条、起点0–60bp、完整引物且最多1个替换；至少20条且匹配90%是工作站推荐阈值，不是官方统一标准。未检出不证明已去除。')
    messages.append('仅去除匹配引物及其前方序列；未知的引物后barcode/spacer未验证，不猜测切除。')
    return {'samples':result,'reliable':reliable,'recommended_primers':{d:sorted(v) for d,v in recommended.items()},
            'warnings':messages,'barcode_identity_confirmed':False,'downstream_additional_sequence_confirmed':False}


def save_report(directory, report):
    directory=Path(directory); directory.mkdir(parents=True,exist_ok=True)
    (directory/'sequence_structure.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    with (directory/'sequence_structure.tsv').open('w',encoding='utf8',newline='') as stream:
        writer=csv.writer(stream,delimiter='\t'); writer.writerow(['SampleID','Group','read','primer','orientation','matched','sampled','fraction','modal_start','prefix_reads','status'])
        for sample in report['samples']:
            for direction,evidence in sample['structure'].items():
                for row in evidence['candidates']:
                    writer.writerow([sample['SampleID'],sample['Group'],direction,row['primer'],row['orientation'],row['matched_reads'],evidence['sampled_reads'],row['fraction'],row['modal_start'],row['prefix_reads'],evidence['status']])
    (directory/'说明.md').write_text('\n'.join(report['warnings'])+'\n',encoding='utf8')


def export_clean_fastq(export_dir, destination, original_samples, config, read_fastq):
    """Read official QIIME MANIFEST, preserve sample IDs and compressed bytes."""
    export_dir=Path(export_dir).resolve(); destination=Path(destination); destination.mkdir(parents=True,exist_ok=True)
    manifest=export_dir/'MANIFEST'
    if not manifest.is_file(): raise ValueError('clean FASTQ 导出缺少 QIIME MANIFEST')
    with manifest.open(encoding='utf-8-sig',newline='') as stream:
        lines=[line for line in stream if line.strip() and not line.startswith('#')]
    if not lines: raise ValueError('clean FASTQ MANIFEST 为空')
    delimiter='\t' if '\t' in lines[0] else ','
    rows=list(csv.DictReader(lines,delimiter=delimiter)); mapped={}
    for row in rows:
        sid=row.get('sample-id'); direction={'forward':'1','reverse':'2'}.get(row.get('direction'))
        filename=(export_dir/row.get('filename','')).resolve()
        if not sid or not direction or not filename.is_relative_to(export_dir) or not filename.is_file(): raise ValueError('clean FASTQ MANIFEST 格式或路径异常')
        if direction in mapped.setdefault(sid,{}): raise ValueError('clean FASTQ 出现重复样本/方向')
        mapped[sid][direction]=filename
    expected={s['SampleID'] for s in original_samples}
    if set(mapped)!=expected: raise ValueError('clean FASTQ 样本缺失或不匹配；可能有样本全部未匹配引物')
    cleaned=[]
    for i,sample in enumerate(original_samples,1):
        needed={'1','2'} if config['mode']=='paired' else {'1'}
        if set(mapped[sample['SampleID']])!=needed: raise ValueError('clean FASTQ 缺少对应方向：'+sample['SampleID'])
        if needed=={'1','2'}:
            ids=[[sid for sid,_ in read_fastq(mapped[sample['SampleID']][d],limit=1000)] for d in ('1','2')]
            if not ids[0] or ids[0]!=ids[1]: raise ValueError('clean R1/R2 抽样配对标识不一致或为空：'+sample['SampleID'])
        files={}
        for direction,source in mapped[sample['SampleID']].items():
            # Safe numbered filenames; mapping keeps original SampleID (incl. 001).
            target=destination/f'sample_{i:04}_R{direction}.fastq.gz'
            if source.suffix!='.gz':
                with source.open('rb') as inp,gzip.open(target,'wb') as out: shutil.copyfileobj(inp,out)
            else: shutil.copyfile(source,target)
            files[direction]=str(target.resolve())
        cleaned.append(dict(sample,files=files))
    with (destination/'clean_sample_sheet.csv').open('w',encoding='utf8',newline='') as stream:
        writer=csv.writer(stream); writer.writerow(['SampleID','Group','forward','reverse'])
        for sample in cleaned: writer.writerow([sample['SampleID'],sample['Group'],sample['files']['1'],sample['files'].get('2','')])
    residual=inspect_samples(cleaned,config,read_fastq)
    residual['reliable']=False
    residual['recommended_primers']={}
    residual['warnings']=['clean FASTQ 抽样残留检查：这里显示的是残留匹配，不是再次切除建议。',
        '引物后未知barcode/spacer不能由此证明已不存在；原始FASTQ保持只读。']
    residual['residual_warning']=any(r['orientation']=='expected' and r['fraction']>.05 for s in residual['samples'] for e in s['structure'].values() for r in e['candidates'])
    save_report(destination/'residual_check',residual)
    return {'sample_sheet':str(destination/'clean_sample_sheet.csv'),'files':[s['files'] for s in cleaned],
            'residual_warning':residual['residual_warning'],'residual_report':str(destination/'residual_check/sequence_structure.json')}
