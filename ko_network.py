"""KO differential results mapped to real KEGG coordinates; never infer flux."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import time
import urllib.request
import urllib.parse
import uuid

from advanced_analysis import prepare_metadata
from external_inputs import read_metadata
from functional_prediction import read_abundance, write_table
from group_style import resolve_group_style, save_group_style
from plot_settings import resolve_plot_settings, save_plot_settings
from ko_map_geometry import normalize_map_id, write_geometry


def _json(path,value):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf8')


def _hash(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''): digest.update(block)
    return digest.hexdigest()


def _number(config,key,default,low,high):
    try:
        raw=config.get(key) if str(config.get(key,'')).strip() else default
        if isinstance(raw,bool): raise ValueError()
        value=float(raw)
        if not math.isfinite(value) or not low<=value<=high: raise ValueError()
        return value
    except (ValueError,TypeError): raise ValueError(f'{key} 必须为 {low}–{high} 的有限数值') from None


def _copy(path,folder,name):
    path=Path(path)
    if not path.is_file(): raise ValueError('输入文件不存在：'+str(path))
    destination=folder/(name+''.join(path.suffixes))
    before=_hash(path);shutil.copyfile(path,destination)
    if _hash(destination)!=before: raise ValueError('复制过程中输入发生变化：'+str(path))
    return destination,{'source':str(path),'copy':destination.name,'sha256':before}


def import_differential(filename,destination,effect_type):
    if effect_type not in ('log2FC','log2_CLR_difference','estimate'):
        raise ValueError('请选择外部差异结果的效应尺度')
    delimiter=',' if Path(filename).suffix.lower()=='.csv' else '\t'
    with Path(filename).open(encoding='utf-8-sig',newline='') as stream:
        reader=csv.DictReader(stream,delimiter=delimiter)
        header=reader.fieldnames or []
        if len(set(header))!=len(header): raise ValueError('差异结果列名重复')
        aliases={'KO':['KO','ko'],'effect':['effect','log2FC','log2FoldChange'],
                 'q':['q','padj','qvalue'],'P':['P','pvalue'],
                 'comparison':['comparison'],'numerator':['numerator'],'denominator':['denominator']}
        keys={key:next((name for name in names if name in header),None) for key,names in aliases.items()}
        if any(not keys[key] for key in ('KO','effect','q','comparison','numerator','denominator')):
            raise ValueError('差异表须含 KO,comparison,numerator,denominator,effect,q；q须为已校正值，不能用原始P冒充')
        if keys['effect']!='effect' and effect_type!='log2FC': raise ValueError('log2FC列与所选效应尺度不一致')
        rows=[];seen=set();comparisons={}
        for line,row in enumerate(reader,2):
            if None in row or any(value is None for value in row.values()): raise ValueError(f'差异表第{line}行列数异常')
            item={key:row[keys[key]] for key in ('KO','effect','q','comparison','numerator','denominator')}
            item['KO']=item['KO'].removeprefix('ko:')
            if not re.fullmatch(r'K\d{5}',item['KO']): raise ValueError('无效KO：'+item['KO'])
            for key in ('comparison','numerator','denominator'):
                value=item[key]
                if not value or value!=value.strip() or any(ord(c)<32 for c in value): raise ValueError('比较名称/方向缺失或含控制字符')
            if item['numerator']==item['denominator']: raise ValueError('比较的两组不能相同')
            pair=(item['numerator'],item['denominator'])
            if item['comparison'] in comparisons and comparisons[item['comparison']]!=pair: raise ValueError('同名比较的方向不一致')
            comparisons[item['comparison']]=pair
            identity=(item['KO'],item['comparison'])
            if identity in seen: raise ValueError('同一比较的KO重复，不能猜测选取某条结果')
            seen.add(identity)
            missing=lambda value:value.strip().lower() in ('','na','nan','null')
            if missing(item['q']) or missing(item['effect']):
                item.update(effect='NA',q='NA',Status='not_estimable',Reason='imported_missing_statistic')
            else:
                effect=float(item['effect']);q=float(item['q'])
                if not math.isfinite(effect) or not math.isfinite(q) or not 0<=q<=1: raise ValueError('效应或q值非有限/越界')
                item.update(effect=effect,q=q,Status='ok',Reason='imported; original model and FDR family supplied by user')
            item['P']=row.get(keys['P'],'NA') if keys['P'] else 'NA'
            rows.append(item)
    if not rows: raise ValueError('差异表为空')
    columns=['KO','comparison','numerator','denominator','effect','P','q','Status','Reason']
    with Path(destination).open('w',encoding='utf8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns,delimiter='\t');writer.writeheader();writer.writerows(rows)
    return [{'id':key,'numerator':pair[0],'denominator':pair[1],'label':pair[0]+' − '+pair[1]} for key,pair in comparisons.items()]


def _seal(run):
    files={p.relative_to(run).as_posix():_hash(p) for p in (run/'02_statistics').rglob('*') if p.is_file()}
    _json(run/'ko_statistics_complete.json',{'files':files,'statistics_refit':False})


def _verify(run):
    marker=run/'ko_statistics_complete.json'
    if not marker.is_file(): raise ValueError('没有已完成的KO差异统计；请先计算或导入差异结果')
    files=json.loads(marker.read_text())['files']
    if not files: raise ValueError('KO统计文件指纹缺失')
    for name,digest in files.items():
        path=run/name
        if not path.resolve().is_relative_to(run) or not path.is_file() or _hash(path)!=digest:
            raise ValueError('KO统计文件已改变，拒绝与原差异结果混用：'+name)


def _map(job,config,run,api):
    map_id=normalize_map_id(config.get('ko_map_id') or 'ko01100')
    folder=run/'01_map';folder.mkdir(exist_ok=True)
    target=folder/(map_id+'.xml')
    if str(config.get('ko_kgml') or '').strip():
        source=api.path(config['ko_kgml']);shutil.copyfile(source,target)
        provenance={'source':str(source),'source_type':'user_local'}
    else:
        url='https://rest.kegg.jp/get/'+map_id+'/kgml'
        try:
            request=urllib.request.Request(url,headers={'User-Agent':'16s-workstation/1'})
            with urllib.request.urlopen(request,timeout=30) as response:
                if urllib.parse.urlparse(response.url).scheme!='https': raise ValueError('KGML下载未保持HTTPS')
                data=response.read(20*1024*1024+1)
            if len(data)>20*1024*1024: raise ValueError('KGML超过20MB限制')
            target.write_bytes(data)
        except Exception as error:
            raise ValueError('KO差异结果已保存；无法下载KEGG KGML。请选择本地对应KGML，再“从已保存差异结果绘图”，不需重跑PICRUSt2或统计。原因：'+str(error)) from error
        provenance={'source':url,'source_type':'KEGG public REST','downloaded_at':time.time()}
    geometry=write_geometry(target,folder/'geometry.json',map_id)
    provenance.update(sha256=_hash(target),map_id=map_id,coordinates='original KGML; no inferred network edges')
    _json(folder/'map_source.json',provenance)
    return geometry,folder/'geometry.json'


def _render(job,config,analysis,run,api):
    _verify(analysis)
    options=json.loads((analysis/'ko_options.json').read_text())
    comparisons=json.loads((analysis/'02_statistics/comparisons.json').read_text())
    if isinstance(comparisons,dict): comparisons=[comparisons]
    comparison=config.get('ko_comparison') or comparisons[0]['id']
    selected=next((item for item in comparisons if item['id']==comparison),None)
    if not selected: raise ValueError('请选择该统计任务中实际存在的比较')
    plotting=resolve_plot_settings(config,fallback=options)
    source=api.approved_run(config.get('ko_map_source_run') or str(analysis))
    cache=source/'ko_map_complete.json'
    desired_map=config.get('ko_map_id') or options.get('ko_map_id','ko01100')
    desired_file=config.get('ko_kgml') or options.get('ko_kgml','')
    cached=json.loads(cache.read_text()) if cache.is_file() else {}
    if cached and normalize_map_id(desired_map)==cached['map_id'] and desired_file==cached.get('original_kgml',''):
        for name,digest in cached['files'].items():
            candidate=source/name
            if not candidate.resolve().is_relative_to(source) or not candidate.is_file() or _hash(candidate)!=digest: raise ValueError('保存的代谢图坐标已改变，不能复用')
        if source!=run:
            shutil.copytree(source/'01_map',run/'01_map')
        path=run/'01_map/geometry.json';geometry=json.loads(path.read_text())
    else:
        geometry,path=_map(job,{**options,**config,'ko_map_id':desired_map,'ko_kgml':desired_file},run,api)
    _json(run/'ko_map_complete.json',{'map_id':geometry['map_id'],'original_kgml':desired_file,
          'files':{p.relative_to(run).as_posix():_hash(p) for p in (run/'01_map').iterdir() if p.is_file()}})
    effective=dict(options,geometry_path=str(path),hits_path=str(analysis/'02_statistics/ko_differential.tsv'),
                   comparison=comparison,comparison_label=selected['label'],plot_settings=plotting,
                   q_threshold=plotting['metabolic_map']['q_threshold'],effect_threshold=plotting['metabolic_map']['effect_threshold'])
    (run/'03_figures').mkdir(exist_ok=True)
    _json(run/'ko_render_options.json',effective)
    runtime=dict(options['runtime'],**{key:config[key] for key in ('conda','r_env','font') if config.get(key)})
    effective['font']=runtime.get('font','DejaVu Serif');_json(run/'ko_render_options.json',effective)
    api.begin_stage(job,'03_figures','差异KO映射到真实代谢图坐标',{'comparison':comparison,'map_id':geometry['map_id'],'statistics_refit':False})
    api.command(job,api.r_command(runtime,str(api.ROOT/'ko_map_render.R'),'--render',str(run/'ko_render_options.json'),str(run)),run,process_env=api.r_execution_env(runtime))
    if not (run/'03_figures/metabolic_map.png').is_file(): raise ValueError('代谢位置图未生成，请保留日志')
    summary_file=run/'03_figures/render_summary.json'
    produced=summary_file if summary_file.is_file() else run/'render_summary.json'
    summary=json.loads(produced.read_text()) if produced.is_file() else {}
    summary.update(source='ko_network',analysis_run=str(analysis),plot_settings=plotting,font=effective['font'],
                   comparison=comparison,map_id=geometry['map_id'],ko_kgml=desired_file)
    _json(summary_file,summary)
    # iPath export is a local text artifact only. No sample/KO data are uploaded.
    output=run/'02_mapping';output.mkdir(exist_ok=True)
    with (analysis/'02_statistics/ko_differential.tsv').open(encoding='utf8') as stream,(output/'ipath_selection.txt').open('w',encoding='utf8') as target:
        for row in csv.DictReader(stream,delimiter='\t'):
            if row['comparison']!=comparison: continue
            try: effect,q=float(row['effect']),float(row['q'])
            except ValueError: continue
            if math.isfinite(effect) and math.isfinite(q) and q<=effective['q_threshold'] and abs(effect)>=effective['effect_threshold'] and effect!=0:
                target.write(row['KO']+' '+plotting['metabolic_map']['up_color' if effect>0 else 'down_color']+' W3\n')
    save_plot_settings(run,plotting)
    api.finish_stage(job,summary,['只改变展示阈值与样式，不重新拟合统计。混合方向的多KO图元单独标示。','代谢图位置不是代谢流、基因表达或实测功能；预测不确定性沿用PICRUSt2/NSTI等报告。'])
    return {'ko_network_run':str(analysis),'ko_comparisons':comparisons,'ko_map_summary':summary,
            'ko_ipath_selection':str(output/'ipath_selection.txt'),'statistics_refit':False}


def run_ko_network(job,config,api):
    source=config.get('ko_source','functional');parent=None;references={}
    if source=='functional':
        parent=api.approved_run(config.get('functional_run',''))
        if not (parent/'functional_complete.json').is_file(): raise ValueError('请选择已完成的PICRUSt2功能预测任务，不是ASV预处理目录')
        source_parameters=json.loads((parent/'parameters.json').read_text());base=source_parameters['config']
        table=parent/'03_tables/KO_abundance.tsv';metadata=config.get('ko_metadata') or ''
        samples=base['samples'];output=parent/'06_KO_network'
    else:
        base=dict(config);output=api.path(config.get('output',''));samples=[]
        if source not in ('table','differential'): raise ValueError('KO输入来源无效')
        table=api.path(config.get('ko_table','')) if source=='table' else None
        metadata=config.get('ko_metadata','')
        if source=='table': samples=read_metadata(api.path(metadata))
    output.mkdir(parents=True,exist_ok=True)
    run=output/('run_ko_network_'+time.strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:8]);run.mkdir();api.marker(run)
    job.update(run=str(run),log=str(run/'workflow.log'))
    inputs=run/'00_inputs';inputs.mkdir()
    runtime=dict(base,**{key:config[key] for key in ('conda','r_env','font') if config.get(key)})
    runtime.setdefault('r_env','16s-workstation-r');runtime.setdefault('font','DejaVu Serif')
    plotting=resolve_plot_settings(config)
    params={'extension':'ko-network','config':dict(config,samples=samples),'extension_parameters':config}
    if parent:
        params.update(functional_run=str(parent),preprocess=source_parameters['preprocess'])
    _json(run/'parameters.json',params)
    options={'runtime':runtime,'plot_settings':plotting,'ko_map_id':config.get('ko_map_id') or 'ko01100',
             'ko_kgml':config.get('ko_kgml',''),'effect_label':'log2 CLR difference (numerator minus denominator)',
             'effect_metadata':{'prediction':'PICRUSt2 inferred abundance; not measured activity',
                                'uncertainty':'P/q conditional on predicted table; phylogenetic imputation uncertainty is not propagated'}}
    _json(run/'ko_options.json',options)
    api.begin_stage(job,'00_inputs','KO差异输入与来源保存',{'source':source})
    if source=='differential':
        original,record=_copy(api.path(config.get('ko_differential_file','')),inputs,'differential_original');references['differential']=record
        stats=run/'02_statistics';stats.mkdir()
        effect_type=config.get('ko_effect_type') or 'log2FC'
        comparisons=import_differential(original,stats/'ko_differential.tsv',effect_type)
        _json(stats/'comparisons.json',comparisons)
        _json(stats/'transformation_report.json',{'method':'imported differential results','effect_type':effect_type,
               'FDR_family':'provided by original analysis; no re-adjustment','statistical_parameters_not_changed':True})
        options['effect_label']=effect_type+' (numerator minus denominator)'
    else:
        original,record=_copy(table,inputs,'KO_abundance_original');references['abundance']=record
        if original.suffix.lower()=='.csv':
            converted=inputs/'converted.tsv'
            with original.open(encoding='utf-8-sig',newline='') as stream,converted.open('w',encoding='utf8',newline='') as dest:
                csv.writer(dest,delimiter='\t').writerows(csv.reader(stream))
            original=converted
        abundances=read_abundance(original,[sample['SampleID'] for sample in samples],ko=True)
        canonical=inputs/'KO_abundance.tsv';write_table(canonical,abundances)
        meta=inputs/'metadata.csv'
        if metadata:
            copied,record=_copy(api.path(metadata),inputs,'metadata_original');references['metadata']=record
            prepare_metadata(copied,samples,abundances['samples'],meta)
        else:
            with meta.open('w',encoding='utf8',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=['SampleID','Group']);writer.writeheader();writer.writerows({key:sample[key] for key in ('SampleID','Group')} for sample in samples)
        covariates=[name.strip() for name in str(config.get('ko_covariates') or '').replace('，',',').split(',') if name.strip()]
        if len(set(covariates))!=len(covariates) or set(covariates)&{'SampleID','Group'}: raise ValueError('协变量不能重复或包含SampleID/Group')
        categorical={name.strip() for name in str(config.get('ko_categorical_covariates') or '').replace('，',',').split(',') if name.strip()}
        if categorical-set(covariates): raise ValueError('分类协变量须同时列在调整协变量中')
        style=resolve_group_style(config,samples)
        model={'counts_path':str(canonical),'metadata_path':str(meta),'covariates':[{'name':name,'type':'categorical'} if name in categorical else name for name in covariates],
               'group_order':style['group_order'],'feature_prevalence':_number(config,'ko_prevalence',10,0,100)/100,
               'pseudocount':_number(config,'ko_pseudocount',1e-6,1e-12,.1)}
        _json(inputs/'model_options.json',model)
        api.finish_stage(job,{'sources':references},['预测KO可含小数，不四舍五入为测序计数；独立样本、CLR线性模型为探索性关联。'])
        api.begin_stage(job,'01_model','预测KO的CLR线性模型与全部两两比较',model)
        api.command(job,api.r_command(runtime,str(api.ROOT/'ko_differential.R'),'--run',str(inputs/'model_options.json'),str(run)),run,process_env=api.r_execution_env(runtime))
        options['effect_metadata']['method']='TSS retained KO + pseudocount + log2 CLR + lm; exploratory; BH all tested KO×contrasts'
    _json(inputs/'source_report.json',references)
    _json(run/'ko_options.json',options)
    _seal(run)
    api.finish_stage(job,{'statistics_complete':True},['所有统计表及方向另存；展示不重新检验。'])
    result=_render(job,config,run,run,api)
    (run/'ko_network_complete.json').write_text('{}')
    return result


def render_saved_ko_network(job,config,api):
    analysis=api.approved_run(config.get('ko_network_run',''));_verify(analysis)
    options=json.loads((analysis/'ko_options.json').read_text())
    folder=analysis/'05_replots';folder.mkdir(exist_ok=True)
    run=folder/('map_'+time.strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:8]);run.mkdir();api.marker(run)
    job.update(run=str(run),log=str(run/'workflow.log'))
    previous=json.loads((analysis/'parameters.json').read_text())
    params={'extension':'ko-network-render','ko_network_run':str(analysis),'config':previous['config'],'statistics_refit':False}
    if previous.get('preprocess'): params['preprocess']=previous['preprocess']
    _json(run/'parameters.json',params)
    return _render(job,config,analysis,run,api)
