"""Read-only history restoration; completed outputs must share their saved parent."""
import json
from datetime import datetime, timezone, timedelta
from pathlib import Path


def _parameters(run):
    return json.loads((run/'parameters.json').read_text(encoding='utf8'))


def list_batches(root, api):
    root=api.path(root)
    if not root.is_dir(): return {'runs':[], 'limited':False}
    runs=[]; visited=0
    for child in sorted(root.iterdir(),key=lambda p:p.name,reverse=True):
        visited+=1
        if visited>2000: break
        if not child.name.startswith(('run_','external_','task_')): continue
        try:
            run=api.approved_run(child); p=_parameters(run); c=p.get('config',p)
            if not isinstance(c.get('samples'),list): continue
            state_file=run/'task_state.json'
            state=json.loads(state_file.read_text(encoding='utf8')) if state_file.is_file() else {}
            runs.append({'run':str(run),'batch_id':c.get('batch_id',str(run)),
                         'batch_name':c.get('batch_name') or run.name,
                         'created':datetime.fromtimestamp(run.stat().st_mtime,timezone(timedelta(hours=8))).isoformat(timespec='seconds'),
                         'sample_count':len(c['samples']),
                         'quality_complete':(run/'quality_complete.json').is_file(),
                         'preprocess_complete':(run/'preprocess_complete.json').is_file(),
                         'external_inputs_complete':(run/'external_inputs_complete.json').is_file(),
                         'task_status':state.get('status','')})
        except (OSError,ValueError,TypeError): continue
        if len(runs)>=100: break
    return {'runs':runs,'limited':visited>2000 or len(runs)>=100}


def restore_batch(selected, api):
    run=api.approved_run(selected); p=_parameters(run)
    if p.get('preprocess'):
        run=api.approved_run(p['preprocess']); p=_parameters(run)
    config=dict(p.get('config',p))
    if not config.get('samples'): raise ValueError('请选择质量检查或预处理批次目录，而非单个图片/检查目录')
    fields={key:'' for key in ('quality_run','run','diversity_run','differential_run','lefse_run','abundance_run','advanced_run','functional_run','ko_network_run','result_run','figure_run')}
    fields['result_run']=str(run); fields['figure_run']=str(run)
    if (run/'quality_complete.json').is_file() and not (run/'preprocess_started.json').exists(): fields['quality_run']=str(run)
    if (run/'preprocess_complete.json').is_file() or (run/'external_inputs_complete.json').is_file(): fields['run']=str(run)
    for directory in (run/'07_diversity',run/'11_extensions'):
        if not directory.is_dir(): continue
        for child in sorted(directory.iterdir(),key=lambda p:p.name):
            try:
                api.approved_run(child); saved=_parameters(child)
                if Path(saved.get('preprocess','')).resolve()!=run: continue
                parent_id=config.get('batch_id'); child_id=saved.get('config',{}).get('batch_id')
                if parent_id and child_id and parent_id!=child_id: continue
                if (child/'qiime_statistics_complete.json').is_file(): fields['diversity_run']=str(child)
                elif (child/'lefse_complete.json').is_file(): fields['lefse_run']=str(child)
                elif (child/'differential_options.json').is_file() and (child/'differential_complete.json').is_file(): fields['differential_run']=str(child)
                elif (child/'10_abundance/manifest.json').is_file() and (child/'abundance_complete.json').is_file(): fields['abundance_run']=str(child)
                elif (child/'advanced_complete.json').is_file(): fields['advanced_run']=str(child)
                elif (child/'functional_complete.json').is_file():
                    fields['functional_run']=str(child)
                    network_folder=child/'06_KO_network'
                    if network_folder.is_dir():
                        for network in sorted(network_folder.iterdir(),key=lambda p:p.name):
                            if ((network/'ko_statistics_complete.json').is_file()
                                    and Path(_parameters(network).get('functional_run','')).resolve()==child):
                                fields['ko_network_run']=str(api.approved_run(network))
            except (OSError,ValueError,TypeError): continue
    config.setdefault('data_state','clean' if config.get('primer_mode')=='removed' else 'raw')
    config.setdefault('input_layout','unmerged')
    for filename,key in [('group_style.json',None),('plot_settings.json','plot_settings')]:
        source=run/filename
        if source.is_file():
            data=json.loads(source.read_text(encoding='utf8'))
            if key: config[key]=data
            else: config.update(data)
    # Preserve the actual saved source configuration, not paths from unrelated jobs.
    return {'config':config,'fields':fields,'batch_id':config.get('batch_id') or str(run),
            'batch_name':config.get('batch_name') or run.name}
