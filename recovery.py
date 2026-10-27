"""Durable task records and conservative, independent retries of failed analyses."""
import hashlib
import json
import os
import shutil
import tempfile
import time
import uuid
from pathlib import Path

from run_layout import artifact_path, initialize_layout

PATH_FIELDS = ('quality_run', 'run', 'diversity_run', 'differential_run', 'lefse_run',
               'abundance_run', 'advanced_run', 'result_run', 'figure_run')
IDENTITY_KEYS = ('input', 'metadata', 'mode', 'data_state', 'primer_mode', 'forward',
                 'reverse', 'env', 'conda', 'input_layout', 'batch_id')
BIOLOGICAL_KEYS = ('trunc_f', 'trunc_r', 'trim_f', 'trim_r', 'amplicon', 'min_frequency',
                   'classifier', 'max_depth', 'parameter_mode')
ABUNDANCE_STAGES = ('04_taxonomy_filter', '05_phylogeny', '06_rarefaction')


def _read(filename, default=None):
    return json.loads(Path(filename).read_text(encoding='utf8')) if Path(filename).is_file() else (default or {})


def _save(filename, data):
    filename = Path(filename)
    temporary = filename.with_name(filename.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf8')
        os.replace(temporary, filename)
    finally:
        temporary.unlink(missing_ok=True)


def _hash(filename):
    digest = hashlib.sha256()
    with Path(filename).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _source_stat(filename):
    stat = Path(filename).stat()
    return {'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns, 'ctime_ns': stat.st_ctime_ns,
            'device': stat.st_dev, 'inode': stat.st_ino}


def _normalized(config):
    config = dict(config)
    config.setdefault('env', 'qiime2-amplicon-2025.7')
    config.setdefault('input_layout', 'unmerged')
    config['data_state'] = config.get('data_state') or ('clean' if config.get('primer_mode') == 'removed' else 'raw')
    config['primer_mode'] = 'removed' if config['data_state'] == 'clean' else ('auto' if config.get('primer_mode') in ('unknown', 'removed', None) else config['primer_mode'])
    return config


def _files(run, paths):
    result = {}
    for path in paths:
        path = Path(path)
        candidates = sorted(path.rglob('*')) if path.is_dir() else [path]
        for filename in candidates:
            if filename.is_symlink() or not filename.resolve().is_relative_to(run):
                raise ValueError('检查点包含链接或批次以外的文件，不能复用')
            if filename.is_file():
                result[str(filename.relative_to(run))] = _hash(filename)
    return result


def _done(run, stage):
    return _read(run / stage / 'report.json').get('status') == 'done'


def seal_checkpoint(run, kind):
    """Seal only completed boundaries, before any dependent stage is attempted."""
    run = Path(run)
    saved = _read(run / 'parameters.json')
    if not saved.get('samples') or not (run / 'input_inventory.json').is_file():
        return
    if kind == 'quality':
        if not _done(run, '01_import'):
            return
        stages = ['01_import']
        logical = ['demux.qza', 'demux.qzv', 'quality-raw']
        if saved.get('primer_mode') != 'removed':
            if not _done(run, '02_primers'):
                return
            stages.append('02_primers')
            logical += ['trimmed.qza', 'trimmed.qzv', 'quality-trimmed']
        paths = [artifact_path(run, name) for name in logical]
        if any(not p.exists() or (p.is_file() and not p.stat().st_size) for p in paths):
            return
    elif kind == 'abundance':
        stages = list(ABUNDANCE_STAGES)
        if not all(_done(run, stage) for stage in stages):
            return
        paths = [run / stage for stage in stages]
        for logical in ('filtered-table.qza', 'taxonomy.qza', 'rooted-tree.qza',
                        'taxa-barplot.qzv', 'alpha-rarefaction.qzv', 'exported-table/feature-table.tsv'):
            filename = artifact_path(run, logical)
            if not filename.is_file() or not filename.stat().st_size:
                return
    else:
        raise ValueError('未知检查点')
    inventory = _read(run / 'input_inventory.json')
    checkpoint = {'format': 1, 'kind': kind, 'config': _normalized(saved),
                  'inventory': inventory, 'source_stats': {entry['path']: _source_stat(entry['path']) for entry in inventory},
                  'stages': stages,
                  'files': _files(run, paths)}
    if kind == 'abundance' and saved.get('classifier'):
        checkpoint['source_stats'][saved['classifier']] = _source_stat(saved['classifier'])
    _save(run / (kind + '_checkpoint.json'), checkpoint)


def ensure_task_run(job, config, api):
    """Keep early validation failures inspectable even before an analysis run exists."""
    if job.get('run'):
        return
    try:
        output = api.path(config['output'])
        raw = api.path(config['input']) if config.get('input') else None
        if raw and (output == raw or output in raw.parents):
            raise ValueError('不能向原始输入目录写入任务记录')
        output.mkdir(parents=True, exist_ok=True)
    except (KeyError, OSError, ValueError):
        output = Path(tempfile.gettempdir()) / '16s-workstation-tasks' / hashlib.sha256(str(api.ROOT).encode()).hexdigest()[:12]
        output.mkdir(parents=True, exist_ok=True)
    run = output / ('task_' + (job.get('id') or uuid.uuid4().hex))
    run.mkdir(exist_ok=True)
    api.marker(run)
    _save(run / 'parameters.json', {'config': config})
    previous = Path(job.get('log', run / 'workflow.log'))
    if previous.is_file():
        shutil.copyfile(previous, run / 'workflow.log')
    job.update(run=str(run), log=str(run / 'workflow.log'))


def persist_task(job):
    if not job.get('run') or '_task' not in job:
        return
    run = Path(job['run'])
    record = dict(job['_task'], id=job.get('id'), run=str(run),
                  status=job.get('status'), step=job.get('step'), error=job.get('error', ''),
                  stage=job.get('stage_folder'), updated_at=time.time())
    saved = _read(run / 'parameters.json')
    if 'source_files' not in record and saved.get('preprocess'):
        parent = Path(saved['preprocess']).resolve()
        paths = [artifact_path(parent, name) for name in
                 ('filtered-table.qza', 'taxonomy.qza', 'filtered-rep-seqs.qza', 'rooted-tree.qza',
                  'exported-table/feature-table.tsv', 'exported-taxonomy/taxonomy.tsv', 'metadata.tsv')]
        record['source_run'] = str(parent)
        record['source_files'] = _files(parent, [p for p in paths if p.is_file()])
        job['_task'].update(source_run=record['source_run'], source_files=record['source_files'])
    if 'cached_files' not in record and (run / 'qiime_statistics_complete.json').is_file():
        record['cached_files'] = _files(run, [run / name for name in
                                             ('00_inputs', '01_core_metrics', '02_exports', '03_statistics')])
        job['_task']['cached_files'] = record['cached_files']
    _save(run / 'task_state.json', record)


def _verify_files(run, files):
    if not files:
        raise ValueError('缺少已完成检查点的文件指纹，不能复用')
    for relative, digest in files.items():
        filename = run / relative
        if not filename.resolve().is_relative_to(run) or filename.is_symlink() or not filename.is_file() or _hash(filename) != digest:
            raise ValueError('检查点文件已改变或缺失，不能复用：' + relative)


def _verify_checkpoint(run, kind, config, api):
    checkpoint = _read(run / (kind + '_checkpoint.json'))
    if checkpoint.get('format') != 1 or checkpoint.get('kind') != kind:
        raise ValueError('没有保存可核验的检查点指纹')
    previous, current = checkpoint['config'], _normalized(config)
    keys = IDENTITY_KEYS + (BIOLOGICAL_KEYS if kind == 'abundance' else ())
    if any(current.get(key, '') != previous.get(key, '') for key in keys):
        raise ValueError('输入、引物、批次或检查点前的参数已改变；请新建任务，不能复用此检查点')
    for entry in checkpoint['inventory']:
        stat = Path(entry['path']).stat()
        if stat.st_size != entry['size'] or stat.st_mtime_ns != entry['mtime_ns']:
            raise ValueError('原始文件已改变，不能复用检查点')
    for filename, signature in checkpoint.get('source_stats', {}).items():
        if _source_stat(filename) != signature:
            raise ValueError('原始文件身份或修改记录已改变，不能复用检查点')
    samples = api.discover(current)
    if samples != previous['samples']:
        raise ValueError('样本、分组或测序文件映射已改变，不能复用检查点')
    if not all(_done(run, stage) for stage in checkpoint['stages']):
        raise ValueError('检查点对应阶段没有成功完成')
    _verify_files(run, checkpoint['files'])
    return checkpoint


def recovery_plan(selected, api, active_jobs=()):
    """Read-only inspection; a running record from another process is interrupted."""
    run = api.approved_run(selected)
    saved = _read(run / 'parameters.json')
    state = _read(run / 'task_state.json') or _read(run / 'job_status.json')
    action = state.get('action') or saved.get('extension') or ('preprocess' if (run / 'preprocess_started.json').exists() else 'quality')
    effective = saved.get('config', saved)
    config = dict(effective, **state.get('config', {}), **saved.get('extension_parameters', {}))
    if action in ('preprocess', 'quality', 'resume-abundance'):
        config.update(effective)
    for key in ('depth', 'permutations'):
        if key in saved:
            config[key] = saved[key]
    status = state.get('status', 'interrupted')
    if status == 'running' and not any(j.get('status') == 'running' and j.get('run') == str(run) for j in active_jobs):
        status = 'interrupted'
    fields = {key: '' for key in PATH_FIELDS}
    stages = [_read(p) for p in sorted(run.glob('*/report.json'))]
    plan = {'run': str(run), 'action': action, 'recoverable': False, 'config': config,
            'fields': fields, 'boundary': '', 'warnings': [], 'status': status,
            'error': state.get('error', ''), 'stages': stages, 'strategy': ''}
    if status not in ('failed', 'cancelled', 'interrupted'):
        plan['boundary'] = '任务仍在运行或已经完成，无需失败恢复'
        return plan
    if action in ('inspect', 'preflight', 'recommend') or (run / '00_checks').exists():
        plan['boundary'] = '检查或参数推荐没有可续跑的分析阶段；请修正设置后重新执行相应检查'
        return plan
    if action == 'external-inputs':
        plan.update(strategy='restart-external', boundary='按原请求的外部产物及分组文件，在独立新目录重新准备分析输入；来源文件保持只读')
    elif action in ('quality', 'preprocess', 'resume-abundance'):
        if (run / 'preprocess_complete.json').exists():
            plan['boundary'] = '预处理已完成；请打开历史批次继续后续分析'
            return plan
        if action != 'quality' and (run / 'abundance_checkpoint.json').is_file():
            plan.update(strategy='abundance', boundary='复制已完成且指纹一致的预处理结果到新目录，仅重新导出丰度；原目录保持不变')
        elif action != 'quality' and (run / 'quality_checkpoint.json').is_file():
            plan.update(strategy='quality', boundary='复制已完成且指纹一致的导入及引物检查点到新目录，从 DADA2 重新运行；不重新导入 FASTQ')
        else:
            plan.update(strategy='restart', boundary='没有可核验的完整检查点；在新目录从 FASTQ 导入重新执行此任务')
            plan['warnings'].append('此前不完整的阶段与输出不会标记为完成，也不会跳过命令')
    else:
        parent = api.approved_run(saved.get('preprocess') or config.get('run', ''))
        parent_parameters = _read(parent / 'parameters.json')
        allowed = parent_parameters.get('supported_modules', [])
        external_complete = (parent / 'external_inputs_complete.json').is_file() and (action.removesuffix('-render') in allowed or action == 'figure-render')
        if not (parent / 'preprocess_complete.json').is_file() and not external_complete:
            plan['boundary'] = '来源预处理尚未完成，不能恢复下游分析'
            return plan
        fields.update(run=str(parent), result_run=str(parent), figure_run=str(parent))
        config['run'] = str(parent)
        plan.update(strategy='analysis', boundary='从已完成预处理结果在独立子目录重新执行本项分析；不重新处理 FASTQ，随机抽平结果可能变化')
        if action == 'diversity' and (run / 'qiime_statistics_complete.json').is_file():
            plan.update(strategy='diversity-render', boundary='从已保存的本次抽平表、指标及检验结果重绘到独立目录；若修改抽平深度或置换次数，则重新执行本项多样性分析')
            config['diversity_run'] = str(run)
            fields['diversity_run'] = str(run)
        elif action.endswith('-render'):
            plan.update(strategy='render', boundary='使用原请求的已保存分析缓存，在独立目录重绘；不重新拟合或执行统计检验')
        if action not in ('diversity', 'abundance', 'differential', 'lefse', 'functional', 'advanced',
                          'diversity-render', 'differential-render', 'advanced-render', 'lefse-render', 'figure-render'):
            plan['boundary'] = '此任务没有已实现的安全恢复边界'
            return plan
    plan['recoverable'] = True
    return plan


def _copy(run, destination, relative):
    source = run / relative
    if not source.exists():
        return
    candidates = sorted(source.rglob('*')) if source.is_dir() else [source]
    for filename in candidates:
        if filename.is_symlink() or not filename.resolve().is_relative_to(run):
            raise ValueError('检查点包含链接或批次以外的文件，不能复用')
        if not filename.is_file():
            continue
        name = filename.relative_to(run)
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(run / name, target)


def _fork(run, config, kind, checkpoint, api):
    output, raw = api.path(config['output']), api.path(config['input'])
    if output == raw or output in raw.parents:
        raise ValueError('恢复结果目录不能是原始目录本身或其父目录')
    if output == run or output.is_relative_to(run):
        raise ValueError('恢复必须写入原失败目录以外的新批次目录')
    output.mkdir(parents=True, exist_ok=True)
    destination = output / ('run_' + time.strftime('%Y%m%d_%H%M%S') + '_' + uuid.uuid4().hex[:8])
    destination.mkdir()
    initialize_layout(destination, 'preprocess')
    api.marker(destination)
    for filename in ('parameters.json', 'input_inventory.json', 'group_style.json', 'plot_settings.json',
                     'parameter_recommendation.json', 'depth_recommendation.json', kind + '_checkpoint.json'):
        _copy(run, destination, filename)
    for name in checkpoint['files']:
        _copy(run, destination, name)
    # An abundance-only fork needs the filtered outputs, metadata and saved lineage.
    # Earlier FASTQ exports and DADA2 partial outputs are deliberately not copied.
    if kind == 'abundance':
        _copy(run, destination, str(artifact_path(run, 'metadata.tsv').relative_to(run)))
        _copy(run, destination, str(artifact_path(run, 'manifest.tsv').relative_to(run)))
    for directory in (ABUNDANCE_STAGES if kind == 'abundance' else ('01_import', '02_primers')):
        for filename in ('report.json', '说明.md'):
            _copy(run, destination, directory + '/' + filename)
    if kind == 'abundance':
        saved = _read(destination / 'parameters.json')
        saved.update({key: config[key] for key in ('output', 'threads', 'r_env', 'font') if key in config})
        _save(destination / 'parameters.json', saved)
    _save(destination / 'recovery_source.json', {'run': str(run), 'boundary': kind, 'files': checkpoint['files']})
    if kind == 'quality':
        _save(destination / 'quality_complete.json', {})
    else:
        _save(destination / 'preprocess_started.json', {})
    return destination


def prepare_retry(job, request, api):
    """Validate identity again immediately before copying or launching any commands."""
    plan = recovery_plan(request.get('recovery_run', ''), api, getattr(api, 'JOBS', {}).values())
    job['recovery_action'] = plan['action']
    job['_task']['action'] = plan['action']
    if not plan['recoverable']:
        raise ValueError(plan['boundary'])
    run = Path(plan['run'])
    config = dict(plan['config'], **{k: v for k, v in request.items() if k != 'recovery_run'})
    # Result path widgets may be cleared by the UI; they cannot redirect recovery.
    for key in PATH_FIELDS:
        config[key] = plan['config'].get(key, '')
    state = _read(run / 'task_state.json')
    if state.get('source_run'):
        _verify_files(Path(state['source_run']), state.get('source_files'))
    strategy = plan['strategy']
    if strategy == 'restart-external':
        # The main form intentionally has no external-file widgets. Its blank
        # values cannot erase or redirect the explicit saved import request.
        for key, value in plan['config'].items():
            if key.startswith('external_') or key == 'metadata':
                config[key] = value
    if strategy in ('quality', 'abundance'):
        checkpoint = _verify_checkpoint(run, strategy, config, api)
        destination = _fork(run, config, strategy, checkpoint, api)
        job.update(run=str(destination), log=str(destination / 'workflow.log'))
        action = 'preprocess' if strategy == 'quality' else 'resume-abundance'
        config.update(quality_run=str(destination) if strategy == 'quality' else '',
                      result_run=str(destination), run='')
    else:
        action = 'diversity-render' if strategy == 'diversity-render' else plan['action']
        if strategy == 'diversity-render':
            if any(str(config.get(key, '')) != str(plan['config'].get(key, '')) for key in ('depth', 'permutations')):
                action, strategy = 'diversity', 'analysis'
            else:
                _verify_files(run, state.get('cached_files'))
        if strategy == 'restart':
            action = 'quality' if plan['action'] == 'quality' else 'preprocess'
            config.update(quality_run='', run='', result_run='')
        if strategy in ('render', 'diversity-render'):
            for key in ('diversity_run', 'differential_run', 'advanced_run', 'lefse_run', 'image', 'plot_target'):
                if key in plan['config']:
                    config[key] = plan['config'][key]
    job['recovery_action'] = action
    job['_task'].update(action=action, config=config, recovery_source=str(run))
    return action, config
