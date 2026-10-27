"""Advanced ecology from saved unrarefied tables; original runs remain read-only."""
import csv
import hashlib
import json
import math
import re
import time
import uuid
from pathlib import Path

from abundance_tables import LEVELS, UNCERTAIN_TAXON, read_table, write_table_pair
from group_style import resolve_group_style, save_group_style
from plot_settings import resolve_plot_settings, save_plot_settings
from differential_analysis import plot_options


DEFAULT_ENV = '16s-workstation-advanced'
MODULES = {'core', 'ordination', 'association'}


def _number(config, name, default, low, high, integer=False):
    raw = config.get(name, default)
    if raw is None or str(raw).strip() == '':
        raw = default
    try:
        if isinstance(raw, bool):
            raise ValueError()
        value = float(raw)
        if not math.isfinite(value) or not low <= value <= high or (integer and value != int(value)):
            raise ValueError()
    except (TypeError, ValueError, OverflowError):
        raise ValueError(f'{name} 必须在 {low}–{high} 之间' + ('且为整数' if integer else '')) from None
    return int(value) if integer else value


def _columns(value):
    if isinstance(value, str):
        values = value.replace('，', ',').split(',')
    elif isinstance(value, list):
        values = value
    else:
        raise ValueError('分析变量须为逗号分隔的列名')
    result = []
    for item in values:
        if not isinstance(item, str):
            raise ValueError('分析变量列名无效')
        name = item.strip()
        if name and name not in result:
            result.append(name)
    if any(any(c in name for c in '\r\n\t\x00') or name == 'SampleID' for name in result):
        raise ValueError('不能将 SampleID 或控制字符作为解释变量')
    return result


def prepare_metadata(filename, samples, observed, destination):
    """Join auxiliary traits by literal SampleID, never by row position."""
    groups = {s['SampleID']: s['Group'] for s in samples}
    if len(groups) != len(samples) or set(observed) != set(groups):
        raise ValueError('高级分析丰度表与当前批次的有效样本不一致')
    header, extra, ignored = ['SampleID', 'Group'], {}, []
    if filename:
        with Path(filename).open(encoding='utf-8-sig', newline='') as stream:
            reader = csv.DictReader(stream)
            names = reader.fieldnames or []
            if ('SampleID' not in names or len(set(names)) != len(names)
                    or any(not n or n != n.strip() or any(c in n for c in '\r\n\t\x00') for n in names)):
                raise ValueError('附加元数据须有唯一 SampleID 表头；列名不能为空、重复或含控制字符')
            for row in reader:
                sid = row.get('SampleID', '')
                if (None in row or not sid or sid != sid.strip() or sid in extra
                        or any(value is None for value in row.values())
                        or any(any(c in str(value) for c in '\r\n\t\x00') for value in row.values())):
                    raise ValueError('附加元数据 SampleID 为空、重复、含空格或列数异常')
                extra[sid] = row
            missing = set(observed) - extra.keys()
            if missing:
                raise ValueError('附加元数据缺少有效样本：' + ', '.join(sorted(missing)))
            ignored = sorted(set(extra) - set(observed))
            header += [name for name in names if name not in header]
    rows = []
    for sid in observed:
        row = dict(extra.get(sid, {}))
        if row.get('Group') and row['Group'] != groups[sid]:
            raise ValueError('附加元数据的 Group 与本批分组不一致：' + sid)
        row.update(SampleID=sid, Group=groups[sid])
        rows.append(row)
    destination = Path(destination)
    with destination.open('w', encoding='utf8', newline='') as stream:
        writer = csv.DictWriter(stream, header, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    repeated = []
    for name in header:
        if name.lower().replace('_', '') in {'subjectid', 'animalid', 'subject', 'animal'}:
            values = [row.get(name, '').strip() for row in rows]
            present = [value for value in values if value and value.upper() != 'NA']
            if len(set(present)) < len(present):
                repeated.append(name)
    return dict(columns=header, samples=len(rows), ignored_samples=ignored,
                repeated_subject_columns=repeated)


def abundance_source(parent, level, api):
    candidates = [parent]
    extensions = parent / '11_extensions'
    if extensions.is_dir():
        for child in sorted(extensions.iterdir(), reverse=True):
            if child.is_dir() and (child / 'abundance_complete.json').is_file():
                saved = json.loads((child / 'parameters.json').read_text(encoding='utf8'))
                if Path(saved.get('preprocess', '/')).resolve() == parent:
                    candidates.append(child)
    for owner in candidates:
        manifest = owner / '10_abundance/manifest.json'
        if not manifest.is_file():
            continue
        data = json.loads(manifest.read_text(encoding='utf8'))
        selected = [item for item in data.get('levels', []) if item.get('level') == level]
        if len(selected) != 1:
            continue
        source = (owner / selected[0]['counts']['tsv']).resolve()
        if not source.is_relative_to(owner.resolve()):
            raise ValueError('丰度清单路径越出任务目录')
        source = api.approved_file(source)
        return source
    raise ValueError('缺少该层级未抽平计数表，请先运行“未抽平 ASV 与各层级丰度导出”')


def options_for(config, saved, parent, run, api):
    level = _number(config, 'advanced_level', 6, 2, 6, True)
    modules = config.get('advanced_modules', ['core', 'ordination'])
    if not isinstance(modules, list) or not modules or set(modules) - MODULES:
        raise ValueError('请选择核心菌群、排序或表型关联中的至少一项')
    inputs = run / '00_inputs'
    inputs.mkdir(exist_ok=True)
    source = abundance_source(parent, level, api)
    total = read_table(source)
    if any(sum(row[i] for row in total[2]) <= 0 for i in range(len(total[0]))):
        raise ValueError('高级分析需要每个有效样本有正的总计数')
    full_paths = write_table_pair(inputs, 'full_counts', total)
    remove_unknown = config.get('advanced_exclude_uncertain', True)
    if not isinstance(remove_unknown, bool):
        raise ValueError('未分类排除选项须为开关')
    keep, omitted = [], []
    for index, name in enumerate(total[1]):
        parts = [part.strip() for part in name.split(';') if part.strip()]
        target = parts[-1] if parts else ''
        label = re.sub(r'^(?:[dkpcofgs]__|D_[0-6]__)', '', target, flags=re.I).strip()
        uncertain = not label or label.lower() in {'others', 'other'} or UNCERTAIN_TAXON.search(target)
        (omitted if remove_unknown and uncertain else keep).append(index)
    if not keep:
        raise ValueError('按当前未分类排除规则没有可分析的分类单元')
    selected = (total[0], [total[1][i] for i in keep], [total[2][i] for i in keep])
    paths = write_table_pair(inputs, 'analysis_counts', selected)
    metadata = config.get('advanced_metadata', '').strip()
    original = api.path(metadata) if metadata else None
    if original:
        original = api.approved_file(original) if (original.parent / 'workstation_run.json').exists() else original
        if not original.is_file():
            raise ValueError('附加元数据 CSV 不存在')
        (inputs / 'metadata_original.csv').write_bytes(original.read_bytes())
    metadata_path = inputs / 'analysis_metadata.csv'
    report = prepare_metadata(original, saved['samples'], total[0], metadata_path)
    aliases = {'predictors': ('advanced_predictors', 'Group'),
               'covariates': ('advanced_covariates', ''),
               'phenotypes': ('advanced_phenotypes', ''),
               'association_covariates': ('advanced_association_covariates', 'Group')}
    options = {key: _columns(config.get(field, default)) for key, (field, default) in aliases.items()}
    if set(options['predictors']) & set(options['covariates']):
        raise ValueError('约束解释变量与控制变量不能重复')
    requested = set().union(*options.values())
    if requested - set(report['columns']):
        raise ValueError('附加元数据缺少所选变量：' + ', '.join(sorted(requested - set(report['columns']))))
    if 'association' in modules and not options['phenotypes']:
        raise ValueError('表型关联须填写附加元数据中的数值型表型列名')
    if set(options['phenotypes']) & set(options['association_covariates']):
        raise ValueError('表型不能同时作为其关联模型的控制变量')
    constrained = config.get('advanced_constrained_method', 'hellinger_rda')
    if constrained not in {'hellinger_rda', 'bray_dbrda', 'none'}:
        raise ValueError('未知约束排序方法')
    association = config.get('advanced_association_method', 'partial_spearman')
    if association not in {'spearman', 'partial_spearman'}:
        raise ValueError('未知关联方法')
    transform = config.get('advanced_association_transform', 'clr')
    if transform not in {'clr', 'relative'}:
        raise ValueError('未知关联数据变换')
    options.update(format=1, counts_path=str(paths['tsv']), total_counts_path=str(full_paths['tsv']),
                   metadata_path=str(metadata_path), level=level, modules=modules,
                   core_min_reads=_number(config, 'advanced_min_reads', 1, 1, 100000000, True),
                   core_min_relative=_number(config, 'advanced_min_relative', .01, 0, 100) / 100,
                   core_prevalence=_number(config, 'advanced_core_prevalence', 80, 1, 100) / 100,
                   feature_prevalence=_number(config, 'advanced_feature_prevalence', 10, 0, 100) / 100,
                   clr_pseudocount=_number(config, 'advanced_pseudocount', .5, .000001, 100),
                   permutations=_number(config, 'advanced_permutations', 999, 99, 99999, True),
                   seed=_number(config, 'advanced_seed', 42, 0, 2147483646, True),
                   top_n=_number(config, 'advanced_top_n', 20, 3, 100, True),
                   constrained_method=constrained, association_method=association, association_transform=transform,
                   repeated_subjects=bool(report['repeated_subject_columns']),
                   group_order=saved['group_order'], group_colors=saved['group_colors'],
                   font=saved.get('font', 'DejaVu Serif'), plot_settings=saved['plot_settings'])
    env = str(config.get('advanced_env', DEFAULT_ENV)).strip()
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', env):
        raise ValueError('高级分析环境名只能包含字母、数字、点、下划线和短横线')
    options['r_env'] = env
    report.update(source_counts=str(source), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  source_metadata_sha256=hashlib.sha256(original.read_bytes()).hexdigest() if original else None,
                  excluded_uncertain=[total[1][i] for i in omitted],
                  excluded_read_fraction={sid: sum(total[2][i][j] for i in omitted) /
                                          sum(row[j] for row in total[2]) for j, sid in enumerate(total[0])},
                  unrarefied=True, original_counts_preserved=True,
                  thresholds_are_heuristics=True)
    (inputs / 'input_report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf8')
    return options


def run_advanced(job, config, saved, parent, api):
    run = Path(job['run'])
    api.begin_stage(job, '00_inputs', '高级分析样本与未抽平输入检查')
    options = options_for(config, saved, parent, run, api)
    filename = run / 'advanced_options.json'
    filename.write_text(json.dumps(options, ensure_ascii=False, indent=2), encoding='utf8')
    api.finish_stage(job, {'samples': len(saved['samples']), 'level': options['level'],
                           'repeated_subjects': options['repeated_subjects']},
                     ['附加元数据按 SampleID 字面匹配；Group 与原批次必须一致。',
                      '未分类单元排除记录于 input_report.json，完整原表及相对丰度分母保留。'])
    runtime = dict(saved, r_env=options['r_env'])
    api.begin_stage(job, '01_analysis', '高级生态、排序与表型关联',
                    {'modules': options['modules'], 'level': options['level'], 'r_env': options['r_env'],
                     'source': 'filtered unrarefied counts', 'seed': options['seed']})
    script = str(api.ROOT / 'advanced_analysis.R')
    api.command(job, api.r_command(runtime, script, '--check', str(filename), str(run)),
                run, process_env=api.r_execution_env(runtime))
    api.command(job, api.r_command(runtime, script, '--run', str(filename), str(run)),
                run, process_env=api.r_execution_env(runtime))
    summary = run / '02_statistics/analysis_summary.json'
    cache = run / '02_statistics/advanced_cache.rds'
    if not summary.is_file() or not cache.is_file():
        raise ValueError('高级分析未生成统计摘要或可复现缓存，请查看日志')
    result = json.loads(summary.read_text(encoding='utf8'))
    api.finish_stage(job, result, [
        '所有分析来自保存的未抽平计数；按模型分别变换，不覆盖原表。',
        '核心菌阈值、伪计数和过滤阈值是可调整规则，并非官方统一最优值。',
        '相关与约束排序用于探索关联，不能证明因果；缺失或不可估计项保存原因。',
        '重复主体被识别时仅保留描述分析，不执行假定独立样本的推断。'])
    (run / 'advanced_complete.json').write_text(json.dumps({'completed_at': time.time(),
                                                            'cache_sha256': hashlib.sha256(cache.read_bytes()).hexdigest()}), encoding='utf8')
    return dict(advanced_run=str(run), advanced_summary=result)


def render_saved_advanced(job, config, api):
    analysis = api.approved_run(config.get('advanced_run', ''))
    marker = analysis / 'advanced_complete.json'
    if not marker.is_file():
        raise ValueError('请选择已完成的高级分析结果')
    original = json.loads((analysis / 'advanced_options.json').read_text(encoding='utf8'))
    cache = analysis / '02_statistics/advanced_cache.rds'
    expected = json.loads(marker.read_text(encoding='utf8')).get('cache_sha256')
    if not expected or hashlib.sha256(cache.read_bytes()).hexdigest() != expected:
        raise ValueError('高级分析缓存已改变，不能与原检验结果混用')
    previous = json.loads((analysis / 'parameters.json').read_text(encoding='utf8'))
    saved = previous['config']
    style = resolve_group_style(config, saved['samples'], fallback=original)
    plotting = resolve_plot_settings(config, fallback=original)
    output = analysis / '03_replots'
    output.mkdir(exist_ok=True)
    run = output / ('figures_' + time.strftime('%Y%m%d_%H%M%S') + '_' + uuid.uuid4().hex[:8])
    run.mkdir(); api.marker(run)
    job.update(run=str(run), log=str(run / 'workflow.log'))
    options = dict(original, cache_run=str(analysis), **style, plot_settings=plotting)
    visual = plot_options(config, fallback=dict(original, diff_plot_top_n=original.get('top_n', 20)))
    options.update({key: visual[key] for key in ('heatmap_low', 'heatmap_mid', 'heatmap_high', 'bubble_low', 'bubble_high')})
    options['top_n'] = visual['diff_plot_top_n']
    options['font'] = config.get('font') or original.get('font', 'DejaVu Serif')
    if config.get('plot_target'):
        target = config['plot_target']
        if not isinstance(target, dict) or not re.fullmatch(r'03_figures/[A-Za-z0-9_-]+', target.get('stem', '')):
            raise ValueError('高级分析单图目标无效')
        options['plot_target'] = target
    filename = run / 'advanced_options.json'
    filename.write_text(json.dumps(options, ensure_ascii=False, indent=2), encoding='utf8')
    save_group_style(run, style); save_plot_settings(run, plotting)
    (run / 'parameters.json').write_text(json.dumps({
        'extension': 'advanced-render', 'preprocess': previous['preprocess'], 'advanced_run': str(analysis),
        'config': dict(saved, **style, plot_settings=plotting), 'statistics_refit': False}, ensure_ascii=False, indent=2), encoding='utf8')
    runtime = dict(saved, r_env=config.get('advanced_env') or original['r_env'])
    if config.get('conda'): runtime['conda'] = config['conda']
    options['r_env'] = runtime['r_env']
    filename.write_text(json.dumps(options, ensure_ascii=False, indent=2), encoding='utf8')
    api.begin_stage(job, '03_figures', '从高级分析缓存重绘', {'statistics_refit': False})
    api.command(job, api.r_command(runtime, str(api.ROOT / 'advanced_analysis.R'), '--render', str(filename), str(run)),
                run, process_env=api.r_execution_env(runtime))
    api.finish_stage(job, {'statistics_refit': False, 'cache_run': str(analysis)})
    (run / 'complete.json').write_text('{"statistics_refit": false}', encoding='utf8')
    job.update(step='高级分析图片重绘完成；未重新拟合统计',
               result={'advanced_run': str(analysis), 'statistics_refit': False})
