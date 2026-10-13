"""Differential taxonomy on saved, unrarefied counts; method names are explicit."""
import csv
import hashlib
import json
import math
import shutil
from collections import Counter
from decimal import Decimal, InvalidOperation
from pathlib import Path
from group_style import resolve_group_style
from plot_settings import resolve_plot_settings

LEVEL_NAMES = {2: 'phylum', 3: 'class', 4: 'order', 5: 'family', 6: 'genus'}
METHODS = {'ancombc2': 'ANCOM-BC2', 'qiime_ancom': 'QIIME 2 ANCOM'}
REFERENCES = {
    'ancombc2': 'https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html',
    'qiime_ancom': 'https://github.com/qiime2/q2-composition/blob/2025.7.0/q2_composition/_ancom.py',
    'original_ancom_ii': 'https://github.com/FrederickHuangLin/ANCOM',
}


def plot_options(config, fallback=None):
    fallback = fallback or {}
    defaults = {'diff_comparisons': 'all_pairs', 'diff_heatmap_scale': 'relative_abundance',
                'heatmap_low': '#2166AC', 'heatmap_mid': '#F7F7F7', 'heatmap_high': '#B2182B',
                'bubble_color_mode': 'abundance', 'bubble_low': '#E8F1FA', 'bubble_high': '#2166AC',
                'diff_plot_top_n': 30}
    values = {key: config.get(key, fallback.get(key, default)) for key, default in defaults.items()}
    for key, allowed in [('diff_comparisons', ('all_pairs', 'reference')),
                         ('diff_heatmap_scale', ('row_zscore', 'relative_abundance')),
                         ('bubble_color_mode', ('abundance', 'group'))]:
        if values[key] not in allowed:
            raise ValueError('差异菌绘图选项无效：' + key)
    import re
    for key in ('heatmap_low', 'heatmap_mid', 'heatmap_high', 'bubble_low', 'bubble_high'):
        if not isinstance(values[key], str) or not re.fullmatch(r'#[0-9a-fA-F]{6}', values[key]):
            raise ValueError('热图和气泡图颜色必须使用 #RRGGBB：' + key)
        values[key] = values[key].upper()
    raw = values['diff_plot_top_n']
    try:
        if isinstance(raw, bool) or str(int(raw)) != str(raw): raise ValueError()
        values['diff_plot_top_n'] = int(raw)
    except (TypeError, ValueError):
        raise ValueError('差异菌展示数量必须是整数') from None
    if not 1 <= values['diff_plot_top_n'] <= 200:
        raise ValueError('差异菌展示数量应为 1 至 200')
    return values


def options(config, samples):
    method = config.get('diff_method', '').strip()
    if method not in METHODS:
        raise ValueError('请选择明确方法：ANCOM-BC2 或 QIIME 2 ANCOM；它们不是原始 ANCOM-II')
    groups = list(dict.fromkeys(x['Group'] for x in samples))
    if len(groups) < 2:
        raise ValueError('差异菌群分析需要至少两个实际分组')
    reference = config.get('reference_group') or ('CON' if 'CON' in groups else groups[0])
    if reference not in groups:
        raise ValueError('参照组不存在于保存的样本分组中')
    levels = config.get('diff_levels', list(LEVEL_NAMES))
    if isinstance(levels, str):
        levels = levels.split(',')
    try:
        levels = sorted({int(level) for level in levels})
    except (TypeError, ValueError):
        raise ValueError('差异菌群分类层级必须为 2 至 6（门至属）') from None
    if not levels or any(level not in LEVEL_NAMES for level in levels):
        raise ValueError('差异菌群分类层级必须为 2 至 6（门至属）')
    values = {'method': method, 'method_name': METHODS[method], 'levels': levels,
              'group_column': 'Group', 'reference_group': reference, 'groups': groups,
              'p_adj_method': 'BH', 'struc_zero': True, 'neg_lb': False,
              'pseudo_sens': True, 'global': len(groups) > 2, 'pairwise': len(groups) > 2,
              'seed': 202507, 'iter_control': {'tol': .01, 'max_iter': 20, 'verbose': True},
              'em_control': {'tol': .00001, 'max_iter': 100},
              'mdfdr_control': {'fwer_ctrl_method': 'holm', 'B': 100},
              'source': 'filtered unrarefied integer counts', 'reference': REFERENCES[method]}
    for key, default, lower, upper in [('alpha', .05, 0, 1), ('prv_cut', .1, 0, 1)]:
        try:
            ui_key = 'diff_alpha' if key == 'alpha' else 'diff_prevalence'
            value = float(config.get(ui_key, config.get(key, default)))
        except (TypeError, ValueError):
            raise ValueError(key + ' 必须为数值') from None
        if not math.isfinite(value) or not lower <= value <= upper or key == 'alpha' and value == 0:
            raise ValueError(key + ' 超出有效范围')
        values[key] = value
    for key, default, maximum in [('lib_cut', 0, 1000000000), ('threads', 1, 256)]:
        try:
            raw = str(config.get(key, default))
            value = int(raw)
        except (TypeError, ValueError):
            raise ValueError(key + ' 必须为整数') from None
        if value < (1 if key == 'threads' else 0) or value > maximum:
            raise ValueError(key + ' 超出有效范围')
        values[key] = value
    values['ancom_env'] = config.get('ancom_env') or '16s-workstation-ancom'
    if method == 'qiime_ancom':
        values['ancombc2_only_options_not_applied'] = {
            key: values[key] for key in ('alpha', 'prv_cut', 'lib_cut', 'reference_group')}
        values.update({'alpha': 'official QIIME/scikit-bio ANCOM default', 'p_adj_method': 'official ANCOM default',
                       'prv_cut': None, 'lib_cut': None, 'struc_zero': False, 'pseudo_sens': False,
                       'global': True, 'pairwise': False,
                       'reference_group_note': 'ANCOM is an omnibus Group test; reference choice is not applied',
                       'pseudocount': 1, 'seed': None, 'iter_control': None,
                       'em_control': None, 'mdfdr_control': None})
    values.update(plot_options(config))
    values['plot_settings'] = resolve_plot_settings(config)
    return values


def read_counts(filename):
    """Read canonical or BIOM TSV without changing literal identifiers."""
    with Path(filename).open(encoding='utf-8-sig', newline='') as stream:
        rows = list(csv.reader(stream, delimiter='\t'))
    headers = [i for i, row in enumerate(rows) if row and row[0] in
               ('FeatureID', '#OTU ID', '#Feature ID', 'feature-id', 'Feature ID')]
    if len(headers) != 1:
        raise ValueError('原始计数表缺少或重复表头')
    header = rows[headers[0]]
    ids = header[1:]
    if not ids or len(ids) != len(set(ids)) or any(not sid for sid in ids):
        raise ValueError('原始计数表样本标识无效')
    features, seen = [], set()
    for row in rows[headers[0] + 1:]:
        if not row:
            continue
        if len(row) != len(header) or not row[0] or row[0] in seen:
            raise ValueError('原始计数表分类标识或列数异常')
        seen.add(row[0])
        try:
            counts = [Decimal(x) for x in row[1:]]
        except InvalidOperation:
            raise ValueError('差异菌群分析必须使用非负整数原始计数') from None
        if any(not n.is_finite() or n < 0 or n != n.to_integral_value() or n > 2 ** 53 for n in counts):
            raise ValueError('差异菌群分析必须使用非负整数原始计数')
        features.append([row[0], *map(int, counts)])
    if not features:
        raise ValueError('原始计数表没有分类特征')
    return ids, features


def aligned_metadata(ids, features, samples):
    saved = {x['SampleID']: x['Group'] for x in samples}
    if len(saved) != len(samples) or any(not sid or not group for sid, group in saved.items()):
        raise ValueError('保存的样本分组重复或缺失')
    unknown = set(ids) - set(saved)
    if unknown:
        raise ValueError('计数表包含没有保存分组的样本：' + ', '.join(sorted(unknown)))
    totals = [sum(row[i + 1] for row in features) for i in range(len(ids))]
    if any(total <= 0 for total in totals):
        raise ValueError('计数表存在零文库样本；请检查过滤结果，不能用伪计数代替真实文库')
    rows = [(sid, saved[sid]) for sid in ids]
    sizes = Counter(group for _, group in rows)
    if len(sizes) < 2 or any(n < 2 for n in sizes.values()):
        raise ValueError('差异检验至少需要两组，且每组至少两个实际计数样本')
    missing = [{'SampleID': sid, 'Group': group, 'reason': 'absent from unrarefied count table'}
               for sid, group in saved.items() if sid not in ids]
    return rows, dict(sizes), missing


def _contained(base, relative):
    target = (base / relative).resolve()
    if not target.is_relative_to(base.resolve()) or not target.is_file():
        raise ValueError('丰度清单中的计数文件不存在或超出结果目录')
    return target


def _manifest_tables(config, parent, selected, api=None):
    source = Path(config.get('abundance_run') or parent)
    filename = source / '10_abundance' / 'manifest.json'
    if not filename.is_file():
        if config.get('abundance_run'):
            raise ValueError('指定的丰度结果没有门至属原始计数清单')
        return {}
    manifest = json.loads(filename.read_text(encoding='utf8'))
    if manifest.get('source') != 'filtered_unrarefied' or Path(manifest.get('source_run', '')).resolve() != Path(parent).resolve():
        raise ValueError('丰度缓存不属于当前预处理批次，请重新导出当前批次丰度')
    if api is None:
        from run_layout import artifact_path
    else:
        artifact_path = api.artifact_path
    source_counts = artifact_path(parent, 'exported-table/feature-table.tsv')
    source_taxonomy = artifact_path(parent, 'exported-taxonomy/taxonomy.tsv')
    for key, target in [('filtered-table.tsv', source_counts), ('taxonomy.tsv', source_taxonomy)]:
        if not Path(target).is_file() or manifest.get('source_sha256', {}).get(key) != hashlib.sha256(Path(target).read_bytes()).hexdigest():
            raise ValueError('丰度缓存来源计数或分类已改变，或缺少来源校验；请重新导出丰度')
    ids, raw_features = read_counts(source_counts)
    from abundance_tables import read_taxonomy
    columns, taxonomy_rows = read_taxonomy(source_taxonomy)
    taxon_column = columns.index('Taxon')
    taxonomy = {feature: row[taxon_column] for feature, row in taxonomy_rows.items()}
    maximum = max(len(taxon.split(';')) for taxon in taxonomy.values())
    entries = {int(x['level']): x for x in manifest['levels']}
    if any(level not in entries for level in selected):
        raise ValueError('丰度结果缺少选定分类层级')
    tables = {}
    for level in selected:
        tsv = _contained(source, entries[level]['counts']['tsv'])
        actual_ids, actual_rows = read_counts(tsv)
        if actual_ids != ids:
            raise ValueError('缓存分类表样本或顺序与当前过滤后表不一致')
        # QIIME 2025.7 collapse strips rank whitespace and pads short lineages with __.
        # This arithmetic check verifies cached count content; QIIME still performs the collapse.
        expected = {}
        if level > maximum:
            raise ValueError('缓存分类层级超出当前实际分类深度')
        for row in raw_features:
            if row[0] not in taxonomy:
                raise ValueError('当前分类表缺少过滤后 ASV 的分类标识')
            ranks = [rank.strip() for rank in taxonomy[row[0]].split(';')]
            label = ';'.join((ranks + ['__'] * (maximum - len(ranks)))[:level])
            counts = expected.setdefault(label, [0] * len(ids))
            for index, value in enumerate(row[1:]):
                counts[index] += value
        if {row[0]: row[1:] for row in actual_rows} != expected:
            raise ValueError('缓存分类计数与当前过滤后 ASV 和分类内容不一致；请重新导出丰度')
        tables[level] = {'tsv': tsv, 'qza': _contained(source, entries[level]['counts_qza']),
                         'source_run': str(source.resolve()), 'source_sha256': manifest['source_sha256']}
    return tables


def run_differential(job, config, saved, parent, api):
    """The app owns isolated run creation, marker, parameters and final inventory."""
    run = Path(job['run'])
    opts = options(config, saved['samples'])
    opts['group_style'] = resolve_group_style(config, saved['samples'], fallback=saved)
    opts['plot_settings'] = resolve_plot_settings(config, fallback=saved)
    tables = _manifest_tables(config, parent, opts['levels'], api)
    api.begin_stage(job, '01_counts', '差异检验输入：未抽平门至属计数', opts)
    entries, expected_ids, metadata_rows, groups, excluded = [], None, None, None, None
    for level in opts['levels']:
        directory = run / '01_counts' / ('L' + str(level) + '_' + LEVEL_NAMES[level])
        directory.mkdir(parents=True, exist_ok=True)
        if level in tables:
            source = tables[level]
        else:
            qza = directory / 'counts.qza'
            api.q(job, saved, 'taxa', 'collapse', '--i-table', api.artifact_path(parent, 'filtered-table.qza'),
                  '--i-taxonomy', api.artifact_path(parent, 'taxonomy.qza'), '--p-level', level, '--o-collapsed-table', qza)
            api.export_artifact(job, saved, qza, directory / 'biom_export')
            tsv = directory / 'biom_export' / 'counts.tsv'
            api.command(job, api.conda_cmd(saved) + ['biom', 'convert', '-i', str(directory / 'biom_export' / 'feature-table.biom'),
                        '-o', str(tsv), '--to-tsv'], run)
            source = {'tsv': tsv, 'qza': qza, 'source_run': str(parent)}
        ids, features = read_counts(source['tsv'])
        rows, sizes, missing = aligned_metadata(ids, features, saved['samples'])
        if expected_ids is not None and ids != expected_ids:
            raise ValueError('不同分类层级的原始计数样本或顺序不一致')
        if opts['reference_group'] not in sizes:
            raise ValueError('参照组在过滤后原始计数表中没有样本')
        expected_ids, metadata_rows, groups, excluded = ids, rows, sizes, missing
        canonical = directory / 'counts.tsv'
        with canonical.open('w', encoding='utf8', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t')
            writer.writerow(['FeatureID', *ids]); writer.writerows(features)
        entries.append({'level': level, 'name': LEVEL_NAMES[level], 'counts_tsv': str(canonical),
                        'counts_sha256': hashlib.sha256(canonical.read_bytes()).hexdigest(),
                        'counts_qza': str(source['qza']), 'source_run': source['source_run'],
                        'source_sha256': source.get('source_sha256'),
                        'feature_count': len(features), 'sample_count': len(ids)})
    metadata = run / '01_counts' / 'metadata.tsv'
    with metadata.open('w', encoding='utf8', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t'); writer.writerow(['sample-id', 'Group'])
        writer.writerow(['#q2:types', 'categorical']); writer.writerows(metadata_rows)
    # All factor levels come from actual saved Groups; lexical labels never imply dosage.
    opts.update({'groups': [g for g in opts['groups'] if g in groups],
                 'global': len(groups) > 2 if opts['method'] == 'ancombc2' else True,
                 'pairwise': len(groups) > 2 if opts['method'] == 'ancombc2' else False,
                 'group_sizes': groups, 'excluded_samples': excluded,
                 'metadata_tsv': str(metadata), 'metadata_sha256': hashlib.sha256(metadata.read_bytes()).hexdigest(), 'tables': entries,
                 'font': config.get('font') or saved.get('font', 'sans')})
    options_file = run / 'differential_options.json'
    options_file.write_text(json.dumps(opts, ensure_ascii=False, indent=2), encoding='utf8')
    notes = ['计数来自过滤后未抽平表；按完整分类谱系合并，未知分类仍保留其谱系标识。',
             'Group 来自本次保存的元数据；按分类表真实样本标识对齐，不推断组间剂量或顺序。']
    if any(n < 5 for n in groups.values()):
        notes.append('某组少于 5 个样本，统计功效和模型稳定性有限，请结合实际结果及警告解释。')
    api.finish_stage(job, {'tables': entries, 'group_sizes': groups, 'excluded_samples': excluded}, notes)
    api.begin_stage(job, '02_differential', METHODS[opts['method']] + ' 差异菌群检验', opts)
    if opts['method'] == 'ancombc2':
        rconfig = dict(saved, r_env=opts['ancom_env'])
        script = api.ROOT / 'differential_ancombc.R'
        api.command(job, api.r_command(rconfig, script, '--check', options_file, run), run,
                    process_env=api.r_execution_env(rconfig))
        try:
            api.command(job, api.r_command(rconfig, script, '--run', options_file, run), run,
                        process_env=api.r_execution_env(rconfig))
        except RuntimeError:
            summary = run / '02_differential' / 'summary.json'
            if summary.is_file():
                job['result'] = json.loads(summary.read_text(encoding='utf8'))
            raise
        summary = run / '02_differential' / 'summary.json'
        if not summary.is_file():
            raise ValueError('ANCOM-BC2 未产生实际结果摘要；请查看依赖和模型运行日志')
        result = json.loads(summary.read_text(encoding='utf8'))
    else:
        # Version probes and the exact official commands are retained in workflow.log.
        api.q(job, saved, 'info')
        statuses = []
        for entry in entries:
            directory = run / '02_differential' / ('L' + str(entry['level']) + '_' + entry['name'])
            directory.mkdir(parents=True, exist_ok=True)
            composition = directory / 'composition.qza'
            visualization = directory / 'ancom.qzv'
            api.q(job, saved, 'composition', 'add-pseudocount', '--i-table', entry['counts_qza'],
                  '--p-pseudocount', 1, '--o-composition-table', composition)
            api.q(job, saved, 'composition', 'ancom', '--i-table', composition, '--m-metadata-file', metadata,
                  '--m-metadata-column', 'Group', '--p-transform-function', 'clr', '--o-visualization', visualization)
            api.export_artifact(job, saved, visualization, directory / 'exported')
            result_file = directory / 'exported' / 'ancom.tsv'
            if not result_file.is_file():
                raise ValueError('QIIME 2 ANCOM 没有导出实际 W 检验结果')
            shutil.copy2(result_file, directory / 'ancom_results.tsv')
            statuses.append({'level': entry['level'], 'method': 'QIIME 2 ANCOM', 'status': 'executed',
                             'results_tsv': str(directory / 'ancom_results.tsv'), 'visualization': str(visualization)})
        result = {'method': 'QIIME 2 ANCOM', 'levels': statuses,
                  'notes': ['W 与拒绝判定由 QIIME 2 官方 ANCOM 计算；本方法不输出 ANCOM-BC2 的 LFC/SE/p/q。',
                            'QIIME 2 ANCOM 是总体 Group 检验，不产生组间两两效应；其检验阈值采用官方动作默认值。']}
        (run / '02_differential' / 'summary.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf8')
        api.command(job, api.r_command(saved, api.ROOT / 'differential_ancombc.R', '--plot-qiime', options_file, run),
                    run, process_env=api.r_execution_env(saved))
    api.finish_stage(job, result, [REFERENCES[opts['method']],
                     '原始 ANCOM-II、ANCOM-BC2 与 QIIME 2 ANCOM 是不同方法；本结果以明确选择的方法命名。'])
    return result
