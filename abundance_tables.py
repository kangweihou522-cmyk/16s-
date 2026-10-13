"""Export filtered, unrarefied ASV and L2–L7 abundance tables.

``run_abundance(job, saved, api)`` receives the caller's app module so this
module can also be used on completed historical runs without importing app.
QIIME performs raw collapse and normalization. Python also exports cleaned
target-rank summaries, retaining the original tables for downstream analysis.
"""
import csv
import hashlib
import json
import re
from decimal import Decimal, InvalidOperation
from pathlib import Path


LEVELS = ((2, 'phylum'), (3, 'class'), (4, 'order'), (5, 'family'),
          (6, 'genus'), (7, 'species'))
TABLE_HEADERS = {'#OTU ID', '#Feature ID', 'Feature ID', 'feature-id', 'FeatureID'}
RANK_PREFIXES = {2: 'p', 3: 'c', 4: 'o', 5: 'f', 6: 'g', 7: 's'}
RANK_PATTERN = re.compile(r'^(?:([dkpcofgs])__|D_([0-6])__)', re.IGNORECASE)
UNCERTAIN_TAXON = re.compile(
    r'incertae[\s_-]*sedis|no[\s_-]*rank|uncultured|unclassified|unassigned|'
    r'unknown|unidentified|unresolved|ambiguous[\s_-]*taxa', re.IGNORECASE)


def read_table(filename, *, relative=False):
    """Read BIOM TSV without treating its provenance comment as the header."""
    with Path(filename).open(encoding='utf-8-sig', newline='') as stream:
        rows = [row for row in csv.reader(stream, delimiter='\t') if row]
    headers = [index for index, row in enumerate(rows) if row[0] in TABLE_HEADERS]
    if len(headers) != 1:
        raise ValueError('丰度表缺少或重复 BIOM TSV 表头：' + str(filename))
    start = headers[0]
    samples = rows[start][1:]
    if not samples or any(not sid for sid in samples) or len(set(samples)) != len(samples):
        raise ValueError('丰度表样本标识为空或重复：' + str(filename))
    features, values, seen = [], [], set()
    for row in rows[start + 1:]:
        if len(row) != len(samples) + 1 or not row[0] or row[0] in seen:
            raise ValueError('丰度表 FeatureID 为空、重复或列数异常：' + str(filename))
        seen.add(row[0]); features.append(row[0])
        numbers = []
        for value in row[1:]:
            try:
                number = Decimal(value)
            except InvalidOperation:
                raise ValueError('丰度表包含非数值：' + str(filename)) from None
            if not number.is_finite() or number < 0 or (relative and number > 1):
                raise ValueError('丰度表包含无效数值：' + str(filename))
            if not relative and number != number.to_integral_value():
                raise ValueError('丰度表序列计数必须为整数：' + str(filename))
            numbers.append(number if relative else int(number))
        values.append(numbers)
    if not features:
        raise ValueError('丰度表没有特征：' + str(filename))
    return samples, features, values


def read_taxonomy(filename):
    """Retain taxon text verbatim; duplicate IDs must never be silently joined."""
    with Path(filename).open(encoding='utf-8-sig', newline='') as stream:
        rows = [row for row in csv.reader(stream, delimiter='\t') if row and not row[0].startswith('#q2:')]
    if not rows or rows[0][0] not in TABLE_HEADERS or 'Taxon' not in rows[0]:
        raise ValueError('分类表缺少 FeatureID 或 Taxon 表头：' + str(filename))
    header = rows[0]
    if len(set(header)) != len(header):
        raise ValueError('分类表列名重复：' + str(filename))
    taxon_index = header.index('Taxon')
    confidence_index = header.index('Confidence') if 'Confidence' in header else None
    annotations = {}
    for row in rows[1:]:
        if len(row) != len(header) or not row[0] or row[0] in annotations:
            raise ValueError('分类表 FeatureID 为空、重复或列数异常：' + str(filename))
        if not row[taxon_index].strip():
            raise ValueError('分类表 Taxon 为空：' + row[0])
        if confidence_index is not None and row[confidence_index]:
            try:
                confidence = Decimal(row[confidence_index])
            except InvalidOperation:
                raise ValueError('分类表 Confidence 非数值：' + row[0]) from None
            # q2-feature-classifier sums class probabilities (roundoff can
            # slightly exceed 1) and writes -1 when confidence is disabled.
            # Preserve exported text; this annotation never changes counts.
            tolerance = Decimal('1e-12')
            if not confidence.is_finite() or not (
                    confidence == -1 or -tolerance <= confidence <= 1 + tolerance):
                raise ValueError('分类表 Confidence 无效：' + row[0] +
                                 '（值：' + row[confidence_index] + '；文件：' + str(filename) + '）')
        annotations[row[0]] = row[1:]
    return header[1:], annotations


def write_rows(filename, header, rows, *, delimiter='\t'):
    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open('w', encoding='utf8', newline='') as stream:
        writer = csv.writer(stream, delimiter=delimiter, lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)


def write_table_pair(directory, stem, table, *, annotations=None):
    samples, features, values = table
    extra_header, extra_rows = annotations if annotations is not None else ([], {})
    missing = [feature for feature in features if annotations is not None and feature not in extra_rows]
    if missing:
        raise ValueError('过滤后 FeatureID 缺少分类注释：' + ', '.join(missing[:10]))
    header = ['FeatureID', *samples, *extra_header]
    rows = [[feature, *numbers, *extra_rows.get(feature, [])] for feature, numbers in zip(features, values)]
    paths = {}
    for extension, delimiter in (('tsv', '\t'), ('csv', ',')):
        target = Path(directory) / (stem + '.' + extension)
        write_rows(target, header, rows, delimiter=delimiter)
        paths[extension] = target
    return paths


def table_totals(table):
    samples, _, values = table
    return {sid: sum(row[index] for row in values) for index, sid in enumerate(samples)}


def _target_rank(lineage, level):
    """Read only the requested rank, including QIIME's D_0__ nomenclature."""
    parts = [part.strip() for part in lineage.split(';')]
    prefixed = [(part, RANK_PATTERN.match(part)) for part in parts]
    targets = [part[match.end():].strip() for part, match in prefixed
               if match is not None and
               ((match.group(1) or '').lower() == RANK_PREFIXES[level] or
                match.group(2) == str(level - 1))]
    if len(targets) > 1:
        raise ValueError('分类谱系的目标层级重复：' + lineage)
    target = targets[0] if targets else None
    if target is None:
        return '', 'Others', 'missing_rank'
    if not target:
        return target, 'Others', 'empty_rank'
    if (UNCERTAIN_TAXON.search(target) or
            target.lower() in {'na', 'n/a', 'nan', 'none', 'null', 'missing',
                               'other', 'others', 'not assigned', 'not_assigned'}):
        return target, 'Others', 'uncertain_rank'
    return target, target, 'retained'


def clean_taxonomy_table(counts, level):
    """Return cleaned counts, relative values, lineage mapping and a report.

    ``counts`` is the (samples, features, integer rows) tuple from ``read_table``.
    ``level`` is 2–7. Features contain taxonomy lineages. Only the requested
    rank determines uncertainty; identical retained target names are summed,
    even across different parents, with every source lineage in the mapping.
    Input objects are untouched. Zero-total samples have all-zero relative
    values and are explicitly listed in the JSON-safe report.
    """
    if level not in RANK_PREFIXES:
        raise ValueError('分类清理仅支持 L2–L7')
    samples, features, values = counts
    if (not samples or len(set(samples)) != len(samples) or
            len(features) != len(values) or not features):
        raise ValueError('分类计数表的样本或特征维度无效')
    grouped, source_counts, mapping = {}, {}, []
    reason_counts = {}
    for lineage, row in zip(features, values):
        if len(row) != len(samples) or any(value < 0 or int(value) != value for value in row):
            raise ValueError('分类清理需要非负整数计数，且列数必须匹配样本')
        target, cleaned, reason = _target_rank(lineage, level)
        grouped.setdefault(cleaned, [0] * len(samples))
        grouped[cleaned] = [total + int(value) for total, value in zip(grouped[cleaned], row)]
        source_counts[cleaned] = source_counts.get(cleaned, 0) + 1
        reason_counts[reason] = reason_counts.get(reason, 0) + 1
        mapping.append({'original_lineage': lineage, 'target_taxon': target,
                        'cleaned_taxon': cleaned, 'reason': reason})
    # Keep the original encounter order; reserve the final row for Others.
    labels = [label for label in grouped if label != 'Others']
    if 'Others' in grouped:
        labels.append('Others')
    cleaned_counts = (list(samples), labels, [grouped[label] for label in labels])
    totals = table_totals(cleaned_counts)
    if totals != table_totals(counts):
        raise ValueError('分类清理后样本计数总和发生变化')
    relative = (list(samples), list(labels),
                [[Decimal(value) / Decimal(totals[sid]) if totals[sid] else Decimal(0)
                  for sid, value in zip(samples, row)] for row in cleaned_counts[2]])
    report = {'level': level, 'name': dict(LEVELS)[level],
              'input_feature_count': len(features), 'cleaned_feature_count': len(labels),
              'reason_counts': reason_counts, 'per_sample_counts': totals,
              'per_sample_relative_sums': {sid: float(value) for sid, value in table_totals(relative).items()},
              'zero_count_samples': [sid for sid in samples if not totals[sid]],
              'others_per_sample_counts': dict(zip(samples, grouped.get('Others', [0] * len(samples)))),
              'merged_taxa': {label: number for label, number in source_counts.items() if number > 1},
              'aggregation': 'sum counts by target-rank label; preserve original lineages in mapping',
              'normalization': 'cleaned counts / per-sample total; zero-total samples remain zero'}
    return cleaned_counts, relative, mapping, report


def export_clean_taxonomy(directory, counts, level):
    """Write additional clean tables and audit files without touching raw files."""
    cleaned, relative, mapping, report = clean_taxonomy_table(counts, level)
    counts_paths = write_table_pair(directory, 'counts_cleaned', cleaned)
    relative_paths = write_table_pair(directory, 'relative_cleaned', relative)
    mapping_paths = {}
    header = ['original_lineage', 'target_taxon', 'cleaned_taxon', 'reason']
    for extension, delimiter in (('tsv', '\t'), ('csv', ',')):
        target = Path(directory) / ('taxonomy_cleaning_mapping.' + extension)
        write_rows(target, header, [[row[key] for key in header] for row in mapping], delimiter=delimiter)
        mapping_paths[extension] = target
    report_path = Path(directory) / 'taxonomy_cleaning_report.json'
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
    return {'counts': counts_paths, 'relative': relative_paths, 'mapping': mapping_paths,
            'report': report_path, 'feature_count': len(cleaned[1])}


def validate_relative(counts, relative):
    if counts[:2] != relative[:2]:
        raise ValueError('相对丰度表的样本或 FeatureID 顺序与计数表不一致')
    totals = table_totals(counts)
    for sid, value in table_totals(relative).items():
        if totals[sid] <= 0 or abs(value - 1) > Decimal('0.000001'):
            raise ValueError('相对丰度样本总和异常：' + sid)
    for count_row, relative_row in zip(counts[2], relative[2]):
        for sid, count, value in zip(counts[0], count_row, relative_row):
            expected = Decimal(count) / Decimal(totals[sid])
            if abs(value - expected) > Decimal('0.000001'):
                raise ValueError('相对丰度与未抽平计数不一致：' + sid)


def _export(job, saved, api, source, destination):
    """A QIIME export directory is created by QIIME, never by this module."""
    destination = Path(destination)
    if not destination.exists():
        destination.parent.mkdir(parents=True, exist_ok=True)
        api.q(job, saved, 'tools', 'export', '--input-path', str(source), '--output-path', str(destination))
    return destination


def _export_table(job, saved, api, source, destination, *, relative=False):
    destination = _export(job, saved, api, source, destination)
    biom = destination / 'feature-table.biom'
    if not biom.is_file():
        raise ValueError('QIIME 特征表导出缺少 feature-table.biom：' + str(destination))
    tsv = destination / 'feature-table.tsv'
    if not tsv.is_file():
        api.command(job, api.conda_cmd(saved) + ['biom', 'convert', '-i', str(biom), '-o', str(tsv), '--to-tsv'], Path(job['run']))
    return biom, read_table(tsv, relative=relative)


def _relative(job, saved, api, count_qza, directory, counts):
    directory = Path(directory)
    output = directory / 'relative-frequency.qza'
    if not output.is_file():
        api.q(job, saved, 'feature-table', 'relative-frequency', '--i-table', str(count_qza), '--o-relative-frequency-table', str(output))
    biom, table = _export_table(job, saved, api, output, directory / 'relative_export', relative=True)
    validate_relative(counts, table)
    paths = write_table_pair(directory, 'relative', table)
    return output, biom, table, paths


def run_abundance(job, saved, api, *, source_run=None):
    """Return a JSON-safe manifest and also save it to 10_abundance/manifest.json.

    `source_run` defaults to `job['run']`; an extension can read an existing
    parent while writing all new files under its own run directory. The source
    may use either the versioned layout or legacy flat paths.
    `saved` contains the existing QIIME conda settings. Existing QIIME outputs
    are reused, permitting completion of an interrupted abundance-only run.
    """
    run = Path(job['run']).resolve()
    source_run = Path(source_run).resolve() if source_run is not None else run
    directory = run / '10_abundance'
    asv_directory = directory / 'asv'
    count_qza = api.artifact_path(source_run, 'filtered-table.qza')
    taxonomy_qza = api.artifact_path(source_run, 'taxonomy.qza')
    barplot_qzv = api.artifact_path(source_run, 'taxa-barplot.qzv')
    for filename in (count_qza, taxonomy_qza, barplot_qzv):
        if not Path(filename).is_file():
            raise ValueError('丰度导出需要已完成预处理的文件：' + str(filename))
    api.begin_stage(job, '10_abundance', '未抽平 ASV 与各分类层级丰度导出',
                    {'source': 'filtered unrarefied table', 'levels': [level for level, _ in LEVELS],
                     'normalization': 'per-sample relative frequency', 'annotation_join': 'FeatureID'})
    asv_directory.mkdir(parents=True, exist_ok=True)
    source_tsv = api.artifact_path(source_run, 'exported-table/feature-table.tsv')
    source_biom = api.artifact_path(source_run, 'exported-table/feature-table.biom')
    if Path(source_tsv).is_file():
        counts = read_table(source_tsv)
    else:
        source_biom, counts = _export_table(job, saved, api, count_qza, asv_directory / 'counts_export')
        source_tsv = source_biom.with_suffix('.tsv')
    declared_samples = saved.get('samples', [])
    declared_ids = [sample['SampleID'] for sample in declared_samples]
    # QIIME taxa filter-table can remove samples left with no reads.
    if declared_ids and (len(set(declared_ids)) != len(declared_ids) or not set(counts[0]).issubset(declared_ids)):
        raise ValueError('过滤后丰度表 SampleID 与保存的样本映射不一致')
    source_taxonomy = api.artifact_path(source_run, 'exported-taxonomy/taxonomy.tsv')
    if not Path(source_taxonomy).is_file():
        source_taxonomy = _export(job, saved, api, taxonomy_qza, directory / 'taxonomy_source_export') / 'taxonomy.tsv'
    annotations = read_taxonomy(source_taxonomy)
    count_paths = write_table_pair(asv_directory, 'counts', counts)
    annotated_count_paths = write_table_pair(asv_directory, 'counts_with_taxonomy', counts, annotations=annotations)
    taxonomy_paths = {}
    taxonomy_rows = [[feature, *annotations[1][feature]] for feature in counts[1]]
    for extension, delimiter in (('tsv', '\t'), ('csv', ',')):
        output = asv_directory / ('filtered_taxonomy.' + extension)
        write_rows(output, ['FeatureID', *annotations[0]], taxonomy_rows, delimiter=delimiter)
        taxonomy_paths[extension] = output
    relative_qza, relative_biom, relative_table, relative_paths = _relative(job, saved, api, count_qza, asv_directory, counts)
    annotated_relative_paths = write_table_pair(asv_directory, 'relative_with_taxonomy', relative_table, annotations=annotations)

    def relative_path(filename):
        filename = Path(filename)
        return str(filename.relative_to(run)) if filename.is_relative_to(run) else str(filename)

    def pair_paths(paths):
        return {key: relative_path(value) for key, value in paths.items()}

    asv = {'counts_qza': relative_path(count_qza), 'counts': pair_paths(count_paths),
           'counts_with_taxonomy': pair_paths(annotated_count_paths),
           'filtered_taxonomy': pair_paths(taxonomy_paths),
           'relative_qza': relative_path(relative_qza), 'relative_biom': relative_path(relative_biom),
           'relative': pair_paths(relative_paths), 'relative_with_taxonomy': pair_paths(annotated_relative_paths)}
    if Path(source_biom).is_file():
        asv['counts_biom'] = relative_path(source_biom)
    levels = []
    for level, name in LEVELS:
        target = directory / ('L' + str(level) + '_' + name)
        target.mkdir(parents=True, exist_ok=True)
        collapsed = target / 'counts.qza'
        visualization = target / 'counts.qzv'
        if not collapsed.is_file():
            api.q(job, saved, 'taxa', 'collapse', '--i-table', str(count_qza), '--i-taxonomy', str(taxonomy_qza), '--p-level', str(level), '--o-collapsed-table', str(collapsed))
        if not visualization.is_file():
            api.q(job, saved, 'feature-table', 'summarize', '--i-table', str(collapsed), '--o-visualization', str(visualization))
        biom, collapsed_counts = _export_table(job, saved, api, collapsed, target / 'counts_export')
        if collapsed_counts[0] != counts[0] or table_totals(collapsed_counts) != table_totals(counts):
            raise ValueError('分类折叠后的样本或序列总数发生变化：L' + str(level))
        paths = write_table_pair(target, 'counts', collapsed_counts)
        relative_qza, relative_biom, _, relative_paths = _relative(job, saved, api, collapsed, target, collapsed_counts)
        cleaned = export_clean_taxonomy(target, collapsed_counts, level)
        levels.append({'level': level, 'name': name, 'feature_count': len(collapsed_counts[1]),
                       'counts_qza': relative_path(collapsed), 'counts_qzv': relative_path(visualization),
                       'counts_biom': relative_path(biom), 'counts': pair_paths(paths),
                       'relative_qza': relative_path(relative_qza), 'relative_biom': relative_path(relative_biom),
                       'relative': pair_paths(relative_paths),
                       'cleaned': {'counts': pair_paths(cleaned['counts']),
                                   'relative': pair_paths(cleaned['relative']),
                                   'mapping': pair_paths(cleaned['mapping']),
                                   'report': relative_path(cleaned['report']),
                                   'feature_count': cleaned['feature_count']}})
    barplot_export = _export(job, saved, api, barplot_qzv, directory / 'taxa_barplot_data')
    barplot_files = [relative_path(filename) for filename in sorted(barplot_export.rglob('*')) if filename.is_file()]
    if not barplot_files:
        raise ValueError('taxa barplot 数据导出为空')
    summary = {'sample_count': len(counts[0]), 'asv_count': len(counts[1]),
               'total_sequences': sum(table_totals(counts).values()), 'per_sample_counts': table_totals(counts),
               'level_feature_counts': {str(item['level']): item['feature_count'] for item in levels}}
    manifest = {'format': 1, 'source': 'filtered_unrarefied',
                'source_run': str(source_run),
                'source_sha256': {'filtered-table.tsv': hashlib.sha256(Path(source_tsv).read_bytes()).hexdigest(),
                                  'taxonomy.tsv': hashlib.sha256(Path(source_taxonomy).read_bytes()).hexdigest()},
                'source_table': relative_path(count_qza), 'source_taxonomy': relative_path(taxonomy_qza),
                'asv': asv, 'levels': levels, 'taxa_barplot': {'qzv': relative_path(barplot_qzv),
                'export_dir': relative_path(barplot_export), 'files': barplot_files}, 'summary': summary}
    (directory / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
    api.finish_stage(job, summary, ['所有计数与相对丰度来自过滤后未抽平表，不改动原始数据或过滤规则。',
                                   '分类注释按 FeatureID 严格匹配；CSV 使用标准引号保留分类名称中的逗号。',
                                   'Confidence 原值保留；容忍 1e-12 的浮点边界误差，-1 表示 QIIME 未计算置信度，不作为低置信度过滤依据。',
                                   '各层级另存清理表：仅按目标层级判断未知分类并合并为 Others，同名分类求和后重新计算相对丰度；原始表保留。',
                                   'QIIME taxa collapse 使用原分类 QZA；官方动作忽略不在过滤后表中的额外分类条目。'])
    return manifest


def resume_abundance(job, config, api):
    """Finish only abundance export after completed preprocessing stages."""
    import time
    run = api.approved_run(config.get('result_run', ''))
    saved = json.loads((run / 'parameters.json').read_text(encoding='utf8'))
    if saved.get('extension') or not saved.get('samples') or (run / 'preprocess_complete.json').exists():
        raise ValueError('请选择丰度导出失败、尚未完成的预处理 run；完整历史任务请使用常规丰度导出')
    if config.get('batch_id') and saved.get('batch_id') and config['batch_id'] != saved['batch_id']:
        raise ValueError('结果目录属于其他批次，请先打开对应历史批次')
    for folder in ('04_taxonomy_filter', '05_phylogeny', '06_rarefaction'):
        report = run / folder / 'report.json'
        if not report.is_file() or json.loads(report.read_text(encoding='utf8')).get('status') != 'done':
            raise ValueError('此前阶段未完成，不能仅补做丰度导出：' + folder)
    for logical in ('filtered-table.qza', 'taxonomy.qza', 'taxa-barplot.qzv',
                    'rooted-tree.qza', 'alpha-rarefaction.qzv', 'exported-table/feature-table.tsv'):
        filename = api.artifact_path(run, logical)
        if not filename.is_file() or not filename.stat().st_size:
            raise ValueError('此前结果缺失，不能补做丰度导出：' + str(filename))
    job.update(run=str(run), log=str(run / 'workflow.log'), step='仅补做丰度导出；复用已完成的预处理')
    abundance = run_abundance(job, saved, api)
    result = {'abundance': abundance, 'resumed_abundance': True}
    recommendation = run / 'depth_recommendation.json'
    if recommendation.is_file():
        result['depth_advice'] = json.loads(recommendation.read_text(encoding='utf8'))
    (run / 'preprocess_complete.json').write_text(json.dumps({
        'completed_at': time.time(), 'completion_mode': 'resumed_abundance'}), encoding='utf8')
    job.update(result=result, step='丰度导出补做完成；可继续选择抽平深度及后续分析')
