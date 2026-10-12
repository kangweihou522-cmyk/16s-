"""Export filtered, unrarefied ASV and L2–L7 abundance tables.

``run_abundance(job, saved, api)`` receives the caller's app module so this
module can also be used on completed historical runs without importing app.
QIIME performs collapse and normalization; Python only validates and exports.
"""
import csv
import hashlib
import json
from decimal import Decimal, InvalidOperation
from pathlib import Path


LEVELS = ((2, 'phylum'), (3, 'class'), (4, 'order'), (5, 'family'),
          (6, 'genus'), (7, 'species'))
TABLE_HEADERS = {'#OTU ID', '#Feature ID', 'Feature ID', 'feature-id', 'FeatureID'}


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
            if not confidence.is_finite() or not 0 <= confidence <= 1:
                raise ValueError('分类表 Confidence 无效：' + row[0])
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
        levels.append({'level': level, 'name': name, 'feature_count': len(collapsed_counts[1]),
                       'counts_qza': relative_path(collapsed), 'counts_qzv': relative_path(visualization),
                       'counts_biom': relative_path(biom), 'counts': pair_paths(paths),
                       'relative_qza': relative_path(relative_qza), 'relative_biom': relative_path(relative_biom),
                       'relative': pair_paths(relative_paths)})
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
                                   'QIIME taxa collapse 使用原分类 QZA；官方动作忽略不在过滤后表中的额外分类条目。'])
    return manifest
