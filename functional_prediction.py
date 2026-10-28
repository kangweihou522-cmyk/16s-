"""Isolated PICRUSt2 predictions and transparent KO annotation summaries.

Only the standard library is used here. Predictions themselves are delegated
to the user's installed PICRUSt2 environment through the app command recorder.
"""
import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
from pathlib import Path


MAPPING_COLUMNS = ['KO', 'L1', 'L2', 'L3', 'pathway_id', 'gene', 'KO_name', 'EC']
OFFICIAL_SOURCES = {
    'pipeline': 'https://raw.githubusercontent.com/picrust/picrust2/v2.6.2/scripts/picrust2_pipeline.py',
    'outputs': 'https://raw.githubusercontent.com/picrust/picrust2/v2.6.2/picrust2/pipeline.py',
    'descriptions': 'https://raw.githubusercontent.com/picrust/picrust2/v2.6.2/scripts/add_descriptions.py',
    'qc': 'https://raw.githubusercontent.com/picrust/picrust2/v2.6.2/picrust2/metagenome_pipeline.py',
}
DEFAULTS = {'picrust_env': 'picrust2', 'min_align': 0.8, 'max_nsti': 2.0}


def _text_open(filename):
    return gzip.open(filename, 'rt', encoding='utf-8-sig', newline='') if str(filename).endswith('.gz') else Path(filename).open(encoding='utf-8-sig', newline='')


def _ko(value):
    value = re.sub(r'^ko:', '', str(value).strip(), flags=re.I)
    if not re.fullmatch(r'K\d{5}', value):
        raise ValueError('KO 必须是 K 加五位数字（允许 ko: 前缀）：' + value)
    return value


def read_mapping(filename):
    """Read TSV without treating commas inside names as separators."""
    with _text_open(filename) as stream:
        first = stream.readline()
        stream.seek(0)
        reader = csv.DictReader(stream, delimiter='\t' if '\t' in first else ',')
        header = reader.fieldnames or []
        if len(header) != len(MAPPING_COLUMNS) or set(header) != set(MAPPING_COLUMNS):
            raise ValueError('KO 注释映射列名必须恰好为 ' + ','.join(MAPPING_COLUMNS))
        result, seen = [], set()
        for line, record in enumerate(reader, 2):
            if None in record or any(record[key] is None for key in MAPPING_COLUMNS):
                raise ValueError(f'KO 映射第 {line} 行列数异常')
            row = {key: record[key].strip() for key in MAPPING_COLUMNS}
            row['KO'] = _ko(row['KO'])
            if any('\n' in value or '\r' in value or '\t' in value for value in row.values()):
                raise ValueError(f'KO 映射第 {line} 行含制表符或换行')
            if row['L2'] and not row['L1'] or row['L3'] and not row['L2']:
                raise ValueError(f'KO 映射第 {line} 行缺少层级父分类')
            key = tuple(row[name] for name in MAPPING_COLUMNS)
            if key not in seen:
                seen.add(key)
                result.append(row)
    if not result:
        raise ValueError('KO 注释映射为空')
    return result


def read_abundance(filename, sample_ids, *, ko=False, integer=False):
    """Read feature rows and sample columns; preserve identifiers as strings."""
    with _text_open(filename) as stream:
        records = csv.reader(stream, delimiter='\t')
        header = next((row for row in records if row and not row[0].startswith('# Constructed')), None)
        if header is None or len(header) < 2:
            raise ValueError('丰度表为空或缺少样本列：' + str(filename))
        actual = header[1:]
        if len(set(actual)) != len(actual) or set(actual) != set(sample_ids) or len(actual) != len(sample_ids):
            raise ValueError('丰度表 SampleID 与父预处理样本不一致：' + str(filename))
        order = [actual.index(sid) for sid in sample_ids]
        values = {}
        for line, row in enumerate(records, 2):
            if not row:
                continue
            if len(row) != len(header) or not row[0] or row[0] != row[0].strip():
                raise ValueError(f'丰度表第 {line} 行标识或列数异常：{filename}')
            feature = _ko(row[0]) if ko else row[0]
            if feature in values:
                raise ValueError('丰度表功能/ASV 标识重复：' + feature)
            try:
                numbers = [float(value) for value in row[1:]]
            except ValueError:
                raise ValueError('丰度表含非数值：' + str(filename)) from None
            if any(not math.isfinite(value) or value < 0 or integer and not value.is_integer() for value in numbers):
                raise ValueError('丰度表含负数、非有限数或无效序列计数：' + str(filename))
            values[feature] = [numbers[index] for index in order]
    if not values:
        raise ValueError('丰度表没有功能/ASV 行：' + str(filename))
    return {'samples': list(sample_ids), 'values': values}


def relative_table(table):
    totals = [sum(values[index] for values in table['values'].values()) for index in range(len(table['samples']))]
    if any(not math.isfinite(total) for total in totals):
        raise ValueError('丰度总和超出有限数值范围')
    return {'samples': table['samples'], 'values': {key: [value / totals[index] if totals[index] else 0.0 for index, value in enumerate(values)] for key, values in table['values'].items()}}, dict(zip(table['samples'], totals))


def write_table(filename, table, delimiter='\t'):
    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open('w', encoding='utf8', newline='') as stream:
        writer = csv.writer(stream, delimiter=delimiter)
        writer.writerow(['feature_id', *table['samples']])
        for feature, values in table['values'].items():
            writer.writerow([feature, *[format(number, '.15g') for number in values]])


def aggregate_categories(table, mapping):
    """Count a KO once per category, even when annotation memberships repeat."""
    memberships = {level: {} for level in ('L1', 'L2', 'L3', 'pathway')}
    annotations = {}
    for row in mapping:
        for level, fields in [('L1', ['L1']), ('L2', ['L1', 'L2']), ('L3', ['L1', 'L2', 'L3']), ('pathway', ['pathway_id'])]:
            if not row[fields[-1]]:
                continue
            parts = [row[field] for field in fields]
            category = parts[0] if len(parts) == 1 else json.dumps(parts, ensure_ascii=False, separators=(',', ':'))
            memberships[level].setdefault(category, set()).add(row['KO'])
            annotations[(level, category)] = {field: row[field] for field in fields}
    result = {}
    for level, groups in memberships.items():
        values = {}
        for category, members in sorted(groups.items()):
            available = sorted(members & table['values'].keys())
            if available:
                values[category] = [sum(table['values'][ko][index] for ko in available) for index in range(len(table['samples']))]
                if any(not math.isfinite(value) for value in values[category]):
                    raise ValueError('分类丰度总和超出有限数值范围')
        result[level] = {'samples': table['samples'], 'values': values}
    matched = {row['KO'] for row in mapping} & table['values'].keys()
    return result, annotations, sorted(table['values'].keys() - matched)


def _fasta_ids(filename):
    ids, sequence, current = set(), '', None
    with _text_open(filename) as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current is not None and not sequence:
                    raise ValueError('FASTA 包含空序列')
                current = line[1:]
                if not current or any(char.isspace() for char in current) or current in ids:
                    raise ValueError('FASTA ASV 标识为空、重复或含空白')
                ids.add(current)
                sequence = ''
            else:
                if current is None or not re.fullmatch(r'[ACGTRYSWKMBDHVNacgtryswkmbdhvn.-]+', line):
                    raise ValueError('FASTA 格式或碱基异常')
                sequence += line
    if not ids or not sequence:
        raise ValueError('FASTA 为空或末条序列为空')
    return ids


def _float_option(config, key, default, maximum=None):
    raw = config.get(key)
    if raw is None or not str(raw).strip():
        return default
    try:
        value = float(raw)
    except (ValueError, TypeError):
        raise ValueError(key + ' 必须是有限非负数') from None
    if isinstance(raw, bool) or not math.isfinite(value) or value < 0 or maximum is not None and value > maximum:
        raise ValueError(key + ' 必须是有限非负数' + (f'，不超过 {maximum}' if maximum is not None else ''))
    return value


def _write_json(filename, value):
    Path(filename).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf8')


def _nsti_rows(filename):
    with _text_open(filename) as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or 'metadata_NSTI' not in header:
            raise ValueError('官方 NSTI 表缺少 metadata_NSTI 列')
        result = {}
        for row in reader:
            sid = row[header[0]]
            if not sid or sid in result or None in row or any(row[key] is None for key in header):
                raise ValueError('官方 NSTI 表 ASV 标识或列数异常')
            value = float(row['metadata_NSTI'])
            if not math.isfinite(value) or value < 0:
                raise ValueError('官方 NSTI 表含无效数值')
            result[sid] = value
    return result


def _weighted_nsti(filename, sample_ids):
    with _text_open(filename) as stream:
        reader = csv.reader(stream, delimiter='\t')
        header = next(reader, [])
        if len(header) != 2 or header[1] != 'weighted_NSTI':
            raise ValueError('官方加权 NSTI 表格式异常')
        result = {}
        for row in reader:
            if len(row) != 2 or row[0] in result:
                raise ValueError('官方加权 NSTI 表样本重复或列数异常')
            value = float(row[1])
            if not math.isfinite(value) or value < 0:
                raise ValueError('官方加权 NSTI 表含无效数值')
            result[row[0]] = value
        if set(result) != set(sample_ids):
            raise ValueError('官方加权 NSTI SampleID 与父预处理不一致')
    return result


def collect_qc(raw, original, folder, max_nsti):
    """Derive exclusion counts only from official retained-ASV output."""
    raw, folder = Path(raw), Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    qc = {'input_asvs': len(original['values']), 'input_reads': dict(zip(original['samples'], [sum(row[i] for row in original['values'].values()) for i in range(len(original['samples']))])), 'weighted_nsti': {}, 'excluded_asvs': {'status': 'unavailable', 'reason': '官方 seqtab_norm.tsv.gz 文件不存在'}}
    for trait in ('KO', 'EC'):
        filename = raw / (trait + '_metagenome_out') / 'weighted_nsti.tsv.gz'
        qc['weighted_nsti'][trait] = {'status': 'available', 'per_sample': _weighted_nsti(filename, original['samples']), 'source': str(filename.relative_to(raw))} if filename.is_file() else {'status': 'unavailable', 'reason': '官方 weighted_nsti.tsv.gz 文件不存在'}
    marker = next((raw / name for name in ('combined_marker_predicted_and_nsti.tsv.gz', 'marker_predicted_and_nsti.tsv.gz') if (raw / name).is_file()), None)
    nsti = _nsti_rows(marker) if marker else {}
    qc['asv_nsti'] = {'status': 'available', 'source': str(marker.relative_to(raw)), 'asv_count': len(nsti), 'above_max_nsti': sum(nsti.get(asv, 0) > max_nsti for asv in original['values'])} if marker else {'status': 'unavailable', 'reason': '官方 ASV NSTI 表不存在'}
    normalized = raw / 'KO_metagenome_out' / 'seqtab_norm.tsv.gz'
    if normalized.is_file():
        retained = read_abundance(normalized, original['samples'])
        unknown = retained['values'].keys() - original['values'].keys()
        if unknown:
            raise ValueError('官方保留 ASV 表包含父输入不存在的 ASV')
        excluded = sorted(original['values'].keys() - retained['values'].keys())
        per_sample = {sid: sum(original['values'][asv][i] for asv in excluded) for i, sid in enumerate(original['samples'])}
        qc['excluded_asvs'] = {'status': 'available', 'source': str(normalized.relative_to(raw)), 'count': len(excluded), 'retained_count': len(retained['values']), 'excluded_input_reads': per_sample, 'fraction_of_input_reads': {sid: per_sample[sid] / qc['input_reads'][sid] if qc['input_reads'][sid] else None for sid in original['samples']}}
        with (folder / 'excluded_asvs.tsv').open('w', encoding='utf8', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t')
            writer.writerow(['ASV', 'metadata_NSTI', 'reason', *original['samples']])
            for asv in excluded:
                reason = 'NSTI_above_max' if asv in nsti and nsti[asv] > max_nsti else 'not_in_official_retained_table; specific_reason_unavailable'
                writer.writerow([asv, nsti.get(asv, ''), reason, *original['values'][asv]])
    _write_json(folder / 'qc.json', qc)
    return qc


def run_functional(job, config, saved, parent, api):
    """The app wrapper owns creation/marking of this fresh extension run."""
    run, parent = Path(job['run']).resolve(), Path(parent).resolve()
    if run == parent:
        raise ValueError('功能预测必须使用新的独立扩展目录')
    samples = [row['SampleID'] for row in saved.get('samples', [])]
    if not samples or len(set(samples)) != len(samples):
        raise ValueError('父预处理样本映射为空或重复')
    env = str(config.get('picrust_env') or '').strip() or DEFAULTS['picrust_env']
    if not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*', env):
        raise ValueError('PICRUSt2 环境名只能包含字母、数字、点、下划线和短横线')
    raw_threads = config.get('threads')
    raw_threads = str(raw_threads).strip() if raw_threads is not None else ''
    if not raw_threads:
        raw_threads = str(min(4, os.cpu_count() or 1))
    if not re.fullmatch(r'\d+', raw_threads) or not 1 <= int(raw_threads) <= max(1, os.cpu_count() or 1):
        raise ValueError('threads 必须是 1 至本机 CPU 数量之间的整数')
    if not str(config.get('ko_mapping') or '').strip():
        raise ValueError('请选择 KO 注释映射 TSV/CSV 文件；此文件不是 PICRUSt2 预测参考数据库')
    mapping_path = api.path(config['ko_mapping'])
    if not mapping_path.is_file():
        raise ValueError('请选择 KO 注释映射 TSV/CSV 文件；此文件不是 PICRUSt2 预测参考数据库')
    inputs = run / '00_inputs'
    inputs.mkdir(parents=True, exist_ok=True)
    mapping_copy = inputs / ('KO_mapping_original' + mapping_path.suffix)
    mapping_copy.write_bytes(mapping_path.read_bytes())
    mapping = read_mapping(mapping_copy)
    effective = {'picrust_env': env, 'threads': int(raw_threads), 'min_align': _float_option(config, 'min_align', 0.8, 1), 'max_nsti': _float_option(config, 'max_nsti', 2.0), 'in_traits': 'EC,KO', 'ko_mapping': str(mapping_path), 'ko_mapping_copy': str(mapping_copy.relative_to(run)), 'mapping_sha256': hashlib.sha256(mapping_copy.read_bytes()).hexdigest(), 'mapping_unique_rows': len(mapping), 'picrust_reference': '', 'picrust_reference2': ''}
    for key in ('picrust_reference', 'picrust_reference2'):
        if str(config.get(key) or '').strip():
            filename = api.path(config[key])
            if not filename.is_dir():
                raise ValueError('PICRUSt2 预测参考路径必须是目录：' + str(filename))
            effective[key] = str(filename)
    prefix = [api.find_conda({**saved, **config}), 'run', '--no-capture-output', '-n', env]
    api.begin_stage(job, '00_inputs', '功能预测输入与环境', effective)
    _write_json(inputs / 'effective_parameters.json', effective)
    _write_json(inputs / 'official_sources.json', {'validated_source_version': '2.6.2', 'urls': OFFICIAL_SOURCES, 'reference_meaning': 'PICRUSt2 phylogenetic references and built-in trait tables perform prediction; KO mapping only annotates predicted KO abundances.'})
    parameters_file = run / 'parameters.json'
    parameters = json.loads(parameters_file.read_text(encoding='utf8')) if parameters_file.is_file() else {}
    parameters['effective_parameters'] = effective
    _write_json(parameters_file, parameters)
    version_start = Path(job['log']).stat().st_size if Path(job['log']).exists() else 0
    api.command(job, prefix + ['picrust2_pipeline.py', '--version'], run)
    version_output = Path(job['log']).read_bytes()[version_start:].decode('utf8', errors='replace') if Path(job['log']).exists() else ''
    versions = re.findall(r'^PICRUSt2\s+(\S+)\s*$', version_output, re.M)
    effective['installed_version'] = versions[-1] if versions else 'unavailable; see recorded --version command output'
    effective['reference_configuration'] = {'ref_dir1': effective['picrust_reference'] or 'installed PICRUSt2 bacterial default', 'ref_dir2': effective['picrust_reference2'] or 'installed PICRUSt2 archaeal default', 'trait_tables': 'installed PICRUSt2 EC/KO defaults'}
    _write_json(inputs / 'effective_parameters.json', effective)
    parameters['effective_parameters'] = effective
    _write_json(parameters_file, parameters)
    help_start = Path(job['log']).stat().st_size if Path(job['log']).exists() else 0
    api.command(job, prefix + ['picrust2_pipeline.py', '--help'], run)
    help_output = Path(job['log']).read_bytes()[help_start:].decode('utf8', errors='replace') if Path(job['log']).exists() else ''
    for key, option in [('picrust_reference', '--ref_dir1'), ('picrust_reference2', '--ref_dir2')]:
        if effective[key] and not re.search(re.escape(option) + r'(?![\w])', help_output):
            raise ValueError('当前 PICRUSt2 CLI 不支持 ' + option + '；请使用 2.6.2 的预测参考选项')
    fasta = inputs / 'filtered-sequences.fasta'
    source = api.artifact_path(parent, 'exported-sequences/dna-sequences.fasta')
    if source.is_file():
        shutil.copyfile(source, fasta)
    else:
        export = inputs / 'sequence_export'
        api.command(job, api.conda_cmd(saved) + ['qiime', 'tools', 'export', '--input-path', str(api.artifact_path(parent, 'filtered-rep-seqs.qza')), '--output-path', str(export)], run)
        source = export / 'dna-sequences.fasta'
        if not source.is_file():
            raise ValueError('QIIME 导出缺少 dna-sequences.fasta')
        shutil.copyfile(source, fasta)
    biom = inputs / 'feature-table.biom'
    source = api.artifact_path(parent, 'exported-table/feature-table.biom')
    if source.is_file():
        shutil.copyfile(source, biom)
    else:
        export = inputs / 'table_export'
        api.command(job, api.conda_cmd(saved) + ['qiime', 'tools', 'export', '--input-path', str(api.artifact_path(parent, 'filtered-table.qza')), '--output-path', str(export)], run)
        source = export / 'feature-table.biom'
        if not source.is_file():
            raise ValueError('QIIME 导出缺少 feature-table.biom')
        shutil.copyfile(source, biom)
    table_tsv = inputs / 'feature-table.tsv'
    # Convert the actual copied BIOM rather than trusting a stale adjacent TSV.
    api.command(job, api.conda_cmd(saved) + ['biom', 'convert', '-i', str(biom), '-o', str(table_tsv), '--to-tsv'], run)
    original = read_abundance(table_tsv, samples, integer=True)
    if _fasta_ids(fasta) != original['values'].keys():
        raise ValueError('过滤后的代表序列 FASTA 与未抽平 ASV 表标识不一致')
    api.finish_stage(job, {'sample_ids': samples, 'input_asvs': len(original['values']), 'mapping_rows': len(mapping)}, ['输入来自父预处理过滤后的未抽平表与对应代表序列。KO 映射文件仅用于预测结果注释；预测参考由 PICRUSt2 自身提供。', 'min_align 默认 0.8，max_nsti 默认 2.0；不根据保留率自动放宽阈值。'])
    api.begin_stage(job, '01_prediction', 'PICRUSt2 KO、EC 与 MetaCyc 预测', effective)
    raw = run / '01_prediction' / 'picrust2_out'
    args = prefix + ['picrust2_pipeline.py', '-s', str(fasta), '-i', str(biom), '-o', str(raw), '-p', str(effective['threads']), '--min_align', str(effective['min_align']), '--max_nsti', str(effective['max_nsti']), '--in_traits', 'EC,KO']
    for key, option in [('picrust_reference', '--ref_dir1'), ('picrust_reference2', '--ref_dir2')]:
        if effective[key]:
            args.extend([option, effective[key]])
    api.command(job, args, run)
    sources = {'KO': raw / 'KO_metagenome_out/pred_metagenome_unstrat.tsv.gz', 'EC': raw / 'EC_metagenome_out/pred_metagenome_unstrat.tsv.gz', 'MetaCyc': raw / 'pathways_out/path_abun_unstrat.tsv.gz'}
    tables = {}
    for name, filename in sources.items():
        if not filename.is_file():
            raise ValueError('PICRUSt2 缺少预期输出：' + str(filename))
        tables[name] = read_abundance(filename, samples, ko=name == 'KO')
    api.finish_stage(job, {'raw_output': str(raw.relative_to(run)), 'features': {name: len(table['values']) for name, table in tables.items()}}, ['完整原始 PICRUSt2 输出及中间文件全部保留。16S 功能预测是推断，不能等同于宏基因组实测或表达量。MetaCyc 通路由 EC 与 MinPath 推断。'])
    api.begin_stage(job, '02_annotations', '官方功能描述与自定义 KO 注释')
    descriptions = run / '02_annotations'
    for name, filename in sources.items():
        api.command(job, prefix + ['add_descriptions.py', '-i', str(filename), '-m', 'METACYC' if name == 'MetaCyc' else name, '-o', str(descriptions / (name + '_descriptions.tsv.gz'))], run)
    with (descriptions / 'KO_annotations.tsv').open('w', encoding='utf8', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow([*MAPPING_COLUMNS, *samples])
        writer.writerows([*[row[key] for key in MAPPING_COLUMNS], *tables['KO']['values'][row['KO']]] for row in mapping if row['KO'] in tables['KO']['values'])
    api.finish_stage(job, {'custom_mapping_rows': len(mapping)}, ['同一 KO 的多条注释保留其分类成员关系；后续分类汇总中同一个 KO 在每个分类仅累加一次。'])
    api.begin_stage(job, '03_tables', 'KO、EC、MetaCyc 丰度与相对丰度')
    written, denominators = {}, {}
    for name, table in tables.items():
        relative, totals = relative_table(table)
        denominators[name] = {'definition': 'sum of all predicted ' + name + ' feature rows within each sample', 'per_sample': totals}
        written[name] = {}
        for kind, data in [('abundance', table), ('relative', relative)]:
            for suffix, delimiter in [('tsv', '\t'), ('csv', ',')]:
                filename = run / '03_tables' / (name + '_' + kind + '.' + suffix)
                write_table(filename, data, delimiter)
                written[name][kind + '_' + suffix] = str(filename.relative_to(run))
    _write_json(run / '03_tables/relative_denominators.json', denominators)
    api.finish_stage(job, {'tables': written, 'relative_denominators': denominators}, ['相对丰度逐样本除以该表所有功能行的丰度总和；总和为 0 的样本保持 0。预测计数允许小数。'])
    api.begin_stage(job, '04_ko_categories', 'KO 对 KEGG 注释分类的求和')
    categories, annotations, unmatched = aggregate_categories(tables['KO'], mapping)
    category_paths, category_denominators = {}, {}
    for level, table in categories.items():
        relative, totals = relative_table(table)
        category_denominators[level] = {'definition': 'sum of mapped KO membership category rows at ' + level + ' within each sample; includes overlapping memberships', 'per_sample': totals}
        category_paths['KEGG_' + level] = {}
        for kind, data in [('abundance', table), ('relative', relative)]:
            for suffix, delimiter in [('tsv', '\t'), ('csv', ',')]:
                filename = run / '04_ko_categories' / ('KEGG_' + level + '_' + kind + '.' + suffix)
                write_table(filename, data, delimiter)
                category_paths['KEGG_' + level][kind + '_' + suffix] = str(filename.relative_to(run))
    with (run / '04_ko_categories/category_annotations.tsv').open('w', encoding='utf8', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['level', 'feature_id', 'L1', 'L2', 'L3', 'pathway_id'])
        for (level, category), row in sorted(annotations.items()):
            if category in categories[level]['values']:
                writer.writerow([level, category, *[row.get(key, '') for key in ('L1', 'L2', 'L3', 'pathway_id')]])
    _write_json(run / '04_ko_categories/relative_denominators.json', category_denominators)
    api.finish_stage(job, {'tables': category_paths, 'unmatched_KO_count': len(unmatched), 'relative_denominators': category_denominators}, ['KEGG L1/L2/L3 与 pathway 是映射到分类的 KO 丰度求和，不能当作 MinPath 推断通路。', '多成员 KO 会计入多个分类，分类总和可能超过 KO 原表总量。相对丰度分母是各分类汇总表逐样本总和。L2/L3 标识包括父级，避免同名分类合并。'])
    api.begin_stage(job, '05_qc', 'NSTI、ASV/read 排除与未映射 KO 记录')
    qc = collect_qc(raw, original, run / '05_qc', effective['max_nsti'])
    write_table(run / '05_qc/unmatched_KO.tsv', {'samples': samples, 'values': {ko: tables['KO']['values'][ko] for ko in unmatched}})
    api.finish_stage(job, qc, ['排除 reads 用输入未抽平计数和官方 seqtab_norm 中的保留 ASV 集合计算。未提供的官方文件明确报告 unavailable，不虚构排除率或 NSTI。', '缺少具体排除原因的 ASV 仅标记未进入官方保留表，不推测为比对失败。官方加权 NSTI 是 PICRUSt2 输出；不自行替代其定义。'])
    job['result'] = {'functional_run':str(run), 'extension': 'functional', 'preprocess': str(parent), 'effective_parameters': effective, 'tables': {**written, **category_paths}, 'qc': qc, 'unmatched_KO_count': len(unmatched), 'prediction_method': 'PICRUSt2 inferred KO/EC; MetaCyc inferred from EC/MinPath; KEGG categories are summed mapped KO memberships', 'official_sources': OFFICIAL_SOURCES}
    return job['result']
