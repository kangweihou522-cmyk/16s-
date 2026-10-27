"""Prepare a separate downstream source without modifying external input files.

QIIME reads/writes artifacts; Python validates literal IDs and count semantics.
No FASTQ, filtering, denoising or preprocessing completion is implied here.
"""
import csv
import hashlib
import json
import re
import shutil
import stat
import time
import uuid
import zipfile
from decimal import Decimal, InvalidOperation
from pathlib import Path, PurePosixPath

from abundance_tables import LEVELS, TABLE_HEADERS, read_taxonomy, write_rows, write_table_pair
from functional_prediction import _fasta_ids


MODULE_REQUIREMENTS = {
    'diversity': ['table', 'metadata', 'rooted_tree'],
    'abundance': ['table', 'metadata', 'taxonomy'],
    'differential': ['table', 'metadata', 'taxonomy'],
    'lefse': ['table', 'metadata', 'taxonomy'],
    'functional': ['table', 'metadata', 'sequences'],
    'advanced': ['table', 'metadata', 'taxonomy'],
}
PATHS = {
    'metadata.tsv': '00_inputs/metadata.tsv',
    'filtered-table.qza': '01_external/filtered-table.qza',
    'exported-table': '01_external/exported-table',
    'taxonomy.qza': '01_external/taxonomy.qza',
    'exported-taxonomy': '01_external/exported-taxonomy',
    'taxa-barplot.qzv': '01_external/taxa-barplot.qzv',
    'filtered-rep-seqs.qza': '01_external/filtered-rep-seqs.qza',
    'exported-sequences': '01_external/exported-sequences',
    'rooted-tree.qza': '01_external/rooted-tree.qza',
    'exports/phylogeny': '01_external/exported-tree',
}
_NONRAW = re.compile(r'relative[ _-]*(?:frequency|abundance)|raref(?:y|ied|action)|'
                     r'normaliz(?:ed|ation)|normalis(?:ed|ation)|proportions?|percent(?:age)?', re.I)


def _json(filename, value):
    Path(filename).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf8')


def _hash(filename):
    digest = hashlib.sha256()
    with Path(filename).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _literal(value, label):
    if (not value or value != value.strip() or any(ord(char) < 32 for char in value)
            or value.startswith('#q2:')):
        raise ValueError(label + ' 不能为空、两端含空格或含控制字符；标识按原文精确匹配')
    return value


def read_counts(filename):
    """Read feature-by-sample TSV/CSV, preserving IDs such as 001 verbatim."""
    delimiter = ',' if Path(filename).suffix.lower() == '.csv' else '\t'
    with Path(filename).open(encoding='utf-8-sig', newline='') as stream:
        rows = [row for row in csv.reader(stream, delimiter=delimiter) if row]
    start = next((index for index, row in enumerate(rows) if row[0] in TABLE_HEADERS), None)
    if start is None:
        raise ValueError('external_table 须为特征×样本原始整数计数表；表头示例：FeatureID,001,002；每行示例：ASV1,12,34（TSV 使用制表符）')
    if any(_NONRAW.search(' '.join(row)) for row in rows[:start]):
        raise ValueError('计数表说明标记为相对丰度、标准化或抽平数据；请提供未抽平原始计数')
    if any(not row[0].startswith('#') for row in rows[:start]):
        raise ValueError('计数表表头前只能有 # 开头的说明行')
    ids = [_literal(value, 'SampleID') for value in rows[start][1:]]
    if not ids or len(ids) != len(set(ids)):
        raise ValueError('计数表 SampleID 为空或重复')
    features, values, seen = [], [], set()
    for row in rows[start + 1:]:
        if len(row) != len(ids) + 1 or row[0] in seen:
            raise ValueError('计数表 FeatureID 重复或列数异常')
        feature = _literal(row[0], 'FeatureID')
        numbers = []
        for raw in row[1:]:
            try:
                number = Decimal(raw)
            except InvalidOperation:
                raise ValueError('计数表只能包含非负整数原始计数') from None
            if not number.is_finite() or number < 0 or number != number.to_integral_value():
                raise ValueError('计数表只能包含非负整数原始计数；不接受相对丰度')
            if number > 2 ** 53:
                raise ValueError('单个计数超出 BIOM 可精确表示范围（2^53），不能无损导入')
            numbers.append(int(number))
        features.append(feature); values.append(numbers); seen.add(feature)
    if not features:
        raise ValueError('计数表没有特征')
    return ids, features, values


def read_metadata(filename):
    with Path(filename).open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ['SampleID', 'Group']:
            raise ValueError('外部输入分组 CSV 表头须为 SampleID,Group；请勿把 FASTQ 路径列作为分组输入')
        samples, seen = [], set()
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError('分组 CSV 每行须恰好两列')
            sid = _literal(row['SampleID'], 'SampleID')
            group = _literal(row['Group'], 'Group')
            if sid in seen:
                raise ValueError('分组 CSV SampleID 重复：' + sid)
            seen.add(sid); samples.append({'SampleID': sid, 'Group': group})
    if not samples:
        raise ValueError('分组 CSV 为空')
    return samples


def artifact_type(filename, expected, *, raw_counts=False):
    """Inspect standard archive metadata without extracting any member."""
    try:
        with zipfile.ZipFile(filename) as archive:
            infos = archive.infolist()
            names = [info.filename for info in infos]
            if len(names) != len(set(names)):
                raise ValueError('QZA 存在重复归档路径')
            for info in infos:
                name = info.filename
                parts = PurePosixPath(name).parts
                if (name.startswith('/') or '\\' in name or '..' in parts or
                        any(part in {'.', ''} for part in name.rstrip('/').split('/')) or
                        any(ord(char) < 32 for char in name) or not parts or ':' in parts[0] or
                        stat.S_ISLNK(info.external_attr >> 16)):
                    raise ValueError('QZA 存在不安全归档路径')
            roots = {PurePosixPath(name).parts[0] for name in names}
            if len(roots) != 1 or not re.fullmatch(r'[0-9a-fA-F-]{36}', next(iter(roots))):
                raise ValueError('QZA 不是标准 QIIME artifact 目录结构')
            root = next(iter(roots))
            metadata = archive.getinfo(root + '/metadata.yaml')
            if metadata.file_size > 65536:
                raise ValueError('QZA metadata.yaml 大小异常')
            content = archive.read(metadata).decode('utf8')
            matches = re.findall(r'^type:\s*[\'"]?([^\r\n\'"]+)[\'"]?\s*$', content, re.M)
            if len(matches) != 1 or matches[0].strip() != expected:
                raise ValueError('QZA 类型须为 ' + expected + '；当前为 ' + (matches[0] if matches else '未知'))
            if raw_counts:
                total = 0
                for info in infos:
                    if not info.filename.endswith('/action/action.yaml'):
                        continue
                    total += info.file_size
                    if info.file_size > 2 * 1024 * 1024 or total > 32 * 1024 * 1024:
                        raise ValueError('QZA provenance 大小异常，无法核验原始计数')
                    action = archive.read(info).decode('utf8')
                    if re.search(r'^\s*action:\s*[\'"]?(?:rarefy|relative[_-]frequency|normalize|normalise)[\'"]?\s*$', action, re.M | re.I):
                        raise ValueError('QZA provenance 显示抽平或相对/标准化处理；请提供原始未抽平计数')
    except (zipfile.BadZipFile, KeyError, UnicodeError) as error:
        raise ValueError('无法读取标准 QZA 元数据：' + str(filename)) from error
    return expected


def _exact(actual, expected, label):
    if set(actual) != set(expected) or len(actual) != len(set(actual)):
        missing, extra = sorted(set(expected) - set(actual)), sorted(set(actual) - set(expected))
        raise ValueError(label + ' 必须精确匹配；缺少：' + ', '.join(missing[:8]) + '；额外：' + ', '.join(extra[:8]))


def _coverage(actual, expected, label):
    if len(actual) != len(set(actual)) or set(expected) - set(actual):
        missing = sorted(set(expected) - set(actual))
        raise ValueError(label + ' 缺少计数表 FeatureID 或存在重复：' + ', '.join(missing[:8]))
    return {'input_feature_count': len(actual), 'count_feature_count': len(expected),
            'unused_feature_ids': sorted(set(actual) - set(expected))}


def _taxonomy(filename, features, alignment):
    columns, rows = read_taxonomy(filename)
    if set(columns) - {'Taxon', 'Confidence'}:
        raise ValueError('分类 TSV 只接受 FeatureID、Taxon 和可选 Confidence 列')
    for feature in rows:
        _literal(feature, '分类 FeatureID')
    alignment.update(_coverage(list(rows), features, '分类 FeatureID 与计数表'))
    canonical = ['Taxon'] + (['Confidence'] if 'Confidence' in columns else [])
    return canonical, {feature: [rows[feature][columns.index(name)] for name in canonical]
                       for feature in features}


def _prepare_fasta(fasta, features, alignment, run):
    alignment.update(_coverage(list(_fasta_ids(fasta)), features, 'FASTA FeatureID 与计数表'))
    if not alignment['unused_feature_ids']:
        return False
    original = run / '00_inputs/sequences_before_subset.fasta'
    shutil.copyfile(fasta, original)
    keep = set(features)
    with original.open(encoding='utf-8-sig') as source, fasta.open('w', encoding='utf8') as target:
        include = False
        for line in source:
            if line.strip().startswith('>'):
                include = line.strip()[1:] in keep
            if include:
                target.write(line)
    _exact(list(_fasta_ids(fasta)), features, '准备后的 FASTA FeatureID 与计数表')
    return True


def _require_file(api, config, key, suffixes):
    if not str(config.get(key) or '').strip():
        raise ValueError('所选模块需要填写 ' + key)
    filename = api.path(config[key])
    if not filename.is_file() or filename.suffix.lower() not in suffixes:
        raise ValueError(key + ' 请选择存在的 ' + '/'.join(sorted(suffixes)) + ' 文件')
    return filename


def _rank_manifest(run, counts, annotations, config, api):
    """Aggregate integer counts by QIIME's lineage labels for advanced analysis."""
    columns, taxonomy = annotations
    taxon_index = columns.index('Taxon')
    lineages = {feature: [part.strip() for part in row[taxon_index].split(';')]
                for feature, row in taxonomy.items()}
    maximum = max(map(len, lineages.values()))
    raw_level = config.get('advanced_level', 6)
    if isinstance(raw_level, bool) or not re.fullmatch(r'[2-6]', str(raw_level)):
        raise ValueError('advanced_level 须为 2 至 6 的整数')
    if int(raw_level) > maximum:
        raise ValueError('分类表深度不足所选高级分析层级 L' + str(raw_level))
    levels = []
    for level, name in LEVELS:
        if level > maximum:
            continue
        aggregated = {}
        for feature, numbers in zip(counts[1], counts[2]):
            lineage = lineages[feature]
            label = ';'.join((lineage + ['__'] * (maximum - len(lineage)))[:level])
            target = aggregated.setdefault(label, [0] * len(counts[0]))
            for index, value in enumerate(numbers):
                target[index] += value
        paths = write_table_pair(run / '10_abundance' / ('L' + str(level) + '_' + name), 'counts',
                                 (counts[0], list(aggregated), list(aggregated.values())))
        levels.append({'level': level, 'name': name, 'counts': {key: str(value.relative_to(run)) for key, value in paths.items()}})
    manifest = {'format': 1, 'source': 'external_confirmed_unrarefied', 'source_run': str(run),
                'aggregation': 'sum integer feature counts by taxonomy lineage; no filtering',
                'source_sha256': {'filtered-table.tsv': _hash(api.artifact_path(run, 'exported-table/feature-table.tsv')),
                                  'taxonomy.tsv': _hash(api.artifact_path(run, 'exported-taxonomy/taxonomy.tsv'))},
                'levels': levels}
    _json(run / '10_abundance/manifest.json', manifest)


def _import_taxonomy(job, saved, api, run, annotations, features):
    columns, rows = annotations
    headerless = run / '00_inputs/taxonomy_import.tsv'
    # Headerless QIIME import retains optional Confidence without rounding it.
    with headerless.open('w', encoding='utf8', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerows([feature, *rows[feature]] for feature in features)
    output = api.artifact_path(run, 'taxonomy.qza')
    if output.exists():
        output.unlink()  # A derived copy in this new run, never an external file.
    api.q(job, saved, 'tools', 'import', '--type', 'FeatureData[Taxonomy]', '--input-format', 'HeaderlessTSVTaxonomyFormat',
          '--input-path', headerless, '--output-path', output)


def prepare_external(job, config, api):
    """Create a fresh marked parent source; return UI-safe paths and capabilities."""
    module = config.get('external_module')
    if module not in MODULE_REQUIREMENTS:
        raise ValueError('请选择独立分析模块 external_module')
    if config.get('external_counts_confirmed') is not True:
        raise ValueError('必须确认 external_counts_confirmed：输入是原始、未抽平的非负整数计数；仅凭数值无法排除抽平表')
    sources = {'table': _require_file(api, config, 'external_table', {'.qza', '.tsv', '.csv'}),
               'metadata': _require_file(api, config, 'metadata', {'.csv'})}
    if 'taxonomy' in MODULE_REQUIREMENTS[module]:
        sources['taxonomy'] = _require_file(api, config, 'external_taxonomy', {'.qza', '.tsv'})
    if module == 'diversity':
        sources['rooted_tree'] = _require_file(api, config, 'external_tree', {'.qza'})
    if module == 'functional':
        sources['sequences'] = _require_file(api, config, 'external_sequences', {'.qza', '.fasta', '.fa', '.fna'})
    types = {'table': 'FeatureTable[Frequency]', 'taxonomy': 'FeatureData[Taxonomy]',
             'rooted_tree': 'Phylogeny[Rooted]', 'sequences': 'FeatureData[Sequence]'}
    samples = read_metadata(sources['metadata'])
    # Validate archive structure and semantic type before any official export.
    for key, source in sources.items():
        if source.suffix.lower() == '.qza':
            artifact_type(source, types[key], raw_counts=key == 'table')
    output = api.path(config.get('output', ''))
    if output.exists() and not output.is_dir():
        raise ValueError('输出位置必须是目录')
    output.mkdir(parents=True, exist_ok=True)
    run = output / ('external_' + module + '_' + time.strftime('%Y%m%d_%H%M%S') + '_' + uuid.uuid4().hex[:8])
    run.mkdir()
    api.marker(run)
    _json(run / 'layout.json', {'version': 1, 'kind': 'external-inputs', 'paths': PATHS})
    job.update(run=str(run), log=str(run / 'workflow.log'))
    (run / '00_inputs').mkdir()
    (run / '01_external').mkdir()
    saved = dict(config, source_kind='external', supported_modules=[module], samples=samples)
    saved.setdefault('env', 'qiime2-amplicon-2025.7')
    saved.setdefault('r_env', '16s-workstation-r')
    _json(run / 'parameters.json', saved)
    api.begin_stage(job, '01_external', '独立模块输入核验与导入', {'module': module, 'raw_unrarefied_counts_confirmed': True})
    reports, copies = {}, {}
    for key, source in sources.items():
        target = run / '00_inputs' / (key + '_original' + source.suffix.lower())
        before = _hash(source)
        shutil.copyfile(source, target)
        if _hash(target) != before:
            raise ValueError('复制过程中输入发生变化，请重新准备：' + str(source))
        reports[key] = {'source': str(source), 'copy': str(target.relative_to(run)),
                        'sha256': before, 'size': target.stat().st_size}
        copies[key] = target
    samples = read_metadata(copies['metadata'])
    saved['samples'] = samples
    for key, copied in copies.items():
        if copied.suffix == '.qza':
            artifact_type(copied, types[key], raw_counts=key == 'table')
    _json(run / 'parameters.json', saved)
    _json(run / '00_inputs/source_report.json', {'sources': reports, 'counts_confirmation': True,
          'limitation': 'Integer values cannot establish absence of prior rarefaction; user confirmation is required.'})
    table_qza = api.artifact_path(run, 'filtered-table.qza')
    table_export = api.artifact_path(run, 'exported-table')
    table_tsv = table_export / 'feature-table.tsv'
    if copies['table'].suffix == '.qza':
        shutil.copyfile(copies['table'], table_qza)
        api.q(job, saved, 'tools', 'export', '--input-path', table_qza, '--output-path', table_export)
        api.command(job, api.conda_cmd(saved) + ['biom', 'convert', '-i', str(table_export / 'feature-table.biom'),
                    '-o', str(table_tsv), '--to-tsv'], run)
        counts = read_counts(table_tsv)
    else:
        counts = read_counts(copies['table'])
        _exact([row['SampleID'] for row in samples], counts[0], '分组 SampleID 与计数表')
        table_export.mkdir()
        write_rows(table_tsv, ['#OTU ID', *counts[0]], [[feature, *numbers] for feature, numbers in zip(counts[1], counts[2])])
        api.command(job, api.conda_cmd(saved) + ['biom', 'convert', '-i', str(table_tsv),
                    '-o', str(table_export / 'feature-table.biom'), '--to-hdf5', '--table-type', 'OTU table'], run)
        api.q(job, saved, 'tools', 'import', '--type', types['table'], '--input-format', 'BIOMV210Format',
              '--input-path', table_export / 'feature-table.biom', '--output-path', table_qza)
    _exact([row['SampleID'] for row in samples], counts[0], '分组 SampleID 与计数表')
    # Keep the source metadata unchanged; save the effective mapping in table order.
    groups = {row['SampleID']: row['Group'] for row in samples}
    saved['samples'] = [{'SampleID': sid, 'Group': groups[sid]} for sid in counts[0]]
    write_rows(api.artifact_path(run, 'metadata.tsv'), ['sample-id', 'Group'],
               [['#q2:types', 'categorical'], *[[sid, groups[sid]] for sid in counts[0]]])
    annotations, alignments = None, {}
    if 'taxonomy' in copies:
        taxonomy_qza = api.artifact_path(run, 'taxonomy.qza')
        taxonomy_export = api.artifact_path(run, 'exported-taxonomy')
        taxonomy_tsv = taxonomy_export / 'taxonomy.tsv'
        if copies['taxonomy'].suffix == '.qza':
            shutil.copyfile(copies['taxonomy'], taxonomy_qza)
            api.q(job, saved, 'tools', 'export', '--input-path', taxonomy_qza, '--output-path', taxonomy_export)
            original_taxonomy = run / '00_inputs/taxonomy_original_export.tsv'
            shutil.copyfile(taxonomy_tsv, original_taxonomy)
        else:
            original_taxonomy = copies['taxonomy']
        alignments['taxonomy'] = {}
        annotations = _taxonomy(original_taxonomy, counts[1], alignments['taxonomy'])
        columns, rows = annotations
        write_rows(taxonomy_tsv, ['Feature ID', *columns], [[feature, *rows[feature]] for feature in counts[1]])
        padding = []
        if module == 'abundance':
            columns, rows = annotations
            for feature, row in rows.items():
                lineage = row[0].split(';')
                if len(lineage) < 7:
                    original = row[0]
                    row[0] = ';'.join(lineage + ['__'] * (7 - len(lineage)))
                    padding.append({'FeatureID': feature, 'original': original, 'derived': row[0],
                                    'unknown_trailing_ranks_added': 7 - len(lineage)})
            _json(run / '00_inputs/taxonomy_padding_report.json', {
                'meaning': 'Only absent trailing ranks receive explicit unknown __ placeholders; no taxonomic classification is inferred.',
                'features': padding})
            if padding:
                shutil.copyfile(taxonomy_tsv, run / '00_inputs/taxonomy_before_padding.tsv')
                write_rows(taxonomy_tsv, ['Feature ID', *columns], [[feature, *rows[feature]] for feature in counts[1]])
        if copies['taxonomy'].suffix != '.qza' or padding or alignments['taxonomy']['unused_feature_ids']:
            _import_taxonomy(job, saved, api, run, annotations, counts[1])
        if module == 'abundance':
            api.q(job, saved, 'taxa', 'barplot', '--i-table', table_qza, '--i-taxonomy', taxonomy_qza,
                  '--m-metadata-file', api.artifact_path(run, 'metadata.tsv'), '--o-visualization', api.artifact_path(run, 'taxa-barplot.qzv'))
        if module == 'advanced':
            _rank_manifest(run, counts, annotations, config, api)
    if 'sequences' in copies:
        sequence_qza = api.artifact_path(run, 'filtered-rep-seqs.qza')
        sequence_export = api.artifact_path(run, 'exported-sequences')
        fasta = sequence_export / 'dna-sequences.fasta'
        if copies['sequences'].suffix == '.qza':
            shutil.copyfile(copies['sequences'], sequence_qza)
            api.q(job, saved, 'tools', 'export', '--input-path', sequence_qza, '--output-path', sequence_export)
        else:
            sequence_export.mkdir()
            shutil.copyfile(copies['sequences'], fasta)
        alignments['sequences'] = {}
        subset = _prepare_fasta(fasta, counts[1], alignments['sequences'], run)
        if copies['sequences'].suffix != '.qza' or subset:
            if sequence_qza.exists():
                sequence_qza.unlink()
            api.q(job, saved, 'tools', 'import', '--type', types['sequences'], '--input-path', fasta, '--output-path', sequence_qza)
        _exact(list(_fasta_ids(fasta)), counts[1], 'FASTA FeatureID 与计数表')
    if 'rooted_tree' in copies:
        tree_qza = api.artifact_path(run, 'rooted-tree.qza')
        shutil.copyfile(copies['rooted_tree'], tree_qza)
        tree_export = api.artifact_path(run, 'exports/phylogeny')
        api.q(job, saved, 'tools', 'export', '--input-path', tree_qza, '--output-path', tree_export)
        tips_json = tree_export / 'tip_ids.json'
        parser = ('import json,sys; from skbio import TreeNode; '
                  'tree=TreeNode.read(sys.argv[1],format="newick",convert_underscores=False); '
                  'json.dump([tip.name for tip in tree.tips()],open(sys.argv[2],"w"))')
        api.command(job, api.conda_cmd(saved) + ['python', '-c', parser, str(tree_export / 'tree.nwk'), str(tips_json)], run)
        tips = json.loads(tips_json.read_text(encoding='utf8'))
        if not isinstance(tips, list) or any(not isinstance(tip, str) for tip in tips):
            raise ValueError('有根树 tip 标识缺失或格式异常')
        alignments['rooted_tree'] = _coverage(tips, counts[1], '有根树 tip FeatureID 与计数表')
    _json(run / '00_inputs/feature_alignment_report.json', alignments)
    for key, source in sources.items():
        if _hash(source) != reports[key]['sha256'] or _hash(copies[key]) != reports[key]['sha256']:
            raise ValueError('准备过程中源文件发生变化，请重新准备：' + str(source))
    required = [table_qza, table_tsv, table_export / 'feature-table.biom', api.artifact_path(run, 'metadata.tsv')]
    if annotations is not None:
        required += [api.artifact_path(run, 'taxonomy.qza'), api.artifact_path(run, 'exported-taxonomy/taxonomy.tsv')]
    if module == 'abundance':
        required.append(api.artifact_path(run, 'taxa-barplot.qzv'))
    if module == 'advanced':
        required.append(run / '10_abundance/manifest.json')
    if module == 'diversity':
        required.append(api.artifact_path(run, 'rooted-tree.qza'))
    if module == 'functional':
        required += [api.artifact_path(run, 'filtered-rep-seqs.qza'), api.artifact_path(run, 'exported-sequences/dna-sequences.fasta')]
    if any(not target.is_file() for target in required):
        raise ValueError('独立输入导入/导出缺少必需文件，任务未完成')
    _json(run / 'parameters.json', saved)
    summary = {'sample_count': len(counts[0]), 'feature_count': len(counts[1]),
               'per_sample_counts': {sid: sum(row[index] for row in counts[2]) for index, sid in enumerate(counts[0])}}
    api.finish_stage(job, summary, ['输入源文件只读；原件副本与 SHA256 保存在 00_inputs。',
                     '本任务仅准备所选独立模块输入，未执行 FASTQ、过滤或 DADA2。',
                     '未抽平原始计数由用户确认；整数值本身不能证明从未抽平。'])
    completion = {'completed_at': time.time(), 'source_kind': 'external', 'supported_modules': [module],
                  'requirements': {module: MODULE_REQUIREMENTS[module]}, 'sources': reports, 'summary': summary,
                  'feature_alignment': alignments,
                  'raw_unrarefied_counts_confirmed': True}
    _json(run / 'external_inputs_complete.json', completion)
    result = {'external_run': str(run), 'external_module': module, 'supported_modules': [module], **summary}
    job['result'] = result
    job['step'] = '独立模块输入准备完成；可执行所选分析'
    return result
