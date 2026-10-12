"""Versioned result paths and factual per-stage reports; legacy runs stay readable."""
import csv
import json
import math
from datetime import datetime, timezone, timedelta
from pathlib import Path

PREPROCESS = {
    'manifest.tsv': '00_inputs/manifest.tsv', 'metadata.tsv': '00_inputs/metadata.tsv',
    **{name: '01_import/' + name for name in ('demux.qza', 'demux.qzv')},
    'quality-raw': '01_import/quality_tables',
    **{name: '02_primers/' + name for name in ('trimmed.qza', 'trimmed.qzv')},
    'quality-trimmed': '02_primers/quality_tables',
    **{name: '03_dada2/' + name for name in ('dada2-table.qza', 'dada2-table.qzv', 'dada2-rep-seqs.qza', 'dada2-stats.qza', 'dada2-stats.qzv')},
    'exports/denoising-stats': '03_dada2/denoising_stats',
    'exported-dada2-table': '03_dada2/asv_table',
    **{name: '04_taxonomy_filter/' + name for name in ('taxonomy.qza', 'taxonomy.qzv', 'filtered-step1.qza', 'filtered-table.qza', 'filtered-table.qzv', 'filtered-rep-seqs.qza', 'taxa-barplot.qzv', 'exported-table', 'exported-taxonomy', 'exported-sequences')},
    **{name: '05_phylogeny/' + name for name in ('aligned.qza', 'masked.qza', 'unrooted-tree.qza', 'rooted-tree.qza')},
    'exports/phylogeny': '05_phylogeny/exported_tree',
    'alpha-rarefaction.qzv': '06_rarefaction/alpha-rarefaction.qzv',
}
DIVERSITY = {
    'metadata.tsv': '00_inputs/metadata.tsv', 'retention.json': '00_inputs/retention.json',
    'core-metrics-results': '01_core_metrics/core-metrics-results', 'exports': '02_exports',
    'exported-rarefied-table': '02_exports/rarefied_table',
    'alpha': '03_statistics/alpha', 'beta': '03_statistics/beta',
    'statistics_status.json': '03_statistics/statistics_status.json',
    'R_original_style': '04_R_figures',
}


def initialize_layout(run, kind):
    run = Path(run)
    mapping = PREPROCESS if kind == 'preprocess' else DIVERSITY
    (run / 'layout.json').write_text(json.dumps({'version': 1, 'kind': kind, 'paths': mapping}, ensure_ascii=False, indent=2), encoding='utf8')
    for destination in mapping.values():
        p = run / destination
        p.parent.mkdir(parents=True, exist_ok=True)


def artifact_path(run, logical):
    run, logical = Path(run), str(logical)
    filename = run / 'layout.json'
    mapping = json.loads(filename.read_text(encoding='utf8'))['paths'] if filename.is_file() else {}
    for prefix in sorted(mapping, key=len, reverse=True):
        if logical == prefix or logical.startswith(prefix + '/'):
            destination = Path(mapping[prefix])
            if destination.is_absolute() or '..' in destination.parts:
                raise ValueError('结果目录布局包含无效路径')
            return run / destination / logical[len(prefix):].lstrip('/') if logical != prefix else run / destination
    return run / logical


def mapped_argument(run, value):
    """Only rewrite known relative output paths, never classifier/input paths."""
    value = str(value)
    if value.startswith('-') or Path(value).is_absolute():
        return value
    resolved = artifact_path(run, value)
    return str(resolved.relative_to(run))


def read_numeric_rows(filename):
    with Path(filename).open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.reader(f, delimiter='\t'))
    rows = [row for row in rows if row and not row[0].startswith('#q2:')]
    if len(rows) < 2:
        raise ValueError('统计表为空：' + str(filename))
    names, result, seen = rows[0], [], set()
    for row in rows[1:]:
        if len(row) != len(names) or not row[0] or row[0] in seen:
            raise ValueError('统计表样本或列数异常：' + str(filename))
        seen.add(row[0])
        item = {'SampleID': row[0]}
        for name, value in zip(names[1:], row[1:]):
            try:
                number = float(value)
            except ValueError:
                raise ValueError('统计表包含非数值计数：' + str(filename)) from None
            if not math.isfinite(number) or number < 0:
                raise ValueError('统计表包含无效数值：' + str(filename))
            item[name] = int(number) if number.is_integer() else number
        result.append(item)
    return result


def quality_counts(folder):
    rows = read_numeric_rows(Path(folder) / 'per-sample-fastq-counts.tsv')
    keys = [key for key in rows[0] if key != 'SampleID']
    if not any('forward' in key.lower() for key in keys):
        raise ValueError('质量统计表缺少正向 reads 数量')
    if any(not float(row[key]).is_integer() for row in rows for key in keys):
        raise ValueError('质量统计表序列数量必须为整数')
    return {'sample_count': len(rows), 'totals': {key: sum(row.get(key, 0) for row in rows) for key in keys}, 'per_sample': rows}


def denoising_counts(folder):
    candidates = list(Path(folder).glob('*.tsv'))
    if len(candidates) != 1:
        raise ValueError('DADA2 统计导出中应有一个 TSV')
    rows = read_numeric_rows(candidates[0])
    keys = [key for key in ('input', 'filtered', 'denoised', 'merged', 'non-chimeric') if key in rows[0]]
    if not all(key in keys for key in ('input', 'filtered', 'non-chimeric')):
        raise ValueError('DADA2 统计表缺少质控或非嵌合体列')
    totals = {key: sum(row[key] for row in rows) for key in keys}
    return {'sample_count': len(rows), 'totals': totals, 'per_sample': rows,
            'final_retention_fraction': totals['non-chimeric'] / totals['input'] if totals['input'] else None}


def table_summary(filename):
    with Path(filename).open(encoding='utf8', newline='') as f:
        rows = list(csv.reader(f, delimiter='\t'))
    headers = [i for i, row in enumerate(rows) if row and row[0] in ('#OTU ID', '#Feature ID', 'Feature ID', 'feature-id')]
    if len(headers) != 1:
        raise ValueError('特征表缺少或重复 BIOM TSV 表头')
    header = rows[headers[0]]
    ids = header[1:]
    if not ids or len(set(ids)) != len(ids):
        raise ValueError('特征表样本列无效')
    counts, seen = dict.fromkeys(ids, 0), set()
    for row in rows[headers[0] + 1:]:
        if not row: continue
        if len(row) != len(header) or row[0] in seen:
            raise ValueError('特征表列数或 ASV 标识异常')
        seen.add(row[0])
        for sid, value in zip(ids, row[1:]):
            number = float(value)
            if not math.isfinite(number) or number < 0 or not number.is_integer():
                raise ValueError('特征表包含无效序列计数')
            counts[sid] += int(number)
    return {'asv_count': len(seen), 'sample_count': len(ids), 'total_sequences': sum(counts.values()), 'per_sample': counts}


def stage_report(run, folder, title, *, status='done', parameters=None, results=None, notes=()):
    run = Path(run)
    directory = run / folder
    directory.mkdir(parents=True, exist_ok=True)
    outputs = [str(p.relative_to(run)) for p in sorted(directory.rglob('*')) if p.is_file() and p.name not in ('说明.md', 'report.json')]
    report = {'stage': folder, 'title': title, 'status': status,
              'updated_at': datetime.now(timezone(timedelta(hours=8))).isoformat(timespec='seconds'),
              'parameters': parameters or {}, 'results': results or {}, 'notes': list(notes), 'outputs': outputs}
    (directory / 'report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf8')
    text = [f'# {title}', '', f'状态：{status}', f'记录时间：{report["updated_at"]}', '', '## 参数', '',
            '```json', json.dumps(report['parameters'], ensure_ascii=False, indent=2), '```', '', '## 处理结果', '',
            '```json', json.dumps(report['results'], ensure_ascii=False, indent=2), '```', '', '## 说明', '']
    text.extend('- ' + str(note) for note in notes)
    text.extend(['', '## 保存文件', ''])
    text.extend('- ' + name for name in outputs)
    (directory / '说明.md').write_text('\n'.join(text) + '\n', encoding='utf8')
    reports = []
    for p in sorted(run.glob('*/report.json')):
        item = json.loads(p.read_text(encoding='utf8'))
        reports.append(f'- {item["title"]}：{item["status"]}，详见 `{item["stage"]}/说明.md`')
    (run / '处理说明.md').write_text('# 本次处理说明\n\n' + '\n'.join(reports) + '\n\n统计来自本批实际输出；未运行或失败的步骤不代表分析完成。\n', encoding='utf8')
    return report
