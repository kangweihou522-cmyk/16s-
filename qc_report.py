"""Report observed DADA2 stage retention; warnings never select samples."""
import csv
import json
import math
from pathlib import Path


LOW_RETENTION = 0.5


def write_dada2_retention(directory, summary, *, paired):
    """Export adjacent-stage fractions from existing counts, without reanalysis."""
    rows = []
    warnings = []
    count_names = ['input', 'filtered', 'denoised'] + (['merged'] if paired else []) + ['non-chimeric']
    for source in summary['per_sample']:
        row = {'SampleID': str(source['SampleID'])}
        for key in count_names:
            value = source[key]
            if isinstance(value, bool) or not math.isfinite(value) or value < 0 or int(value) != value:
                raise ValueError('DADA2 保留率报告需要非负整数阶段计数')
            row[key] = int(value)
        stages = [('filter_retention_fraction', 'filtered', 'input', '质量过滤')]
        if paired:
            stages.append(('merge_retention_fraction', 'merged', 'denoised', '双端拼接'))
        stages.append(('nonchimera_retention_fraction', 'non-chimeric', 'merged' if paired else 'denoised', '去嵌合体'))
        alerts = []
        for key, numerator, denominator, label in stages:
            fraction = row[numerator] / row[denominator] if row[denominator] else None
            row[key] = fraction
            if fraction is not None and fraction < LOW_RETENTION:
                alerts.append(f'{label}保留率 {fraction:.1%} 低于排查提示阈值 {LOW_RETENTION:.0%}')
        row['overall_retention_fraction'] = row['non-chimeric'] / row['input'] if row['input'] else None
        if row['non-chimeric'] == 0:
            alerts.append('去嵌合体后为零，请检查该样本的各阶段损失')
        row['warnings'] = alerts
        warnings.extend(f'{row["SampleID"]}：{message}' for message in alerts)
        rows.append(row)
    report = {
        'paired': paired, 'count_unit': 'read_pairs_or_merged_sequences' if paired else 'reads',
        'criteria': {'low_adjacent_stage_retention_fraction': LOW_RETENTION,
                     'type': 'workstation_heuristic_warning', 'automatically_excludes_samples': False},
        'fraction_definitions': {
            'filter_retention_fraction': 'filtered / input',
            'merge_retention_fraction': 'merged / denoised (paired only)',
            'nonchimera_retention_fraction': 'non-chimeric / ' + ('merged' if paired else 'denoised'),
            'overall_retention_fraction': 'non-chimeric / input',
        },
        'notes': ['比例范围0–1；分母为零时留空/null，不虚构0%或100%。',
                  '50%仅为工作站排查提示，不是官方统一合格线；不自动删除任何样本。',
                  '质量过滤低时核对长度、maxEE和质量；拼接低时核对方向、重叠与残留附加序列；非嵌合体保留低时核对PCR与文库结构。'],
        'per_sample': rows, 'warnings': warnings,
    }
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    json_path = directory / 'dada2_retention.json'
    csv_path = directory / 'dada2_retention.csv'
    json_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
    fractions = ['filter_retention_fraction'] + (['merge_retention_fraction'] if paired else []) + ['nonchimera_retention_fraction', 'overall_retention_fraction']
    with csv_path.open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['SampleID', *count_names, *fractions, 'warnings'])
        writer.writeheader()
        for row in rows:
            writer.writerow(dict(row, warnings='；'.join(row['warnings'])))
    return {'csv': str(csv_path), 'json': str(json_path), 'sample_count': len(rows),
            'warnings': warnings, 'criteria': report['criteria']}
