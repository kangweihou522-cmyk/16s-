"""Minimum-read depth advice and separate saved-curve diagnostics."""
import csv
import io
import math
import re
import zipfile
from pathlib import Path
from statistics import mean

THRESHOLDS = {'observed_features': 0.05, 'shannon': 0.02, 'faith_pd': 0.05}
TARGET = 0.90


def read_curves(folder=None, visualization=None):
    curves = {}
    for metric in THRESHOLDS:
        filename = Path(folder) / (metric + '.csv') if folder else None
        content = None
        if filename and filename.is_file():
            content = filename.read_text(encoding='utf-8-sig')
        elif visualization and Path(visualization).is_file():
            with zipfile.ZipFile(visualization) as archive:
                matches = [name for name in archive.namelist() if name.endswith('/data/' + metric + '.csv')]
                if len(matches) == 1:
                    content = archive.read(matches[0]).decode('utf-8-sig')
        if content is None: continue
        reader = csv.DictReader(io.StringIO(content))
        columns = [(key, int(match.group(1))) for key in reader.fieldnames or []
                   if (match := re.fullmatch(r'depth-(\d+)_iter-\d+', key))]
        if not columns or not reader.fieldnames or reader.fieldnames[0] != 'sample-id':
            raise ValueError(f'{metric} 抽平曲线 CSV 格式无法识别')
        samples = {}
        for row in reader:
            sid = row['sample-id']
            if not sid or sid in samples or None in row:
                raise ValueError(f'{metric} 抽平曲线样本或列数异常')
            values = {}
            for key, depth in columns:
                raw = row[key]
                if raw is None: raise ValueError(f'{metric} 抽平曲线缺少列值')
                if raw.strip().lower() in ('', 'nan', 'na'): continue
                value = float(raw)
                if not math.isfinite(value) or value < 0 or depth < 1:
                    raise ValueError(f'{metric} 抽平曲线包含无效数值')
                values.setdefault(depth, []).append(value)
            samples[sid] = {depth: mean(items) for depth, items in values.items()}
        curves[metric] = samples
    return curves


def recommend_depth(samples, counts, curves=None):
    """Preserve every original sample; curves are diagnostic only."""
    curves = curves or {}
    groups = {}
    for sample in samples:
        sid, group = sample['SampleID'], sample['Group']
        count = counts.get(sid, 0)
        if isinstance(count, bool) or not math.isfinite(count) or count < 0 or int(count) != count:
            raise ValueError('抽平推荐需要非负整数样本计数')
        groups.setdefault(group, []).append((sid, int(count)))
    known = {sid for members in groups.values() for sid, _ in members}
    if any(set(data) - known for data in curves.values()):
        raise ValueError('抽平曲线中的样本与本批分组不一致')
    warnings = []
    required = {group: len(members) for group, members in groups.items()}
    diagnostic_required = {group: max(2, math.ceil(TARGET * len(members))) for group, members in groups.items()}
    minimum_count = min((count for members in groups.values() for _, count in members), default=0)
    usable_groups = len(groups) >= 2 and all(len(members) >= 2 for members in groups.values())
    selected = ((minimum_count - 1) // 100) * 100 if minimum_count > 100 and usable_groups else None
    if not usable_groups:
        status = 'insufficient_samples'
        warnings.append('无法推荐正式抽平深度：需要至少2个组，且每组至少2个原始样本。请人工检查分组设计后再决定分析方案。')
    elif minimum_count <= 100:
        status = 'insufficient_depth'
        warnings.append(f'无法推荐正式抽平深度：原始样本中过滤后最小reads为{minimum_count}，不存在严格小于它的正100倍数；不能自动丢弃低计数或零计数样本。请人工检查过滤结果和低测序量样本。')
    else:
        status = 'minimum_read_rule'
    caps = []
    for group, members in groups.items():
        ordered = sorted((count for _, count in members), reverse=True)
        caps.append(ordered[diagnostic_required[group] - 1] if len(ordered) >= diagnostic_required[group] else 0)
    cap = min(caps) if caps and len(groups) >= 2 else 0
    curve_candidate = None
    evidence = None
    if cap >= 1 and set(THRESHOLDS).issubset(curves):
        available = [set(depth for values in curves[metric].values() for depth in values) for metric in THRESHOLDS]
        depths = sorted(set.intersection(*available))
        for index in range(len(depths) - 2):
            window = depths[index:index + 3]
            base = window[0]
            if base > cap or window[-1] < base * 1.2: continue
            supported, evaluated = {}, {}
            for group, members in groups.items():
                stable = tested = 0
                for sid, count in members:
                    if count < window[-1]: continue
                    vectors = [curves[metric].get(sid, {}) for metric in THRESHOLDS]
                    if any(any(depth not in vector for depth in window) for vector in vectors): continue
                    tested += 1
                    if all(max(abs(vector[depth] - vector[base]) / max(abs(vector[base]), 1e-12)
                               for depth in window[1:]) <= THRESHOLDS[metric]
                           for metric, vector in zip(THRESHOLDS, vectors)):
                        stable += 1
                supported[group], evaluated[group] = stable, tested
            if all(supported[group] >= diagnostic_required[group] for group in groups):
                curve_candidate = base
                evidence = {'depths': window, 'stable_samples_by_group': supported,
                            'evaluated_samples_by_group': evaluated}
                break
    if not set(THRESHOLDS).issubset(curves):
        warnings.append('缺少完整的三项抽平曲线数据，无法完成曲线平台诊断；请人工查看已有曲线。曲线诊断不改变最小reads规则的推荐值。')
    elif curve_candidate is None:
        warnings.append('现有曲线未找到同时满足三点平台启发式及各组90%留样的诊断候选深度；请人工查看曲线。曲线诊断不改变最小reads规则的推荐值。')
    retention = [{'SampleID': sid, 'Group': group, 'count': count,
                  'retained': selected is not None and count >= selected}
                 for group, members in groups.items() for sid, count in members]
    summary = [{'Group': group, 'original': len(members),
                'retained': sum(row['retained'] for row in retention if row['Group'] == group),
                'required': required[group]} for group, members in groups.items()]
    return {'recommended_depth': selected, 'status': status, 'groups': summary,
            'retention': retention, 'retention_depth_cap': minimum_count or None, 'curve_evidence': evidence,
            'curve_plateau_candidate_depth': curve_candidate,
            'criteria': {'target_group_retention': 1.0, 'minimum_samples_per_group': 2,
                         'minimum_groups': 2, 'depth_multiple': 100, 'strictly_below_minimum_reads': True,
                         'curve_target_group_retention': TARGET,
                         'curve_relative_change_limits': THRESHOLDS, 'curve_points': 3,
                         'minimum_depth_span_ratio': 1.2},
            'method': '正式推荐取所有原始样本中过滤后最小reads严格以下的最大正100倍数，保留全部样本。三点曲线平台启发式仅作独立诊断，不改变推荐深度，也不证明推荐深度这个具体点已通过曲线验证；诊断阈值不是官方统一最佳深度标准。',
            'warnings': warnings}
