"""Cautious DADA2 candidates from QIIME demux seven-number summaries.

These quality and read-length heuristics do not identify universal optimal
parameters. Summaries must describe reads after confirmed primer handling.
"""
import csv
import math
import os
from pathlib import Path
from statistics import median


_STATS = ("count", "2%", "9%", "25%", "50%", "75%", "91%", "98%")
_COVERAGE = 0.99
_Q25 = 30
_WINDOW = 10
_MIN_CANDIDATE = 50
_MIN_OVERLAP = 12


def _number(value, description):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{description}: expected a numeric value.") from None
    if not math.isfinite(number) or number < 0:
        raise ValueError(f"{description}: values must be finite and nonnegative.")
    return number


def _summary(path, direction):
    try:
        with Path(path).open(encoding="utf-8-sig", newline="") as stream:
            rows = [row for row in csv.reader(stream, delimiter="\t") if row]
    except (OSError, TypeError, UnicodeError, csv.Error) as exc:
        raise ValueError(f"{direction}: cannot read the seven-number TSV: {exc}") from exc
    if len(rows) < 2 or len(rows[0]) < 2:
        raise ValueError(f"{direction}: the seven-number TSV is empty or incomplete.")
    rows = [[cell.strip() for cell in row] for row in rows]
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError(f"{direction}: TSV rows have different lengths.")

    # QIIME's exported demux summaries use positions as columns. Accept a
    # transposed copy as well, but never guess statistics from an unknown label.
    if set(rows[0][1:]) == set(_STATS):
        if len(rows[0][1:]) != len(_STATS):
            raise ValueError(f"{direction}: duplicate statistic labels.")
        positions = [row[0] for row in rows[1:]]
        columns = {
            label: [row[index] for row in rows[1:]]
            for index, label in enumerate(rows[0][1:], 1)
        }
    else:
        positions = rows[0][1:]
        labels = [row[0] for row in rows[1:]]
        if len(labels) != len(set(labels)):
            raise ValueError(f"{direction}: duplicate statistic labels.")
        if set(labels) != set(_STATS):
            raise ValueError(
                f"{direction}: expected count and 2%, 9%, 25%, 50%, 75%, 91%, 98% rows."
            )
        columns = {row[0]: row[1:] for row in rows[1:]}
    if positions != [str(position) for position in range(1, len(positions) + 1)]:
        raise ValueError(f"{direction}: positions must be consecutive, starting at 1.")
    columns = {
        label: [_number(value, f"{direction} {label} at position {i}")
                for i, value in enumerate(values, 1)]
        for label, values in columns.items()
    }
    counts = columns["count"]
    if any(count <= 0 or not count.is_integer() for count in counts):
        raise ValueError(f"{direction}: sampled-read counts must be positive integers.")
    if any(later > earlier for earlier, later in zip(counts, counts[1:])):
        raise ValueError(f"{direction}: sampled-read counts must not increase with position.")
    for position in range(len(positions)):
        quantiles = [columns[label][position] for label in _STATS[1:]]
        if quantiles != sorted(quantiles):
            raise ValueError(f"{direction}: unordered quality percentiles at position {position + 1}.")

    coverage_limit = max(i for i, count in enumerate(counts, 1)
                         if count >= _COVERAGE * counts[0])
    deterioration_start = None
    window_median = None
    for index in range(coverage_limit - _WINDOW + 1):
        quality = median(columns["25%"][index:index + _WINDOW])
        if quality < _Q25:
            deterioration_start, window_median = index + 1, quality
            break
    candidate = min(coverage_limit, deterioration_start - 1
                    if deterioration_start is not None else coverage_limit)
    if candidate < _MIN_CANDIDATE:
        raise ValueError(
            f"{direction}: the cautious candidate is only {candidate} bases; "
            "inspect the quality profile and choose parameters manually. "
            "This heuristic requires at least 50 bases, not a universal DADA2 limit."
        )
    return candidate, {
        "summary_positions": len(positions),
        "sampled_reads": int(counts[0]),
        "coverage_limit": coverage_limit,
        "quality_window_start": deterioration_start,
        "quality_window_median_q25": window_median,
        "selected_truncation": candidate,
        "retained_read_fraction": counts[candidate - 1] / counts[0],
    }


def recommend_parameters(forward_tsv, reverse_tsv=None, *, amplicon=None, cpu_count=None):
    """Return JSON-compatible candidates, evidence and actionable warnings.

    ``amplicon`` is the known positive length after primer removal, in bases;
    it is never inferred from quality summaries. Invalid summaries raise
    ``ValueError``. The 50-base floor only limits this heuristic's scope.
    """
    length = None
    if amplicon not in (None, ""):
        if isinstance(amplicon, bool) or not isinstance(amplicon, (str, int, float)):
            raise ValueError("Amplicon length must be a known positive integer or blank.")
        length = _number(amplicon, "Amplicon length")
        if length <= 0 or not length.is_integer():
            raise ValueError("Amplicon length must be a known positive integer or blank.")
        length = int(length)
    forward, forward_evidence = _summary(forward_tsv, "Forward")
    reverse, reverse_evidence = "", None
    if reverse_tsv is not None:
        reverse, reverse_evidence = _summary(reverse_tsv, "Reverse")
    try:
        threads = max(1, min(4, int(cpu_count if cpu_count is not None else os.cpu_count() or 1)))
    except (ValueError, TypeError, OverflowError):
        raise ValueError("CPU count must be an integer.") from None

    warnings = []
    overlap = None
    sufficient = None
    requires_review = False
    if reverse_evidence is not None:
        if length is None:
            warnings.append(
                "扩增子长度未知，质量统计不能推断拼接重叠。可留空运行，但须核对 DADA2 实际拼接保留率；有可靠的去引物后长度时请填写。"
            )
            requires_review = True
        else:
            overlap = forward + reverse - length
            sufficient = overlap >= _MIN_OVERLAP
            if not sufficient:
                warnings.append(
                    f"预计重叠为 {overlap} bp，低于 {_MIN_OVERLAP} bp。请核对实际扩增子长度和质量图并手动调整；推荐没有自动延伸到低质量区域。"
                )
                requires_review = True
    return {
        "parameters": {
            "trunc_f": forward, "trunc_r": reverse, "trim_f": 0, "trim_r": 0,
            "threads": threads, "min_frequency": 1, "max_depth": "",
            "amplicon": amplicon if length is not None else "",
        },
        "evidence": {
            "forward": forward_evidence, "reverse": reverse_evidence,
            "overlap": {"amplicon_length": length, "predicted_bases": overlap,
                        "required_bases": _MIN_OVERLAP, "sufficient": sufficient},
            "criteria": {"read_coverage": _COVERAGE, "q25_threshold": _Q25,
                         "window_bases": _WINDOW, "minimum_candidate_bases": _MIN_CANDIDATE},
            "primer_handling": "统计必须对应已确认的引物处理结果，额外 trim 默认 0；不从质量值判断是否去引物。",
            "depth": "曲线最大深度将在过滤后的特征表生成后按实际样本计数确定。",
        },
        "warnings": warnings,
        "method": (
            "保守启发式，不是官方固定阈值或保证最优参数：保留至少 99% 抽样 reads 覆盖的位置，"
            "在首个连续 10 个位置的 Q25 中位数低于 Q30 的窗口前截断；实际效果须用 DADA2 质控和拼接保留率验证。"
        ),
        "requires_review": requires_review,
    }
