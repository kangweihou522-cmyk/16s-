"""Named group colors and display order, independent of statistical contrasts."""
import json
import re
from pathlib import Path

DEFAULT_COLORS = {'CON': '#D55E00', 'LDC': '#0072B2', 'HDC': '#CC79A7'}
PALETTE = ('#56B4E9', '#009E73', '#E69F00', '#882255', '#44AA99',
           '#332288', '#999933', '#AA4499', '#88CCEE', '#661100')


def resolve_group_style(config, samples, fallback=None):
    groups = list(dict.fromkeys(sample['Group'] for sample in samples))
    if any(not isinstance(group, str) or not group or group != group.strip() for group in groups):
        raise ValueError('绘图分组名称缺失或无效')
    fallback = fallback or {}
    order = config.get('group_order', fallback.get('group_order', []))
    colors = config.get('group_colors', fallback.get('group_colors', {}))
    if not isinstance(order, list) or any(not isinstance(g, str) or not g for g in order) or len(set(order)) != len(order):
        raise ValueError('绘图分组顺序必须是没有重复名称的列表')
    if not isinstance(colors, dict) or any(not isinstance(g, str) or not isinstance(c, str)
                                         or not re.fullmatch(r'#[0-9A-Fa-f]{6}', c) for g, c in colors.items()):
        raise ValueError('分组配色必须使用 #RRGGBB 十六进制颜色')
    order = [group for group in order if group in groups] + [group for group in groups if group not in order]
    defaults = dict(DEFAULT_COLORS)
    extra = [group for group in groups if group not in defaults]
    defaults.update({group: PALETTE[i % len(PALETTE)] for i, group in enumerate(extra)})
    return {'group_order': order,
            'group_colors': {group: colors.get(group, defaults[group]).upper() for group in order}}


def save_group_style(run, style):
    filename = Path(run) / 'group_style.json'
    filename.write_text(json.dumps(style, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
    return filename
