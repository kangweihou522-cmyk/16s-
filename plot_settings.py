"""Presentation-only options; default None retains each plotter's prior style."""
from copy import deepcopy
import json
import math
from pathlib import Path

GENERAL = {'dpi': None, 'font_family': 'original', 'base_size': None, 'title_size': None, 'axis_text_size': None,
           'legend_text_size': None, 'legend_position': 'original', 'theme': 'original',
           'background': 'white', 'show_caption': True, 'formats': 'original'}
COMMON = {'width': None, 'height': None, 'title': '', 'x_label': '', 'y_label': ''}
DEFAULTS = {
    'general': GENERAL,
    'alpha': {**COMMON, 'point_size': None, 'box_width': None, 'jitter_width': None,
              'show_points': True, 'show_letters': True, 'show_outliers': False, 'plot_type': 'boxplot'},
    'beta': {**COMMON, 'point_size': None, 'point_alpha': None, 'show_ellipse': True,
             'ellipse_level': .95, 'ellipse_alpha': .15, 'show_sample_labels': False,
             'show_stats': True, 'panel_columns': 2},
    'taxonomy': {**COMMON, 'top_n': 15, 'levels': [2, 6], 'view': 'both',
                 'bar_width': .85, 'palette': 'original'},
    'heatmap': {**COMMON, 'label_size': None, 'row_order': 'abundance', 'show_values': False},
    'bubble': {**COMMON, 'max_size': None, 'point_alpha': .9, 'show_values': False},
    'effect': {**COMMON, 'point_size': None, 'show_errorbars': True},
    'lefse_bar': {**COMMON, 'feature_font_size': None, 'class_font_size': None,
                  'label_levels': 2, 'top_n': 0, 'orientation': 'horizontal'},
    'lefse_cladogram': {**COMMON, 'label_font_size': None, 'max_level': None,
                        'min_point_size': None, 'max_point_size': None},
}
RANGES = {
    'dpi': (72, 600, True), 'base_size': (6, 30, False), 'title_size': (6, 40, False),
    'axis_text_size': (5, 30, False), 'legend_text_size': (5, 30, False),
    'width': (2, 30, False), 'height': (2, 40, False),
    'point_size': (.1, 12, False), 'box_width': (.1, 1, False), 'jitter_width': (0, .5, False),
    'point_alpha': (.05, 1, False), 'ellipse_level': (.5, .99, False),
    'ellipse_alpha': (0, 1, False), 'panel_columns': (1, 4, True),
    'top_n': (1, 100, True), 'bar_width': (.1, 1, False), 'label_size': (5, 24, False),
    'max_size': (1, 20, False), 'feature_font_size': (5, 24, False),
    'class_font_size': (5, 24, False), 'label_levels': (1, 7, True),
    'label_font_size': (5, 24, False), 'max_level': (1, 7, True),
    'min_point_size': (.1, 20, False), 'max_point_size': (.1, 30, False),
}
ENUMS = {
    'font_family': ('original', 'DejaVu Serif', 'DejaVu Sans', 'Times New Roman'),
    'legend_position': ('original', 'right', 'bottom', 'none'),
    'theme': ('original', 'classic', 'minimal', 'bw'), 'background': ('white', 'transparent'),
    'formats': ('original', 'png_svg', 'png_svg_pdf'), 'plot_type': ('boxplot', 'violin'),
    'view': ('both', 'sample', 'group'), 'palette': ('original', 'Dark 3', 'Set 2', 'Viridis'),
    'row_order': ('abundance', 'alphabetical'), 'orientation': ('horizontal', 'vertical'),
}


def resolve_plot_settings(config, fallback=None):
    fallback = fallback or {}
    previous = fallback.get('plot_settings', fallback if 'general' in fallback else {})
    requested = config.get('plot_settings', {})
    if not isinstance(previous, dict) or not isinstance(requested, dict):
        raise ValueError('绘图参数必须是分图种保存的对象')
    if set(requested) - set(DEFAULTS): raise ValueError('包含未知绘图类别')
    result = deepcopy(DEFAULTS)
    for profile, defaults in DEFAULTS.items():
        for values in (previous.get(profile, {}), requested.get(profile, {})):
            if not isinstance(values, dict) or set(values) - set(defaults):
                raise ValueError('绘图参数无效或包含未知字段：' + profile)
            result[profile].update(deepcopy(values))
        for key, value in result[profile].items():
            label = profile + '.' + key
            default = defaults[key]
            if value in (None, '') and default is None:
                result[profile][key] = None
            elif isinstance(default, bool):
                if not isinstance(value, bool): raise ValueError(label + ' 必须是开关值')
            elif key in RANGES:
                low, high, whole = RANGES[key]
                if profile == 'lefse_bar' and key == 'top_n': low, high = 0, 300
                try:
                    if isinstance(value, bool): raise ValueError()
                    number = float(value)
                    if not math.isfinite(number) or not low <= number <= high or whole and number != int(number):
                        raise ValueError()
                except (TypeError, ValueError):
                    raise ValueError(label + f' 必须为 {low}–{high} 范围内的有效数值') from None
                result[profile][key] = int(number) if whole else number
            elif key in ENUMS:
                if value not in ENUMS[key]: raise ValueError('无效绘图选项：' + label)
            elif key == 'levels':
                if not isinstance(value, list) or not value or any(isinstance(v, bool) or not isinstance(v, int) or v not in range(2, 8) for v in value) or len(set(value)) != len(value):
                    raise ValueError('分类柱图层级请选择 L2–L7，不能重复或全不选')
            elif not isinstance(value, str) or len(value) > 300 or '\x00' in value:
                raise ValueError(label + ' 必须是最多300字的标签文本')
    points = result['lefse_cladogram']
    minimum = points['min_point_size'] if points['min_point_size'] is not None else 1
    maximum = points['max_point_size'] if points['max_point_size'] is not None else 6
    if minimum >= maximum:
        raise ValueError('LEfSe 分支图最小节点必须小于最大节点；留空时分别使用1和6')
    return result


def save_plot_settings(run, settings):
    filename = Path(run) / 'plot_settings.json'
    filename.write_text(json.dumps(settings, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
    return filename


def ui_schema():
    titles = {'general': '通用导出与文字', 'alpha': 'α 多样性图', 'beta': 'β 多样性 / PCoA',
              'taxonomy': '分类丰度柱图', 'heatmap': '差异菌丰度热图', 'bubble': '差异菌丰度气泡图',
              'effect': '差异效应 / ANCOM W 图', 'lefse_bar': 'LEfSe LDA 柱图',
              'lefse_cladogram': 'LEfSe 分类分支图'}
    labels = {'dpi': 'PNG 分辨率 DPI（留空原版）', 'font_family': '字体（原版沿用各图字体）',
              'base_size': '基础字号（留空原版）', 'title_size': '标题字号', 'axis_text_size': '坐标刻度字号',
              'legend_text_size': '图例字号', 'legend_position': '图例位置', 'theme': '绘图主题',
              'background': '背景', 'show_caption': '显示说明文字', 'formats': '导出格式',
              'width': '图宽 / 英寸（留空自动）', 'height': '图高 / 英寸（留空自动）',
              'title': '自定义标题（留空原版）', 'x_label': '自定义横轴标题', 'y_label': '自定义纵轴标题',
              'point_size': '散点大小', 'box_width': '箱体宽度', 'jitter_width': '散点横向抖动宽度',
              'show_points': '显示样本散点', 'show_letters': '显示已保存检验的显著性字母',
              'show_outliers': '显示箱线图离群点', 'plot_type': '图形类型', 'point_alpha': '散点透明度',
              'show_ellipse': '显示置信椭圆', 'ellipse_level': '椭圆置信水平', 'ellipse_alpha': '椭圆填充透明度',
              'show_sample_labels': '显示样本 ID', 'show_stats': '显示已保存的检验注释',
              'panel_columns': 'PCoA 合图列数', 'top_n': '最多展示数量', 'levels': '分类层级（可多选）',
              'view': '展示单位', 'bar_width': '柱体宽度', 'palette': '分类单元配色',
              'label_size': '菌名标签字号', 'row_order': '菌名行顺序', 'show_values': '显示平均相对丰度数值（%）',
              'max_size': '最大气泡大小', 'show_errorbars': '显示原始模型 SE 误差线',
              'feature_font_size': '菌名标签字号', 'class_font_size': '分组图例字号',
              'label_levels': '标签显示的谱系层数', 'orientation': '柱图方向',
              'label_font_size': '分支标签字号', 'max_level': '最深显示层级（留空全部）',
              'min_point_size': '最小节点大小', 'max_point_size': '最大节点大小'}
    choices = {'original': '沿用原版', 'right': '右侧', 'bottom': '下方', 'none': '隐藏',
               'classic': '简洁坐标轴', 'minimal': '简约', 'bw': '黑白网格', 'white': '白色',
               'transparent': '透明', 'png_svg': 'PNG + SVG', 'png_svg_pdf': 'PNG + SVG + PDF',
               'boxplot': '箱线图', 'violin': '小提琴图', 'both': '样本图及分组图',
               'sample': '每个样本', 'group': '分组均值', 'abundance': '平均丰度从高到低',
               'alphabetical': '按菌名排序', 'horizontal': '横向', 'vertical': '纵向'}
    rows=[]
    for profile, defaults in DEFAULTS.items():
        fields=[]
        for key, default in defaults.items():
            field={'key':key,'label':labels[key],'default':default}
            if isinstance(default,bool): field['type']='checkbox'
            elif key in ENUMS:
                field.update(type='select',options=[{'value':v,'label':choices.get(v,v)} for v in ENUMS[key]])
            elif key in RANGES:
                lo,hi,whole=RANGES[key]
                if profile=='lefse_bar' and key=='top_n': lo,hi=0,300;field['label']='最多显示菌数（0 表示全部）'
                field.update(type='number',min=lo,max=hi,step=1 if whole else 'any',nullable=default is None)
            elif key=='levels':
                field.update(type='levels',options=[{'value':v,'label':'L'+str(v)+' '+name} for v,name in zip(range(2,8),('门','纲','目','科','属','种'))])
            else: field.update(type='text',maxlength=300)
            fields.append(field)
        rows.append({'key':profile,'title':titles[profile],'fields':fields})
    return rows
