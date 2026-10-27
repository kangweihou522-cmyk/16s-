"""Describe saved figures and redraw one image without rerunning an analysis."""
import hashlib
import json
import os
import re
import time
import uuid
from pathlib import Path

from differential_analysis import plot_options
from group_style import resolve_group_style, save_group_style
from plot_settings import resolve_plot_settings, save_plot_settings

PREVIEW_FOLDERS = {'04_previews', '05_replots', '03_replots', '12_figure_previews'}
ALPHA = {'observed_features', 'shannon', 'evenness', 'faith_pd', 'chao1', 'simpson', 'goods_coverage'}
BETA = {'bray_curtis', 'jaccard', 'unweighted_unifrac', 'weighted_unifrac', 'four_metrics'}
RANKS = ('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
LABELS = {'alpha': 'α 多样性', 'beta': 'PCoA', 'taxonomy': '分类丰度',
          'heatmap': '差异菌丰度热图', 'bubble': '差异菌丰度气泡图', 'effect': '差异效应',
          'lefse_bar': 'LEfSe LDA 柱图', 'lefse_cladogram': 'LEfSe 分类分支图'}


def _json(filename, default=None):
    return json.loads(Path(filename).read_text(encoding='utf8')) if Path(filename).is_file() else (default or {})


def _save(filename, value):
    Path(filename).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf8')


def _owner(image, api):
    for parent in image.parents:
        if (parent / 'workstation_run.json').is_file():
            return api.approved_run(str(parent))
    raise ValueError('图片不在工作站结果目录内')


def describe_figure(image, api):
    image = api.approved_file(image)
    if image.suffix.lower() != '.png':
        raise ValueError('请选择工作站生成的 PNG 原图')
    owner = _owner(image, api)
    previous = _json(owner / 'parameters.json')
    config = previous.get('config', previous)
    record = _json(owner / 'figure_request.json')
    profile = stem = kind = None
    settings = _json(owner / 'plot_settings.json') or config.get('plot_settings', {})
    summary = {}
    for directory in image.parents:
        if directory == owner.parent:
            break
        report = directory / 'render_summary.json'
        if not report.is_file():
            continue
        value = _json(report)
        relative = image.relative_to(directory).with_suffix('').as_posix()
        if value.get('source') == 'saved_QIIME2_outputs':
            if re.fullmatch(r'alpha/[^/]+', relative) and image.stem in ALPHA:
                profile = 'alpha'
            elif relative.startswith('beta/PCoA_') and relative.removeprefix('beta/PCoA_') in BETA:
                profile = 'beta'
            elif re.fullmatch(r'taxonomy/(phylum|class|order|family|genus|species)_(sample|group_mean)_unrarefied', relative):
                profile = 'taxonomy'
            if profile:
                kind, stem, summary = 'diversity', relative, value
                analysis = record.get('analysis_run') or value.get('input_run') or previous.get('diversity_run') or owner
                settings = value.get('plot_settings') or settings
                break
        elif value.get('source') == 'advanced_analysis':
            relative = image.relative_to(owner).with_suffix('').as_posix()
            entries = value.get('figures', [])
            if isinstance(entries, dict): entries = [entries]
            matching = [item for item in entries if item.get('stem') == relative]
            if len(matching) == 1 and matching[0].get('profile') in {'beta', 'heatmap', 'bubble'}:
                kind, profile, stem, summary = 'advanced', matching[0]['profile'], relative, value
                analysis = record.get('analysis_run') or value.get('advanced_run') or previous.get('advanced_run') or owner
                settings = value.get('plot_settings') or settings
                break
        elif directory.name == '03_figures' and image.stem in ('bar', 'cladogram') and value.get('lefse_run'):
            kind, profile, stem, summary = 'lefse', 'lefse_' + image.stem, '03_figures/' + image.stem, value
            analysis = record.get('analysis_run') or value['lefse_run']
            settings = value.get('plot_settings') or settings
            break
    if profile is None:
        relative = image.relative_to(owner).with_suffix('').as_posix()
        if re.fullmatch(r'02_differential/L[2-6]_[^/]+/abundance/abundance_(heatmap|bubble)', relative):
            profile = image.stem.removeprefix('abundance_')
        elif re.fullmatch(r'02_differential/L[2-6]_[^/]+/(comparisons/pair_\d+/effect|ancom_W)', relative):
            profile = 'effect'
        if profile:
            kind, stem = 'differential', relative
            options = _json(owner / 'differential_options.json')
            if options.get('method') not in ('ancombc2', 'qiime_ancom'):
                raise ValueError('差异图片缺少原始方法记录')
            analysis = record.get('analysis_run') or options.get('analysis_run') or previous.get('differential_run') or owner
            settings = options.get('plot_settings') or settings
    if profile is None:
        raise ValueError('此图片不是支持单图调整的 R 或 LEfSe 分析图')
    analysis = api.approved_run(str(analysis))
    saved = _json(analysis / 'parameters.json')
    original_config = saved.get('config', saved)
    style = _json(owner / 'group_style.json') or summary.get('group_style') or config
    samples = config.get('samples') or original_config.get('samples')
    if not samples and kind == 'lefse':
        samples = [{'Group': name} for name in _json(analysis / 'lefse_complete.json').get('groups', {})]
    style = resolve_group_style({}, samples or [], fallback=style)
    target = {'profile': profile, 'stem': stem}
    if profile == 'effect' and image.stem == 'effect':
        comparison = _json(image.parent / 'comparison.json')
        if comparison.get('contrast'):
            target['comparison_key'] = comparison['contrast']
    result = {'image': str(image), 'source_image': record.get('source_image') or str(image),
              'kind': kind, 'profile': profile, 'stem': stem, 'plot_target': target,
              'title': LABELS[profile] + ' · ' + (target.get('comparison_key') or image.stem),
              'analysis_run': str(analysis), 'render_run': str(owner),
              'font': summary.get('font') or (options.get('font') if kind == 'differential' else None) or config.get('font') or original_config.get('font', 'sans'),
              'plot_settings': resolve_plot_settings({'plot_settings': settings}),
              'group_style': style, 'visual_options': {},
              **{ext: str(image.with_suffix('.' + ext)) for ext in ('png', 'svg', 'pdf') if image.with_suffix('.' + ext).is_file()}}
    if kind == 'differential':
        result['visual_options'] = plot_options({}, fallback=options)
    if kind == 'advanced':
        original = _json(analysis / 'advanced_options.json')
        visual = plot_options({}, fallback={**original, **summary.get('visual_options', {}),
                                             'diff_plot_top_n': summary.get('top_n', original.get('top_n', 20))})
        result['visual_options'] = {key: visual[key] for key in
            ('diff_plot_top_n', 'heatmap_low', 'heatmap_mid', 'heatmap_high', 'bubble_low', 'bubble_high')}
        result['title'] = '高级生态分析 · ' + image.stem
    if kind == 'lefse':
        result['lda_threshold'] = summary.get('lda_threshold', 2)
    return result


def figure_catalog(run, api):
    run = api.approved_run(run)
    figures, examined = [], 0
    for directory, folders, names in os.walk(run, followlinks=False):
        folders[:] = sorted(name for name in folders if name not in PREVIEW_FOLDERS and name != '09_execution_code')
        for name in sorted(names):
            if not name.endswith('.png'):
                continue
            image = Path(directory) / name
            if not image.resolve().is_relative_to(run):
                continue
            examined += 1
            try:
                figures.append(describe_figure(str(image), api))
            except (ValueError, KeyError, OSError, json.JSONDecodeError):
                pass
            if len(figures) >= 300 or examined >= 3000:
                return {'run': str(run), 'figures': figures, 'limited': True}
    return {'run': str(run), 'figures': figures, 'limited': len(figures) >= 300 or examined >= 3000}


def render_figure(job, request, api):
    figure = describe_figure(request.get('image', ''), api)
    profile = figure['profile']
    requested = request.get('plot_settings', {})
    if not isinstance(requested, dict) or set(requested) - {'general', profile}:
        raise ValueError('单图预览只允许修改通用参数和当前图的参数')
    settings = resolve_plot_settings({'plot_settings': requested}, fallback=figure['plot_settings'])
    analysis = Path(figure['analysis_run'])
    previous = _json(analysis / 'parameters.json')
    saved = previous.get('config', previous)
    samples = saved.get('samples', [])
    if not samples and figure['kind'] == 'lefse':
        samples = [{'Group': g} for g in _json(analysis / 'lefse_complete.json').get('groups', {})]
    style = resolve_group_style(request, samples, fallback=figure['group_style'])
    config = {**style, 'plot_settings': settings, 'plot_target': figure['plot_target'], 'font': figure['font']}
    config.update({key: request[key] for key in ('env','conda','r_env','advanced_env','lefse_env','ancom_env','font')
                   if request.get(key) and str(request[key]).strip()})
    if figure['kind'] == 'diversity':
        api.render_saved_diversity(job, dict(config, diversity_run=str(analysis)))
        output = Path(job['run']) / '04_R_figures' / 'R_original_style' / (figure['stem'] + '.png')
    elif figure['kind'] == 'differential':
        visual = {key: request[key] for key in figure['visual_options'] if key in request}
        config.update(plot_options(visual, fallback=figure['visual_options']))
        api.render_saved_differential(job, dict(config, differential_run=str(analysis)))
        output = Path(job['run']) / (figure['stem'] + '.png')
    elif figure['kind'] == 'advanced':
        from advanced_analysis import render_saved_advanced
        visual = {key: request[key] for key in figure['visual_options'] if key in request}
        render_saved_advanced(job, dict(config, advanced_run=str(analysis), **visual), api)
        output = Path(job['run']) / (figure['stem'] + '.png')
    else:
        from lefse_analysis import render_lefse, _number
        config['lda_threshold'] = _number(request.get('lda_threshold', figure['lda_threshold']), 'LDA 阈值')
        destination = analysis / '04_previews'
        destination.mkdir(exist_ok=True)
        run = destination / ('figure_' + time.strftime('%Y%m%d_%H%M%S') + '_' + uuid.uuid4().hex[:8])
        run.mkdir(); api.marker(run)
        job.update(run=str(run), log=str(run / 'workflow.log'))
        save_plot_settings(run, settings); save_group_style(run, style)
        _save(run / 'parameters.json', {'extension': 'figure-render', 'preprocess': previous['preprocess'],
              'lefse_run': str(analysis), 'config': {**saved, **config}, 'extension_parameters': config,
              'statistics_refit': False})
        render_lefse(job, config, analysis, api)
        output = run / (figure['stem'] + '.png')
    if not output.is_file() or output.stat().st_size == 0:
        raise ValueError('单图预览未生成所选图片；原图保留')
    run = Path(job['run'])
    produced = list(run.rglob('*.png'))
    if produced != [output]:
        raise ValueError('单图预览产生了额外图片；原图保留，请检查日志')
    _save(run / 'figure_request.json', {'source_image': figure['source_image'], 'previous_image': figure['image'],
          'previous_image_sha256': hashlib.sha256(Path(figure['image']).read_bytes()).hexdigest(),
          'analysis_run': str(analysis), 'plot_target': figure['plot_target'], 'statistics_refit': False})
    _save(run / 'complete.json', {'statistics_refit': False, 'analysis_run': str(analysis), 'single_figure': True})
    job['result'] = {**job.get('result', {}), 'figure': describe_figure(str(output), api), 'statistics_refit': False}
    job['step'] = '所选图片预览完成；原图与统计结果保留'
