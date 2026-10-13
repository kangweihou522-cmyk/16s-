"""Archive observed commands and source snapshots without executing replay code."""
import hashlib
import json
import shlex
import uuid
from datetime import datetime, timezone
from pathlib import Path


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def capture_source(job, filename, category, explanation):
    record = job.get('_execution')
    if record is None:
        return
    filename = Path(filename).resolve()
    try:
        raw = filename.read_bytes()
        content = raw.decode('utf-8-sig')
        digest = hashlib.sha256(raw).hexdigest()
    except (OSError, UnicodeError) as error:
        record['warnings'].append(f'未能保存源码 {filename}：{error}')
        return
    key = str(filename) + ':' + digest
    record['sources'][key] = {'path': str(filename), 'sha256': digest,
                              'category': category, 'explanation': explanation,
                              'content': content, 'captured_at': now()}


def prepare_execution(job, root, action, config, runtime):
    job['_execution'] = {'id': job.get('id') or uuid.uuid4().hex,
                         'action': action, 'started_at': now(),
                         'requested_parameters': dict(config), 'runtime': dict(runtime),
                         'sources': {}, 'commands': [], 'warnings': []}
    for name, explanation in [
        ('app.py', '任务编排、输入检查、QIIME 2 命令与统计执行、R 调用。'),
        ('run_layout.py', '结果路径、计数汇总和每一步处理说明。'),
        ('quality_advice.py', '读取逐位置质量统计，生成 DADA2 参数建议。'),
        ('depth_advice.py', '依据过滤后样本计数推荐抽平深度，另附曲线诊断。'),
        ('group_style.py', '统一的分组名称、显示顺序和十六进制配色；不改变统计设计。'),
        ('plot_settings.py', '各图种的尺寸、字号、显示和导出参数；不重算统计结果。'),
        ('taxonomy_palettes.json', '分类柱图内置配色，包括用户 R 配色和常用分类色。'),
        ('figure_editor.py', '读取已出图的实际参数，选择单张图片从缓存重绘、保存预览与原图关系。'),
        ('abundance_tables.py', 'ASV 与 L2–L7 计数、相对丰度和分类注释导出。'),
        ('differential_analysis.py', '门到属的 ANCOM-BC2 / QIIME ANCOM 差异分析编排。'),
        ('lefse_analysis.py', 'LEfSe 分类谱系预处理、检验和基于缓存的 LDA 阈值重绘。'),
        ('functional_prediction.py', 'PICRUSt2、NSTI 与用户 KO 数据库注释和层级汇总。'),
        ('execution_record.py', '记录实际执行命令、参数及源码文本。'),
    ]:
        capture_source(job, Path(root) / name, 'JSON' if name.endswith('.json') else 'Python', explanation)


def capture_r_sources(job, args, cwd=None):
    """Snapshot the invoked script and only the assets render.R actually sources."""
    values = list(map(str, args))
    if 'Rscript' not in values:
        return
    index = values.index('Rscript') + 1
    if values[index:index + 1] == ['--vanilla']:
        index += 1
    if index >= len(values) or values[index].startswith('-'):
        return
    script = Path(values[index])
    if not script.is_absolute():
        script = Path(cwd or Path.cwd()) / script
    capture_source(job, script, 'R', '本次 Rscript 实际调用的完整脚本；参数见命令记录。')
    if script.name != 'render.R':
        return
    capture_source(job, script.with_name('taxonomy_plot.R'), 'R', '按目标层级清洗与合并 Others，选择 Top N 和分类配色。')
    capture_source(job, script.with_name('taxonomy_palettes.json'), 'JSON', 'render.R 本次读取的分类调色板定义。')
    custom = values[index + 4] if len(values) > index + 4 else ''
    assets = Path(custom) if custom else script.parent / 'plots'
    if not assets.is_absolute():
        assets = Path(cwd or Path.cwd()) / assets
    names = ['qiime2_io.R', 'plot_alpha_original_style.R', 'plot_beta_original_style.R']
    if all((assets / name).is_file() for name in names):
        for name in names:
            capture_source(job, assets / name, 'R', 'render.R 本次加载的绘图/结果读取函数（含用户自定义版本）。')


def command_explanation(args):
    values = list(map(str, args))
    if 'qiime' in values:
        tail = values[values.index('qiime') + 1:]
        if '--help' in tail:
            return '检查 QIIME 2 插件是否可用，读取其帮助信息。'
        action = tuple(tail[:2])
        explanations = {
            ('tools', 'import'): '按 manifest 将 FASTQ 导入为 QIIME 2 序列产物。',
            ('tools', 'export'): '导出该 QZA/QZV 的标准数据文件。',
            ('tools', 'peek'): '检查产物类型和格式。',
            ('demux', 'summarize'): '生成质量图、样本序列数和逐位置质量统计。',
            ('diversity', 'alpha-rarefaction'): '在多个深度重复抽样，计算 α 指标曲线；不覆盖未抽平表。',
            ('diversity', 'core-metrics-phylogenetic'): '按本次深度抽平一次，保存抽平表、核心 α/β 指标和 PCoA。',
            ('diversity', 'alpha'): '从本次已保存的抽平表计算额外 α 指标。',
            ('diversity', 'alpha-group-significance'): '执行 α 总体及两两 Kruskal–Wallis 检验，保存校正结果。',
            ('diversity', 'beta-group-significance'): '按所示方法、置换次数和分组执行 β 检验。',
            ('feature-classifier', 'classify-sklearn'): '使用用户提供的训练器进行分类注释。',
            ('phylogeny', 'align-to-tree-mafft-fasttree'): '进行 MAFFT 比对、掩蔽和 FastTree 建树。',
        }
        if action in explanations:
            return explanations[action]
        if tail[:1] == ['dada2']:
            return 'DADA2 质量过滤、降噪、双端拼接（双端模式）和 consensus 去嵌合体。'
        if tail[:1] == ['cutadapt']:
            return '按用户确认的引物及本次参数去除引物。'
        return '执行 QIIME 2 ' + ' '.join(tail[:2]) + '，输入/输出及参数如下。'
    if 'Rscript' in values:
        return '在专用 R 环境运行版本检查、依赖检查或读取已保存结果绘图；完整源码附后。'
    if 'biom' in values:
        return '将 BIOM 特征表转换为 TSV，保留计数供复核。'
    return '执行本次任务调用的外部程序。'


def source_text(sources, category=None):
    parts = []
    for item in sources.values():
        if category and item['category'] != category:
            continue
        parts.extend([f"\n# 文件：{item['path']}", f"# 说明：{item['explanation']}",
                      f"# SHA256：{item['sha256']}；保存时间：{item['captured_at']}",
                      item['content']])
    return '\n'.join(parts) + '\n'


def finalize_execution(job):
    record = job.get('_execution')
    if record is None or not job.get('run'):
        return
    run = Path(job['run'])
    folder = run / '09_execution_code' / (record['action'] + '_' + record['id'][:12])
    folder.mkdir(parents=True, exist_ok=True)
    saved = run / 'parameters.json'
    effective = json.loads(saved.read_text(encoding='utf8')) if saved.is_file() else None
    parameters = {'requested': record['requested_parameters'], 'saved_effective': effective}
    (folder / '本次参数.txt').write_text(json.dumps(parameters, ensure_ascii=False, indent=2), encoding='utf8')
    lines = ['# 本次完整运行代码与参数（UTF-8 文本）',
             f"# 任务：{record['action']}；状态：{job.get('status')}",
             f"# 开始：{record['started_at']}；结束：{now()}", f'# 结果目录：{run}',
             '# 此文件用于阅读和复核；原路径已存在时直接重跑命令可能报错。',
             '# 随机抽平/置换重跑可能不同；以本次保存的表、指标和检验结果为准。',
             '# 命令列出实际显式参数；未显式指定的默认值请参照 QZA/QZV provenance 和对应软件版本。',
             '# 完整程序包含未走到的条件分支；只有下面的命令记录代表实际执行的外部命令。',
             '# 工作站源码版本：' + json.dumps(record['runtime'], ensure_ascii=False),
             '\n# 本次请求参数及保存的最终生效参数',
             json.dumps(parameters, ensure_ascii=False, indent=2)]
    if job.get('error'):
        lines.append('# 任务错误：' + str(job['error']))
    for number, item in enumerate(record['commands'], start=1):
        lines.extend([f"\n# 命令 {number}；阶段：{item['stage']}",
                      '# 说明：' + command_explanation(item['args']),
                      '# 阶段参数：' + json.dumps(item['stage_parameters'], ensure_ascii=False),
                      f"# 状态：{item['status']}；退出码：{item['exit_code']}",
                      f"# 开始：{item['started_at']}；结束：{item['finished_at']}"])
        if item.get('launch_error'):
            lines.append('# 程序启动错误：' + item['launch_error'])
        lines.extend(['cd -- ' + shlex.quote(item['cwd']), shlex.join(item['args'])])
    if not record['commands']:
        lines.append('\n# 本任务未通过命令执行器运行外部命令；Python 处理及检查结果见阶段报告。')
    for logical in ['manifest.tsv', 'metadata.tsv', 'layout.json', 'differential_options.json',
                    'lefse_parameters.json', 'group_style.json', 'plot_settings.json', 'plot_target.json', 'figure_request.json', '00_inputs/effective_parameters.json',
                    '03_figures/render_summary.json', '10_abundance/manifest.json']:
        # Use the saved layout, preserving both organized and older flat runs.
        from run_layout import artifact_path
        filename = artifact_path(run, logical)
        if filename.is_file():
            lines.extend(['\n# 本次文件：' + str(filename), filename.read_text(encoding='utf-8-sig')])
    for category in ['R', 'Python', 'JSON']:
        text = source_text(record['sources'], category)
        if category == 'R' and not text.strip():
            text = '# 本任务没有实际调用 R 脚本。\n'
        (folder / (category + '完整代码.txt')).write_text(text, encoding='utf8')
        lines.extend(['\n# ' + category + ' 完整源码快照及说明', text])
    if effective and 'preprocess' in effective:
        parent = Path(effective['preprocess'])
        prior = sorted((parent / '09_execution_code').glob('*/完整运行代码.txt'))
        if prior:
            historical = '\n'.join('# 关联预处理记录：' + str(name) + '\n' + name.read_text(encoding='utf8') for name in prior)
            (folder / '关联预处理代码.txt').write_text(historical, encoding='utf8')
            lines.extend(['\n# 本次所用预处理的历史运行记录（各次状态以原记录为准）', historical])
        else:
            lines.append('\n# 关联预处理由旧版生成，未保存当时的完整源码；不能用当前源码冒充旧版实际代码。')
            for name in ['commands.jsonl', 'parameters.json']:
                filename = parent / name
                if filename.is_file():
                    lines.extend(['# 旧版原始记录：' + str(filename), filename.read_text(encoding='utf8')])
        dependencies=effective.get('extension_parameters',{})
        for key,label in [('abundance_run','复用的丰度导出'),('lefse_run','用于重绘的 LEfSe 原始检验'),('diversity_run','用于重绘的多样性原始指标及检验'),('differential_run','用于重绘的差异分析原始模型及检验')]:
            if key=='abundance_run' and record['action']!='differential': continue
            if key=='lefse_run' and record['action']!='lefse-render': continue
            if key=='diversity_run' and record['action']!='diversity-render': continue
            if key=='differential_run' and record['action']!='differential-render': continue
            value=dependencies.get(key)
            if not value: continue
            dependency=Path(value)
            if re_match := __import__('re').match(r'^([A-Za-z]):[\\/](.*)$',str(value)):
                dependency=Path('/mnt/'+re_match[1].lower())/re_match[2].replace('\\','/')
            dependency=dependency.expanduser().resolve()
            if dependency in (parent,run) or not (dependency/'workstation_run.json').is_file(): continue
            for filename in sorted((dependency/'09_execution_code').glob('*/完整运行代码.txt')):
                lines.extend(['\n# 关联历史任务：'+label+'；'+str(filename),filename.read_text(encoding='utf8')])
    lines.extend('# 源码记录提示：' + warning for warning in record['warnings'])
    (folder / '完整运行代码.txt').write_text('\n'.join(lines) + '\n', encoding='utf8')
    job['execution_code_directory'] = str(folder)
