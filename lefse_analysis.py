"""LEfSe on saved, unrarefied counts; plot previews reuse the original fit.

CLI contracts follow https://github.com/SegataLab/lefse (format/run/plot modules).
Only newly generated formatter pickles are opened, inside the LEfSe environment.
"""
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path

from execution_record import capture_source
from group_style import resolve_group_style
from plot_settings import resolve_plot_settings
from run_layout import artifact_path

RANKS = ('k', 'p', 'c', 'o', 'f', 'g', 's')
CLI_NAMES = {
    'format': ('lefse_format_input.py', 'format_input.py', 'lefse_format_input', 'format_input'),
    'run': ('lefse_run.py', 'run_lefse.py', 'lefse_run', 'run_lefse'),
    'bar': ('lefse_plot_res.py', 'plot_res.py', 'lefse_plot_res', 'plot_res'),
    'cladogram': ('lefse_plot_cladogram.py', 'plot_cladogram.py', 'lefse_plot_cladogram', 'plot_cladogram'),
}
UNKNOWN = re.compile(r'(?:^|_)(?:incertae_?sedis|norank|no_rank|uncultured|unclassified|unidentified|unknown|unassigned|ambiguous|not_assigned)(?:_|$)')
FIT_RESULT = '02_statistics/lefse_results.res'
PLOT_TARGETS = {'lefse_bar': '03_figures/bar', 'lefse_cladogram': '03_figures/cladogram'}


def _json(filename, value):
    Path(filename).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf8')


def _tsv(filename, header, rows):
    with Path(filename).open('w', encoding='utf8', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(header)
        writer.writerows(rows)


def _number(value, label, minimum=0, maximum=None):
    if isinstance(value, bool):
        raise ValueError(label + ' 必须为有效数值')
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(label + ' 必须为有效数值') from None
    if not math.isfinite(number) or number < minimum or (maximum is not None and number > maximum):
        raise ValueError(label + ' 超出允许范围')
    return number


def parameters(config):
    env = str(config.get('lefse_env') or 'lefse').strip()
    if not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*', env):
        raise ValueError('LEfSe 环境名只能包含字母、数字、点、下划线和短横线')
    strategy = config.get('multiclass_strategy', 0)
    strategies = {'0': 0, '1': 1, 'one-against-all': 0, 'all-against-all': 1}
    if isinstance(strategy, bool) or str(strategy) not in strategies:
        raise ValueError('LEfSe 多组策略必须是 one-against-all (0) 或 all-against-all (1)')
    kw = _number(config.get('kw_alpha', 0.05), 'KW alpha', 0, 1)
    wilcoxon = _number(config.get('wilcoxon_alpha', 0.05), 'Wilcoxon alpha', 0, 1)
    if kw == 0 or wilcoxon == 0:
        raise ValueError('KW / Wilcoxon alpha 必须大于 0')
    return {'lefse_env': env, 'lda_threshold': _number(config.get('lda_threshold', 2), 'LDA 阈值'),
            'kw_alpha': kw, 'wilcoxon_alpha': wilcoxon,
            'multiclass_strategy': strategies[str(strategy)], 'fit_lda_threshold': 0,
            'wilcoxon_min_samples': 10, 'subclass': None,
            'subject': 'unique SampleID (formatter metadata only)', 'paired_design': False,
            'normalization': 'CPM from each original filtered ASV sample total; formatter -o -1'}


def _samples(samples):
    rows, seen = [], set()
    for item in samples:
        sid, group = item.get('SampleID'), item.get('Group')
        if not isinstance(sid, str) or not sid or sid != sid.strip() or sid in seen:
            raise ValueError('保存的 SampleID 缺失、重复或无效')
        if not isinstance(group, str) or not group or group != group.strip():
            raise ValueError('保存的 Group 缺失或无效')
        if any(char in sid + group for char in '\t\r\n\x00'):
            raise ValueError('SampleID / Group 包含无效控制字符')
        seen.add(sid)
        rows.append({'SampleID': sid, 'Group': group})
    groups = Counter(item['Group'] for item in rows)
    if len(groups) < 2 or any(count < 2 for count in groups.values()):
        raise ValueError('LEfSe 至少需要两组，每组至少 2 个原始样本')
    return rows


def read_counts(filename, samples):
    """Match by exact sample IDs, never by column order or numeric conversion."""
    ids = [sample['SampleID'] for sample in samples]
    with Path(filename).open(encoding='utf-8-sig', newline='') as stream:
        rows = list(csv.reader(stream, delimiter='\t'))
    heads = [i for i, row in enumerate(rows) if row and row[0] in ('#OTU ID', '#Feature ID', 'Feature ID', 'feature-id')]
    if len(heads) != 1:
        raise ValueError('LEfSe 特征表缺少或重复 BIOM TSV 表头')
    header = rows[heads[0]][1:]
    if len(header) != len(set(header)) or set(header) != set(ids):
        raise ValueError('LEfSe 特征表样本与保存的 SampleID 分组不一致')
    indices = [header.index(sid) + 1 for sid in ids]
    counts = {}
    for row in rows[heads[0] + 1:]:
        if not row:
            continue
        if len(row) != len(header) + 1 or not row[0] or row[0] in counts:
            raise ValueError('LEfSe 特征表 ASV 标识重复或列数异常')
        values = []
        for index in indices:
            try:
                value = Decimal(row[index])
            except InvalidOperation:
                raise ValueError('LEfSe 特征表包含无效序列计数') from None
            if not value.is_finite() or value < 0 or value != value.to_integral_value():
                raise ValueError('LEfSe 特征表计数必须为非负整数')
            values.append(int(value))
        counts[row[0]] = values
    totals = [sum(values[i] for values in counts.values()) for i in range(len(ids))]
    if any(value == 0 for value in totals):
        raise ValueError('LEfSe 原始样本总计数为 0；不能自动丢弃样本')
    return counts, totals


def parse_taxon(taxon):
    result = [''] * 7
    previous = -1
    for position, raw in enumerate(taxon.split(';')):
        raw = raw.strip()
        match = re.match(r'^(?:([dkpcofgs])__|D_([0-6])__)(.*)$', raw, re.I)
        if match:
            rank = int(match[2]) if match[2] is not None else RANKS.index('k' if match[1].lower() == 'd' else match[1].lower())
            label = match[3].strip()
        else:
            rank, label = position, raw
        if rank > 6:
            raise ValueError('分类注释包含超过 L7 的层级')
        if rank <= previous:
            raise ValueError('分类注释层级重复或顺序异常')
        previous = rank
        result[rank] = label
    return tuple(result)


def read_taxonomy(filename, features):
    with Path(filename).open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        names = reader.fieldnames or []
        taxon = next((name for name in names if name.lower() in ('taxon', 'taxonomy')), None)
        if not taxon or not names:
            raise ValueError('分类 TSV 缺少 Taxon 列')
        values = {}
        for row in reader:
            sid = row.get(names[0])
            if not sid or sid in values or row.get(taxon) is None or None in row:
                raise ValueError('分类 TSV 的 ASV 标识重复或列数异常')
            values[sid] = row[taxon]
    missing = set(features) - values.keys()
    if missing:
        raise ValueError('分类 TSV 缺少特征表 ASV：' + ', '.join(sorted(missing)[:5]))
    return {sid: parse_taxon(values[sid]) for sid in features}


def unknown(label):
    token = re.sub(r'[^a-z0-9]+', '_', label.casefold()).strip('_')
    return not token or bool(UNKNOWN.search(token)) or token in {
        'other', 'others', 'na', 'none', 'nan', 'bacterium', 'archaeon', 'microorganism',
        'environmental_sample', 'sp', 'spp',
    }


def safe_names(labels, prefix=''):
    """Whitelist characters and disambiguate all collisions independent of order."""
    bases = defaultdict(list)
    for label in sorted(set(labels)):
        base = prefix + re.sub(r'[^A-Za-z0-9_]', '_', label)
        if not base or not base[0].isalpha():
            base = 'f_' + base
        bases[base].append(label)
    mapping, collisions = {}, set()
    for base, originals in sorted(bases.items()):
        for original in originals:
            candidate = base
            if len(originals) > 1:
                collisions.add(original)
                candidate += '_h' + hashlib.sha256(original.encode('utf8')).hexdigest()[:16]
            mapping[original] = candidate
    if len(set(mapping.values())) != len(mapping):
        raise ValueError('名称清理仍存在冲突；请检查分类名称')
    return mapping, collisions


def build_abundances(counts, taxonomy, samples, totals):
    labels = {RANKS[i] + '__' + label for lineage in taxonomy.values() for i, label in enumerate(lineage) if not unknown(label)}
    names, collisions = safe_names(labels)
    group_names, group_collisions = safe_names([sample['Group'] for sample in samples])
    unknown_branches = {}
    for lineage in taxonomy.values():
        for i, label in enumerate(lineage):
            if unknown(label):
                context = tuple(RANKS[j] + '__' + value for j, value in enumerate(lineage[:i + 1]))
                digest = hashlib.sha256(json.dumps(context, ensure_ascii=False).encode('utf8')).hexdigest()[:16]
                unknown_branches[context] = RANKS[i] + '__Others_at_L' + str(i + 1) + '_h' + digest
    if len(set(unknown_branches.values())) != len(unknown_branches) or set(unknown_branches.values()) & set(names.values()):
        raise ValueError('未知分类上下文编码冲突')
    aggregates, others = {}, {}
    for sid in sorted(counts):
        lineage, values = taxonomy[sid], counts[sid]
        for level in range(2, 8):
            original = tuple(RANKS[i] + '__' + label for i, label in enumerate(lineage[:level]))
            if unknown(lineage[level - 1]):
                # Keep each original parent context separate. Never analyze a pooled Others.
                context = (level, original[:-1])
                entry = others.setdefault(context, {'asv_count': 0, 'counts': [0] * len(totals)})
                entry['asv_count'] += 1
                entry['counts'] = [a + b for a, b in zip(entry['counts'], values)]
                continue
            entry = aggregates.setdefault(original, [0] * len(totals))
            aggregates[original] = [a + b for a, b in zip(entry, values)]
    features = []
    for original, values in sorted(aggregates.items()):
        # Unknown intermediate nodes keep their exact parent context, but have no tested row.
        key = '|'.join(names.get(label, unknown_branches.get(original[:i + 1], '')) for i, label in enumerate(original))
        if any(not branch for branch in key.split('|')):
            raise ValueError('分类分支缺少名称映射')
        features.append({'input_feature': key, 'result_feature': key.replace('|', '.'),
                         'rank': 'L' + str(len(original)), 'taxon': original[-1],
                         'lineage': '; '.join(original), 'counts': values,
                         'cpm': [count / total * 1e6 for count, total in zip(values, totals)]})
    exclusions = []
    for (level, parent), value in sorted(others.items()):
        exclusions.append({'rank': 'L' + str(level), 'parent_lineage': '; '.join(parent),
                           'display_taxon': 'Others', **value})
    return {'features': features, 'others': exclusions, 'names': names, 'name_collisions': collisions,
            'groups': group_names, 'group_collisions': group_collisions,
            'unknown_branches': unknown_branches}


def _prefix(config, api):
    return [api.find_conda(config), 'run', '--no-capture-output', '-n', config['lefse_env']]


def _runtime(job, config, api, operations):
    # Check wrappers on the PATH of the selected env, never execute user shell text.
    code = """import importlib.metadata,json,os,shutil,sys
names=json.loads(sys.argv[1]);paths={}
for operation,choices in names.items():
    paths[operation]=None
    for choice in choices:
        path=shutil.which(choice)
        if path and os.path.commonpath([os.path.abspath(path),os.path.abspath(sys.prefix)])==os.path.abspath(sys.prefix):
            paths[operation]=path;break
try: version=importlib.metadata.version('lefse')
except importlib.metadata.PackageNotFoundError: version='unknown'
print('LEFSE_RUNTIME='+json.dumps({'prefix':sys.prefix,'python':sys.executable,'version':version,'paths':paths}))
"""
    proc = api.probe(job, _prefix(config, api) + ['python', '-c', code, json.dumps({op: CLI_NAMES[op] for op in operations})], 90)
    install_hint = ('请在 Ubuntu 的工作站源码目录执行：WORKSTATION_LEFSE_ENV='
                    + config['lefse_env'] + ' bash install_lefse.sh')
    lines = [line[len('LEFSE_RUNTIME='):] for line in proc.stdout.splitlines() if line.startswith('LEFSE_RUNTIME=')]
    if proc.returncode or len(lines) != 1:
        raise ValueError('无法检查已安装的 LEfSe 环境：' + (proc.stderr or proc.stdout)[-1500:]
                         + '\n' + install_hint)
    runtime = json.loads(lines[0])
    missing = [op for op in operations if not runtime['paths'].get(op)]
    if missing:
        raise ValueError('LEfSe 环境缺少官方命令：' + ', '.join(CLI_NAMES[op][0] for op in missing)
                         + '\n检测环境：' + str(runtime.get('prefix', 'unknown'))
                         + '；LEfSe版本：' + str(runtime.get('version', 'unknown'))
                         + '\n' + install_hint)
    return runtime


PRUNE_SCRIPT = """# Only open the formatter pickle generated in this same isolated run.
import json,pickle,sys
from pathlib import Path
filename,mapping,report=map(Path,sys.argv[1:])
allowed=set(json.loads(mapping.read_text(encoding='utf8')))
with filename.open('rb') as stream: data=pickle.load(stream)
if not allowed.issubset(data['feats']): raise ValueError('Formatter lost or renamed requested features')
extra=sorted(set(data['feats'])-allowed)
constant=sorted(key for key in allowed if len(set(data['feats'][key]))<=1)
for key in extra+constant: del data['feats'][key]
with filename.open('wb') as stream: pickle.dump(data,stream)
report.write_text(json.dumps({'removed_auto_ancestors':extra,'constant_features':constant,'tested_features':len(data['feats'])},indent=2),encoding='utf8')
"""

EMPTY_PLOT_SCRIPT = """# An explicit no-hit image, without fabricated biomarkers or tree branches.
import json,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
folder=Path(sys.argv[1]);threshold=float(sys.argv[2]);status=sys.argv[3]
settings=json.loads(Path(sys.argv[4]).read_text(encoding='utf8')) if len(sys.argv)>4 else {}
target=sys.argv[5] if len(sys.argv)>5 else ''
if target not in ('','bar','cladogram'): raise ValueError('Unknown LEfSe plot target')
general=settings.get('general',{}); formats=('png','svg','pdf') if general.get('formats')=='png_svg_pdf' else ('png','svg')
if general.get('base_size') is not None: matplotlib.rcParams['font.size']=general['base_size']
if general.get('font_family','original')!='original':
    from matplotlib.font_manager import findfont,FontProperties
    findfont(FontProperties(family=[general['font_family']]),fallback_to_default=False)
    matplotlib.rcParams['font.family']=[general['font_family']]
matplotlib.rcParams['svg.fonttype']='none'
for kind,title in [('bar','LEfSe LDA scores'),('cladogram','LEfSe taxonomy cladogram')]:
    if target and kind!=target: continue
    profile=settings.get('lefse_'+kind,{})
    figure,axis=plt.subplots(figsize=(profile.get('width') or 8,profile.get('height') or 4.5));axis.axis('off')
    title=(profile.get('title') or title).replace('{threshold}',str(threshold))
    axis.set_title(title,fontsize=general.get('title_size') or general.get('base_size'))
    axis.text(.5,.55,'No biomarkers above LDA threshold '+str(threshold),ha='center',va='center',transform=axis.transAxes)
    note='No variable known taxa; statistical analysis skipped' if status.startswith('skipped') else 'Saved analysis; no statistical refit'
    if general.get('show_caption',True): axis.text(.5,.4,note,ha='center',va='center',transform=axis.transAxes,color='#555555')
    for extension in formats:
        figure.savefig(folder/(kind+'.'+extension),format=extension,dpi=general.get('dpi') or 160,bbox_inches='tight',
                       transparent=general.get('background')=='transparent')
    plt.close(figure)
"""

OFFICIAL_PLOT_SCRIPT = """# Presentation adapter for the installed official LEfSe plotter, not a model fit.
# Apply display controls, named group colors/order, and Python 3 vertical-plot compatibility.
# Saved LDA scores and the official taxonomy/plotting functions remain in use.
import ast,hashlib,importlib.util,json,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

kind,style_file=sys.argv[1:3]
style=json.loads(Path(style_file).read_text(encoding='utf8'))
settings=style.get('plot_settings',{}); general=settings.get('general',{}); profile=settings.get('lefse_'+kind,{})
orientation=profile.get('orientation','horizontal') if kind=='bar' else None
original_order=style['group_order']; groups=style['groups']
safe_order=[groups[name] for name in original_order]
positions={name:i for i,name in enumerate(safe_order)}
original_names={safe:original for original,safe in groups.items()}
palette=[style['group_colors'][name] for name in original_order]
name='lefse.lefse_plot_res' if kind=='bar' else 'lefse.lefse_plot_cladogram'
spec=importlib.util.find_spec(name)
if spec is None or not spec.origin: raise ValueError('Official LEfSe plot module is unavailable: '+name)
source_path=Path(spec.origin); source=source_path.read_bytes()

def group_order(values):
    values=list(values)
    if any(value not in positions for value in values):
        raise ValueError('Official LEfSe plot contains an unknown group')
    return sorted(values,key=positions.__getitem__)

# Change explicit class sorts only; leave taxonomic node/feature sorting intact.
targets={ast.dump(ast.parse(expr,mode='eval').body) for expr in
         ("params['all_feats'].split(':')", "data['cls']", "tree['classes']")}
class GroupSortAdapter(ast.NodeTransformer):
    count=0
    vertical=False
    def visit_FunctionDef(self,node):
        previous=self.vertical; self.vertical=node.name=='plot_histo_ver'
        self.generic_visit(node); self.vertical=previous
        return node
    def visit_Subscript(self,node):
        self.generic_visit(node)
        # Official 1.1.2's vertical plot uses zip(...)[1], invalid in Python 3.
        if self.vertical and isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Name) and node.value.func.id=='zip':
            node.value=ast.copy_location(ast.Call(func=ast.Name(id='list',ctx=ast.Load()),args=[node.value],keywords=[]),node.value)
        return node
    def visit_Assign(self,node):
        self.generic_visit(node)
        if self.vertical and any(isinstance(target,ast.Name) and target.id=='nam' for target in node.targets):
            node.value=ast.copy_location(ast.parse("['.'.join(d[0].split('.')[-params['n_scl']:]) for d in data['rows']]",mode='eval').body,node.value)
        return node
    def visit_Call(self,node):
        self.generic_visit(node)
        if (isinstance(node.func,ast.Name) and node.func.id=='sorted' and
            len(node.args)==1 and not node.keywords and ast.dump(node.args[0]) in targets):
            node.func=ast.copy_location(ast.Name(id='_workstation_group_order',ctx=ast.Load()),node.func)
            self.count+=1
        if self.vertical and isinstance(node.func,ast.Attribute) and node.func.attr=='sort' and any(k.arg=='key' for k in node.keywords):
            # Keep group blocks ordered even when an LDA score exceeds log abundance.
            for item in node.keywords:
                if item.arg=='key': item.value=ast.copy_location(ast.parse("lambda ab: (_workstation_positions[ab[2]], abs(float(ab[3])))",mode='eval').body,item.value)
        return node
adapter=GroupSortAdapter(); tree=adapter.visit(ast.parse(source,filename=str(source_path)))
if adapter.count<2: raise ValueError('This official LEfSe plot version cannot apply the requested group order safely')
namespace={'__name__':'workstation_official_lefse_plot','__file__':str(source_path)}
exec(compile(ast.fix_missing_locations(tree),str(source_path),'exec'),namespace)
namespace['_workstation_group_order']=group_order
namespace['_workstation_positions']=positions
namespace['colors']=palette
if kind=='cladogram': namespace['dark_colors']=palette
if kind=='bar' and orientation=='vertical':
    original_read=namespace['read_data']
    def ordered_data(*args,**kwargs):
        data=original_read(*args,**kwargs); data['cls']=list(safe_order)
        return data
    namespace['read_data']=ordered_data
runtime_parameters={}; effective_canvas={}
original_params=namespace['read_params']
def visual_params(*args,**kwargs):
    params=original_params(*args,**kwargs)
    if kind=='bar':
        size=profile.get('feature_font_size') or general.get('base_size')
        if size is not None: params['feature_font_size']=size
        # Official vertical function interprets rs as a coordinate, while its
        # CLI (and horizontal function) use the right margin proportion.
        if params.get('orientation')=='v': params['rs']=1-params['rs']
    runtime_parameters.update(params)
    return params
namespace['read_params']=visual_params
if general.get('base_size') is not None: matplotlib.rcParams['font.size']=general['base_size']
if general.get('font_family','original')!='original':
    from matplotlib.font_manager import findfont,FontProperties
    findfont(FontProperties(family=[general['font_family']]),fallback_to_default=False)
    matplotlib.rcParams['font.family']=[general['font_family']]

original_legend=Axes.legend
def ordered_legend(axis,*args,**kwargs):
    if len(args)>=2: handles,labels=args[:2]; extra=args[2:]
    elif not args: handles,labels=axis.get_legend_handles_labels(); extra=()
    else: return original_legend(axis,*args,**kwargs)
    grouped=bool(labels) and all(label in original_names for label in labels)
    if grouped:
        pairs=sorted(zip(handles,labels),key=lambda pair:positions[pair[1]])
        handles=[pair[0] for pair in pairs]
        labels=[original_names[pair[1]] for pair in pairs]
    legend_size=(profile.get('class_font_size') if grouped else None) or general.get('legend_text_size') or general.get('base_size')
    if legend_size is not None: kwargs['prop']={'size':legend_size}
    position=general.get('legend_position','original')
    if position=='right':
        kwargs.update(loc='upper left',bbox_to_anchor=(1.02,1 if grouped else max(.2,.9-.055*len(safe_order))))
        kwargs['ncol']=1
    elif position=='bottom' and grouped:
        kwargs.update(loc='upper center',bbox_to_anchor=(.5,-.16),ncol=max(1,len(labels)))
    legend=original_legend(axis,handles,labels,*extra,**kwargs)
    if position=='none': legend.set_visible(False)
    return legend
Axes.legend=ordered_legend

# Horizontal bar groups follow the chosen order from top to bottom.
original_savefig=plt.savefig
def ordered_savefig(*args,**kwargs):
    figure=plt.gcf()
    width,height=figure.get_size_inches()
    figure.set_size_inches(profile.get('width') or width,profile.get('height') or height,forward=True)
    actual_width,actual_height=figure.get_size_inches()
    effective_canvas.update(width=float(actual_width),height=float(actual_height),dpi=kwargs.get('dpi'),
                            font_family=list(matplotlib.rcParams['font.family']),base_size=matplotlib.rcParams['font.size'])
    for axis in figure.axes:
        if kind=='bar' and orientation=='horizontal': axis.invert_yaxis()
        if kind=='bar' and orientation=='vertical' and axis.get_legend() is None: axis.legend()
        if profile.get('x_label'): axis.set_xlabel(profile['x_label'])
        if profile.get('y_label'): axis.set_ylabel(profile['y_label'])
        tick_size=general.get('axis_text_size') or general.get('base_size')
        if tick_size is not None: axis.tick_params(labelsize=tick_size)
        if kind=='bar' and orientation=='vertical' and profile.get('feature_font_size') is not None:
            axis.tick_params(axis='x',labelsize=profile['feature_font_size'])
        theme=general.get('theme','original')
        if theme!='original':
            axis.set_frame_on(theme!='minimal')
            for name,spine in axis.spines.items(): spine.set_visible(theme=='bw' or (theme=='classic' and name in ('left','bottom')))
            if theme=='classic': axis.grid(False)
            else: axis.grid(True,color='#DDDDDD',linewidth=.5)
        if kind=='bar' and profile.get('top_n',0) and general.get('show_caption',True):
            figure.text(.5,.01,'Display up to top '+str(profile['top_n'])+' biomarkers by |LDA|; saved analysis unchanged',ha='center',fontsize=general.get('axis_text_size') or 8)
    if general.get('background')=='transparent': kwargs.update(transparent=True,facecolor='none',edgecolor='none')
    if general.get('legend_position','original')!='original': kwargs['bbox_inches']='tight'
    return original_savefig(*args,**kwargs)
plt.savefig=ordered_savefig
matplotlib.rcParams['svg.fonttype']='none'
sys.argv=[name,*sys.argv[3:]]
display_empty=False
if kind=='cladogram' and profile.get('max_level') is not None:
    lines=Path(sys.argv[1]).read_text(encoding='utf8').splitlines()
    display_empty=not any(line.strip() and line.split()[0].count('.')<profile['max_level'] for line in lines)
if display_empty:
    params=namespace['read_params'](sys.argv)
    figure,axis=plt.subplots(); axis.axis('off')
    axis.set_title(params['title'],fontsize=general.get('title_size') or general.get('base_size'))
    axis.text(.5,.5,'No taxa available at the selected display level',ha='center',va='center',transform=axis.transAxes)
    plt.savefig(params['output_file'],format=params['format'],dpi=params['dpi']); plt.close(figure)
else: namespace['plot_res' if kind=='bar' else 'plot_cladogram']()
Path(sys.argv[2]+'.style.json').write_text(json.dumps({
    'official_module':name,'official_source':str(source_path),
    'official_source_sha256':hashlib.sha256(source).hexdigest(),
    'presentation_group_sort_adaptations':adapter.count,
    'group_order':original_order,'group_colors':style['group_colors'],
    'plot_settings':settings,'orientation':orientation,
    'official_parameters':runtime_parameters,'effective_canvas':effective_canvas,
    'empty_display_level':display_empty,
    'statistics_refit':False},ensure_ascii=False,indent=2)+'\\n',encoding='utf8')
"""


def _script(filename, source, job):
    Path(filename).write_text(source, encoding='utf8')
    capture_source(job, filename, 'Python', '本次实际执行的 LEfSe 输入整理、官方绘图配色与顺序适配或无命中提示脚本。')
    return str(filename)


def read_results(filename):
    result, seen = [], set()
    with Path(filename).open(encoding='utf8', newline='') as stream:
        for row in csv.reader(stream, delimiter='\t'):
            if not row:
                continue
            if len(row) != 5 or not row[0] or row[0] in seen:
                raise ValueError('LEfSe 结果格式无效或特征重复')
            seen.add(row[0])
            abundance = _number(row[1], 'LEfSe abundance', minimum=-math.inf)
            score = _number(row[3], 'LEfSe LDA score', minimum=-math.inf) if row[3] else None
            pvalue = _number(row[4], 'LEfSe KW P', 0, 1) if row[4] not in ('', '-') else None
            if bool(row[2]) != (score is not None):
                raise ValueError('LEfSe 结果的分组与 LDA 分数不一致')
            result.append({'feature': row[0], 'log_abundance': abundance, 'class': row[2],
                           'lda_score': score, 'kw_p': pvalue, 'raw': row})
    return result


def filter_results(results, threshold):
    """Use the official strict abs(LDA)>threshold rule; retain tree context rows."""
    threshold = _number(threshold, 'LDA 阈值')
    rows, selected = [], []
    for item in results:
        row = list(item['raw'])
        if item['lda_score'] is not None and abs(item['lda_score']) > threshold:
            selected.append(item)
        else:
            row[2:4] = ['', '']
        rows.append(row)
    return rows, selected


def _plot_options(kind, settings, threshold, class_order):
    """Official supported CLI parameters; figure/layout overrides live in the adapter."""
    general, profile = settings['general'], settings['lefse_' + kind]
    title = profile['title'] or ('LEfSe LDA scores' if kind == 'bar' else 'LEfSe taxonomy cladogram') + ' (LDA > ' + str(threshold) + ')'
    options = ['--dpi', str(general['dpi'] or 160), '--title', title.replace('{threshold}', str(threshold)),
               '--all_feats', class_order]
    title_size = general['title_size'] or general['base_size']
    if title_size is not None:
        options += ['--title_font_size', str(title_size)]
    class_size = profile.get('class_font_size') or general['legend_text_size'] or general['base_size']
    if class_size is not None:
        options += ['--class_legend_font_size', str(class_size)]
    if kind == 'bar':
        options += ['--subclades', str(profile['label_levels']), '--max_feature_len', '10000',
                    '--orientation', 'v' if profile['orientation'] == 'vertical' else 'h']
        fields = {'width': '--width', 'height': '--height'}
        feature_size = profile['feature_font_size'] or general['base_size']
        if feature_size is not None:
            # The official CLI declares this flag as int; adapter retains any fractional size.
            options += ['--feature_font_size', str(round(feature_size))]
    else:
        fields = {'max_level': '--max_lev', 'min_point_size': '--min_point_size', 'max_point_size': '--max_point_size'}
        feature_size = profile['label_font_size'] or general['base_size']
        if feature_size is not None:
            options += ['--label_font_size', str(feature_size)]
    for field, flag in fields.items():
        if profile[field] is not None:
            options += [flag, str(profile[field])]
    return options


def _save_input(run, counts, totals, samples, data):
    folder = run / '01_inputs'
    folder.mkdir(exist_ok=True)
    ids = [sample['SampleID'] for sample in samples]
    _tsv(folder / 'sample_mapping.tsv', ['InputColumn', 'SampleID', 'Group', 'LEfSeClass', 'OriginalTotalReads'],
         [(i + 1, row['SampleID'], row['Group'], data['groups'][row['Group']], totals[i]) for i, row in enumerate(samples)])
    _tsv(folder / 'feature_mapping.tsv', ['InputFeature', 'ResultFeature', 'Rank', 'Taxon', 'OriginalLineage'],
         [(row['input_feature'], row['result_feature'], row['rank'], row['taxon'], row['lineage']) for row in data['features']])
    _tsv(folder / 'name_mapping.tsv', ['Kind', 'OriginalName', 'SafeName', 'CollisionDisambiguated'],
         [('taxon', original, safe, original in data['name_collisions']) for original, safe in sorted(data['names'].items())]
         + [('group', original, safe, original in data['group_collisions']) for original, safe in sorted(data['groups'].items())]
         + [('unknown_parent_context', '; '.join(original), safe, False) for original, safe in sorted(data['unknown_branches'].items())])
    _tsv(folder / 'others_counts.tsv', ['Rank', 'OriginalParentLineage', 'DisplayTaxon', 'ExcludedASVCount', *ids],
         [(row['rank'], row['parent_lineage'], 'Others', row['asv_count'], *row['counts']) for row in data['others']])
    _tsv(folder / 'taxon_counts.tsv', ['Feature', *ids], [(row['input_feature'], *row['counts']) for row in data['features']])
    # LEfSe 1.1.2 requires subject metadata even for independent samples.
    # Each SampleID is unique; this labels samples without introducing pairs.
    with (folder / 'lefse_input_cpm.tsv').open('w', encoding='utf8', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['Class', *(data['groups'][sample['Group']] for sample in samples)])
        writer.writerow(['SampleID', *ids])
        writer.writerows([row['input_feature'], *(format(value, '.17g') for value in row['cpm'])] for row in data['features'])
    _json(folder / 'allowed_features.json', [row['result_feature'] for row in data['features']])
    ranks = {}
    for level in range(2, 8):
        rank = 'L' + str(level)
        excluded = [row for row in data['others'] if row['rank'] == rank]
        ranks[rank] = {'candidate_taxa': sum(row['rank'] == rank for row in data['features']),
                       'excluded_asv_count': sum(row['asv_count'] for row in excluded),
                       'excluded_reads': dict(zip(ids, [sum(row['counts'][i] for row in excluded) for i in range(len(ids))]))}
    summary = {'asv_count': len(counts), 'sample_count': len(samples), 'group_counts': dict(Counter(row['Group'] for row in samples)),
               'original_total_reads': dict(zip(ids, totals)), 'ranks': ranks,
               'normalization': 'Each taxon count / original sample ASV total * 1e6, exactly once',
               'unknown_policy': 'Only unknown terminal taxa are excluded at a rank. Known ancestors and descendants remain eligible. Unknown intermediate branches use reversible, context-specific Others codes and have no tested row.'}
    _json(folder / 'input_summary.json', summary)
    return summary


def run_lefse(job, config, saved, parent, api):
    """The caller creates/marks job.run and validates the completed preprocess parent."""
    run, parent = Path(job['run']), Path(parent)
    settings = parameters(config)
    runtime_config = {**saved, **settings}
    if config.get('conda'):
        runtime_config['conda'] = config['conda']
    samples = _samples(saved.get('samples', []))
    style = resolve_group_style(config, samples, fallback=saved)
    _json(run / 'group_style.json', style)
    _json(run / 'plot_settings.json', resolve_plot_settings(config, fallback=saved))
    api.begin_stage(job, '01_inputs', 'LEfSe：未抽平计数、分组与分类输入', settings)
    counts, totals = read_counts(artifact_path(parent, 'exported-table/feature-table.tsv'), samples)
    taxonomy = read_taxonomy(artifact_path(parent, 'exported-taxonomy/taxonomy.tsv'), counts)
    data = build_abundances(counts, taxonomy, samples, totals)
    summary = _save_input(run, counts, totals, samples, data)
    runtime = _runtime(job, runtime_config, api, ('format', 'run', 'bar', 'cladogram'))
    _json(run / '01_inputs/runtime.json', runtime)
    _json(run / 'lefse_parameters.json', settings)
    api.finish_stage(job, summary, ['CPM 分母是每个原始样本的过滤后 ASV 总计数；不按层级行之和重复归一化。',
                                   '唯一 SampleID 行用于满足官方 formatter 的 subject 元数据要求；不是配对或重复测量设计。',
                                   '未知末级分类显示为 Others；不把无关未知 ASV 合并作为候选 biomarker。已注释的下级分类保留；未知中间分支使用保存映射的上下文编码，且其本身不检验。'])
    api.begin_stage(job, '02_statistics', '官方 LEfSe：KW、组间一致性与 LDA', settings)
    (run / '02_statistics').mkdir(exist_ok=True)
    pruning = {'tested_features': 0, 'removed_auto_ancestors': [], 'constant_features': []}
    prefix = _prefix(runtime_config, api)
    if data['features']:
        api.command(job, prefix + [runtime['paths']['format'], str(run / '01_inputs/lefse_input_cpm.tsv'),
                                  str(run / '01_inputs/lefse_input.pkl'), '-c', '1', '-u', '2', '-f', 'r', '-o', '-1'], run)
        script = _script(run / '01_inputs/prune_formatted.py', PRUNE_SCRIPT, job)
        api.command(job, prefix + ['python', script, str(run / '01_inputs/lefse_input.pkl'),
                                  str(run / '01_inputs/allowed_features.json'), str(run / '01_inputs/pruning.json')], run)
        pruning = json.loads((run / '01_inputs/pruning.json').read_text(encoding='utf8'))
    else:
        _json(run / '01_inputs/pruning.json', pruning)
    result_file = run / FIT_RESULT
    if pruning['tested_features']:
        api.command(job, prefix + [runtime['paths']['run'], str(run / '01_inputs/lefse_input.pkl'), str(result_file),
                                  '-a', str(settings['kw_alpha']), '-w', str(settings['wilcoxon_alpha']),
                                  '-l', '0', '-y', str(settings['multiclass_strategy']), '-r', 'lda',
                                  '--wilc', '1', '--min_c', '10'], run)
        if not result_file.is_file():
            raise ValueError('LEfSe 命令结束但没有生成结果文件')
    else:
        result_file.write_text('', encoding='utf8')
    results = read_results(result_file)
    features = {row['result_feature']: row for row in data['features']}
    inverse_groups = {safe: original for original, safe in data['groups'].items()}
    if any(row['feature'] not in features or (row['class'] and row['class'] not in inverse_groups) for row in results):
        raise ValueError('LEfSe 结果与本次特征 / 分组输入不一致')
    _tsv(run / '02_statistics/lefse_results.tsv', ['Feature', 'Rank', 'Taxon', 'OriginalLineage', 'Group', 'LDA_score', 'KW_P'],
         [(row['feature'], features[row['feature']]['rank'], features[row['feature']]['taxon'], features[row['feature']]['lineage'],
           inverse_groups.get(row['class'], ''), row['lda_score'] if row['lda_score'] is not None else '',
           row['kw_p'] if row['kw_p'] is not None else '') for row in results])
    status = {'statistics': 'executed' if pruning['tested_features'] else 'skipped_no_variable_known_taxa',
              **pruning, 'positive_lda_features': sum(row['lda_score'] is not None for row in results),
              'result_file': FIT_RESULT, 'result_sha256': hashlib.sha256(result_file.read_bytes()).hexdigest(),
              'settings': settings, 'groups': data['groups'], 'group_style': style, 'preprocess': str(parent)}
    _json(run / 'lefse_complete.json', status)
    api.finish_stage(job, status, ['L2–L7 分类层级联合分析；格式器自动补出的 L1 / 未知中间节点与常数分类不参加检验。',
                                  '原始标准 LEfSe 多组策略为 one-against-all (0)，all-against-all (1) 更严格。',
                                  '只有 Group，没有 subclass / 配对 subject。官方默认一致性步骤在每组不足 10 个样本时使用中位数方向规则，不能称为已进行 Wilcoxon P 检验。',
                                  '官方 .res 最后一列是通过一致性步骤的 KW P，不是 Wilcoxon P；LEfSe 的此流程不提供 BH 校正 P。',
                                  'LDA 以阈值 0 拟合一次并保存正分数；预览按用户阈值严格筛选 abs(LDA) > 阈值，不重新拟合。'])
    result = render_lefse(job, config, run, api)
    job['result'] = result
    return result


def render_lefse(job, config, analysis, api):
    """Render only into the current job's isolated run, never mutate cached scores."""
    analysis, output = Path(analysis), Path(job['run'])
    target = config.get('plot_target')
    if target is not None and (not isinstance(target, dict) or set(target) != {'profile', 'stem'}
            or not isinstance(target.get('profile'), str) or target.get('profile') not in PLOT_TARGETS
            or target.get('stem') != PLOT_TARGETS.get(target.get('profile'))):
        raise ValueError('LEfSe 单图预览目标无效；请选择原图旁的预览按钮')
    kinds = (target['profile'][len('lefse_'):],) if target else ('bar', 'cladogram')
    marker = analysis / 'lefse_complete.json'
    if not marker.is_file():
        raise ValueError('所选 LEfSe 分析尚未完成统计')
    status = json.loads(marker.read_text(encoding='utf8'))
    viewed = analysis / '03_figures/render_summary.json'
    inherited_threshold = (json.loads(viewed.read_text(encoding='utf8')).get('lda_threshold')
                           if viewed.is_file() else None)
    if inherited_threshold is None:
        inherited_threshold = status.get('settings', {}).get('lda_threshold', 2)
    threshold = _number(config.get('lda_threshold', inherited_threshold), 'LDA 阈值')
    result_file = analysis / FIT_RESULT
    if hashlib.sha256(result_file.read_bytes()).hexdigest() != status['result_sha256']:
        raise ValueError('保存的 LEfSe 结果已被改动，不能使用旧缓存预览')
    results = read_results(result_file)
    filtered, selected = filter_results(results, threshold)
    saved_file = analysis / 'parameters.json'
    saved = json.loads(saved_file.read_text(encoding='utf8')) if saved_file.is_file() else {}
    style_file = analysis / 'group_style.json'
    fallback = (json.loads(style_file.read_text(encoding='utf8')) if style_file.is_file()
                else status.get('group_style') or saved.get('config', {}))
    style = resolve_group_style(config, [{'Group': group} for group in status['groups']], fallback=fallback)
    _json(output / 'group_style.json', style)
    visual_file = analysis / 'plot_settings.json'
    visual_fallback = (json.loads(visual_file.read_text(encoding='utf8')) if visual_file.is_file()
                       else saved.get('config', {}))
    visual = resolve_plot_settings(config, fallback=visual_fallback)
    _json(output / 'plot_settings.json', visual)
    runtime_config = {**saved.get('config', {}), **status['settings']}
    if config.get('conda'):
        runtime_config['conda'] = config['conda']
    api.begin_stage(job, '03_figures', 'LEfSe：从保存分数筛选并绘图',
                    {'lda_threshold': threshold, 'statistics_refit': False, 'plot_settings': visual,
                     'plot_target': target, **style})
    folder = output / '03_figures'
    folder.mkdir(exist_ok=True)
    plot_result = folder / 'threshold_results.res'
    with plot_result.open('w', encoding='utf8', newline='') as stream:
        csv.writer(stream, delimiter='\t').writerows(filtered)
    runtime = _runtime(job, runtime_config, api, kinds)
    _json(folder / 'runtime.json', runtime)
    prefix = _prefix(runtime_config, api)
    formats = ('png', 'svg', 'pdf') if visual['general']['formats'] == 'png_svg_pdf' else ('png', 'svg')
    figures = {kind + '_' + ext: str(folder / (kind + '.' + ext)) for kind in kinds for ext in formats}
    top_n = visual['lefse_bar']['top_n']
    bar_selected = (sorted(selected, key=lambda row: (-abs(row['lda_score']), row['feature']))[:top_n]
                    if top_n else selected)
    bar_result = plot_result
    if selected and top_n and 'bar' in kinds:
        bar_result = folder / 'bar_display_results.res'
        with bar_result.open('w', encoding='utf8', newline='') as stream:
            csv.writer(stream, delimiter='\t').writerows(row['raw'] for row in bar_selected)
    if selected:
        class_order = ':'.join(status['groups'][group] for group in style['group_order'])
        _json(folder / 'plot_style.json', {**style, 'groups': status['groups'], 'plot_settings': visual})
        script = _script(folder / 'render_official.py', OFFICIAL_PLOT_SCRIPT, job)
        for kind in kinds:
            for ext in formats:
                args = prefix + ['python', script, kind, str(folder / 'plot_style.json'),
                                 str(bar_result if kind == 'bar' else plot_result), figures[kind + '_' + ext],
                                 '--format', ext, *_plot_options(kind, visual, threshold, class_order)]
                api.command(job, args, output)
    else:
        script = _script(folder / 'render_empty.py', EMPTY_PLOT_SCRIPT, job)
        api.command(job, prefix + ['python', script, str(folder), str(threshold), status['statistics'],
                                  str(output / 'plot_settings.json'), kinds[0] if target else ''], output)
    for filename in figures.values():
        path = Path(filename)
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError('LEfSe 绘图未生成非空 PNG / SVG：' + path.name)
    summary = {'lefse_run': str(analysis), 'render_run': str(output), 'lda_threshold': threshold,
               'biomarker_count': len(selected), 'statistics_refit': False, 'figures': figures,
               'result_sha256': status['result_sha256'], 'no_hits': not selected,
               'statistics_status': status['statistics'], 'plot_settings': visual,
               'bar_biomarker_count': len(bar_selected), 'bar_display_top_n': top_n,
               'cladogram_uses_full_threshold_selection': True, 'plot_target': target,
               'figure_manifest': [{'profile': 'lefse_' + kind, 'stem': '03_figures/' + kind,
                                    **{ext: figures[kind + '_' + ext] for ext in formats}} for kind in kinds], **style}
    _json(folder / 'render_summary.json', summary)
    api.finish_stage(job, summary, ['柱图与分类层级 cladogram 来自同一份原始 LDA 缓存与当前筛选阈值。',
                                   'cladogram 表示分类注释层级，不是序列系统发育树。',
                                   '阈值调整只改变显示，不重新计算 P 值或优化显著性；没有命中时生成明确提示图。',
                                   '分组配色与图例顺序使用当前统一设置；LDA 柱图按指定分组从上到下排列。改变配色或顺序仅重绘，保留原始统计结果。',
                                   '柱图的 Top N 仅按当前阈值命中的绝对 LDA 分数筛选展示；cladogram 仍使用完整阈值结果。字号、方向、主题、尺寸和显示层级只改变图形。',
                                   '名称与原始 Group / 分类的对应关系见原分析 01_inputs/*_mapping.tsv。'])
    job['result'] = summary
    return summary
