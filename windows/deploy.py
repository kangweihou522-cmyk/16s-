"""Deploy immutable 16S source versions, reusing the existing WSL Conda environments."""
from pathlib import Path, PureWindowsPath
from contextlib import contextmanager
from datetime import datetime, timezone
import argparse
import base64
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import uuid

sys.dont_write_bytecode = True

try:
    from . import station_runtime
except ImportError:
    import station_runtime


SERVICE = '16s-workstation'
DIRECTORIES = ['程序', '原始数据', '参考文件', '分析结果', '日志', '配置']
DEFAULT_ENVS = dict(qiime='qiime2-amplicon-2025.7', r='16s-workstation-r',
                    ancom='16s-workstation-ancom', lefse='lefse', picrust2='picrust2')
SOURCE_DIRS = {'plots', 'windows', 'assets', 'static', 'fonts', 'vendor', 'docs'}
SOURCE_FILES = {'README.md', 'VALIDATION.md', 'SHA256SUMS.txt', '启动.sh', '安装依赖.sh',
                '分组模板.csv', '分组模板_两列.csv', '部署工作站.bat', '启动工作站.bat',
                '检查环境.bat', '更新工作站.bat', '停止工作站.bat', '安装高级分析环境.bat'}
IGNORE = {'.git', '.venv', 'runs', 'updates', '__pycache__', '.cache', '.pytest_cache',
          '原始数据', '参考文件', '分析结果', '日志', '配置', '程序', 'node_modules'}
REQUIRED = ['app.py', 'index.html', 'launch.py', 'run_layout.py', 'batch_context.py', 'deployment_context.py',
            'advanced_analysis.py', 'advanced_analysis.R', 'qc_report.py', 'recovery.py', 'external_inputs.py',
            'abundance_tables.py', 'depth_advice.py', 'differential_analysis.py',
            'execution_record.py', 'figure_editor.py', 'functional_prediction.py',
            'ko_network.py', 'ko_map_geometry.py', 'ko_differential.R', 'ko_map_render.R',
            'group_style.py', 'lefse_analysis.py', 'plot_settings.py', 'quality_advice.py',
            'sequence_structure.py', 'render.R', 'taxonomy_plot.R', 'differential_ancombc.R',
            'taxonomy_palettes.json', 'plots/qiime2_io.R', 'plots/plot_alpha_original_style.R',
            'plots/plot_beta_original_style.R', 'windows/deploy.py', 'windows/station_runtime.py',
            'windows/Workstation.ps1', 'windows/Station.ps1', 'docs/USER_GUIDE.md',
            'docs/USER_GUIDE.docx', 'static/guide/index.html', 'docs/guide-manifest.json']


def digest(data):
    return hashlib.sha256(data).hexdigest()


def _json_bytes(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2).encode('utf-8')


def _fingerprint(value):
    return digest(_json_bytes(value))


def _safe_path(path):
    path = Path(path)
    if any(c in str(path) for c in '\r\n\x00'):
        raise ValueError('路径包含无效字符')
    # Resolve only after checking ancestors: never write through a user symlink.
    if any(p.is_symlink() for p in [path, *path.parents]):
        raise ValueError('部署路径不能包含符号链接：'+str(path))
    return path.resolve()


def _root(root):
    root = _safe_path(root)
    if root == Path(root.anchor):
        raise ValueError('请选择工作站子目录，不能使用磁盘或文件系统根目录')
    for name in DIRECTORIES:
        path = root/name
        if path.is_symlink() or (path.exists() and not path.is_dir()):
            raise ValueError('工作站目录无效，未覆盖：'+str(path))
    return root


def read_config(root):
    root = _root(root)
    path = root/'配置/workstation.json'
    if path.is_symlink():
        raise ValueError('配置不能是符号链接，未覆盖：'+str(path))
    if not path.exists():
        return {}
    try:
        value = json.loads(path.read_text(encoding='utf-8-sig'))
    except (ValueError, OSError) as error:
        raise ValueError('工作站配置损坏，未覆盖：'+str(path)) from error
    if (not isinstance(value, dict) or value.get('service_id') != SERVICE
            or value.get('schema') != 1 or value.get('root_wsl') != str(root)):
        raise ValueError('目标配置不属于本工作站，未覆盖：'+str(path))
    try:
        station_runtime._port(value['port'])
        if (not isinstance(value['program_wsl'], str)
                or Path(value['program_wsl']).parent != root/'程序'
                or value['runs_wsl'] != str(root/'分析结果')
                or not isinstance(value['envs'], dict)
                or any(not isinstance(value['envs'].get(key), str) or not value['envs'][key].strip()
                       for key in DEFAULT_ENVS)):
            raise ValueError('配置字段无效')
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError('工作站配置损坏，未覆盖：'+str(path)) from error
    return value


def atomic_write(path, data):
    path = Path(path)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _source_file(name):
    return (name in SOURCE_FILES or (not name.startswith(('.', 'test_'))
            and Path(name).suffix.lower() in {'.py', '.r', '.json', '.html', '.css', '.js'}))


def _ignored(name):
    return (name.lower() in IGNORE or name.lower().endswith(
        ('.fastq', '.fastq.gz', '.fq', '.fq.gz', '.bundle', '.pyc')))


def inventory(source):
    source = _safe_path(source)
    if not source.is_dir():
        raise ValueError('源码目录不存在：'+str(source))
    for name in REQUIRED:
        path = source/name
        if not path.is_file() or path.is_symlink():
            raise ValueError('源码不完整，请使用新版完整 ZIP：'+name)
    files = {}
    candidates = [p for p in source.iterdir() if _source_file(p.name) and not _ignored(p.name)]
    # Walk only known source directories, so a large raw-data directory is never scanned.
    for name in sorted(SOURCE_DIRS):
        folder = source/name
        if folder.is_symlink():
            raise ValueError('源码包含符号链接，未复制：'+name)
        if not folder.exists():
            continue
        for directory, folders, names in os.walk(folder, followlinks=False):
            for item in folders[:]:
                path = Path(directory)/item
                if _ignored(item):
                    folders.remove(item)
                elif path.is_symlink():
                    raise ValueError('源码包含符号链接，未复制：'+str(path.relative_to(source)))
            candidates.extend(Path(directory)/item for item in names if not _ignored(item))
    for path in sorted(candidates):
        name = path.relative_to(source).as_posix()
        if path.is_symlink():
            raise ValueError('源码包含符号链接，未复制：'+name)
        if path.is_file():
            content = path.read_bytes()
            if path.suffix == '.py':
                try:
                    compile(content, name, 'exec')
                except (SyntaxError, ValueError) as error:
                    raise ValueError('Python 源码格式错误，请重新下载完整源码：'+name) from error
            files[name] = digest(content)
    identifier = digest(json.dumps(files, sort_keys=True).encode('utf-8'))
    return identifier, files


def _manifest(release, identifier):
    if release.is_symlink() or not release.is_dir():
        raise ValueError('程序版本目录无效')
    marker = release/'.source-manifest.json'
    if marker.is_symlink():
        raise ValueError('程序清单无效')
    try:
        files = json.loads(marker.read_text(encoding='utf-8'))
        if (not isinstance(files, dict) or not files
                or digest(json.dumps(files, sort_keys=True).encode('utf-8')) != identifier):
            raise ValueError('程序清单与部署计划不一致')
        for name, expected in files.items():
            relative = Path(name)
            if (relative.is_absolute() or '..' in relative.parts or '\\' in name
                    or not re.fullmatch('[0-9a-f]{64}', expected)):
                raise ValueError('程序清单包含无效路径或校验值')
            path = release/relative
            if (any(p.is_symlink() for p in [path, *path.parents] if p != release.parent)
                    or not path.is_file() or digest(path.read_bytes()) != expected):
                raise ValueError('程序文件已改变：'+name)
        expected_names = set(files)|{'.source-manifest.json'}
        for directory, folders, names in os.walk(release, followlinks=False):
            if any((Path(directory)/name).is_symlink() for name in folders+names):
                raise ValueError('程序版本包含符号链接')
            if any((Path(directory)/name).relative_to(release).as_posix() not in expected_names
                   for name in names):
                raise ValueError('程序版本包含额外文件')
        if any(name not in files for name in REQUIRED):
            raise ValueError('程序版本不完整')
    except (OSError, TypeError, json.JSONDecodeError) as error:
        raise ValueError('程序清单损坏') from error
    return files


def find_conda(preferred=None):
    candidates = ([preferred] if preferred else []) + [
        str(Path.home()/name/'bin/conda') for name in ('.miniconda3', 'miniconda3', 'anaconda3')]
    candidates.append(shutil.which('conda'))
    for candidate in candidates:
        if candidate and Path(candidate).is_file() and os.access(candidate, os.X_OK):
            return str(Path(candidate).resolve())
    raise ValueError('找不到已有的 Conda。请在所选 WSL 用户中确认 ~/.miniconda3、~/miniconda3、~/anaconda3 或 PATH；本部署不会安装或迁移环境。')


def _windows_path(value):
    if not isinstance(value, str) or any(ord(c) < 32 for c in value):
        raise ValueError('Windows 路径无效')
    path = PureWindowsPath(value)
    if (not re.fullmatch('[A-Za-z]:', path.drive) or not path.is_absolute()
            or len(path.parts) < 2 or '..' in path.parts
            or any(re.search(r'[<>:"|?*]', part) or part.endswith((' ', '.'))
                   or re.match(r'(?i)^(CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])(?:\.|$)', part)
                   for part in path.parts[1:])):
        raise ValueError('请选择 Windows 磁盘中的完整工作站子目录')
    return str(path)


def ensure_idle(root, configs):
    for port in {c['port'] for c in configs if c and 'port' in c}:
        state = station_runtime.probe(root, port)
        if state['state'] == 'occupied':
            raise ValueError('端口由其他程序或其他工作站使用，未更改配置。')
        if state['state'] == 'running' and state['health']['busy']:
            raise ValueError('还有未结束的分析，请先完成或在网页中取消任务。')


def _version(source, identifier):
    try:
        result = subprocess.run(['git', '-C', str(source), 'rev-parse', '--short=12', 'HEAD'],
                                capture_output=True, text=True, timeout=5, check=True)
        revision = result.stdout.strip()
        if re.fullmatch('[0-9a-f]{7,40}', revision):
            return 'git:'+revision+'-'+identifier[:12]
    except (OSError, subprocess.SubprocessError):
        pass
    return 'sha256:'+identifier[:16]


def prepare(source, root, windows_root, distro, port, user=''):
    source, root = _safe_path(source), _root(root)
    if root == source or root.is_relative_to(source):
        raise ValueError('工作站目录不能位于待复制的源码目录内')
    station_runtime._port(port)
    windows_root = _windows_path(windows_root)
    if any(c in distro+user for c in '\r\n\x00') or not distro.strip():
        raise ValueError('WSL 发行版或用户无效')
    old = read_config(root)
    conda = find_conda(old.get('conda_wsl'))
    ensure_idle(root, [old, dict(port=port)])
    identifier, files = inventory(source)
    for name in DIRECTORIES:
        (root/name).mkdir(parents=True, exist_ok=True)
    release = root/'程序'/identifier[:16]
    suffix = 0
    while release.exists() or release.is_symlink():
        try:
            if _manifest(release, identifier) == files:
                break
        except ValueError:
            pass
        suffix += 1
        release = root/'程序'/(identifier[:16]+'-'+str(suffix))
    if not release.exists():
        stage = root/'程序'/('.staging-'+uuid.uuid4().hex)
        stage.mkdir()
        try:
            for name, expected in files.items():
                target = stage/name
                target.parent.mkdir(parents=True, exist_ok=True)
                origin = source/name
                if any(p.is_symlink() for p in [origin, *origin.parents]):
                    raise ValueError('复制时源码路径改变，请重试：'+name)
                shutil.copyfile(origin, target)
                if digest(target.read_bytes()) != expected:
                    raise ValueError('复制时源码改变，请重试：'+name)
            atomic_write(stage/'.source-manifest.json', _json_bytes(files))
            _manifest(stage, identifier)
            stage.rename(release)
        finally:
            if stage.exists():
                shutil.rmtree(stage)
    now = datetime.now(timezone.utc).isoformat()
    envs = dict(DEFAULT_ENVS, **old.get('envs', {}))
    config = dict(old, schema=1, service_id=SERVICE, root_windows=windows_root,
                  root_wsl=str(root), program_wsl=str(release),
                  program_windows=windows_root+'\\程序\\'+release.name,
                  runs_wsl=str(root/'分析结果'), runs_windows=windows_root+'\\分析结果',
                  distro=distro, wsl_user=user, port=port, conda_wsl=conda, envs=envs,
                  source_digest=identifier, version=_version(source, identifier),
                  created_at=old.get('created_at', now), updated_at=now)
    plan = dict(config=config, previous_config=_fingerprint(old))
    plan_file = root/'配置'/('pending-'+uuid.uuid4().hex+'.json')
    atomic_write(plan_file, _json_bytes(plan))
    return dict(config=config, plan_file=str(plan_file))


def _plan(root, plan_file, verify_release=True):
    root, plan_file = _root(root), _safe_path(plan_file)
    if plan_file.parent != root/'配置' or not re.fullmatch('pending-[0-9a-f]{32}\\.json', plan_file.name):
        raise ValueError('部署计划路径无效')
    try:
        plan = json.loads(plan_file.read_text(encoding='utf-8'))
        config = plan['config']
        if (not isinstance(config, dict) or config.get('schema') != 1 or config.get('service_id') != SERVICE
                or config['root_wsl'] != str(root)
                or config['runs_wsl'] != str(root/'分析结果')):
            raise ValueError('部署计划与目标目录不一致')
        if not isinstance(plan.get('previous_config'), str) or not re.fullmatch('[0-9a-f]{64}', plan['previous_config']):
            raise ValueError('部署计划缺少原始配置校验值')
        release = _safe_path(config['program_wsl'])
        identifier = config['source_digest']
        if (release.parent != root/'程序' or not re.fullmatch('[0-9a-f]{64}', identifier)
                or not re.fullmatch(identifier[:16]+r'(-[1-9][0-9]*)?', release.name)):
            raise ValueError('部署程序路径无效')
        windows_root = _windows_path(config['root_windows'])
        if (config['program_windows'] != windows_root+'\\程序\\'+release.name
                or config['runs_windows'] != windows_root+'\\分析结果'):
            raise ValueError('Windows 部署路径与计划不一致')
        station_runtime._port(config['port'])
        if (not isinstance(config.get('envs'), dict)
                or any(not isinstance(config['envs'].get(key), str) or not config['envs'][key].strip()
                       for key in DEFAULT_ENVS)):
            raise ValueError('Conda 环境配置无效')
        if verify_release:
            _manifest(release, identifier)
    except (KeyError, TypeError, OSError, json.JSONDecodeError) as error:
        raise ValueError('部署计划损坏，未更改配置') from error
    return plan


def launchers(config):
    helpers = {'Station.ps1': (Path(config['program_wsl'])/'windows/Station.ps1').read_bytes()}
    for name, action in [('启动工作站.bat', 'Start'), ('检查环境.bat', 'Check'),
                         ('更新工作站.bat', 'Update'), ('停止工作站.bat', 'Stop'),
                         ('部署工作站.bat', 'Install'), ('安装高级分析环境.bat', 'InstallAdvanced')]:
        pause = 'pause\r\n' if action == 'Check' else 'if not "%STATION_EXIT%"=="0" pause\r\n'
        helpers[name] = ('@echo off\r\nsetlocal\r\npowershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass '
                         '-File "%~dp0Station.ps1" -Action '+action+' %*\r\n'
                         'set "STATION_EXIT=%ERRORLEVEL%"\r\n'+pause+'exit /b %STATION_EXIT%\r\n').encode('ascii')
    return helpers


@contextmanager
def _activation_lock(root):
    lock = root/'配置/.deployment-lock'
    try:
        lock.mkdir()
    except FileExistsError as error:
        raise ValueError('另一部署正在进行；请等待它完成，未覆盖配置。') from error
    try:
        yield
    finally:
        lock.rmdir()


def activate(root, plan_file):
    root = _root(root)
    with _activation_lock(root):
        plan = _plan(root, plan_file)
        config, old = plan['config'], read_config(root)
        if _fingerprint(old) != plan['previous_config']:
            raise ValueError('另一部署已改变当前版本，请重试，未覆盖配置')
        find_conda(config.get('conda_wsl'))
        ensure_idle(root, [old, config])
        helpers = launchers(config)
        originals = {}
        for name, content in helpers.items():
            target = root/name
            if target.is_symlink() or (target.exists() and not target.is_file()):
                raise ValueError('启动文件不是普通文件，未覆盖：'+str(target))
            previous = target.read_bytes() if target.exists() else None
            if (previous is not None and digest(previous) not in
                    {digest(content), old.get('managed_files', {}).get(name)}):
                raise ValueError('启动文件已由你修改，未覆盖，请先保留或重命名：'+str(target))
            originals[name] = previous
        # The server checks busy again under its job lock when accepting this request.
        for port in {c['port'] for c in [old, config] if c and 'port' in c}:
            station_runtime.stop(root, port)
        config['managed_files'] = {name: digest(content) for name, content in helpers.items()}
        plan['rollback'] = dict(config=old, active_config=_fingerprint(config),
                                helpers={name: base64.b64encode(content).decode('ascii')
                                         if content is not None else None
                                         for name, content in originals.items()})
        atomic_write(Path(plan_file), _json_bytes(plan))
        written = []
        try:
            for name, content in helpers.items():
                atomic_write(root/name, content)
                written.append(name)
            # The pointer is written last. A failed deployment keeps the prior configuration.
            atomic_write(root/'配置/workstation.json', _json_bytes(config))
        except BaseException:
            for name in reversed(written):
                if originals[name] is None:
                    (root/name).unlink(missing_ok=True)
                else:
                    atomic_write(root/name, originals[name])
            raise
        return config


def rollback(root, plan_file):
    root = _root(root)
    with _activation_lock(root):
        # Recovery must also work when the failed candidate has damaged its files.
        plan = _plan(root, plan_file, verify_release=False)
        state, active = plan.get('rollback'), read_config(root)
        if not isinstance(state, dict) or _fingerprint(active) != state.get('active_config'):
            raise ValueError('当前版本已改变，不能回滚此部署，未覆盖配置。')
        expected = {}
        managed = active.get('managed_files', {})
        if set(managed) != {'Station.ps1', '启动工作站.bat', '检查环境.bat', '更新工作站.bat', '停止工作站.bat', '部署工作站.bat', '安装高级分析环境.bat'}:
            raise ValueError('回滚配置缺少受管启动文件，未覆盖配置。')
        for name, expected_digest in managed.items():
            path = root/name
            if path.is_symlink() or not path.is_file():
                raise ValueError('启动文件已由你修改，不能回滚，未覆盖：'+str(path))
            content = path.read_bytes()
            if digest(content) != expected_digest:
                raise ValueError('启动文件已由你修改，不能回滚，未覆盖：'+str(path))
            expected[name] = content
        ensure_idle(root, [active])
        station_runtime.stop(root, active['port'])
        originals = state['helpers']
        if set(originals) != set(expected):
            raise ValueError('回滚计划损坏，未更改配置')
        previous = state['config']
        if previous and (previous.get('service_id') != SERVICE or previous.get('root_wsl') != str(root)):
            raise ValueError('回滚配置不属于本工作站')
        decoded = {name: base64.b64decode(content, validate=True) if content is not None else None
                   for name, content in originals.items()}
        restored = []
        try:
            for name, content in decoded.items():
                if content is None:
                    (root/name).unlink(missing_ok=True)
                else:
                    atomic_write(root/name, content)
                restored.append(name)
            path = root/'配置/workstation.json'
            if previous:
                atomic_write(path, _json_bytes(previous))
            else:
                path.unlink(missing_ok=True)
        except BaseException:
            for name in restored:
                atomic_write(root/name, expected[name])
            raise
        Path(plan_file).unlink(missing_ok=True)
        return previous


def complete(root, plan_file):
    root = _root(root)
    with _activation_lock(root):
        plan = _plan(root, plan_file)
        if _fingerprint(read_config(root)) != plan.get('rollback', {}).get('active_config'):
            raise ValueError('当前版本已改变，未删除部署记录。')
        Path(plan_file).unlink(missing_ok=True)
    return dict(completed=True)


def check(root, plan_file=None):
    root = _root(root)
    config = _plan(root, plan_file)['config'] if plan_file else read_config(root)
    if not config:
        raise ValueError('尚未部署工作站，请先部署。')
    conda = find_conda(config.get('conda_wsl'))
    warnings, available = [], set()
    try:
        result = subprocess.run([conda, 'env', 'list', '--json'], capture_output=True,
                                text=True, timeout=30, check=True)
        values = json.loads(result.stdout)
        available = {Path(item).name for item in values['envs'] if isinstance(item, str)}
    except (OSError, subprocess.SubprocessError, ValueError, KeyError, TypeError):
        warnings.append('无法列出已有 Conda 环境；请在 WSL 中运行 conda env list 核对。')
    envs = {key: dict(name=name, available=name in available)
            for key, name in config['envs'].items()}
    missing = [value['name'] for value in envs.values() if not value['available']]
    if missing:
        warnings.append('以下环境未找到：'+', '.join(missing)+'。仍可打开网页；对应分析须先使用已有可用环境。')
    optional = {'advanced': {'name': '16s-workstation-advanced',
                              'available': '16s-workstation-advanced' in available}}
    return dict(conda_wsl=conda, conda_available=True, envs=envs, optional_envs=optional,
                missing_envs=missing, ready=not missing and not warnings, warnings=warnings)


def install_advanced(root):
    """Prepare a separate optional R environment, never mutate QIIME/plot R."""
    root = _root(root)
    config = read_config(root)
    if not config:
        raise ValueError('请先完成工作站部署')
    ensure_idle(root, [config])
    name = '16s-workstation-advanced'
    if name in config['envs'].values():
        raise ValueError('高级环境名与现有分析环境冲突，未修改任何环境')
    conda = find_conda(config.get('conda_wsl'))
    try:
        listed = subprocess.run([conda, 'env', 'list', '--json'], capture_output=True,
                                text=True, timeout=30, check=True)
        envs = json.loads(listed.stdout)['envs']
        exists = any(Path(item).name == name for item in envs)
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
        raise ValueError('无法确认已有 Conda 环境，未安装；请查看 conda env list') from error
    packages = ['r-base=4.4', 'r-vegan>=2.6,<3', 'r-permute', 'r-ppcor',
                'r-ggplot2', 'r-jsonlite', 'r-svglite', 'r-systemfonts', 'r-mass']
    command = [conda, 'install' if exists else 'create', '-y', '-n', name,
               '--override-channels', '-c', 'conda-forge', '--strict-channel-priority', *packages]
    print('仅准备独立高级分析环境 '+name+'；已有 QIIME、R绘图、ANCOM、LEfSe、PICRUSt2 环境保留。', flush=True)
    result = subprocess.run(command)
    if result.returncode:
        raise ValueError('高级环境安装失败（退出码 '+str(result.returncode)+'）；现有分析环境未被此命令修改')
    env = dict(os.environ)
    for key in ('R_HOME', 'LD_LIBRARY_PATH', 'LD_PRELOAD'):
        env.pop(key, None)
    for key in ('R_LIBS', 'R_LIBS_USER', 'R_LIBS_SITE', 'R_PROFILE', 'R_PROFILE_USER', 'R_ENVIRON', 'R_ENVIRON_USER'):
        env[key] = '/dev/null'
    check_code = 'stopifnot(all(vapply(c("vegan","permute","ppcor","ggplot2","jsonlite","svglite","systemfonts"),requireNamespace,logical(1),quietly=TRUE)))'
    checked = subprocess.run([conda, 'run', '--no-capture-output', '-n', name,
                              'Rscript', '--vanilla', '-e', check_code], env=env)
    if checked.returncode:
        raise ValueError('高级环境已准备，但 R 依赖检查未通过，请保留终端错误')
    return {'environment': name, 'available': True}


def start(root):
    root = _root(root)
    config = read_config(root)
    if not config:
        raise ValueError('尚未部署工作站，请先部署。')
    find_conda(config.get('conda_wsl'))
    ensure_idle(root, [config])
    if station_runtime.probe(root, config['port'])['state'] == 'running':
        raise ValueError('此工作站已运行，请打开对应网页。')
    release = _safe_path(config['program_wsl'])
    if release.parent != root/'程序':
        raise ValueError('程序目录不属于本工作站')
    _manifest(release, config['source_digest'])
    # Imports must not create __pycache__ inside an immutable source version.
    environment = dict(os.environ, STATION_ROOT=str(root), PYTHONDONTWRITEBYTECODE='1',
                       PYTHONUNBUFFERED='1')
    command = [sys.executable, '-B', '-u', str(release/'app.py'), '--port', str(config['port']), '--no-browser']
    log_file = root/'日志/server.log'
    if log_file.is_symlink():
        raise ValueError('服务日志不能是符号链接')
    with log_file.open('ab', buffering=0) as log:
        process = subprocess.Popen(command, cwd=release, env=environment,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        try:
            while True:
                try:
                    chunk = os.read(process.stdout.fileno(), 65536)
                except KeyboardInterrupt:
                    if process.poll() is None:
                        process.send_signal(signal.SIGINT)
                    continue
                if not chunk:
                    break
                log.write(chunk)
                output = getattr(sys.stdout, 'buffer', sys.stdout)
                try:
                    output.write(chunk if hasattr(sys.stdout, 'buffer') else chunk.decode('utf-8', errors='replace'))
                    output.flush()
                except BrokenPipeError:
                    pass
            return process.wait()
        except KeyboardInterrupt:
            if process.poll() is None:
                process.send_signal(signal.SIGINT)
            return process.wait()
        finally:
            process.stdout.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    stage = sub.add_parser('prepare')
    for name in ['source', 'root', 'windows_root', 'distro']:
        stage.add_argument(name)
    stage.add_argument('port', type=int)
    stage.add_argument('user', nargs='?', default='')
    for action in ['activate', 'rollback', 'complete']:
        commit = sub.add_parser(action)
        commit.add_argument('root'); commit.add_argument('plan_file')
    verify = sub.add_parser('check')
    verify.add_argument('root'); verify.add_argument('--plan', dest='plan_file')
    run = sub.add_parser('start'); run.add_argument('root')
    advanced = sub.add_parser('install-advanced'); advanced.add_argument('root')
    args = parser.parse_args()
    try:
        if args.action == 'prepare':
            result = prepare(args.source, args.root, args.windows_root, args.distro, args.port, args.user)
        elif args.action == 'activate':
            result = activate(args.root, args.plan_file)
        elif args.action == 'rollback':
            result = rollback(args.root, args.plan_file)
        elif args.action == 'complete':
            result = complete(args.root, args.plan_file)
        elif args.action == 'check':
            result = check(args.root, args.plan_file)
        elif args.action == 'install-advanced':
            result = install_advanced(args.root)
        else:
            sys.exit(start(args.root))
        print(json.dumps(result, ensure_ascii=False))
    except (ValueError, OSError) as error:
        parser.exit(1, str(error)+'\n')
