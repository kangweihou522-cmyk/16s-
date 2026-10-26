"""Read scheme-A deployment settings without changing existing analysis paths."""
import json
import os
from pathlib import Path


SERVICE_ID = '16s-workstation'


def load_context(program, root=None):
    selected = root if root is not None else os.environ.get('STATION_ROOT', '')
    if not selected:
        return {}
    directory = Path(selected).expanduser().resolve()
    config = json.loads((directory / '配置/workstation.json').read_text(encoding='utf-8-sig'))
    if config.get('schema') != 1 or config.get('service_id') != SERVICE_ID:
        raise ValueError('部署配置不属于本 16S 工作站')
    if Path(config.get('root_wsl', '/')).resolve() != directory:
        raise ValueError('部署配置的安装目录不匹配')
    if Path(config.get('program_wsl', '/')).resolve() != Path(program).resolve():
        raise ValueError('当前程序不是已激活版本，请使用安装目录中的启动工作站.bat')
    if Path(config.get('runs_wsl', '/')).resolve() != directory / '分析结果':
        raise ValueError('部署配置的默认结果目录不匹配')
    return config


def frontend_defaults(config):
    if not config:
        return {}
    envs = config.get('envs', {})
    names = {'env': 'qiime', 'r_env': 'r', 'ancom_env': 'ancom',
             'lefse_env': 'lefse', 'picrust_env': 'picrust2'}
    defaults = {field: envs[key] for field, key in names.items() if envs.get(key)}
    defaults.update(output=config['runs_wsl'], conda=config.get('conda_wsl', ''))
    return defaults


def public_context(config):
    keys = ('root_windows', 'root_wsl', 'program_windows', 'program_wsl',
            'runs_windows', 'runs_wsl', 'distro', 'wsl_user', 'port', 'version')
    return {'deployed': bool(config), **{key: config[key] for key in keys if key in config},
            'defaults': frontend_defaults(config)}
