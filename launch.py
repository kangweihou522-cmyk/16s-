"""Update a clean main checkout without prompting, then start the local app."""
from pathlib import Path
import os, shutil, subprocess, sys

ROOT = Path(__file__).resolve().parent

def auto_update():
    if not shutil.which('git') or not (ROOT/'.git').exists():
        print('当前目录没有 Git 信息，使用本机版本。',flush=True)
        return False
    env = os.environ.copy()
    env['GIT_TERMINAL_PROMPT'] = '0'
    def git(*args):
        return subprocess.run(['git','-C',str(ROOT),*args],env=env,text=True,capture_output=True,timeout=30)
    try:
        status = git('status','--porcelain')
        branch = git('symbolic-ref','--quiet','--short','HEAD')
        if status.returncode or branch.returncode:
            print('无法检查 Git 状态，使用本机版本。',flush=True)
            return False
        if status.stdout.strip():
            print('发现本地文件改动，保留改动并跳过自动更新。',flush=True)
            return False
        if branch.stdout.strip() != 'main':
            print('当前不是 main 分支，使用本机版本。',flush=True)
            return False
        before = git('rev-parse','HEAD').stdout.strip()
        print('正在检查工作站源码更新…',flush=True)
        result = git('pull','--ff-only','origin','main')
        if result.returncode:
            print('自动更新未完成（网络、认证或分支状态）；使用本机版本。',flush=True)
            return False
        after = git('rev-parse','HEAD').stdout.strip()
        changed = before != after
        print('源码已更新。' if changed else '源码已经是最新版本。',flush=True)
        return changed
    except (OSError, subprocess.TimeoutExpired):
        print('更新检查不可用或超时，使用本机版本。',flush=True)
        return False

def main():
    args = sys.argv[1:]
    skip = '--no-update' in args
    args = [a for a in args if a != '--no-update']
    if not skip and os.environ.get('WORKSTATION_UPDATED_ONCE') != '1' and auto_update():
        os.environ['WORKSTATION_UPDATED_ONCE'] = '1'
        os.execv(sys.executable,[sys.executable,str(ROOT/'launch.py'),*args])
    os.execv(sys.executable,[sys.executable,str(ROOT/'app.py'),*args])

if __name__ == '__main__': main()
