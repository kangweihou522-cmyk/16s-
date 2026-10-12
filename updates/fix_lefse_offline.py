"""Run from the workstation directory; repair the known old file without network."""
import hashlib
import os
import subprocess
import tempfile
from pathlib import Path

p = Path('lefse_analysis.py')
old = p.read_bytes()
digest = hashlib.sha256(old).hexdigest()
fixed = 'd50dfd4dc1809af59f548a743ace5cc1ed87811f062a28a3f90d4e9d6b2d3700'
original = '8d74390c0b5d9f747e2a416f0f9965926ecdd960b7111e88ee69134fed78c746'
index = subprocess.run(['git', 'show', ':lefse_analysis.py'], capture_output=True, check=True).stdout
if hashlib.sha256(index).hexdigest() not in (original, fixed):
    raise SystemExit('Git 暂存区包含其他修改，已保留，未执行修复。')
if digest == fixed:
    print('LEfSe 源码已修复；请重启工作站加载它。', flush=True)
else:
    if digest != original:
        raise SystemExit('源码与已确认的旧版不一致，未修改。请提供此信息：' + digest)
    text = old.decode('utf8')
    changes = [
        ("'subclass': None, 'subject': None,", "'subclass': None,\n            'subject': 'unique SampleID (formatter metadata only)', 'paired_design': False,"),
        ('# No SampleID row: -u means repeated subject, which this metadata does not contain.',
         '# LEfSe 1.1.2 requires subject metadata even for independent samples.\n    # Each SampleID is unique; this labels samples without introducing pairs.'),
        ("        writer.writerows([row['input_feature'],", "        writer.writerow(['SampleID', *ids])\n        writer.writerows([row['input_feature'],"),
        ("                                   '未知末级分类显示为 Others；",
         "                                   '唯一 SampleID 行用于满足官方 formatter 的 subject 元数据要求；不是配对或重复测量设计。',\n                                   '未知末级分类显示为 Others；"),
        ("'-c', '1', '-f', 'r', '-o', '-1'", "'-c', '1', '-u', '2', '-f', 'r', '-o', '-1'"),
    ]
    for before, after in changes:
        if text.count(before) != 1:
            raise SystemExit('修复位置不唯一，未修改文件。')
        text = text.replace(before, after, 1)
    new = text.encode('utf8')
    if hashlib.sha256(new).hexdigest() != fixed:
        raise SystemExit('修复内容校验失败，未修改文件。')
    compile(new, str(p), 'exec')
    backup = Path(tempfile.mkdtemp(prefix='16s-lefse-backup-', dir=p.resolve().parent.parent)) / p.name
    backup.write_bytes(old)
    with tempfile.NamedTemporaryFile(dir=p.resolve().parent, delete=False) as f:
        temporary = Path(f.name)
        f.write(new)
    try:
        temporary.chmod(p.stat().st_mode & 0o777)
        os.replace(temporary, p)
    finally:
        temporary.unlink(missing_ok=True)
    print('LEfSe 修复及校验成功；原文件备份：' + str(backup), flush=True)
# Match the published file in the index, allowing a later ff-only pull to absorb it.
subprocess.run(['git', 'add', '--', str(p)], check=True)
