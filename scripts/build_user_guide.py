"""Build offline HTML/Word from one manuscript; --check needs no Pandoc."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE = 'docs/USER_GUIDE.md'
MANIFEST = 'docs/guide-manifest.json'
OUTPUTS = ('static/guide/index.html', 'docs/USER_GUIDE.docx')
RUNTIME = ('app.py', 'index.html', 'sequence_structure.py', 'quality_advice.py', 'qc_report.py',
           'depth_advice.py', 'abundance_tables.py', 'differential_analysis.py',
           'differential_ancombc.R', 'lefse_analysis.py', 'functional_prediction.py',
           'advanced_analysis.py', 'advanced_analysis.R', 'render.R',
           'plots/qiime2_io.R', 'plot_settings.py', 'figure_editor.py')
CSS = '''<style>
:root{color-scheme:light}body{max-width:1080px;margin:0 auto;padding:24px;
font:16px/1.8 system-ui,"Microsoft Yahei",sans-serif;color:#183531;background:white}
h1,h2,h3{line-height:1.4;scroll-margin-top:16px;color:#125a4c}h2{border-top:1px solid #d4e0d9;padding-top:22px}
a{color:#126657;overflow-wrap:anywhere}p,li{overflow-wrap:anywhere}table{border-collapse:collapse;width:100%;font-size:14px}
td,th{border:1px solid #d4e0d9;padding:8px;vertical-align:top}th{background:#eef5f1}
pre{white-space:pre-wrap;overflow:auto;background:#eef5f1;padding:14px;border-radius:8px}
code{font-size:.9em}blockquote{border-left:4px solid #64a28b;margin-left:0;padding:4px 16px;background:#f3f8f5}
#TOC{background:#f3f8f5;padding:16px;border-radius:8px}::selection{background:#ffe07a;color:#172b24}
@media print{body{font-size:11pt;padding:0}h2{break-after:avoid}table,pre{break-inside:avoid}}
</style>'''


def digest(data):
    return hashlib.sha256(data).hexdigest()


def runtime_hashes():
    return {name: digest((ROOT / name).read_bytes().replace(b'\r\n', b'\n')) for name in RUNTIME}


def check():
    value = json.loads((ROOT / MANIFEST).read_text(encoding='utf8'))
    if value.get('source_sha256') != digest((ROOT / SOURCE).read_bytes()):
        raise ValueError('使用说明文稿已修改，请重新构建HTML/Word')
    if value.get('runtime') != runtime_hashes():
        raise ValueError('分析代码已修改，请复核说明并重新构建HTML/Word')
    for name in OUTPUTS:
        if value.get('outputs', {}).get(name) != digest((ROOT / name).read_bytes()):
            raise ValueError('使用说明构建产物已修改：' + name)
    return value


def build():
    pandoc = shutil.which('pandoc')
    if not pandoc:
        raise ValueError('构建说明需要Pandoc；工作站运行及说明阅读不需要安装它')
    raw = (ROOT / SOURCE).read_bytes()
    text = raw.decode('utf8')
    if re.search(r'<script\b|javascript\s*:', text, re.I):
        raise ValueError('说明文稿不允许脚本')
    # Turn GitHub-compatible HTML anchors into native heading identifiers so
    # the same table of contents also has real bookmarks in Word.
    text = re.sub(r'<a id="(guide-\d+)"></a>\s*\n(## [^\n]+)',
                  r'\2 {#\1}', text)
    # Downloaded Word/HTML stays useful outside the repository checkout.
    text = text.replace('](WINDOWS_DEPLOYMENT.md)',
                        '](https://github.com/kangweihou522-cmyk/16s-/blob/main/docs/WINDOWS_DEPLOYMENT.md)')
    source_hash = digest(raw)
    text += '\n\n---\n\n文稿 SHA256：`' + source_hash + '`。网页与 Word 由同一 Markdown 生成。\n'
    environment = dict(os.environ, SOURCE_DATE_EPOCH='946684800')
    with tempfile.TemporaryDirectory(prefix='16s-guide-') as folder:
        folder = Path(folder)
        header = folder / 'header.html'; header.write_text(CSS)
        html, word = folder / 'index.html', folder / 'USER_GUIDE.docx'
        common = [pandoc, '--from=markdown', '--standalone', '--metadata=lang:zh-CN']
        for command in ([*common, '--to=html5', '--section-divs', '--toc', '--toc-depth=2',
                         '--metadata=pagetitle:16S工作站使用说明', '--include-in-header=' + str(header), '-o', str(html)],
                        [*common, '--to=docx', '-o', str(word)]):
            subprocess.run(command, input=text.encode(), env=environment,
                           check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        html_text = html.read_text()
        html_text = re.sub(r'<a href="(https?://[^"]+)"', r'<a target="_blank" rel="noopener noreferrer" href="\1"', html_text)
        if '<script' in html_text.lower():
            raise ValueError('说明HTML不能包含脚本')
        if len(re.findall(r'<h2>\d+\.', html_text)) != 19:
            raise ValueError('说明的19个编号章节未正确生成标题，请检查Markdown分段')
        html.write_text(html_text)
        payloads = {OUTPUTS[0]: html.read_bytes(), OUTPUTS[1]: word.read_bytes()}
        manifest = {'schema': 1, 'source': SOURCE, 'source_sha256': source_hash,
                    'runtime': runtime_hashes(),
                    'outputs': {name: digest(data) for name, data in payloads.items()}}
        for name, data in payloads.items():
            path = ROOT / name; path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        (ROOT / MANIFEST).write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
    return check()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = check() if args.check else build()
    print('使用说明同步检查通过：' + value['source_sha256'][:12])
