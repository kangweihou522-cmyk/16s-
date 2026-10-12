"""Focused real-R verification of legacy and organized saved-result layouts."""
import os
import shutil
import subprocess

from run_layout import artifact_path, initialize_layout
from test_render_fixture import fixture, validate_positive, RSCRIPT, SCRIPT, ASSETS, ROOT


def main():
    outputs = []
    for organized in (False, True):
        name = 'layout_organized' if organized else 'layout_legacy'
        run, parent = fixture(name)
        if organized:
            for base, kind in ((parent, 'preprocess'), (run, 'diversity')):
                original_files = [p for p in base.rglob('*') if p.is_file()]
                initialize_layout(base, kind)
                for source in original_files:
                    target = artifact_path(base, source.relative_to(base))
                    if source != target:
                        target.parent.mkdir(parents=True, exist_ok=True)
                        shutil.move(str(source), str(target))
        env = os.environ.copy(); env['XDG_CACHE_HOME'] = str(ROOT / 'cache')
        result = subprocess.run([str(RSCRIPT), str(SCRIPT), str(run), 'DejaVu Serif', str(parent), str(ASSETS)], text=True, capture_output=True, env=env)
        if result.returncode: raise AssertionError(result.stderr)
        out = artifact_path(run, 'R_original_style')
        count = validate_positive(run, True, out=out)
        outputs.append(out)
        print(f'{name}: {count} genuine SVG/PDF/PNG outputs validated', flush=True)
    legacy, organized = outputs
    for old in legacy.rglob('*.csv'):
        new = organized / old.relative_to(legacy)
        assert new.read_bytes() == old.read_bytes(), f'Layout changed saved-result interpretation: {old.name}'
    print('PASS: old/new layouts preserve alpha inputs, PCoA coordinates and unrarefied taxonomy tables.')


if __name__ == '__main__': main()
