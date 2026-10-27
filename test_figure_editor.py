"""Image-level requests preserve cached analyses and other images' settings."""
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import app
from figure_editor import describe_figure, figure_catalog


class FigureEditorTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(); self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.parent = self.root / 'preprocess'; self.parent.mkdir(); app.marker(self.parent)
        self.saved = {'samples': [{'SampleID': '001', 'Group': 'CON'}, {'SampleID': '002', 'Group': 'LDC'},
                                  {'SampleID': '003', 'Group': 'HDC'}], 'font': 'DejaVu Serif',
                      'plot_settings': {'general': {'dpi': 144}, 'alpha': {'show_points': False}}}
        self.write(self.parent / 'parameters.json', self.saved)
        self.analysis = self.parent / '07_diversity' / 'diversity'; self.analysis.mkdir(parents=True); app.marker(self.analysis)
        self.write(self.analysis / 'parameters.json', {'preprocess': str(self.parent), 'config': self.saved})
        self.write(self.analysis / 'qiime_statistics_complete.json', {})
        self.plotroot = self.analysis / '04_R_figures'
        self.figure(self.plotroot / 'alpha/shannon.png')
        self.figure(self.plotroot / 'alpha/simpson.png')
        self.summary(self.plotroot)

    def write(self, filename, value):
        filename.parent.mkdir(parents=True, exist_ok=True)
        filename.write_text(json.dumps(value), encoding='utf8')

    def figure(self, filename):
        filename.parent.mkdir(parents=True, exist_ok=True)
        filename.write_bytes(b'FIXTURE IMAGE - not rendered by R')
        filename.with_suffix('.svg').write_text('FIXTURE SVG')
        return filename

    def summary(self, directory, settings=None):
        self.write(directory / 'render_summary.json', {'source': 'saved_QIIME2_outputs',
                   'input_run': str(self.analysis), 'plot_settings': settings or self.saved['plot_settings']})

    def snapshot(self):
        return {file: file.read_bytes() for file in self.parent.rglob('*') if file.is_file()}

    def test_catalog_recognizes_images_and_skips_nested_preview_versions(self):
        preview = self.analysis / '05_replots' / 'preview'; preview.mkdir(parents=True); app.marker(preview)
        self.write(preview / 'parameters.json', {'diversity_run': str(self.analysis), 'config': self.saved})
        self.figure(preview / '04_R_figures/R_original_style/alpha/shannon.png')
        self.summary(preview / '04_R_figures/R_original_style')
        self.figure(self.analysis / 'native_qzv_export/random.png')
        base = figure_catalog(str(self.parent), app.module_api())
        self.assertEqual(len(base['figures']), 2)
        self.assertEqual({f['stem'] for f in base['figures']}, {'alpha/shannon', 'alpha/simpson'})
        self.assertEqual(len(figure_catalog(str(preview), app.module_api())['figures']), 1)

    def test_unapproved_or_unsupported_images_and_other_profiles_rejected(self):
        external = self.figure(self.root / 'external.png')
        with self.assertRaises(ValueError): describe_figure(str(external), app.module_api())
        unsupported = self.figure(self.parent / 'unknown.png')
        with self.assertRaises(ValueError): describe_figure(str(unsupported), app.module_api())
        for settings in ({'heatmap': {'width': 10}}, {'alpha': {'width': 0}}):
            job = {'status': 'running'}
            with patch.object(app, 'command') as command:
                app.worker(job, 'figure-render', {'image': str(self.plotroot / 'alpha/shannon.png'), 'plot_settings': settings})
            self.assertEqual(job['status'], 'failed'); command.assert_not_called()
            # Early failures now get a separate durable task record; no
            # rejected preview or source analysis is written or executed.
            failure = Path(job['run'])
            self.assertFalse(failure.is_relative_to(self.parent))
            self.assertEqual(json.loads((failure/'task_state.json').read_text())['status'], 'failed')
            self.assertFalse(list(failure.rglob('*.png')))

    def test_single_diversity_preview_and_repeat_inherit_viewed_settings(self):
        old = self.snapshot()
        calls = []
        def command(job, args, cwd, **kwargs):
            self.assertEqual(args[2], str(self.analysis))
            settings = json.loads(Path(args[-2]).read_text())
            target = json.loads(Path(args[-1]).read_text())
            self.assertEqual(target, {'profile': 'alpha', 'stem': 'alpha/shannon'})
            self.assertFalse(settings['alpha']['show_points'])
            self.assertEqual(settings['general']['dpi'], 96)
            calls.append(settings)
            output = Path(args[-3]); self.figure(output / 'alpha/shannon.png'); self.summary(output, settings)
        with patch.object(app, 'q') as qiime, patch.object(app, 'command', side_effect=command), \
             patch.object(app, 'r_command', side_effect=lambda config, *args: ['Rscript', *args]), \
             patch.object(app, 'r_execution_env', return_value={}):
            job = {'status': 'running'}
            app.worker(job, 'figure-render', {'image': str(self.plotroot / 'alpha/shannon.png'),
                       'plot_settings': {'general': {'dpi': 96}, 'alpha': {'width': 7}}})
            self.assertEqual(job['status'], 'done', job.get('error'))
            new = job['result']['figure']
            self.assertEqual(new['source_image'], str(self.plotroot / 'alpha/shannon.png'))
            self.assertEqual(new['plot_settings']['alpha']['width'], 7)
            next_job = {'status': 'running'}
            app.worker(next_job, 'figure-render', {'image': new['image'], 'plot_settings': {'alpha': {'title': 'Changed title'}}})
            self.assertEqual(next_job['status'], 'done', next_job.get('error'))
            self.assertEqual(calls[-1]['alpha']['width'], 7)
            self.assertEqual(calls[-1]['alpha']['title'], 'Changed title')
            qiime.assert_not_called()
        for filename, data in old.items(): self.assertEqual(filename.read_bytes(), data)
        sibling = describe_figure(str(self.plotroot / 'alpha/simpson.png'), app.module_api())
        self.assertEqual(sibling['plot_settings']['general']['dpi'], 144)
        self.assertIsNone(sibling['plot_settings']['alpha']['width'])

    def test_failed_preview_keeps_old_image_and_parameters(self):
        old = self.snapshot()
        job = {'status': 'running'}
        with patch.object(app, 'command', side_effect=RuntimeError('synthetic R failure')), \
             patch.object(app, 'r_command', side_effect=lambda config, *args: ['Rscript', *args]), \
             patch.object(app, 'r_execution_env', return_value={}):
            app.worker(job, 'figure-render', {'image': str(self.plotroot / 'alpha/shannon.png'), 'plot_settings': {}})
        self.assertEqual(job['status'], 'failed')
        self.assertNotIn('figure', job.get('result', {}))
        for filename, data in old.items(): self.assertEqual(filename.read_bytes(), data)

    def test_picture_font_inherits_viewed_version_instead_of_analysis_default(self):
        preview = self.analysis / '05_replots' / 'old_preview'; preview.mkdir(parents=True); app.marker(preview)
        config = dict(self.saved, font='DejaVu Sans')
        self.write(preview / 'parameters.json', {'diversity_run': str(self.analysis), 'config': config})
        image = self.figure(preview / '04_R_figures/R_original_style/alpha/shannon.png')
        self.summary(image.parent.parent)
        with patch.object(app, 'render_saved_diversity', side_effect=RuntimeError('Stop before rendering')) as renderer:
            job = {'status': 'running'}
            app.worker(job, 'figure-render', {'image': str(image), 'plot_settings': {}})
            self.assertEqual(renderer.call_args.args[1]['font'], 'DejaVu Sans')

    def test_differential_preview_preserves_statistical_controls_and_comparison(self):
        source = self.parent / 'differential'; source.mkdir(); app.marker(source)
        self.write(source / 'parameters.json', {'preprocess': str(self.parent), 'config': self.saved})
        counts = source / 'counts.tsv'; counts.write_text('FeatureID\t001\t002\t003\na\t1\t2\t3\n')
        meta = source / 'metadata.tsv'; meta.write_text('sample-id\tGroup\n001\tCON\n002\tLDC\n003\tHDC\n')
        folder = source / '02_differential/L6_genus'; folder.mkdir(parents=True)
        (folder / 'ancombc2_fit.rds').write_bytes(b'SAVED FIT FIXTURE - not read by R')
        image = self.figure(folder / 'comparisons/pair_002/effect.png')
        self.write(image.parent / 'comparison.json', {'contrast': 'GroupHDC_GroupLDC'})
        self.write(source / 'differential_options.json', {'method': 'ancombc2', 'alpha': .0123,
              'reference_group': 'CON', 'groups': ['CON', 'LDC', 'HDC'], 'metadata_tsv': str(meta),
              'metadata_sha256': hashlib.sha256(meta.read_bytes()).hexdigest(),
              'tables': [{'level': 2, 'name': 'phylum', 'counts_tsv': str(source / 'unrelated_missing_counts.tsv')},
                         {'level': 6, 'name': 'genus', 'counts_tsv': str(counts), 'counts_sha256': hashlib.sha256(counts.read_bytes()).hexdigest()}]})
        old = self.snapshot()
        def command(job, args, cwd, **kwargs):
            opts = json.loads(Path(args[-2]).read_text())
            self.assertEqual(opts['alpha'], .0123); self.assertEqual(opts['reference_group'], 'CON')
            self.assertEqual(len(opts['tables']), 1)
            self.assertEqual(opts['plot_target']['comparison_key'], 'GroupHDC_GroupLDC')
            output = self.figure(Path(cwd) / (opts['plot_target']['stem'] + '.png'))
            self.write(output.parent / 'comparison.json', {'contrast': opts['plot_target']['comparison_key']})
            self.write(Path(cwd) / '02_differential/summary.json', {'statistics_refit': False})
        with patch.object(app, 'q') as qiime, patch('differential_analysis.run_differential') as fit, \
             patch.object(app, 'command', side_effect=command), \
             patch.object(app, 'r_command', side_effect=lambda config, *args: ['Rscript', *args]), \
             patch.object(app, 'r_execution_env', return_value={}):
            job = {'status': 'running'}
            app.worker(job, 'figure-render', {'image': str(image), 'diff_alpha': .99,
                       'group_order': ['HDC', 'CON', 'LDC'], 'plot_settings': {'effect': {'width': 9}}})
            self.assertEqual(job['status'], 'done', job.get('error'))
            qiime.assert_not_called(); fit.assert_not_called()
        for filename, data in old.items(): self.assertEqual(filename.read_bytes(), data)

    def test_lefse_card_threshold_and_saved_style_inherit_without_refit(self):
        source = self.parent / 'lefse'; source.mkdir(); app.marker(source)
        self.write(source / 'parameters.json', {'preprocess': str(self.parent), 'config': self.saved})
        self.write(source / 'lefse_complete.json', {'groups': {'CON': 'g1', 'LDC': 'g2', 'HDC': 'g3'}})
        image = self.figure(source / '03_figures/bar.png')
        self.write(image.parent / 'render_summary.json', {'lefse_run': str(source), 'lda_threshold': 4.2})
        old = self.snapshot()
        def renderer(job, config, analysis, api):
            self.assertEqual(config['lda_threshold'], 4.2)
            self.assertEqual(config['plot_target'], {'profile': 'lefse_bar', 'stem': '03_figures/bar'})
            output = self.figure(Path(job['run']) / '03_figures/bar.png')
            self.write(output.parent / 'render_summary.json', {'lefse_run': str(source), 'lda_threshold': 4.2,
                      'plot_settings': config['plot_settings']})
            job['result'] = {'statistics_refit': False}
        with patch('lefse_analysis.render_lefse', side_effect=renderer), patch('lefse_analysis.run_lefse') as fit:
            job = {'status': 'running'}
            app.worker(job, 'figure-render', {'image': str(image), 'plot_settings': {'lefse_bar': {'width': 8}}})
            self.assertEqual(job['status'], 'done', job.get('error')); fit.assert_not_called()
            self.assertEqual(job['result']['figure']['lda_threshold'], 4.2)
        for filename, data in old.items(): self.assertEqual(filename.read_bytes(), data)


if __name__ == '__main__': unittest.main()
