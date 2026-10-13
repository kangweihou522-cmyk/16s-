"""Group presentation never changes saved sample mapping or statistical reference."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import app
from differential_analysis import options, plot_options
from group_style import resolve_group_style


class GroupStyleTests(unittest.TestCase):
    samples = [{'SampleID': '001', 'Group': 'CON'}, {'SampleID': '002', 'Group': 'LDC'},
               {'SampleID': '003', 'Group': 'HDC'}, {'SampleID': '004', 'Group': 'CON'}]
    style = {'group_order': ['HDC', 'CON', 'LDC'],
             'group_colors': {'HDC': '#123456', 'CON': '#abcdef', 'LDC': '#654321'}}

    def test_named_style_order_and_case_do_not_change_sample_mapping_or_reference(self):
        before = json.dumps(self.samples)
        result = resolve_group_style(self.style, self.samples)
        self.assertEqual(result['group_order'], ['HDC', 'CON', 'LDC'])
        self.assertEqual(result['group_colors']['CON'], '#ABCDEF')
        self.assertEqual(json.dumps(self.samples), before)
        settings = options(dict(self.style, diff_method='ancombc2'), self.samples)
        self.assertEqual(settings['reference_group'], 'CON')
        self.assertEqual(settings['groups'], ['CON', 'LDC', 'HDC'])

    def test_inheritance_new_groups_and_inactive_groups_keep_named_colors(self):
        result = resolve_group_style({}, self.samples, fallback=self.style)
        self.assertEqual(result['group_order'], self.style['group_order'])
        config = {'group_order': ['unused', 'LDC'], 'group_colors': {'unused': '#000000', 'LDC': '#ffffff'}}
        result = resolve_group_style(config, self.samples)
        self.assertEqual(result['group_order'], ['LDC', 'CON', 'HDC'])
        self.assertEqual(result['group_colors'], {'LDC': '#FFFFFF', 'CON': '#D55E00', 'HDC': '#CC79A7'})

    def test_invalid_colors_or_duplicate_order_rejected(self):
        for config in ({'group_order': ['CON', 'CON']}, {'group_order': 'CON,HDC'},
                       {'group_colors': {'CON': 'red'}}, {'group_colors': []},
                       {'group_colors': {'CON': '#123456; command'}}):
            with self.subTest(config=config), self.assertRaises(ValueError):
                resolve_group_style(config, self.samples)

    def setUp(self):
        temp = tempfile.TemporaryDirectory(); self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.parent = self.root / 'parent'; self.parent.mkdir(); app.marker(self.parent)
        self.saved = dict(self.style, samples=self.samples, font='DejaVu Serif', env='fixture')
        (self.parent / 'parameters.json').write_text(json.dumps(self.saved))
        (self.parent / 'preprocess_complete.json').write_text('{}')

    def test_csv_and_historical_group_loading_preserve_literal_names(self):
        filename = self.root / 'groups.csv'
        filename.write_text('SampleID,Group\n001,组 B\n002,CON\n003,组 B\n', encoding='utf8')
        result = app.group_style_input(metadata=str(filename))
        self.assertEqual(result['group_order'], ['组 B', 'CON'])
        self.assertEqual(result['sample_counts'], {'组 B': 2, 'CON': 1})
        result = app.group_style_input(run=str(self.parent))
        self.assertEqual(result['group_order'], self.style['group_order'])

    def test_extension_inherits_style_and_preserves_parent(self):
        old = (self.parent / 'parameters.json').read_bytes()
        job = {'status': 'running'}
        def model(job, config, saved, parent, api):
            self.assertEqual(config['group_order'], self.style['group_order'])
            self.assertEqual(saved['group_colors']['CON'], '#ABCDEF')
            self.assertTrue(Path(job['run'], 'group_style.json').is_file())
            return {'test_only': True}
        with patch.object(app, 'totals', return_value={s['SampleID']: 100 for s in self.samples}), \
             patch('differential_analysis.run_differential', side_effect=model):
            app.worker(job, 'differential', {'run': str(self.parent), 'diff_method': 'ancombc2'})
        self.assertEqual(job['status'], 'done', job.get('error'))
        self.assertEqual((self.parent / 'parameters.json').read_bytes(), old)

    def test_cached_diversity_replot_keeps_source_and_does_not_run_qiime(self):
        analysis = self.parent / 'diversity'; analysis.mkdir(); app.marker(analysis)
        (analysis / 'parameters.json').write_text(json.dumps({'preprocess': str(self.parent), 'config': self.saved}))
        (analysis / 'qiime_statistics_complete.json').write_text('{}')
        (analysis / 'saved-vector.tsv').write_text('sample-id\tshannon\n001\t3.5\n')
        original = {f: f.read_bytes() for f in analysis.iterdir() if f.is_file()}
        style = {'group_order': ['LDC', 'HDC', 'CON'], 'group_colors': self.style['group_colors']}
        calls = []
        def command(job, args, cwd, **kwargs):
            calls.append(args)
            self.assertEqual(args[2], str(analysis))
            actual = json.loads(Path(args[-3]).read_text())
            self.assertEqual(actual['group_order'], style['group_order'])
            plotting = json.loads(Path(args[-1]).read_text())
            self.assertEqual(plotting['alpha']['plot_type'], 'violin')
            self.assertEqual(plotting['general']['dpi'], 144)
            output = Path(args[-2]); output.mkdir(parents=True)
            (output / 'render_summary.json').write_text(json.dumps({'fixture': True}))
        job = {'status': 'running'}
        with patch.object(app, 'q') as qiime, patch.object(app, 'command', side_effect=command), \
             patch.object(app, 'r_command', side_effect=lambda config, *args: ['Rscript', *args]), \
             patch.object(app, 'r_execution_env', return_value={}):
            app.worker(job, 'diversity-render', dict(style, diversity_run=str(analysis),
                       plot_settings={'alpha': {'plot_type': 'violin'}, 'general': {'dpi': 144}}))
        self.assertEqual(job['status'], 'done', job.get('error'))
        qiime.assert_not_called(); self.assertEqual(len(calls), 1)
        for file, data in original.items(): self.assertEqual(file.read_bytes(), data)
        self.assertNotEqual(Path(job['run']), analysis)
        self.assertFalse(job['result']['statistics_refit'])

    def test_differential_plot_options_are_separate_from_saved_statistics(self):
        config = {'diff_comparisons': 'all_pairs', 'heatmap_low': '#abcdef',
                  'bubble_high': '#123456', 'diff_plot_top_n': '18'}
        settings = plot_options(config)
        self.assertEqual(settings['heatmap_low'], '#ABCDEF')
        self.assertEqual(settings['diff_plot_top_n'], 18)
        self.assertEqual(settings['diff_heatmap_scale'], 'relative_abundance')
        self.assertNotIn('alpha', settings)
        for bad in ({'heatmap_low': 'red'}, {'bubble_high': '#12345'},
                    {'diff_plot_top_n': 0}, {'diff_plot_top_n': 3.5},
                    {'bubble_color_mode': 'q_value'}):
            with self.subTest(config=bad), self.assertRaises(ValueError): plot_options(bad)

    def differential_fixture(self):
        analysis = self.parent / 'differential'; analysis.mkdir(); app.marker(analysis)
        counts = analysis / '01_counts' / 'counts.tsv'; counts.parent.mkdir()
        counts.write_text('FeatureID\t001\t002\t003\t004\na\t10\t20\t30\t40\n')
        meta = analysis / '01_counts' / 'metadata.tsv'
        meta.write_text('sample-id\tGroup\n001\tCON\n002\tLDC\n003\tHDC\n004\tCON\n')
        cached = analysis / '02_differential' / 'L2_phylum'; cached.mkdir(parents=True)
        (cached / 'ancombc2_fit.rds').write_bytes(b'BACKEND FIXTURE ONLY - not read by real R')
        (analysis / 'parameters.json').write_text(json.dumps({'preprocess': str(self.parent), 'config': self.saved}))
        import hashlib
        opts = {'method': 'ancombc2', 'alpha': .0123, 'reference_group': 'CON', 'group_style': self.style,
                'metadata_tsv': str(meta), 'metadata_sha256': hashlib.sha256(meta.read_bytes()).hexdigest(),
                'tables': [{'counts_tsv': str(counts), 'counts_sha256': hashlib.sha256(counts.read_bytes()).hexdigest(), 'level': 2, 'name': 'phylum'}]}
        (analysis / 'differential_options.json').write_text(json.dumps(opts))
        return analysis, counts

    def test_cached_differential_replot_keeps_model_statistics_and_input(self):
        analysis, counts = self.differential_fixture()
        original = {f: f.read_bytes() for f in analysis.rglob('*') if f.is_file()}
        calls = []
        def command(job, args, cwd, **kwargs):
            calls.append(args)
            self.assertIn('--plot-ancom', args); self.assertNotIn('--run', args)
            opts = json.loads(Path(args[-2]).read_text())
            self.assertEqual(opts['alpha'], .0123)
            self.assertEqual(opts['reference_group'], 'CON')
            self.assertEqual(opts['bubble_high'], '#123456')
            self.assertEqual(opts['group_style']['group_order'], ['LDC', 'HDC', 'CON'])
            self.assertEqual(opts['plot_settings']['heatmap']['width'], 11)
            self.assertEqual(opts['plot_settings']['bubble']['max_size'], 8)
            self.assertEqual(json.loads(Path(cwd, 'plot_settings.json').read_text()), opts['plot_settings'])
            directory = Path(cwd) / '02_differential'; directory.mkdir(exist_ok=True)
            (directory / 'summary.json').write_text(json.dumps({'statistics_refit': False, 'fixture_only': True}))
            output = directory / 'L2_phylum' / 'abundance'; output.mkdir(parents=True)
            for kind in ('heatmap', 'bubble'):
                for ext in ('png', 'svg'): (output / ('abundance_' + kind + '.' + ext)).write_text('SIMULATED IMAGE')
        job = {'status': 'running'}
        with patch.object(app, 'q') as qiime, patch('differential_analysis.run_differential') as model, \
             patch.object(app, 'command', side_effect=command), \
             patch.object(app, 'r_command', side_effect=lambda config, *args: ['Rscript', *args]), \
             patch.object(app, 'r_execution_env', return_value={}):
            app.worker(job, 'differential-render', {'differential_run': str(analysis), 'group_order': ['LDC', 'HDC', 'CON'],
                                                 'bubble_high': '#123456', 'diff_alpha': '.99',
                                                 'plot_settings': {'heatmap': {'width': 11}, 'bubble': {'max_size': 8}}})
        self.assertEqual(job['status'], 'done', job.get('error'))
        qiime.assert_not_called(); model.assert_not_called(); self.assertEqual(len(calls), 1)
        for filename, data in original.items(): self.assertEqual(filename.read_bytes(), data)
        self.assertEqual(len(job['result']['differential_figures']), 1)
        self.assertEqual(job['result']['differential_run'], str(analysis))

    def test_cached_differential_modified_counts_rejected_before_rendering(self):
        analysis, counts = self.differential_fixture()
        counts.write_text(counts.read_text() + '\nchanged\n')
        job = {'status': 'running'}
        with patch.object(app, 'command') as command:
            app.worker(job, 'differential-render', {'differential_run': str(analysis)})
        self.assertEqual(job['status'], 'failed')
        self.assertIn('计数已改变', job['error']); command.assert_not_called()


if __name__ == '__main__': unittest.main()
