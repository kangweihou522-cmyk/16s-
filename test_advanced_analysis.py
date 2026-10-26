"""Check metadata identity, full denominators and cache-only orchestration."""
import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import app
from advanced_analysis import options_for, prepare_metadata
from figure_editor import describe_figure
from plot_settings import resolve_plot_settings


class AdvancedAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.parent = self.root / 'parent'; self.parent.mkdir(); app.marker(self.parent)
        self.samples = [{'SampleID': sid, 'Group': group} for sid, group in
                        [('001', 'CON'), ('NA', 'CON'), ('010,case', 'LDC'), ('x4', 'LDC')]]
        self.saved = {'samples': self.samples, 'env': 'qiime-fixture', 'font': 'DejaVu Serif',
                      'group_order': ['CON', 'LDC'], 'group_colors': {'CON': '#AA0000', 'LDC': '#0000AA'},
                      'plot_settings': resolve_plot_settings({})}
        (self.parent / 'parameters.json').write_text(json.dumps(self.saved))
        (self.parent / 'preprocess_complete.json').write_text('{}')
        abundance = self.parent / '10_abundance'; (abundance / 'L6_genus').mkdir(parents=True)
        self.source = abundance / 'L6_genus/counts.tsv'
        self.source.write_text('FeatureID\t001\tNA\t010,case\tx4\n'
                               'd__Bacteria;p__P;g__Good\t10\t20\t30\t40\n'
                               'd__Bacteria;p__P;g__unclassified\t90\t80\t70\t60\n')
        (abundance / 'manifest.json').write_text(json.dumps({'levels': [
            {'level': 6, 'counts': {'tsv': '10_abundance/L6_genus/counts.tsv'}}]}))
        table = self.parent / 'exported-table'; table.mkdir()
        (table / 'feature-table.tsv').write_text(self.source.read_text().replace('FeatureID', '#OTU ID'))

    def new_run(self, name='advanced'):
        run = self.parent / name; run.mkdir(); app.marker(run)
        return run

    def test_literal_sample_join_original_denominator_and_uncertain_report(self):
        auxiliary = self.root / 'traits.csv'
        with auxiliary.open('w', newline='') as stream:
            writer = csv.writer(stream); writer.writerow(['SampleID', 'IL6', 'AnimalID'])
            writer.writerows([['x4', 40, 'animal4'], ['010,case', 30, 'shared'],
                              ['NA', 20, 'shared'], ['001', 10, 'animal1'], ['extra', 99, 'other']])
        before = {p: p.read_bytes() for p in self.parent.rglob('*') if p.is_file()}
        run = self.new_run()
        options = options_for({'advanced_metadata': str(auxiliary)}, self.saved, self.parent, run, app.module_api())
        with Path(options['metadata_path']).open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual([r['SampleID'] for r in rows], ['001', 'NA', '010,case', 'x4'])
        self.assertEqual([r['IL6'] for r in rows], ['10', '20', '30', '40'])
        self.assertTrue(options['repeated_subjects'])
        self.assertEqual(options['association_transform'], 'clr')
        self.assertEqual(options['core_min_relative'], .0001)
        self.assertIn('unclassified', Path(options['total_counts_path']).read_text())
        self.assertNotIn('unclassified', Path(options['counts_path']).read_text())
        report = json.loads((run / '00_inputs/input_report.json').read_text())
        self.assertEqual(report['ignored_samples'], ['extra'])
        self.assertEqual(report['excluded_read_fraction']['001'], .9)
        for filename, raw in before.items():
            self.assertEqual(filename.read_bytes(), raw)

    def test_missing_duplicate_or_conflicting_metadata_is_rejected(self):
        filename = self.root / 'traits.csv'; destination = self.root / 'joined.csv'
        for content in ['SampleID,IL6\n001,1\n',
                        'SampleID,IL6\n001,1\n001,2\n',
                        'SampleID,Group\n001,LDC\nNA,CON\n"010,case",LDC\nx4,LDC\n']:
            filename.write_text(content)
            with self.subTest(content=content), self.assertRaises(ValueError):
                prepare_metadata(filename, self.samples, [s['SampleID'] for s in self.samples], destination)

    def test_invalid_thresholds_modules_and_manifest_paths_fail_before_external_commands(self):
        for index, config in enumerate([{'advanced_level': 9}, {'advanced_modules': ['magic']},
                                        {'advanced_core_prevalence': 101}, {'advanced_pseudocount': 0},
                                        {'advanced_modules': ['association']},
                                        {'advanced_predictors': 'Group', 'advanced_covariates': 'Group'}]):
            with self.subTest(config=config), self.assertRaises(ValueError):
                options_for(config, self.saved, self.parent, self.new_run('bad' + str(index)), app.module_api())
        manifest = self.parent / '10_abundance/manifest.json'
        manifest.write_text(json.dumps({'levels': [{'level': 6, 'counts': {'tsv': '../../elsewhere.tsv'}}]}))
        with self.assertRaisesRegex(ValueError, '越出'):
            options_for({}, self.saved, self.parent, self.new_run('escape'), app.module_api())

    def prepare_cache(self):
        analysis = self.new_run('cache')
        (analysis / '02_statistics').mkdir(); (analysis / '03_figures').mkdir()
        cache = analysis / '02_statistics/advanced_cache.rds'; cache.write_bytes(b'synthetic serialized cache')
        (analysis / 'advanced_complete.json').write_text(json.dumps({'cache_sha256': hashlib.sha256(cache.read_bytes()).hexdigest()}))
        options = {**self.saved, 'r_env': 'advanced-fixture', 'top_n': 20}
        (analysis / 'advanced_options.json').write_text(json.dumps(options))
        (analysis / 'parameters.json').write_text(json.dumps({'extension': 'advanced', 'preprocess': str(self.parent), 'config': self.saved}))
        png = analysis / '03_figures/core_heatmap.png'; png.write_bytes(b'synthetic png')
        summary = {'source': 'advanced_analysis', 'advanced_run': str(analysis), 'font': 'DejaVu Serif',
                   'top_n': 20, 'plot_settings': self.saved['plot_settings'],
                   'figures': [{'stem': '03_figures/core_heatmap', 'profile': 'heatmap'}]}
        (png.parent / 'render_summary.json').write_text(json.dumps(summary))
        return analysis, png, summary

    def test_figure_preview_uses_only_cache_and_preserves_source(self):
        analysis, png, summary = self.prepare_cache()
        self.assertEqual(describe_figure(png, app.module_api())['kind'], 'advanced')
        before = {p: p.read_bytes() for p in analysis.rglob('*') if p.is_file()}
        calls = []
        def render(job, args, cwd=None, **kwargs):
            self.assertIn('--render', args); self.assertNotIn('--run', args)
            options = json.loads(Path(args[-2]).read_text())
            self.assertEqual(options['cache_run'], str(analysis))
            self.assertEqual(options['top_n'], 7)
            self.assertEqual(options['heatmap_low'], '#123456')
            calls.append(args)
            out = Path(cwd) / '03_figures'; out.mkdir(exist_ok=True)
            (out / 'core_heatmap.png').write_bytes(b'preview png')
            (out / 'core_heatmap.svg').write_text('<svg/>')
            (out / 'render_summary.json').write_text(json.dumps(dict(summary, top_n=7,
                plot_settings=options['plot_settings'], visual_options={'heatmap_low': '#123456'})))
        job = {'status': 'running'}
        with patch.object(app, 'q') as qiime, patch.object(app, 'command', side_effect=render), \
                patch.object(app, 'r_command', side_effect=lambda c, *args: ['Rscript', *args]), \
                patch.object(app, 'r_execution_env', return_value={}):
            app.worker(job, 'figure-render', {'image': str(png), 'diff_plot_top_n': 7,
                                             'heatmap_low': '#123456', 'plot_settings': {'heatmap': {'width': 9}}})
        self.assertEqual(job['status'], 'done', job.get('error'))
        self.assertEqual(len(calls), 1); qiime.assert_not_called()
        self.assertEqual(job['result']['figure']['visual_options']['diff_plot_top_n'], 7)
        for path, raw in before.items():
            self.assertEqual(path.read_bytes(), raw)

    def test_changed_cache_is_refused(self):
        analysis, png, _ = self.prepare_cache()
        (analysis / '02_statistics/advanced_cache.rds').write_bytes(b'changed')
        with patch.object(app, 'command') as command:
            job = {'status': 'running'}
            app.worker(job, 'figure-render', {'image': str(png)})
        self.assertEqual(job['status'], 'failed'); command.assert_not_called()
        self.assertIn('缓存已改变', job['error'])

    def test_complete_worker_reuses_stage_folder_and_finishes_history_marker(self):
        before = {p: p.read_bytes() for p in self.parent.rglob('*') if p.is_file()}
        calls = []
        def external(job, args, cwd=None, **kwargs):
            calls.append(args)
            run = Path(job['run'])
            self.assertTrue((run/'00_inputs/report.json').is_file())
            if '--run' in args:
                folder = run/'02_statistics'; folder.mkdir()
                (folder/'advanced_cache.rds').write_bytes(b'synthetic cache')
                (folder/'analysis_summary.json').write_text(json.dumps({'diagnostics': [
                    {'Module': 'core', 'Status': 'complete', 'Reason': 'synthetic test only'}]}))
        job = {'status': 'running'}
        with patch.object(app, 'command', side_effect=external), patch.object(app, 'q') as qiime, \
                patch.object(app, 'r_command', side_effect=lambda c, *args: ['Rscript', *args]), \
                patch.object(app, 'r_execution_env', return_value={}):
            app.worker(job, 'advanced', {'run': str(self.parent)})
        self.assertEqual(job['status'], 'done', job.get('error'))
        run = Path(job['result']['advanced_run'])
        self.assertTrue((run/'advanced_complete.json').is_file())
        self.assertEqual(len(calls), 2)
        self.assertIn('--check', calls[0]); self.assertIn('--run', calls[1])
        qiime.assert_not_called()
        from batch_context import restore_batch
        self.assertEqual(restore_batch(self.parent, app.module_api())['fields']['advanced_run'], str(run))
        for filename, raw in before.items():
            self.assertEqual(filename.read_bytes(), raw)


if __name__ == '__main__':
    unittest.main()
