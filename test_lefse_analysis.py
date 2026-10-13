"""Input/statistic-cache contracts; the simulated CLI is explicitly not LEfSe."""
import csv
import hashlib
import importlib.util
import json
import os
import pickle
import re
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from types import SimpleNamespace

import lefse_analysis as lefse
from run_layout import stage_report


class SimulatedCLI:
    """Exercise generated data and prune script without claiming an actual fit."""
    def __init__(self):
        self.commands = []
        self.probes = []
        self.fitted_features = []

    def find_conda(self, config):
        return 'simulated-conda'

    def probe(self, job, args, timeout):
        self.probes.append(args)
        operations = json.loads(args[-1])
        runtime = {'prefix': '/simulated/env', 'python': sys.executable, 'version': 'SIMULATED',
                   'paths': {op: '/simulated/env/bin/' + names[0] for op, names in operations.items()}}
        return SimpleNamespace(returncode=0, stdout='LEFSE_RUNTIME=' + json.dumps(runtime) + '\n', stderr='')

    def begin_stage(self, job, folder, title, parameters=None):
        job.update(stage_folder=folder, stage_title=title, stage_parameters=parameters or {})
        stage_report(job['run'], folder, title, status='running', parameters=parameters)

    def finish_stage(self, job, results=None, notes=()):
        stage_report(job['run'], job['stage_folder'], job['stage_title'], parameters=job['stage_parameters'], results=results, notes=notes)

    def command(self, job, args, cwd):
        self.commands.append(list(args))
        executable, tail = Path(args[5]).name, args[6:]
        if executable == 'lefse_format_input.py':
            with Path(tail[0]).open(newline='') as stream:
                rows = list(csv.reader(stream, delimiter='\t'))
            class_row = int(tail[tail.index('-c') + 1]) - 1
            subject_row = int(tail[tail.index('-u') + 1]) - 1
            features = {row[0].replace('|', '.'): list(map(float, row[1:]))
                        for i, row in enumerate(rows) if i not in (class_row, subject_row)}
            for key, values in list(features.items()):
                for level in range(1, len(key.split('.'))):
                    features.setdefault('.'.join(key.split('.')[:level]), values)
            with Path(tail[1]).open('wb') as stream:
                pickle.dump({'feats': features, 'cls': {'class': rows[class_row][1:],
                                                     'subject': rows[subject_row][1:]}}, stream)
        elif executable == 'python':
            if Path(tail[0]).name == 'prune_formatted.py':
                # Execute the actual stdlib prune script against a newly generated test pickle.
                source = Path(tail[0]).read_text()
                previous = sys.argv
                try:
                    sys.argv = tail
                    exec(compile(source, tail[0], 'exec'), {'__name__': '__main__'})
                finally:
                    sys.argv = previous
            elif Path(tail[0]).name == 'render_empty.py':
                settings = json.loads(Path(tail[4]).read_text()) if len(tail) > 4 else {}
                formats = ('png', 'svg', 'pdf') if settings.get('general', {}).get('formats') == 'png_svg_pdf' else ('png', 'svg')
                for kind in ('bar', 'cladogram'):
                    for ext in formats:
                        (Path(tail[1]) / (kind + '.' + ext)).write_text('SIMULATED no hits image')
            elif Path(tail[0]).name == 'render_official.py':
                Path(tail[4]).write_text('SIMULATED official plot presentation adapter output')
        elif executable == 'lefse_run.py':
            with Path(tail[0]).open('rb') as stream:
                data = pickle.load(stream)
            self.fitted_features = sorted(data['feats'])
            with Path(tail[1]).open('w', newline='') as stream:
                writer = csv.writer(stream, delimiter='\t')
                for i, key in enumerate(self.fitted_features):
                    # Deterministic SIMULATED scores for filtering, not computed statistics.
                    writer.writerow([key, '5.0', data['cls']['class'][0], '3.5' if i % 2 else '1.5', '.01'])
        elif executable in ('lefse_plot_res.py', 'lefse_plot_cladogram.py'):
            Path(tail[1]).write_text('SIMULATED official plot output')
        else:
            raise AssertionError('Unexpected simulated command: ' + executable)


class LefseInputTests(unittest.TestCase):
    def setUp(self):
        self.samples = [{'SampleID': sid, 'Group': group} for sid, group in [('001', 'A'), ('002', 'A'), ('003', 'B'), ('004', 'B')]]

    def test_exact_sample_mapping_reorders_columns_and_retains_leading_zero(self):
        with tempfile.TemporaryDirectory() as temporary:
            table = Path(temporary) / 'counts.tsv'
            table.write_text('# Constructed from biom\n#OTU ID\t004\t001\t003\t002\nasv1\t40.0\t10\t30\t20\n')
            counts, totals = lefse.read_counts(table, self.samples)
            self.assertEqual(counts, {'asv1': [10, 20, 30, 40]})
            self.assertEqual(totals, [10, 20, 30, 40])

    def test_sample_mismatch_duplicate_asv_and_invalid_counts_rejected(self):
        invalid = [
            '#OTU ID\t001\t002\t003\t003\na\t1\t1\t1\t1\n',
            '#OTU ID\t1\t002\t003\t004\na\t1\t1\t1\t1\n',
            '#OTU ID\t001\t002\t003\t004\na\t1\t1\t1\t1\na\t2\t2\t2\t2\n',
        ]
        invalid += ['#OTU ID\t001\t002\t003\t004\na\t' + value + '\t1\t1\t1\n' for value in ('nan', 'inf', '-1', '1.5', 'oops', '0')]
        with tempfile.TemporaryDirectory() as temporary:
            table = Path(temporary) / 'table.tsv'
            for text in invalid:
                with self.subTest(text=text):
                    table.write_text(text)
                    with self.assertRaises(ValueError):
                        lefse.read_counts(table, self.samples)

    def test_taxonomy_missing_and_duplicate_asv_are_errors(self):
        with tempfile.TemporaryDirectory() as temporary:
            table = Path(temporary) / 'taxonomy.tsv'
            for text in ('Feature ID\tTaxon\na\tk__Bacteria\n',
                         'Feature ID\tTaxon\na\tk__Bacteria\na\tk__Bacteria\n'):
                table.write_text(text)
                with self.assertRaises(ValueError):
                    lefse.read_taxonomy(table, ['a', 'b'])

    def test_rank_formats_and_invalid_order(self):
        standard = lefse.parse_taxon('d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__S')
        qiime = lefse.parse_taxon('D_0__Bacteria; D_1__P; D_2__C; D_3__O; D_4__F; D_5__G; D_6__S')
        self.assertEqual(standard, qiime)
        self.assertEqual(lefse.parse_taxon('k__Bacteria; p__P; g__G'), ('Bacteria', 'P', '', '', '', 'G', ''))
        with self.assertRaisesRegex(ValueError, '顺序'):
            lefse.parse_taxon('k__Bacteria; g__G; f__F')

    def test_unknown_taxa_markers(self):
        for label in ('', 'Incertae_Sedis', 'incertae sedis', 'norank', 'uncultured bacterium', 'unclassified_Foo', 'unknown', 'unassigned', 'Others'):
            with self.subTest(label=label):
                self.assertTrue(lefse.unknown(label))
        self.assertFalse(lefse.unknown('Lactobacillus'))

    def test_cpm_uses_original_asv_total_once_and_parent_sums(self):
        counts = {'a': [10, 20, 30, 40], 'b': [30, 20, 10, 0], 'c': [60, 60, 60, 60]}
        taxonomy = {'a': lefse.parse_taxon('k__Bacteria; p__P; c__C; o__O; f__F; g__G; s__S1'),
                    'b': lefse.parse_taxon('k__Bacteria; p__P; c__C; o__O; f__F; g__G; s__S2'),
                    'c': lefse.parse_taxon('k__Bacteria; p__uncultured')}
        data = lefse.build_abundances(counts, taxonomy, self.samples, [100] * 4)
        genus = next(row for row in data['features'] if row['rank'] == 'L6')
        self.assertEqual(genus['counts'], [40] * 4)
        self.assertEqual(genus['cpm'], [400000] * 4)
        species = [row for row in data['features'] if row['rank'] == 'L7']
        self.assertEqual([sum(row['cpm'][i] for row in species) for i in range(4)], genus['cpm'])
        self.assertFalse(any('uncultured' in row['input_feature'] or 'Others' in row['input_feature'] for row in data['features']))

    def test_unknown_children_do_not_create_false_pooled_others_biomarker(self):
        counts = {'a': [90, 90, 1, 1], 'b': [1, 1, 90, 90]}
        taxonomy = {'a': lefse.parse_taxon('k__Bacteria; p__P; c__C; o__O; f__F1; g__uncultured'),
                    'b': lefse.parse_taxon('k__Bacteria; p__P; c__C; o__O; f__F2; g__norank')}
        data = lefse.build_abundances(counts, taxonomy, self.samples, [91] * 4)
        self.assertEqual(len([row for row in data['features'] if row['rank'] == 'L5']), 2)
        self.assertEqual(len([row for row in data['features'] if row['rank'] == 'L6']), 0)
        others = [row for row in data['others'] if row['rank'] == 'L6']
        self.assertEqual(len(others), 2)
        self.assertNotEqual(others[0]['parent_lineage'], others[1]['parent_lineage'])
        self.assertEqual(others[0]['counts'], [90, 90, 1, 1])
        self.assertEqual(others[1]['counts'], [1, 1, 90, 90])

    def test_same_genus_under_two_families_has_two_full_lineage_keys(self):
        taxonomy = {sid: lefse.parse_taxon('k__Bacteria; p__P; c__C; o__O; f__' + family + '; g__Same; s__Species')
                    for sid, family in [('a', 'Family1'), ('b', 'Family2')]}
        data = lefse.build_abundances({'a': [1] * 4, 'b': [2] * 4}, taxonomy, self.samples, [3] * 4)
        genera = [row for row in data['features'] if row['rank'] == 'L6']
        self.assertEqual(len(genera), 2)
        self.assertNotEqual(genera[0]['result_feature'], genera[1]['result_feature'])
        self.assertTrue(all(row['input_feature'].count('|') == 5 for row in genera))

    def test_name_collisions_are_detected_and_disambiguated_independent_of_order(self):
        labels = ['g__a-b', 'g__a b', 'g__a|b', 'g__a.b', 'g__中文', 'g__另一个']
        mapping, collisions = lefse.safe_names(labels)
        self.assertEqual(lefse.safe_names(list(reversed(labels))), (mapping, collisions))
        self.assertEqual(len(mapping), len(set(mapping.values())))
        self.assertIn('g__a-b', collisions)
        self.assertTrue(all(__import__('re').fullmatch('[A-Za-z][A-Za-z0-9_]*', safe) for safe in mapping.values()))

    def test_known_descendants_below_unknown_intermediate_rank_are_preserved(self):
        taxonomy = {'a': lefse.parse_taxon('k__Bacteria; p__P; c__Incertae_Sedis; o__O; f__F; g__Named')}
        data = lefse.build_abundances({'a': [1] * 4}, taxonomy, self.samples, [1] * 4)
        self.assertEqual([row['rank'] for row in data['features']], ['L2', 'L4', 'L5', 'L6'])
        self.assertEqual([row['rank'] for row in data['others']], ['L3', 'L7'])
        genus = next(row for row in data['features'] if row['rank'] == 'L6')
        self.assertIn('c__Others_at_L3_h', genus['input_feature'])
        self.assertTrue(genus['input_feature'].endswith('g__Named'))
        original, safe = next((original, safe) for original, safe in data['unknown_branches'].items() if len(original) == 3)
        self.assertEqual(original[-1], 'c__Incertae_Sedis')
        self.assertIn(safe, genus['input_feature'])

    def test_unknown_intermediate_contexts_never_pool_unrelated_named_descendants(self):
        taxonomy = {sid: lefse.parse_taxon('k__Bacteria; p__' + phylum + '; c__Incertae_Sedis; o__O; f__F; g__Named')
                    for sid, phylum in [('a', 'P1'), ('b', 'P2')]}
        data = lefse.build_abundances({'a': [9, 9, 1, 1], 'b': [1, 1, 9, 9]}, taxonomy, self.samples, [10] * 4)
        genera = [row for row in data['features'] if row['rank'] == 'L6']
        self.assertEqual(len(genera), 2)
        self.assertNotEqual(genera[0]['input_feature'].split('|')[2], genera[1]['input_feature'].split('|')[2])
        self.assertEqual({tuple(row['counts']) for row in genera}, {(9, 9, 1, 1), (1, 1, 9, 9)})
        self.assertFalse(any(row['rank'] == 'L3' for row in data['features']))

    def test_settings_defaults_and_numeric_validation(self):
        params = lefse.parameters({})
        self.assertEqual((params['lefse_env'], params['lda_threshold'], params['multiclass_strategy']), ('lefse', 2, 0))
        self.assertEqual(lefse.parameters({'multiclass_strategy': 'all-against-all'})['multiclass_strategy'], 1)
        for config in ({'lda_threshold': -1}, {'lda_threshold': 'nan'}, {'kw_alpha': 0}, {'wilcoxon_alpha': 1.1}, {'lefse_env': 'foo; echo bad'}, {'multiclass_strategy': 2}):
            with self.subTest(config=config), self.assertRaises(ValueError):
                lefse.parameters(config)


class LefseCacheTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.parent = self.root / 'preprocess'
        self.run = self.root / 'lefse'
        self.parent.mkdir()
        self.run.mkdir()
        table, taxonomy = self.parent / 'exported-table', self.parent / 'exported-taxonomy'
        table.mkdir()
        taxonomy.mkdir()
        table.joinpath('feature-table.tsv').write_text('#OTU ID\t004\t001\t003\t002\na\t90\t10\t80\t20\nb\t10\t90\t20\t80\n')
        taxonomy.joinpath('taxonomy.tsv').write_text('Feature ID\tTaxon\tConfidence\na\tk__Bacteria; p__P; c__C; o__O; f__F1; g__G1; s__S1\t.99\nb\tk__Bacteria; p__P; c__C; o__O; f__F2; g__G2; s__S2\t.99\n')
        self.saved = {'samples': [{'SampleID': sid, 'Group': group} for sid, group in [('001', 'A'), ('002', 'A'), ('003', 'B'), ('004', 'B')]]}
        self.original_parameters = {'extension': 'lefse', 'preprocess': str(self.parent), 'config': self.saved, 'extension_parameters': {'lda_threshold': 2}}
        (self.run / 'parameters.json').write_text(json.dumps(self.original_parameters))
        self.job = {'run': str(self.run), 'log': str(self.run / 'workflow.log')}
        self.api = SimulatedCLI()

    def fit(self):
        return lefse.run_lefse(self.job, {'lda_threshold': 2}, self.saved, self.parent, self.api)

    @unittest.skipUnless(importlib.util.find_spec('lefse') and importlib.util.find_spec('biom'),
                         'Official LEfSe formatter dependencies are not installed')
    def test_real_formatter_keeps_unique_sample_ids_groups_and_cpm_aligned(self):
        samples = self.saved['samples']
        counts, totals = lefse.read_counts(self.parent / 'exported-table/feature-table.tsv', samples)
        taxonomy = lefse.read_taxonomy(self.parent / 'exported-taxonomy/taxonomy.tsv', counts)
        data = lefse.build_abundances(counts, taxonomy, samples, totals)
        lefse._save_input(self.run, counts, totals, samples, data)
        output = self.run / '01_inputs/official_formatter.pkl'
        proc = subprocess.run([sys.executable, '-c',
                               'from lefse.lefse_format_input import format_input; format_input()',
                               str(self.run / '01_inputs/lefse_input_cpm.tsv'), str(output),
                               '-c', '1', '-u', '2', '-f', 'r', '-o', '-1'],
                              capture_output=True, text=True, timeout=30)
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        # Only deserialize this test's freshly generated official formatter output.
        with output.open('rb') as stream:
            formatted = pickle.load(stream)
        subjects = formatted['cls']['subject']
        self.assertEqual(set(subjects), {sample['SampleID'] for sample in samples})
        self.assertEqual(len(subjects), len(set(subjects)))
        expected_groups = {sample['SampleID']: data['groups'][sample['Group']] for sample in samples}
        self.assertEqual(dict(zip(subjects, formatted['cls']['class'])), expected_groups)
        self.assertNotIn('SampleID', formatted['feats'])
        for feature in data['features']:
            expected = dict(zip((sample['SampleID'] for sample in samples), feature['cpm']))
            self.assertEqual(formatted['feats'][feature['result_feature']], [expected[sid] for sid in subjects])

    def test_fit_uses_official_cli_flags_with_single_normalization_and_preserves_wrapper_parameters(self):
        result = self.fit()
        commands = self.api.commands
        formatter = next(args for args in commands if Path(args[5]).name == 'lefse_format_input.py')
        model = next(args for args in commands if Path(args[5]).name == 'lefse_run.py')
        self.assertEqual(formatter[-8:], ['-c', '1', '-u', '2', '-f', 'r', '-o', '-1'])
        self.assertNotIn('-s', formatter)
        with (self.run / '01_inputs/lefse_input_cpm.tsv').open(newline='') as stream:
            input_rows = list(csv.reader(stream, delimiter='\t'))
        self.assertEqual(input_rows[1], ['SampleID', '001', '002', '003', '004'])
        self.assertNotIn('SampleID', self.api.fitted_features)
        self.assertFalse(lefse.parameters({})['paired_design'])
        self.assertEqual(model[model.index('-l') + 1], '0')
        self.assertEqual(model[model.index('-y') + 1], '0')
        self.assertEqual(json.loads((self.run / 'parameters.json').read_text()), self.original_parameters)
        self.assertTrue(all(key.count('.') >= 1 for key in self.api.fitted_features))
        self.assertTrue(all('f__' in key for key in self.api.fitted_features))
        self.assertEqual(result['lefse_run'], str(self.run))
        self.assertEqual(set(result['figures']), {'bar_png', 'bar_svg', 'cladogram_png', 'cladogram_svg'})
        self.assertEqual(len([args for args in commands if Path(args[6]).name == 'render_official.py']), 4)
        self.assertTrue((self.run / 'lefse_complete.json').is_file())

    def test_replot_changes_threshold_without_model_or_formatter_and_preserves_fit(self):
        original = self.fit()
        cached = (self.run / lefse.FIT_RESULT).read_bytes()
        initial_figures = {key: Path(path).read_bytes() for key, path in original['figures'].items()}
        command_start = len(self.api.commands)
        output = self.run / 'renders' / 'threshold_1'
        output.mkdir(parents=True)
        job = {'run': str(output), 'log': str(output / 'workflow.log')}
        preview = lefse.render_lefse(job, {'lda_threshold': 1}, self.run, self.api)
        later = self.api.commands[command_start:]
        self.assertFalse(any(Path(args[5]).name in ('lefse_run.py', 'lefse_format_input.py') for args in later))
        self.assertEqual((self.run / lefse.FIT_RESULT).read_bytes(), cached)
        self.assertGreater(preview['biomarker_count'], original['biomarker_count'])
        self.assertEqual(preview['lefse_run'], str(self.run))
        self.assertFalse(preview['statistics_refit'])
        self.assertTrue(all(Path(path).is_relative_to(output) for path in preview['figures'].values()))
        self.assertEqual(initial_figures, {key: Path(path).read_bytes() for key, path in original['figures'].items()})

    def test_no_hit_renders_explicit_placeholder_without_running_official_empty_plot(self):
        self.fit()
        output = self.run / 'renders' / 'threshold_99'
        output.mkdir(parents=True)
        command_start = len(self.api.commands)
        result = lefse.render_lefse({'run': str(output)}, {'lda_threshold': 99}, self.run, self.api)
        self.assertTrue(result['no_hits'])
        self.assertEqual(result['biomarker_count'], 0)
        later = self.api.commands[command_start:]
        self.assertEqual(len(later), 1)
        self.assertEqual(Path(later[0][6]).name, 'render_empty.py')
        self.assertTrue(all(Path(path).stat().st_size for path in result['figures'].values()))

    def test_style_only_replot_uses_original_names_and_preserves_saved_statistics(self):
        self.fit()
        cached = (self.run / lefse.FIT_RESULT).read_bytes()
        marker = (self.run / 'lefse_complete.json').read_bytes()
        output = self.run / 'renders' / 'style_only'
        output.mkdir(parents=True)
        start = len(self.api.commands)
        style = {'group_order': ['B', 'A'], 'group_colors': {'A': '#112233', 'B': '#AABBCC'}}
        result = lefse.render_lefse({'run': str(output)}, {'lda_threshold': 2, **style}, self.run, self.api)
        self.assertEqual(result['group_order'], ['B', 'A'])
        self.assertEqual(result['group_colors'], style['group_colors'])
        self.assertEqual((self.run / lefse.FIT_RESULT).read_bytes(), cached)
        self.assertEqual((self.run / 'lefse_complete.json').read_bytes(), marker)
        later = self.api.commands[start:]
        self.assertEqual(len(later), 4)
        self.assertTrue(all(Path(args[6]).name == 'render_official.py' for args in later))
        self.assertTrue(all(args[args.index('--all_feats') + 1] == 'B:A' for args in later))
        plot_style = json.loads((output / '03_figures/plot_style.json').read_text())
        self.assertEqual(plot_style['groups'], {'A': 'A', 'B': 'B'})
        self.assertEqual(plot_style['group_colors'], style['group_colors'])
        inherited = self.run / 'renders' / 'inherited_style'
        inherited.mkdir(parents=True)
        original_style = json.loads((self.run / 'group_style.json').read_text())
        preview = lefse.render_lefse({'run': str(inherited)}, {'lda_threshold': 2}, self.run, self.api)
        self.assertEqual(preview['group_order'], original_style['group_order'])
        self.assertEqual(preview['group_colors'], original_style['group_colors'])

    def test_changed_cache_rejected_without_recomputation(self):
        self.fit()
        with (self.run / lefse.FIT_RESULT).open('a') as stream:
            stream.write('tampered\n')
        start = len(self.api.commands)
        with self.assertRaisesRegex(ValueError, '改动'):
            lefse.render_lefse(self.job, {'lda_threshold': 3}, self.run, self.api)
        self.assertEqual(len(self.api.commands), start)

    def test_plot_options_top_n_and_pdf_only_change_display_not_saved_model(self):
        self.fit()
        cached = (self.run / lefse.FIT_RESULT).read_bytes()
        marker = (self.run / 'lefse_complete.json').read_bytes()
        output = self.run / 'renders' / 'plot_settings'
        output.mkdir(parents=True)
        start = len(self.api.commands)
        config = {'lda_threshold': 1, 'plot_settings': {
            'general': {'dpi': 210, 'formats': 'png_svg_pdf', 'background': 'transparent', 'legend_position': 'bottom'},
            'lefse_bar': {'orientation': 'vertical', 'top_n': 1, 'label_levels': 3, 'feature_font_size': 12,
                          'class_font_size': 11, 'width': 9, 'height': 6, 'title': 'Custom LDA {threshold}'},
            'lefse_cladogram': {'max_level': 5, 'min_point_size': 2, 'max_point_size': 8, 'label_font_size': 10}}}
        result = lefse.render_lefse({'run': str(output)}, config, self.run, self.api)
        self.assertEqual((self.run / lefse.FIT_RESULT).read_bytes(), cached)
        self.assertEqual((self.run / 'lefse_complete.json').read_bytes(), marker)
        self.assertEqual(result['bar_biomarker_count'], 1)
        self.assertGreater(result['biomarker_count'], 1)
        self.assertEqual(len(result['figures']), 6)
        later = self.api.commands[start:]
        self.assertEqual(len(later), 6)
        self.assertTrue(all(Path(args[6]).name == 'render_official.py' for args in later))
        bars = [args for args in later if args[7] == 'bar']
        clades = [args for args in later if args[7] == 'cladogram']
        self.assertTrue(all(args[9].endswith('bar_display_results.res') for args in bars))
        self.assertTrue(all(args[9].endswith('threshold_results.res') for args in clades))
        bar_rows = lefse.read_results(output / '03_figures/bar_display_results.res')
        self.assertEqual(len(bar_rows), 1)
        self.assertEqual(abs(bar_rows[0]['lda_score']), max(abs(row['lda_score']) for row in lefse.read_results(self.run / lefse.FIT_RESULT)))
        self.assertTrue(all(args[args.index('--orientation') + 1] == 'v' for args in bars))
        self.assertTrue(all(args[args.index('--max_lev') + 1] == '5' for args in clades))
        self.assertTrue(all(args[args.index('--dpi') + 1] == '210' for args in later))
        self.assertTrue(all(args[args.index('--title') + 1] == 'Custom LDA 1.0' for args in bars))

    def test_no_hit_preserves_truthful_placeholder_with_optional_pdf(self):
        self.fit()
        output = self.run / 'renders' / 'no_hits_pdf'
        output.mkdir(parents=True)
        result = lefse.render_lefse({'run': str(output)}, {'lda_threshold': 99,
            'plot_settings': {'general': {'formats': 'png_svg_pdf', 'show_caption': False}}}, self.run, self.api)
        self.assertTrue(result['no_hits'])
        self.assertEqual(result['bar_biomarker_count'], 0)
        self.assertEqual(len(result['figures']), 6)
        self.assertIn('No biomarkers above LDA threshold', lefse.EMPTY_PLOT_SCRIPT)

    def test_fit_prunes_formatter_added_unknown_ancestors_but_fits_named_descendants(self):
        (self.parent / 'exported-taxonomy/taxonomy.tsv').write_text(
            'Feature ID\tTaxon\na\tk__Bacteria; p__P; c__Incertae_Sedis; o__O; f__F1; g__Named1\n'
            'b\tk__Bacteria; p__P; c__Incertae_Sedis; o__O; f__F2; g__Named2\n')
        self.fit()
        status = json.loads((self.run / 'lefse_complete.json').read_text())
        self.assertTrue(any(key.endswith('g__Named1') for key in self.api.fitted_features))
        self.assertTrue(any(key.endswith('g__Named2') for key in self.api.fitted_features))
        self.assertFalse(any('Others' in key.split('.')[-1] for key in self.api.fitted_features))
        self.assertTrue(any('Others_at_L3_h' in key.split('.')[-1] for key in status['removed_auto_ancestors']))
        mapping = (self.run / '01_inputs/name_mapping.tsv').read_text()
        self.assertIn('unknown_parent_context', mapping)
        self.assertIn('c__Incertae_Sedis', mapping)

    def test_constant_features_are_skipped_with_honest_status(self):
        (self.parent / 'exported-table/feature-table.tsv').write_text('#OTU ID\t001\t002\t003\t004\na\t50\t50\t50\t50\nb\t50\t50\t50\t50\n')
        result = self.fit()
        status = json.loads((self.run / 'lefse_complete.json').read_text())
        self.assertEqual(status['statistics'], 'skipped_no_variable_known_taxa')
        self.assertEqual(status['tested_features'], 0)
        self.assertTrue(result['no_hits'])
        self.assertFalse(any(Path(args[5]).name == 'lefse_run.py' for args in self.api.commands))

    def test_read_result_blank_lda_and_strict_threshold_boundary(self):
        result_file = self.root / 'scores.res'
        result_file.write_text('a\t3\tG_A\t2\t.01\nb\t2\t\t\t-\nc\t4\tG_B\t-3\t.02\n')
        results = lefse.read_results(result_file)
        rows, selected = lefse.filter_results(results, 2)
        self.assertEqual([row['feature'] for row in selected], ['c'])
        self.assertEqual(rows[0][2:4], ['', ''])
        self.assertEqual(rows[1][2:4], ['', ''])
        self.assertEqual(rows[2][2:4], ['G_B', '-3'])
        self.assertEqual(results[0]['raw'][2:4], ['G_A', '2'])


@unittest.skipUnless(os.environ.get('LEFSE_OFFICIAL_PLOT_ROOT'),
                     'Set LEFSE_OFFICIAL_PLOT_ROOT to the previously downloaded official LEfSe source')
class OfficialPlotStyleTests(unittest.TestCase):
    def test_actual_no_hit_images_keep_mandatory_status_when_caption_is_hidden(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            script = root / 'render_empty.py'
            script.write_text(lefse.EMPTY_PLOT_SCRIPT)
            settings = lefse.resolve_plot_settings({'plot_settings': {
                'general': {'dpi': 90, 'show_caption': False, 'formats': 'png_svg_pdf',
                            'background': 'transparent', 'font_family': 'DejaVu Serif'},
                'lefse_bar': {'width': 7, 'height': 4, 'title': 'No hits at {threshold}'}}})
            setting_file = root / 'settings.json'
            setting_file.write_text(json.dumps(settings))
            proc = subprocess.run([sys.executable, str(script), str(root), '99', 'executed', str(setting_file)],
                                  capture_output=True, text=True, timeout=30,
                                  env={**os.environ, 'MPLCONFIGDIR': str(root / 'matplotlib')})
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            for kind in ('bar', 'cladogram'):
                for ext in ('png', 'svg', 'pdf'):
                    self.assertGreater((root / (kind + '.' + ext)).stat().st_size, 0)
                svg = (root / (kind + '.svg')).read_text()
                self.assertIn('No biomarkers above LDA threshold 99.0', svg)
                self.assertNotIn('Saved analysis; no statistical refit', svg)
            self.assertIn('No hits at 99.0', (root / 'bar.svg').read_text())

    def test_custom_vertical_pdf_fonts_transparency_and_display_level(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            script = root / 'render_official.py'
            script.write_text(lefse.OFFICIAL_PLOT_SCRIPT)
            order = ['Z/group', 'A group', 'B group']
            settings = lefse.resolve_plot_settings({'plot_settings': {
                'general': {'dpi': 90, 'base_size': 11.5, 'title_size': 16, 'axis_text_size': 9, 'font_family': 'DejaVu Serif',
                            'legend_text_size': 12, 'background': 'transparent', 'theme': 'bw',
                            'legend_position': 'bottom', 'formats': 'png_svg_pdf', 'show_caption': False},
                'lefse_bar': {'orientation': 'vertical', 'width': 8, 'height': 5, 'feature_font_size': 12.5,
                              'title': 'Custom title {threshold}', 'x_label': 'Named taxa', 'y_label': 'Saved LDA'},
                'lefse_cladogram': {'width': 8, 'height': 5, 'max_level': 7, 'min_point_size': 1, 'max_point_size': 8}}})
            style = {'group_order': order, 'group_colors': {'Z/group': '#AA1122', 'A group': '#228844', 'B group': '#1144BB'},
                     'groups': {'Z/group': 'Z_group', 'A group': 'A_group', 'B group': 'B_group'}, 'plot_settings': settings}
            style_file = root / 'style.json'
            style_file.write_text(json.dumps(style))
            fixture = root / 'fixture.res'
            # Missing middle class tests that vertical colors do not shift; large
            # LDA vs log abundance tests that groups remain in requested order.
            fixture.write_text('k__Bacteria.p__P1.c__C1.o__O1.f__F1.g__G1\t1\tZ_group\t10\t.01\n'
                               'k__Bacteria.p__P3.c__C3.o__O3.f__F3.g__G3\t1\tB_group\t2.5\t.03\n')
            original = fixture.read_bytes()
            runner = """import os,sys,types,runpy
sys.path.insert(0,os.environ['LEFSE_OFFICIAL_PLOT_ROOT'])
rpy2=types.ModuleType('rpy2'); rpy2.__path__=[]
sys.modules['rpy2']=rpy2; sys.modules['rpy2.robjects']=types.ModuleType('rpy2.robjects')
sys.argv=sys.argv[1:]; runpy.run_path(sys.argv[0],run_name='__main__')
"""
            env = {**os.environ, 'MPLCONFIGDIR': str(root / 'matplotlib')}

            def render(kind, extension, output=None):
                output = output or root / (kind + '.' + extension)
                args = [sys.executable, '-c', runner, str(script), kind, str(style_file), str(fixture), str(output),
                        '--format', extension, *lefse._plot_options(kind, settings, 2, 'Z_group:A_group:B_group')]
                proc = subprocess.run(args, capture_output=True, text=True, env=env, timeout=30)
                self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
                self.assertGreater(output.stat().st_size, 0)
                return output

            for kind in ('bar', 'cladogram'):
                render(kind, 'png')
                self.assertTrue(render(kind, 'pdf').read_bytes().startswith(b'%PDF'))
                svg = render(kind, 'svg')
                texts = [(''.join(node.itertext()).strip(), node.attrib) for node in ET.parse(svg).iter()
                         if node.tag.endswith('}text')]
                self.assertEqual([label for label, _ in texts if label in order], ['Z/group', 'B group'])
                for color in ('#aa1122', '#1144bb'):
                    self.assertIn(color, svg.read_text().lower())
                audit = json.loads(Path(str(svg) + '.style.json').read_text())
                self.assertEqual(audit['plot_settings'], settings)
                if kind == 'bar':
                    self.assertIn('Custom title 2', [label for label, _ in texts])
                    self.assertIn('Named taxa', [label for label, _ in texts])
                    self.assertIn('Saved LDA', [label for label, _ in texts])
                    features = {label: float(attrs['x']) if 'x' in attrs else
                                float(re.search(r'translate\(([-+0-9.e]+)', attrs['transform'])[1])
                                for label, attrs in texts if label.startswith('f__')}
                    self.assertLess(features['f__F1.g__G1'], features['f__F3.g__G3'])
                    self.assertTrue(all('font-size: 12.5px' in attrs['style'] for label, attrs in texts if label.startswith('f__')))
                    self.assertTrue(all('font-size: 12px' in attrs['style'] for label, attrs in texts if label in order))
                    self.assertTrue(all('DejaVu Serif' in attrs['style'] for label, attrs in texts))
            # Hiding a legend still leaves valid official bar drawing functions.
            settings['general']['legend_position'] = 'none'
            style_file.write_text(json.dumps(style))
            hidden = render('bar', 'svg', root / 'hidden.svg')
            text = hidden.read_text()
            self.assertNotIn('Z/group', text)
            self.assertNotIn('B group', text)
            # No fabricated ancestor effects when a chosen display level has no rows.
            settings['lefse_cladogram']['max_level'] = 2
            style_file.write_text(json.dumps(style))
            empty = render('cladogram', 'svg', root / 'display_empty.svg')
            self.assertIn('No taxa available at the selected display level', empty.read_text())
            self.assertTrue(json.loads(Path(str(empty) + '.style.json').read_text())['empty_display_level'])
            self.assertEqual(fixture.read_bytes(), original)

    def test_official_renderers_follow_named_colors_original_labels_and_order(self):
        """Run real official plotting functions, with no statistical fit or R calls."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            script = root / 'render_official.py'
            script.write_text(lefse.OFFICIAL_PLOT_SCRIPT)
            order = ['Z/group', 'A group', 'B group']
            colors = {'Z/group': '#AA1122', 'A group': '#228844', 'B group': '#1144BB'}
            style = {'group_order': order, 'group_colors': colors,
                     'groups': {'Z/group': 'Z_group', 'A group': 'A_group', 'B group': 'B_group'}}
            (root / 'style.json').write_text(json.dumps(style))
            # Explicit rendering fixture, not scores claimed to come from a LEfSe fit.
            fixture = root / 'fixture.res'
            fixture.write_text('k__Bacteria.p__P1.c__C1.o__O1.f__F1.g__G1\t5\tZ_group\t3.1\t.01\n'
                               'k__Bacteria.p__P2.c__C2.o__O2.f__F2.g__G2\t4\tA_group\t4\t.02\n'
                               'k__Bacteria.p__P3.c__C3.o__O3.f__F3.g__G3\t3\tB_group\t2.5\t.03\n')
            original = fixture.read_bytes()
            # The official plot modules import lefse.lefse, whose unused rpy2 import
            # otherwise requires installing an R bridge. No R routine is provided.
            runner = """import os,sys,types,runpy
sys.path.insert(0,os.environ['LEFSE_OFFICIAL_PLOT_ROOT'])
rpy2=types.ModuleType('rpy2'); rpy2.__path__=[]
robjects=types.ModuleType('rpy2.robjects')
sys.modules['rpy2']=rpy2; sys.modules['rpy2.robjects']=robjects
sys.argv=sys.argv[1:]; runpy.run_path(sys.argv[0],run_name='__main__')
"""
            for kind in ('bar', 'cladogram'):
                for ext in ('png', 'svg'):
                    output = root / (kind + '.' + ext)
                    args = [sys.executable, '-c', runner, str(script), kind, str(root / 'style.json'),
                            str(fixture), str(output), '--format', ext, '--all_feats', 'Z_group:A_group:B_group']
                    if kind == 'bar':
                        args += ['--subclades', '2', '--max_feature_len', '10000']
                    env = {**os.environ, 'MPLCONFIGDIR': str(root / 'matplotlib')}
                    proc = subprocess.run(args, capture_output=True, text=True, env=env, timeout=30)
                    self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
                    self.assertGreater(output.stat().st_size, 0)
                    audit = json.loads(Path(str(output) + '.style.json').read_text())
                    self.assertEqual(audit['group_order'], order)
                    self.assertFalse(audit['statistics_refit'])
                svg = root / (kind + '.svg')
                tree = ET.parse(svg)
                texts = [(''.join(node.itertext()).strip(), node.attrib) for node in tree.iter()
                         if node.tag.endswith('}text')]
                self.assertEqual([label for label, _ in texts if label in order], order)
                for color in colors.values():
                    self.assertIn(color.lower(), svg.read_text().lower())
                if kind == 'bar':
                    features = {label: float(attrs['y']) for label, attrs in texts if label.startswith('f__')}
                    self.assertLess(features['f__F1.g__G1'], features['f__F2.g__G2'])
                    self.assertLess(features['f__F2.g__G2'], features['f__F3.g__G3'])
            self.assertEqual(fixture.read_bytes(), original)


if __name__ == '__main__':
    unittest.main()
