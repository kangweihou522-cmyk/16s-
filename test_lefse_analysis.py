"""Input/statistic-cache contracts; the simulated CLI is explicitly not LEfSe."""
import csv
import hashlib
import json
import pickle
import sys
import tempfile
import unittest
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
            features = {row[0].replace('|', '.'): list(map(float, row[1:])) for row in rows[1:]}
            for key, values in list(features.items()):
                for level in range(1, len(key.split('.'))):
                    features.setdefault('.'.join(key.split('.')[:level]), values)
            with Path(tail[1]).open('wb') as stream:
                pickle.dump({'feats': features, 'cls': {'class': rows[0][1:]}}, stream)
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
                for kind in ('bar', 'cladogram'):
                    for ext in ('png', 'svg'):
                        (Path(tail[1]) / (kind + '.' + ext)).write_text('SIMULATED no hits image')
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

    def test_fit_uses_official_cli_flags_with_single_normalization_and_preserves_wrapper_parameters(self):
        result = self.fit()
        commands = self.api.commands
        formatter = next(args for args in commands if Path(args[5]).name == 'lefse_format_input.py')
        model = next(args for args in commands if Path(args[5]).name == 'lefse_run.py')
        self.assertEqual(formatter[-6:], ['-c', '1', '-f', 'r', '-o', '-1'])
        self.assertNotIn('-u', formatter)
        self.assertNotIn('-s', formatter)
        self.assertEqual(model[model.index('-l') + 1], '0')
        self.assertEqual(model[model.index('-y') + 1], '0')
        self.assertEqual(json.loads((self.run / 'parameters.json').read_text()), self.original_parameters)
        self.assertTrue(all(key.count('.') >= 1 for key in self.api.fitted_features))
        self.assertTrue(all('f__' in key for key in self.api.fitted_features))
        self.assertEqual(result['lefse_run'], str(self.run))
        self.assertEqual(set(result['figures']), {'bar_png', 'bar_svg', 'cladogram_png', 'cladogram_svg'})
        self.assertEqual(len([args for args in commands if Path(args[5]).name.startswith('lefse_plot_')]), 4)
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

    def test_changed_cache_rejected_without_recomputation(self):
        self.fit()
        with (self.run / lefse.FIT_RESULT).open('a') as stream:
            stream.write('tampered\n')
        start = len(self.api.commands)
        with self.assertRaisesRegex(ValueError, '改动'):
            lefse.render_lefse(self.job, {'lda_threshold': 3}, self.run, self.api)
        self.assertEqual(len(self.api.commands), start)

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


if __name__ == '__main__':
    unittest.main()
