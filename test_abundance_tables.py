"""Synthetic CLI outputs exercise abundance orchestration, not QIIME analysis."""
import csv
import copy
import json
import tempfile
import unittest
from pathlib import Path

from abundance_tables import (LEVELS, clean_taxonomy_table, export_clean_taxonomy,
                              read_table, read_taxonomy, run_abundance,
                              table_totals, validate_relative, write_table_pair)
from run_layout import artifact_path, initialize_layout, stage_report


class FakeCLI:
    artifact_path = staticmethod(artifact_path)

    def __init__(self, source, *, bad_relative=False, lost_counts=False):
        self.calls = []
        self.tables = {str(artifact_path(source, 'filtered-table.qza')):
                       (['001', '010,case'], ['ASV1', 'ASV2'], [[2, 6], [6, 2]])}
        self.bad_relative = bad_relative
        self.lost_counts = lost_counts

    def conda_cmd(self, saved):
        return ['mock-conda', 'run', '-n', saved['env']]

    def q(self, job, saved, *args):
        self.calls.append(args)
        if args[:2] == ('tools', 'export'):
            source = args[args.index('--input-path') + 1]
            destination = Path(args[args.index('--output-path') + 1])
            assert destination.is_absolute(), 'New export path must avoid app logical mapping'
            assert not destination.exists(), 'QIIME creates its own export directory'
            assert destination.parent.is_dir(), 'Create parent before calling CLI'
            destination.mkdir()
            if source.endswith('taxa-barplot.qzv'):
                (destination / 'level-2.csv').write_text('index,k__Bacteria;p__Firmicutes\n001,8\n')
            elif source.endswith('taxonomy.qza'):
                (destination / 'taxonomy.tsv').write_text('Feature ID\tTaxon\nASV1\tk__Bacteria;p__P1;c__C;o__O;f__F;g__G;s__S\nASV2\tk__Bacteria;p__P2;c__C;o__O;f__F;g__G;s__S\n')
            else:
                (destination / 'feature-table.biom').write_text(json.dumps(self.tables[source]))
            return
        output_index = next(index for index, value in enumerate(args) if value.startswith('--o-'))
        output = Path(args[output_index + 1])
        assert output.is_absolute(), 'New artifact path must avoid app logical mapping'
        assert output.parent.is_dir()
        output.write_bytes(b'synthetic artifact, not real QIIME')
        if args[:2] == ('taxa', 'collapse'):
            table = self.tables[args[args.index('--i-table') + 1]]
            level = int(args[args.index('--p-level') + 1])
            taxon = ';'.join(['k__Bacteria', 'p__Firmicutes', 'c__C', 'o__O', 'f__F', 'g__G', 's__S'][:level])
            self.tables[str(output)] = (table[0], [taxon], [[7 if self.lost_counts else 8, 8]])
        elif args[:2] == ('feature-table', 'relative-frequency'):
            samples, features, values = self.tables[args[args.index('--i-table') + 1]]
            totals = [sum(row[index] for row in values) for index in range(len(samples))]
            numbers = [[value / totals[index] for index, value in enumerate(row)] for row in values]
            if self.bad_relative and len(numbers) == 2:
                numbers = [[0.5, 0.5], [0.5, 0.5]]
            self.tables[str(output)] = (samples, features, numbers)

    def command(self, job, args, cwd=None):
        assert 'biom' in args
        source = Path(args[args.index('-i') + 1])
        target = Path(args[args.index('-o') + 1])
        samples, features, numbers = json.loads(source.read_text())
        with target.open('w', newline='') as stream:
            stream.write('# Constructed from biom file\n')
            writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
            writer.writerow(['#OTU ID', *samples])
            writer.writerows([feature, *row] for feature, row in zip(features, numbers))

    def begin_stage(self, job, folder, title, parameters=None):
        job['stage_folder'] = folder
        job['stage_title'] = title
        job['stage_parameters'] = parameters
        stage_report(job['run'], folder, title, status='running', parameters=parameters)

    def finish_stage(self, job, results=None, notes=()):
        stage_report(job['run'], job.pop('stage_folder'), job.pop('stage_title'),
                     parameters=job.pop('stage_parameters'), results=results, notes=notes)


class AbundanceTablesTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.saved = {'env': 'mock-qiime', 'samples': [{'SampleID': '001'}, {'SampleID': '010,case'}]}

    def source(self, *, layout=True):
        source = self.root / ('parent-layout' if layout else 'parent-legacy')
        source.mkdir()
        if layout:
            initialize_layout(source, 'preprocess')
        for name in ('filtered-table.qza', 'taxonomy.qza', 'taxa-barplot.qzv'):
            target = artifact_path(source, name)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(b'original artifact')
        target = artifact_path(source, 'exported-table/feature-table.tsv')
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text('# Constructed from biom file\n#OTU ID\t001\t010,case\nASV1\t2.0\t6\nASV2\t6\t2\n')
        target = artifact_path(source, 'exported-taxonomy/taxonomy.tsv')
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text('Feature ID\tTaxon\tConfidence\nASV2\tk__Bacteria; p__Other; c__C; o__O; f__F; g__G; s__S\t0.99\nASV1\tk__Bacteria; p__Firmicutes, group; c__C; o__O; f__F; g__G; s__S\t0.95\nremoved\tk__Archaea; p__Other\t0.8\n')
        return source

    def test_annotation_csv_joins_by_feature_id_and_preserves_sample_ids(self):
        source = self.source()
        table = read_table(artifact_path(source, 'exported-table/feature-table.tsv'))
        annotations = read_taxonomy(artifact_path(source, 'exported-taxonomy/taxonomy.tsv'))
        pair = write_table_pair(self.root / 'out', 'annotated', table, annotations=annotations)
        with pair['csv'].open(newline='') as stream:
            rows = list(csv.reader(stream))
        self.assertEqual(rows[0], ['FeatureID', '001', '010,case', 'Taxon', 'Confidence'])
        self.assertEqual(rows[1][:3], ['ASV1', '2', '6'])
        self.assertIn('p__Firmicutes, group', rows[1][3])
        self.assertIn('p__Other', rows[2][3])
        self.assertIn('"010,case"', pair['csv'].read_text())
        self.assertIn('"k__Bacteria; p__Firmicutes, group;', pair['csv'].read_text())

    def test_numeric_values_and_duplicate_ids_are_strict(self):
        filename = self.root / 'table.tsv'
        for value in ('NaN', 'Infinity', '-1', 'n/a', '0.5'):
            with self.subTest(value=value):
                filename.write_text('#OTU ID\t001\nASV1\t' + value + '\n')
                with self.assertRaises(ValueError):
                    read_table(filename)
        for text in ('#OTU ID\t001\t001\nASV1\t1\t2\n',
                     '#OTU ID\t001\nASV1\t1\nASV1\t2\n',
                     '#OTU ID\t\nASV1\t1\n', '#OTU ID\t001\n\t1\n'):
            filename.write_text(text)
            with self.assertRaises(ValueError):
                read_table(filename)
        filename.write_text('#OTU ID\t001\nASV1\t1.01\n')
        with self.assertRaises(ValueError):
            read_table(filename, relative=True)

    def test_annotation_rejects_missing_and_duplicate_feature_ids(self):
        filename = self.root / 'taxonomy.tsv'
        filename.write_text('Feature ID\tTaxon\nASV1\tk__Bacteria\nASV1\tk__Archaea\n')
        with self.assertRaisesRegex(ValueError, '重复'):
            read_taxonomy(filename)
        filename.write_text('Feature ID\tTaxon\nASV2\tk__Bacteria\n')
        with self.assertRaisesRegex(ValueError, '缺少分类'):
            write_table_pair(self.root / 'out', 'annotated', (['001'], ['ASV1'], [[1]]),
                             annotations=read_taxonomy(filename))

    def test_cleaning_aggregates_labels_conserves_counts_and_preserves_input(self):
        counts = (['001', '010,case', '007'],
                  ['k__Bacteria;p__P;f__F1;g__Shared;s__Species',
                   'k__Bacteria;p__P;f__F2;g__Shared;s__Species2',
                   'k__Bacteria;p__P;f__F;g__uncultured',
                   'k__Bacteria;p__P;f__F;g__Incertae_Sedis',
                   'k__Bacteria;p__P;f__F;s__NamedSpecies',
                   'k__Bacteria;p__P;f__F;g__  ',
                   'k__Bacteria;p__P;f__F;g__Known;s__unknown',
                   'D_0__Bacteria;D_1__P;D_4__F3;D_5__Shared;D_6__uncultured'],
                  [[1, 2, 0], [3, 4, 0], [5, 6, 0], [7, 8, 0],
                   [9, 10, 0], [11, 12, 0], [13, 14, 0], [15, 16, 0]])
        before = copy.deepcopy(counts)
        cleaned, relative, mapping, report = clean_taxonomy_table(counts, 6)
        self.assertEqual(counts, before)
        self.assertEqual(cleaned[0], ['001', '010,case', '007'])
        self.assertEqual(cleaned[1], ['Shared', 'Known', 'Others'])
        self.assertEqual(cleaned[2], [[19, 22, 0], [13, 14, 0], [32, 36, 0]])
        self.assertEqual(table_totals(cleaned), table_totals(counts))
        for sid, total in table_totals(relative).items():
            self.assertAlmostEqual(float(total), 0 if sid == '007' else 1)
        self.assertAlmostEqual(float(relative[2][0][0]), 19 / 64)
        self.assertTrue(all(row[2] == 0 for row in relative[2]))
        self.assertEqual(report['zero_count_samples'], ['007'])
        self.assertEqual(report['others_per_sample_counts'], {'001': 32, '010,case': 36, '007': 0})
        self.assertEqual(report['merged_taxa'], {'Shared': 3, 'Others': 4})
        self.assertEqual([row['original_lineage'] for row in mapping], counts[1])
        self.assertEqual([row['reason'] for row in mapping[4:6]], ['missing_rank', 'empty_rank'])
        self.assertEqual(len([row for row in mapping if row['cleaned_taxon'] == 'Shared']), 3)

    def test_cleaning_is_target_specific_at_all_levels_and_supports_qiime_prefixes(self):
        prefixes = ['p', 'c', 'o', 'f', 'g', 's']
        for (level, name), prefix in zip(LEVELS, prefixes):
            with self.subTest(level=level):
                counts = (['001'],
                          ['k__Bacteria;' + prefix + '__' + name + ';x__uncultured',
                           'D_0__Bacteria;D_' + str(level - 1) + '__' + name,
                           'k__unclassified;' + prefix + '__Known',
                           'k__Bacteria;' + prefix + '__unclassified'],
                          [[1], [2], [3], [4]])
                cleaned, _, _, report = clean_taxonomy_table(counts, level)
                self.assertEqual(cleaned[1:], ([name, 'Known', 'Others'], [[3], [3], [4]]))
                self.assertEqual(report['reason_counts'], {'retained': 3, 'uncertain_rank': 1})
        counts = (['001'], ['k__Bacteria;p__Firmicutes;f__Family;g__uncultured;s__Species'], [[5]])
        for level, expected in [(2, 'Firmicutes'), (5, 'Family'), (6, 'Others'), (7, 'Species')]:
            self.assertEqual(clean_taxonomy_table(counts, level)[0][1], [expected])

    def test_cleaning_missing_empty_and_uncertain_target_names(self):
        uncertain = ['Incertae_Sedis', 'INCERTAE SEDIS', 'incertae-sedis', 'norank',
                     'uncultured bacterium', 'Unclassified', 'unassigned', 'unknown',
                     'unidentified', 'unresolved', 'ambiguous_taxa', 'NA', 'N/A',
                     'nan', 'None', 'null', 'missing', 'Others', 'not_assigned']
        lineages = ['p__P;g__' + name for name in uncertain]
        lineages += ['p__P;g__', 'p__P;f__Family;s__Species',
                     'Bacteria;P;C;O;F;PlainGenus', 'p__P;g__Candidatus_Nanobacterium']
        counts = (['001'], lineages, [[1] for _ in lineages])
        cleaned, _, mapping, report = clean_taxonomy_table(counts, 6)
        self.assertEqual(cleaned[1:], (['Candidatus_Nanobacterium', 'Others'], [[1], [len(lineages) - 1]]))
        self.assertEqual(report['reason_counts'], {'uncertain_rank': len(uncertain),
                                                 'empty_rank': 1, 'missing_rank': 2, 'retained': 1})
        self.assertEqual(mapping[-2]['reason'], 'missing_rank')

    def test_clean_exports_keep_original_tables_and_roundtrip_csv(self):
        counts = (['001', '010,case'], ['p__P;g__Named, group', 'p__P;g__uncultured'], [[2, 0], [6, 8]])
        raw = write_table_pair(self.root, 'counts', counts)
        raw_relative = write_table_pair(self.root, 'relative',
                                        (counts[0], counts[1], [[0.25, 0], [0.75, 1]]))
        before = {path: path.read_bytes() for path in [*raw.values(), *raw_relative.values()]}
        outputs = export_clean_taxonomy(self.root, counts, 6)
        self.assertEqual(before, {path: path.read_bytes() for path in before})
        cleaned = read_table(outputs['counts']['tsv'])
        relative = read_table(outputs['relative']['tsv'], relative=True)
        validate_relative(cleaned, relative)
        with outputs['counts']['csv'].open(newline='') as stream:
            rows = list(csv.reader(stream))
        self.assertEqual(rows, [['FeatureID', '001', '010,case'], ['Named, group', '2', '0'], ['Others', '6', '8']])
        with outputs['mapping']['csv'].open(newline='') as stream:
            mapping = list(csv.DictReader(stream))
        self.assertEqual(mapping[0]['original_lineage'], counts[1][0])
        self.assertEqual(mapping[0]['cleaned_taxon'], 'Named, group')
        self.assertEqual(json.loads(outputs['report'].read_text())['per_sample_counts'], {'001': 8, '010,case': 8})

    def test_all_six_levels_use_filtered_counts_and_actual_relative_outputs(self):
        source = self.source()
        api = FakeCLI(source)
        manifest = run_abundance({'run': str(source)}, self.saved, api)
        self.assertEqual(manifest, json.loads((source / '10_abundance/manifest.json').read_text()))
        self.assertEqual(manifest['summary']['per_sample_counts'], {'001': 8, '010,case': 8})
        collapses = [args for args in api.calls if args[:2] == ('taxa', 'collapse')]
        self.assertEqual([args[args.index('--p-level') + 1] for args in collapses], [str(level) for level, _ in LEVELS])
        self.assertEqual(len([args for args in api.calls if args[:2] == ('feature-table', 'summarize')]), 6)
        self.assertEqual(len([args for args in api.calls if args[:2] == ('feature-table', 'relative-frequency')]), 7)
        for item, (level, name) in zip(manifest['levels'], LEVELS):
            self.assertEqual((item['level'], item['name']), (level, name))
            counts = read_table(source / item['counts']['tsv'])
            relative = read_table(source / item['relative']['tsv'], relative=True)
            validate_relative(counts, relative)
            self.assertEqual(counts[2], [[8, 8]])
            cleaned = item['cleaned']
            clean_counts = read_table(source / cleaned['counts']['tsv'])
            clean_relative = read_table(source / cleaned['relative']['tsv'], relative=True)
            validate_relative(clean_counts, clean_relative)
            self.assertEqual(clean_counts[2], counts[2])
            self.assertNotIn('__', clean_counts[1][0])
            self.assertIn('__', counts[1][0])
            self.assertEqual(cleaned['feature_count'], 1)
            for kind in ('counts', 'relative', 'mapping'):
                for extension in ('tsv', 'csv'):
                    self.assertTrue((source / cleaned[kind][extension]).is_file())
            self.assertEqual(json.loads((source / cleaned['report']).read_text())['level'], level)
            for key in ('counts_qza', 'counts_qzv', 'counts_biom', 'relative_qza', 'relative_biom'):
                self.assertTrue((source / item[key]).is_file())
        for args in collapses:
            self.assertEqual(args[args.index('--i-table') + 1], str(artifact_path(source, 'filtered-table.qza')))
            self.assertEqual(args[args.index('--i-taxonomy') + 1], str(artifact_path(source, 'taxonomy.qza')))
        with (source / manifest['asv']['filtered_taxonomy']['csv']).open(newline='') as stream:
            rows = list(csv.reader(stream))
        self.assertEqual([row[0] for row in rows[1:]], ['ASV1', 'ASV2'])
        self.assertEqual(json.loads((source / '10_abundance/report.json').read_text())['status'], 'done')
        self.assertTrue(manifest['taxa_barplot']['files'])
        previous_calls = len(api.calls)
        run_abundance({'run': str(source)}, self.saved, api)
        self.assertEqual(len(api.calls), previous_calls)

    def test_historical_parent_remains_unchanged_for_legacy_and_versioned_layout(self):
        for layout in (False, True):
            with self.subTest(layout=layout):
                source = self.source(layout=layout)
                before = {filename.relative_to(source): filename.read_bytes() for filename in source.rglob('*') if filename.is_file()}
                output = self.root / ('extension-' + str(layout))
                output.mkdir()
                api = FakeCLI(source)
                manifest = run_abundance({'run': str(output)}, self.saved, api, source_run=source)
                after = {filename.relative_to(source): filename.read_bytes() for filename in source.rglob('*') if filename.is_file()}
                self.assertEqual(before, after)
                self.assertEqual(manifest['source_table'], str(artifact_path(source, 'filtered-table.qza')))
                self.assertTrue((output / manifest['levels'][0]['counts']['tsv']).is_file())
                self.assertFalse(any(args[0] == 'dada2' for args in api.calls))
                barplot = next(args for args in api.calls if args[:2] == ('tools', 'export') and args[args.index('--input-path') + 1].endswith('taxa-barplot.qzv'))
                self.assertEqual(barplot[barplot.index('--input-path') + 1], str(artifact_path(source, 'taxa-barplot.qzv')))

    def test_inconsistent_normalization_or_collapsed_totals_fail(self):
        for option, message in (('bad_relative', '相对丰度与未抽平计数不一致'), ('lost_counts', '序列总数发生变化')):
            with self.subTest(option=option):
                source = self.source(layout=option == 'bad_relative')
                with self.assertRaisesRegex(ValueError, message):
                    run_abundance({'run': str(source)}, self.saved, FakeCLI(source, **{option: True}))
                self.assertFalse((source / '10_abundance/manifest.json').exists())


if __name__ == '__main__':
    unittest.main()
