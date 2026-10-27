"""Synthetic input adapters verify validation and isolation, not real QIIME fits."""
import csv
import hashlib
import json
import tempfile
import unittest
import zipfile
from pathlib import Path
from types import SimpleNamespace

from external_inputs import MODULE_REQUIREMENTS, artifact_type, prepare_external, read_counts
from run_layout import artifact_path


UUID = '12345678-1234-1234-1234-123456789abc'
TABLE = '#OTU ID\t001\t002\nASV1\t12\t34\nASV2\t5\t6\n'
TAXONOMY = ('FeatureID\tTaxon\tConfidence\n'
            'ASV1\td__Bacteria; p__P; c__C; o__O; f__F; g__G; s__S1\t0.9999999999999999\n'
            'ASV2\td__Bacteria; p__P; c__C; o__O; f__F; g__G; s__S2\t-1\n')


def qza(filename, semantic_type, files=None, action=None):
    with zipfile.ZipFile(filename, 'w') as archive:
        archive.writestr(UUID + '/metadata.yaml', 'uuid: ' + UUID + '\ntype: ' + semantic_type + '\nformat: Synthetic\n')
        archive.writestr(UUID + '/VERSION', 'QIIME 2\narchive: 5\nframework: 2025.7.0\n')
        for name, content in (files or {}).items():
            archive.writestr(UUID + '/data/' + name, content)
        if action:
            archive.writestr(UUID + '/provenance/action/action.yaml', 'action:\n    type: method\n    action: ' + action + '\n')


class ExternalInputTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / 'external original; $(untouched)'
        self.source.mkdir()
        self.output = self.root / 'output with spaces'
        self.table = self.source / 'counts.tsv'
        self.table.write_text(TABLE, encoding='utf8')
        self.taxonomy = self.source / 'taxonomy.tsv'
        self.taxonomy.write_text(TAXONOMY, encoding='utf8')
        self.metadata = self.source / 'groups.csv'
        self.metadata.write_text('SampleID,Group\n002,B\n001,A\n', encoding='utf8')
        self.fasta = self.source / 'representative.fasta'
        self.fasta.write_text('>ASV1\nACGT\n>ASV2\nAGGT\n', encoding='utf8')
        self.tree = self.source / 'rooted.qza'
        qza(self.tree, 'Phylogeny[Rooted]', {'tree.nwk': '(ASV1:0.1,ASV2:0.2);'})
        self.calls = []
        self.tree_ids = ['ASV1', 'ASV2']
        self.api = SimpleNamespace(path=lambda value: Path(value).resolve(), artifact_path=artifact_path,
                                   marker=lambda run: (run / 'workstation_run.json').write_text('{"application":"16S-workstation","format":2}'),
                                   q=self.q, command=self.command, conda_cmd=lambda config: ['conda', 'run', '-n', config['env']],
                                   begin_stage=lambda *args, **kwargs: None, finish_stage=lambda *args, **kwargs: None)

    def config(self, module='lefse', **changes):
        return dict(external_module=module, external_counts_confirmed=True, external_table=str(self.table),
                    metadata=str(self.metadata), external_taxonomy=str(self.taxonomy), external_tree=str(self.tree),
                    external_sequences=str(self.fasta), output=str(self.output), **changes)

    def q(self, job, config, *raw):
        args = list(map(str, raw)); self.calls.append(args)
        if args[:2] == ['tools', 'export']:
            source = Path(args[args.index('--input-path') + 1])
            destination = Path(args[args.index('--output-path') + 1]); destination.mkdir(parents=True)
            with zipfile.ZipFile(source) as archive:
                for name in archive.namelist():
                    if name.startswith(UUID + '/data/') and not name.endswith('/'):
                        (destination / Path(name).name).write_bytes(archive.read(name))
        elif args[:2] == ['tools', 'import']:
            source = Path(args[args.index('--input-path') + 1])
            target = Path(args[args.index('--output-path') + 1]); target.parent.mkdir(parents=True, exist_ok=True)
            semantic_type = args[args.index('--type') + 1]
            files = {}
            if semantic_type == 'FeatureTable[Frequency]':
                files['feature-table.biom'] = source.read_bytes()
            elif semantic_type == 'FeatureData[Sequence]':
                files['dna-sequences.fasta'] = source.read_bytes()
            qza(target, semantic_type, files)
        elif args[:2] == ['taxa', 'barplot']:
            Path(args[args.index('--o-visualization') + 1]).write_bytes(b'synthetic visualization')
        else:
            self.fail('Unexpected QIIME action: ' + repr(args))

    def command(self, job, raw, cwd=None):
        args = list(map(str, raw)); self.calls.append(args)
        if 'biom' in args:
            source = Path(args[args.index('-i') + 1]); target = Path(args[args.index('-o') + 1])
            target.write_bytes(source.read_bytes())
        elif 'python' in args:
            Path(args[-1]).write_text(json.dumps(self.tree_ids), encoding='utf8')
        else:
            self.fail('Unexpected command: ' + repr(args))

    def prepare(self, config=None):
        job = {}
        result = prepare_external(job, config or self.config(), self.api)
        return job, result, Path(result['external_run'])

    def test_all_modules_prepare_only_own_prerequisites_without_changing_sources(self):
        original = {str(filename): filename.read_bytes() for filename in self.source.iterdir()}
        for module in MODULE_REQUIREMENTS:
            with self.subTest(module=module):
                self.calls.clear()
                job, result, run = self.prepare(self.config(module))
                completion = json.loads((run / 'external_inputs_complete.json').read_text())
                saved = json.loads((run / 'parameters.json').read_text())
                self.assertEqual(result['supported_modules'], [module])
                self.assertEqual(saved['supported_modules'], [module])
                self.assertEqual(completion['supported_modules'], [module])
                self.assertEqual(saved['source_kind'], 'external')
                self.assertEqual(saved['samples'], [{'SampleID': '001', 'Group': 'A'}, {'SampleID': '002', 'Group': 'B'}])
                self.assertEqual(set(completion['sources']), set(MODULE_REQUIREMENTS[module]))
                self.assertTrue((run / 'workstation_run.json').is_file())
                self.assertFalse((run / 'preprocess_complete.json').exists())
                self.assertFalse(any('dada2' in str(call) for call in self.calls))
                self.assertEqual(read_counts(artifact_path(run, 'exported-table/feature-table.tsv')), (['001', '002'], ['ASV1', 'ASV2'], [[12, 34], [5, 6]]))
                self.assertEqual(artifact_path(run, 'rooted-tree.qza').exists(), module == 'diversity')
                self.assertEqual(artifact_path(run, 'filtered-rep-seqs.qza').exists(), module == 'functional')
                self.assertEqual(artifact_path(run, 'taxa-barplot.qzv').exists(), module == 'abundance')
                for source in completion['sources'].values():
                    self.assertEqual(hashlib.sha256((run / source['copy']).read_bytes()).hexdigest(), source['sha256'])
                self.assertEqual(job['result'], result)
        for filename, content in original.items():
            self.assertEqual(Path(filename).read_bytes(), content)

    def test_repeated_preparation_uses_new_parent(self):
        _, _, first = self.prepare()
        original = {filename.relative_to(first): filename.read_bytes() for filename in first.rglob('*') if filename.is_file()}
        _, _, second = self.prepare()
        self.assertNotEqual(first, second)
        for filename, content in original.items():
            self.assertEqual((first / filename).read_bytes(), content)

    def test_module_guards_and_history_accept_only_prepared_capabilities(self):
        import app
        from batch_context import restore_batch, list_batches
        for module in MODULE_REQUIREMENTS:
            _, _, run = self.prepare(self.config(module))
            app.require_analysis_source(run, module)
            other = next(name for name in MODULE_REQUIREMENTS if name != module)
            with self.assertRaisesRegex(ValueError, '所需文件'):
                app.require_analysis_source(run, other)
            self.assertEqual(restore_batch(run, app.module_api())['fields']['run'], str(run))
        history = list_batches(self.output, app.module_api())['runs']
        self.assertEqual(len(history), len(MODULE_REQUIREMENTS))
        self.assertTrue(all(row['external_inputs_complete'] and not row['preprocess_complete'] for row in history))

    def test_worker_prepares_external_source_without_fastq_processing(self):
        import app
        from unittest.mock import patch
        before = {file: file.read_bytes() for file in self.source.iterdir()}
        job = {'status': 'running', 'action': 'external-inputs'}
        with patch.object(app, 'q', side_effect=self.q), \
                patch.object(app, 'command', side_effect=lambda j, a, cwd=None, **kw: self.command(j, a, cwd)), \
                patch.object(app, 'conda_cmd', return_value=['synthetic-conda']):
            app.worker(job, 'external-inputs', self.config('advanced'))
        self.assertEqual(job['status'], 'done', job.get('error'))
        run = Path(job['result']['external_run'])
        self.assertEqual(json.loads((run/'task_state.json').read_text())['action'], 'external-inputs')
        self.assertFalse((run/'preprocess_complete.json').exists())
        for file, content in before.items():
            self.assertEqual(file.read_bytes(), content)

    def test_frequency_qza_uses_official_export_and_validates_actual_counts(self):
        table = self.source / 'input.qza'
        qza(table, 'FeatureTable[Frequency]', {'feature-table.biom': TABLE})
        config = self.config(); config['external_table'] = str(table)
        _, _, run = self.prepare(config)
        self.assertEqual(artifact_path(run, 'filtered-table.qza').read_bytes(), table.read_bytes())
        self.assertTrue(any(call[:2] == ['tools', 'export'] for call in self.calls))
        qza(table, 'FeatureTable[Frequency]', {'feature-table.biom': TABLE.replace('\t12\t', '\t0.12\t')})
        with self.assertRaisesRegex(ValueError, '整数'):
            self.prepare(config)

    def test_qza_taxonomy_and_sequences_export(self):
        taxonomy = self.source / 'taxonomy.qza'
        qza(taxonomy, 'FeatureData[Taxonomy]', {'taxonomy.tsv': TAXONOMY})
        config = self.config(); config['external_taxonomy'] = str(taxonomy)
        _, _, run = self.prepare(config)
        self.assertEqual(artifact_path(run, 'exported-taxonomy/taxonomy.tsv').read_text().replace('Feature ID', 'FeatureID'), TAXONOMY)
        seqs = self.source / 'sequences.qza'
        qza(seqs, 'FeatureData[Sequence]', {'dna-sequences.fasta': self.fasta.read_bytes()})
        config = self.config('functional'); config['external_sequences'] = str(seqs)
        self.prepare(config)

    def test_csv_quoted_ids_counts_and_taxonomy_confidence_are_preserved(self):
        table = self.source / 'counts.csv'
        with table.open('w', newline='', encoding='utf8') as stream:
            csv.writer(stream).writerows([['FeatureID', '001', 'a,b'], ['ASV1', 12, 34], ['ASV2', 5, 6]])
        self.metadata.write_text('SampleID,Group\n001,A\n"a,b",B\n')
        config = self.config(); config['external_table'] = str(table)
        _, _, run = self.prepare(config)
        self.assertEqual(read_counts(artifact_path(run, 'exported-table/feature-table.tsv'))[0], ['001', 'a,b'])
        self.assertIn('0.9999999999999999', artifact_path(run, 'exported-taxonomy/taxonomy.tsv').read_text())
        self.assertTrue(any(str(self.output) in argument for call in self.calls for argument in call))
        self.assertFalse(any(isinstance(argument, list) for call in self.calls for argument in call))

    def test_advanced_aggregation_preserves_totals_and_lineages(self):
        _, _, run = self.prepare(self.config('advanced'))
        manifest = json.loads((run / '10_abundance/manifest.json').read_text())
        genus = next(entry for entry in manifest['levels'] if entry['level'] == 6)
        ids, features, numbers = read_counts(run / genus['counts']['tsv'])
        self.assertEqual(ids, ['001', '002'])
        self.assertEqual(numbers, [[17, 40]])
        self.assertEqual(features, ['d__Bacteria;p__P;c__C;o__O;f__F;g__G'])

    def test_genus_only_taxonomy_for_abundance_marks_absent_species_as_unknown(self):
        original = TAXONOMY.replace('; s__S1', '').replace('; s__S2', '')
        self.taxonomy.write_text(original)
        for use_qza in (False, True):
            config = self.config('abundance')
            if use_qza:
                artifact = self.source / 'genus-taxonomy.qza'
                qza(artifact, 'FeatureData[Taxonomy]', {'taxonomy.tsv': original})
                config['external_taxonomy'] = str(artifact)
            _, _, run = self.prepare(config)
            derived = artifact_path(run, 'exported-taxonomy/taxonomy.tsv').read_text()
            self.assertIn('g__G;__\t0.9999999999999999', derived)
            self.assertEqual((run / '00_inputs/taxonomy_before_padding.tsv').read_text().replace('Feature ID', 'FeatureID'), original)
            report = json.loads((run / '00_inputs/taxonomy_padding_report.json').read_text())
            self.assertEqual([row['unknown_trailing_ranks_added'] for row in report['features']], [1, 1])
            self.assertEqual(self.taxonomy.read_text(), original)

    def test_source_change_during_preparation_is_detected(self):
        original_command = self.api.command
        def modify_source(job, args, cwd=None):
            original_command(job, args, cwd)
            self.table.write_text(TABLE + 'ASV3\t1\t2\n')
        self.api.command = modify_source
        job = {}
        with self.assertRaisesRegex(ValueError, '源文件发生变化'):
            prepare_external(job, self.config(), self.api)
        self.assertFalse(Path(job['run'], 'external_inputs_complete.json').exists())

    def test_required_fields_and_raw_count_confirmation(self):
        for value in [False, 'true', 1, None]:
            config = self.config(); config['external_counts_confirmed'] = value
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, '确认'):
                self.prepare(config)
        for module, key in [('diversity', 'external_tree'), ('functional', 'external_sequences'),
                            ('abundance', 'external_taxonomy'), ('differential', 'external_taxonomy'),
                            ('lefse', 'external_taxonomy'), ('advanced', 'external_taxonomy')]:
            config = self.config(module); config[key] = ''
            with self.subTest(module=module), self.assertRaisesRegex(ValueError, key):
                self.prepare(config)

    def test_metadata_requires_literal_exact_sample_sets(self):
        for contents in ['SampleID,Group\n1,A\n002,B\n', 'SampleID,Group\n001,A\n001,B\n',
                         'SampleID,Group\n001,A\n', 'SampleID,Group\n001,A\n002,B\n003,C\n',
                         'SampleID,Group\n 001,A\n002,B\n']:
            self.metadata.write_text(contents)
            with self.subTest(contents=contents), self.assertRaisesRegex(ValueError, 'SampleID'):
                self.prepare()

    def test_invalid_counts_and_semantic_comments_are_rejected(self):
        for contents in [TABLE.replace('\t12\t', '\t0.12\t'), TABLE.replace('\t12\t', '\t-1\t'),
                         TABLE.replace('\t12\t', '\tNaN\t'), TABLE.replace('ASV2', 'ASV1'),
                         TABLE.replace('002', '001'), TABLE.replace('\t12\t', '\t9007199254740993\t'),
                         '# relative abundance\n' + TABLE, '# rarefied at 1000 reads\n' + TABLE,
                         'sample,A,B\nASV1,1,2\n']:
            self.table.write_text(contents)
            with self.subTest(contents=contents), self.assertRaises(ValueError):
                self.prepare()

    def test_taxonomy_requires_all_count_features_without_duplicates(self):
        for contents in [TAXONOMY.replace('ASV2', 'ASV3'), TAXONOMY.replace('ASV2', 'ASV1')]:
            self.taxonomy.write_text(contents)
            with self.subTest(contents=contents), self.assertRaisesRegex(ValueError, 'FeatureID'):
                self.prepare()

    def test_auxiliary_supersets_are_preserved_and_unused_features_reported(self):
        original = TAXONOMY + 'extra\td__Bacteria; p__Other\t1\n'
        self.taxonomy.write_text(original)
        for use_qza in (False, True):
            config = self.config()
            if use_qza:
                artifact = self.source / 'all-taxonomy.qza'
                qza(artifact, 'FeatureData[Taxonomy]', {'taxonomy.tsv': original})
                config['external_taxonomy'] = str(artifact)
            _, _, run = self.prepare(config)
            report = json.loads((run / '00_inputs/feature_alignment_report.json').read_text())
            self.assertEqual(report['taxonomy']['unused_feature_ids'], ['extra'])
            self.assertNotIn('extra', artifact_path(run, 'exported-taxonomy/taxonomy.tsv').read_text())
            self.assertEqual(self.taxonomy.read_text(), original)
        sequence_original = '>ASV1\nACGT\n>ASV2\nAGGT\n>extra\nTTTT\n'
        self.fasta.write_text(sequence_original)
        for use_qza in (False, True):
            config = self.config('functional')
            if use_qza:
                artifact = self.source / 'all-sequences.qza'
                qza(artifact, 'FeatureData[Sequence]', {'dna-sequences.fasta': sequence_original})
                config['external_sequences'] = str(artifact)
            _, _, run = self.prepare(config)
            report = json.loads((run / '00_inputs/feature_alignment_report.json').read_text())
            self.assertEqual(report['sequences']['unused_feature_ids'], ['extra'])
            self.assertNotIn('extra', artifact_path(run, 'exported-sequences/dna-sequences.fasta').read_text())
            self.assertEqual(self.fasta.read_text(), sequence_original)
        self.tree_ids = ['ASV1', 'ASV2', 'extra']
        qza(self.tree, 'Phylogeny[Rooted]', {'tree.nwk': '(ASV1:0.1,ASV2:0.2,extra:0.3);'})
        tree_original = self.tree.read_bytes()
        _, _, run = self.prepare(self.config('diversity'))
        report = json.loads((run / '00_inputs/feature_alignment_report.json').read_text())
        self.assertEqual(report['rooted_tree']['unused_feature_ids'], ['extra'])
        self.assertEqual(artifact_path(run, 'rooted-tree.qza').read_bytes(), tree_original)

    def test_sequences_and_tree_require_exact_feature_sets_and_no_duplicates(self):
        for content in ['>ASV1\nACGT\n', '>ASV1\nACGT\n>ASV1\nAGGT\n', '>ASV1\nACGT\n>ASV3\nAGGT\n']:
            self.fasta.write_text(content)
            with self.subTest(content=content), self.assertRaises(ValueError):
                self.prepare(self.config('functional'))
        for tips in [['ASV1'], ['ASV1', 'ASV3'], ['ASV1', 'ASV1']]:
            self.tree_ids = tips
            with self.subTest(tips=tips), self.assertRaisesRegex(ValueError, 'FeatureID'):
                self.prepare(self.config('diversity'))
        qza(self.tree, 'Phylogeny[Unrooted]')
        with self.assertRaisesRegex(ValueError, 'Phylogeny\\[Rooted\\]'):
            self.prepare(self.config('diversity'))

    def test_qza_types_provenance_and_archive_paths_are_checked(self):
        artifact = self.source / 'unsafe.qza'
        for kind in ['FeatureTable[RelativeFrequency]', 'FeatureTable[PresenceAbsence]']:
            qza(artifact, kind)
            with self.subTest(kind=kind), self.assertRaisesRegex(ValueError, 'FeatureTable\\[Frequency\\]'):
                artifact_type(artifact, 'FeatureTable[Frequency]', raw_counts=True)
        for action in ['rarefy', 'relative_frequency', 'normalize']:
            qza(artifact, 'FeatureTable[Frequency]', action=action)
            with self.subTest(action=action), self.assertRaisesRegex(ValueError, 'provenance'):
                artifact_type(artifact, 'FeatureTable[Frequency]', raw_counts=True)
        for path in [UUID + '/../escape', '/tmp/escape', UUID + '/data/..\\escape']:
            qza(artifact, 'FeatureTable[Frequency]')
            with zipfile.ZipFile(artifact, 'a') as archive:
                archive.writestr(path, 'never extracted')
            with self.subTest(path=path), self.assertRaisesRegex(ValueError, '路径'):
                artifact_type(artifact, 'FeatureTable[Frequency]', raw_counts=True)
        self.assertFalse((self.root / 'escape').exists())

    def test_failed_import_never_gets_completion_marker(self):
        config = self.config(); self.taxonomy.write_text(TAXONOMY.replace('ASV2', 'wrong'))
        job = {}
        with self.assertRaises(ValueError):
            prepare_external(job, config, self.api)
        self.assertFalse(Path(job['run'], 'external_inputs_complete.json').exists())


if __name__ == '__main__':
    unittest.main()
