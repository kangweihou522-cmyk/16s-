"""Small synthetic tables verify orchestration; PICRUSt2 is not run here."""
import csv
import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import functional_prediction as fp
from run_layout import artifact_path, initialize_layout, stage_report


def compressed(filename, content):
    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(filename, 'wt', encoding='utf8') as stream:
        stream.write(content)


class FunctionalTableTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)

    def mapping(self, filename='map.tsv'):
        filename = self.root / filename
        with filename.open('w', encoding='utf8', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t' if filename.suffix == '.tsv' else ',')
            writer.writerow(fp.MAPPING_COLUMNS)
            row = ['ko:K00001', 'Metabolism', 'Shared', 'Leaf', 'map00010', 'geneA', 'name, with commas', '1.1.1.1, 2.2.2.2']
            writer.writerow(row)
            writer.writerow(row)
            writer.writerow(['K00001', 'Metabolism', 'Shared', 'Leaf', 'map00010', 'geneB', 'other annotation', ''])
            writer.writerow(['K00001', 'Metabolism', 'Shared', 'Leaf', 'map00020', 'geneA', 'second pathway', ''])
            writer.writerow(['K00002', 'Other parent', 'Shared', 'Leaf', 'map00020', 'geneC', 'other KO', ''])
        return filename

    def test_tsv_commas_csv_duplicates_and_ko_prefix(self):
        for filename in ('map.tsv', 'map.csv'):
            rows = fp.read_mapping(self.mapping(filename))
            self.assertEqual(len(rows), 4)
            self.assertEqual(rows[0]['KO'], 'K00001')
            self.assertEqual(rows[0]['KO_name'], 'name, with commas')
            self.assertEqual(rows[0]['EC'], '1.1.1.1, 2.2.2.2')

    def test_mapping_rejects_nonexact_headers_and_bad_ko(self):
        filename = self.mapping()
        original = filename.read_text()
        for invalid in (original.replace('KO\tL1', 'ko\tL1', 1), original.replace('ko:K00001', 'K1', 1), original.replace('KO\tL1', 'KO\tKO', 1)):
            filename.write_text(invalid)
            with self.assertRaises(ValueError):
                fp.read_mapping(filename)

    def test_gzip_fractional_counts_reordered_string_sample_ids(self):
        filename = self.root / 'ko.tsv.gz'
        compressed(filename, 'function\t010\t001\nko:K00001\t1.5\t2.5\nK00002\t2\t0\n')
        table = fp.read_abundance(filename, ['001', '010'], ko=True)
        self.assertEqual(table['values']['K00001'], [2.5, 1.5])
        relative, totals = fp.relative_table(table)
        self.assertEqual(totals, {'001': 2.5, '010': 3.5})
        self.assertEqual(relative['values']['K00001'][0], 1)
        with self.assertRaisesRegex(ValueError, '序列计数'):
            fp.read_abundance(filename, ['001', '010'], ko=True, integer=True)

    def test_invalid_counts_samples_and_ko_ids(self):
        filename = self.root / 'bad.tsv.gz'
        for text in ('function\t001\t999\nK00001\t1\t2\n', 'function\t001\t001\nK00001\t1\t2\n', 'function\t001\t010\nK00001\t-1\t2\n', 'function\t001\t010\nK00001\tnan\t2\n', 'function\t001\t010\nK00001\tinf\t2\n', 'function\t001\t010\nKO00001\t1\t2\n', 'function\t001\t010\nK00001\t1\t2\nko:K00001\t1\t2\n'):
            compressed(filename, text)
            with self.assertRaises(ValueError):
                fp.read_abundance(filename, ['001', '010'], ko=True)

    def test_membership_dedup_hierarchies_and_relative_denominators(self):
        table = {'samples': ['001', '010'], 'values': {'K00001': [2.5, 0], 'K00002': [1.5, 4], 'K99999': [8, 8]}}
        categories, annotations, unmatched = fp.aggregate_categories(table, fp.read_mapping(self.mapping()))
        self.assertEqual(categories['pathway']['values']['map00010'], [2.5, 0])
        self.assertEqual(categories['pathway']['values']['map00020'], [4, 4])
        self.assertEqual(len(categories['L2']['values']), 2)
        self.assertEqual(len(categories['L3']['values']), 2)
        self.assertEqual(unmatched, ['K99999'])
        self.assertEqual(annotations[('L2', '["Metabolism","Shared"]')]['L1'], 'Metabolism')
        relative, totals = fp.relative_table(categories['pathway'])
        self.assertEqual(totals, {'001': 6.5, '010': 4})
        self.assertAlmostEqual(sum(values[0] for values in relative['values'].values()), 1)

    def test_qc_missing_outputs_explicitly_unavailable(self):
        qc = fp.collect_qc(self.root, {'samples': ['001'], 'values': {'ASV1': [10]}}, self.root / 'qc', 2)
        self.assertEqual(qc['excluded_asvs']['status'], 'unavailable')
        self.assertEqual(qc['weighted_nsti']['KO']['status'], 'unavailable')
        self.assertEqual(qc['asv_nsti']['status'], 'unavailable')
        self.assertNotIn('count', qc['excluded_asvs'])
        self.assertFalse((self.root / 'qc/excluded_asvs.tsv').exists())


class FunctionalOrchestrationTests(unittest.TestCase):
    mapping = FunctionalTableTests.mapping

    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        self.parent = self.root / 'parent'
        self.parent.mkdir()
        initialize_layout(self.parent, 'preprocess')
        seq = artifact_path(self.parent, 'exported-sequences/dna-sequences.fasta')
        seq.parent.mkdir(parents=True, exist_ok=True)
        seq.write_text('>ASV1\nACGT\n>ASV2\nTGCA\n')
        biom = artifact_path(self.parent, 'exported-table/feature-table.biom')
        biom.parent.mkdir(parents=True, exist_ok=True)
        biom.write_bytes(b'placeholder: synthetic BIOM')
        self.run = self.parent / '11_extensions' / 'new_functional'
        self.run.mkdir(parents=True)
        (self.run / 'parameters.json').write_text(json.dumps({'extension': 'functional', 'preprocess': str(self.parent), 'config': {'env': 'qiime'}}))
        self.saved = {'env': 'qiime', 'samples': [{'SampleID': '001', 'Group': 'A'}, {'SampleID': '010', 'Group': 'B'}]}
        self.config = {'ko_mapping': str(self.mapping()), 'threads': '1'}
        self.job = {'run': str(self.run), 'log': str(self.run / 'workflow.log')}
        self.calls = []
        self.raw_table = '# Constructed from biom file\n#OTU ID\t001\t010\nASV1\t10\t20\nASV2\t5\t3\n'
        self.help_text = 'usage --ref_dir1 PATH --ref_dir2 PATH --in_traits EC,KO\n'
        self.api = SimpleNamespace(path=lambda value: Path(value).resolve(), artifact_path=artifact_path, find_conda=lambda config: '/chosen/conda', conda_cmd=lambda config: ['/chosen/conda', 'run', '--no-capture-output', '-n', config['env']], command=self.command, begin_stage=self.begin, finish_stage=self.finish)

    def begin(self, job, folder, title, parameters=None):
        job.update(stage_folder=folder, stage_title=title, stage_parameters=parameters or {})
        stage_report(job['run'], folder, title, status='running', parameters=parameters)

    def finish(self, job, results=None, notes=()):
        stage_report(job['run'], job['stage_folder'], job['stage_title'], parameters=job['stage_parameters'], results=results, notes=notes)
        for key in ('stage_folder', 'stage_title', 'stage_parameters'):
            job.pop(key)

    def command(self, job, args, cwd=None):
        args = list(map(str, args))
        self.calls.append(args)
        with Path(job['log']).open('a', encoding='utf8') as log:
            log.write(json.dumps(args) + '\n')
            if '--version' in args:
                log.write('PICRUSt2 2.6.2\n')
            if '--help' in args:
                log.write(self.help_text)
        if 'biom' in args:
            Path(args[args.index('-o') + 1]).write_text(self.raw_table)
        if '--output-path' in args:
            destination = Path(args[args.index('--output-path') + 1])
            destination.mkdir(parents=True)
            source = args[args.index('--input-path') + 1]
            if source.endswith('filtered-table.qza'):
                (destination / 'feature-table.biom').write_bytes(b'placeholder: exported BIOM')
            else:
                (destination / 'dna-sequences.fasta').write_text('>ASV1\nACGT\n>ASV2\nTGCA\n')
        if 'picrust2_pipeline.py' in args and '-o' in args:
            raw = Path(args[args.index('-o') + 1])
            self.assertFalse(raw.exists(), 'PICRUSt2 requires a fresh output directory')
            compressed(raw / 'KO_metagenome_out/pred_metagenome_unstrat.tsv.gz', 'function\t001\t010\nK00001\t2.5\t0\nK00002\t1.5\t4\n')
            compressed(raw / 'EC_metagenome_out/pred_metagenome_unstrat.tsv.gz', 'function\t001\t010\n1.1.1.1\t3.5\t4\n')
            compressed(raw / 'pathways_out/path_abun_unstrat.tsv.gz', 'pathway\t001\t010\nPWY-1\t1.5\t2\n')
            compressed(raw / 'combined_marker_predicted_and_nsti.tsv.gz', 'sequence\t16S_rRNA_Count\tmetadata_NSTI\nASV1\t2\t0.2\nASV2\t2\t3.0\n')
            for trait in ('KO', 'EC'):
                compressed(raw / (trait + '_metagenome_out/seqtab_norm.tsv.gz'), 'normalized\t001\t010\nASV1\t5\t10\n')
                compressed(raw / (trait + '_metagenome_out/weighted_nsti.tsv.gz'), '\tweighted_NSTI\n001\t0.2\n010\t0.2\n')
        if 'add_descriptions.py' in args:
            compressed(args[args.index('-o') + 1], 'synthetic annotation output\n')

    def test_isolated_recorded_prediction_defaults_outputs_and_qc(self):
        original_mapping = Path(self.config['ko_mapping'])
        mapping_bytes = original_mapping.read_bytes()
        before = {str(path): path.read_bytes() for path in self.parent.rglob('*') if path.is_file() and not path.is_relative_to(self.run)}
        result = fp.run_functional(self.job, self.config, self.saved, self.parent, self.api)
        prediction = next(args for args in self.calls if 'picrust2_pipeline.py' in args and '-o' in args)
        self.assertEqual(prediction[:5], ['/chosen/conda', 'run', '--no-capture-output', '-n', 'picrust2'])
        for option, value in [('--min_align', '0.8'), ('--max_nsti', '2.0'), ('--in_traits', 'EC,KO')]:
            self.assertEqual(prediction[prediction.index(option) + 1], value)
        self.assertNotIn('--ref_dir', prediction)
        self.assertEqual(result['effective_parameters']['installed_version'], '2.6.2')
        self.assertEqual(result['qc']['excluded_asvs']['excluded_input_reads'], {'001': 5, '010': 3})
        self.assertEqual(result['qc']['weighted_nsti']['KO']['per_sample'], {'001': 0.2, '010': 0.2})
        for table in result['tables'].values():
            for filename in table.values():
                self.assertTrue((self.run / filename).is_file())
        self.assertEqual(len([args for args in self.calls if 'add_descriptions.py' in args]), 3)
        parameters = json.loads((self.run / 'parameters.json').read_text())
        self.assertEqual(parameters['extension'], 'functional')
        self.assertEqual(parameters['config']['env'], 'qiime')
        for filename, contents in before.items():
            self.assertEqual(Path(filename).read_bytes(), contents)
        copied_mapping = self.run / result['effective_parameters']['ko_mapping_copy']
        original_mapping.write_text('modified after successful prediction\n')
        self.assertEqual(copied_mapping.read_bytes(), mapping_bytes)
        self.assertEqual(result['effective_parameters']['mapping_sha256'], hashlib.sha256(mapping_bytes).hexdigest())
        self.assertEqual(len(fp.read_mapping(copied_mapping)), 4)

    def test_exports_missing_parent_exports_into_child(self):
        artifact_path(self.parent, 'exported-table/feature-table.biom').unlink()
        artifact_path(self.parent, 'exported-sequences/dna-sequences.fasta').unlink()
        fp.run_functional(self.job, self.config, self.saved, self.parent, self.api)
        exports = [args for args in self.calls if '--output-path' in args]
        self.assertEqual(len(exports), 2)
        for args in exports:
            self.assertTrue(Path(args[args.index('--output-path') + 1]).is_relative_to(self.run))
        self.assertFalse(artifact_path(self.parent, 'exported-table/feature-table.biom').exists())

    def test_verified_dual_reference_flags(self):
        reference = self.root / 'prediction_reference'
        reference.mkdir()
        self.config.update(picrust_reference=str(reference), picrust_reference2=str(reference), min_align='0.7', max_nsti='1.5')
        fp.run_functional(self.job, self.config, self.saved, self.parent, self.api)
        args = next(args for args in self.calls if 'picrust2_pipeline.py' in args and '-o' in args)
        self.assertEqual(args[args.index('--ref_dir1') + 1], str(reference))
        self.assertEqual(args[args.index('--ref_dir2') + 1], str(reference))
        self.assertEqual(args[args.index('--min_align') + 1], '0.7')

    def test_missing_predicted_sample_fails_without_completion_claim(self):
        original_command = self.api.command
        def missing_sample(job, args, cwd=None):
            original_command(job, args, cwd)
            if 'picrust2_pipeline.py' in args and '-o' in args:
                raw = Path(args[args.index('-o') + 1])
                compressed(raw / 'KO_metagenome_out/pred_metagenome_unstrat.tsv.gz', 'function\t001\nK00001\t2.5\n')
        self.api.command = missing_sample
        with self.assertRaisesRegex(ValueError, 'SampleID'):
            fp.run_functional(self.job, self.config, self.saved, self.parent, self.api)
        self.assertNotIn('result', self.job)
        self.assertFalse((self.run / 'functional_complete.json').exists())
        self.assertTrue(artifact_path(self.parent, 'exported-table/feature-table.biom').exists())
        self.assertEqual(json.loads((self.run / '01_prediction/report.json').read_text())['status'], 'running')

    def test_unsupported_reference_and_mismatched_input_stop_before_prediction(self):
        reference = self.root / 'reference'
        reference.mkdir()
        self.config['picrust_reference'] = str(reference)
        self.help_text = 'usage --ref_dir PATH\n'
        with self.assertRaisesRegex(ValueError, '不支持 --ref_dir1'):
            fp.run_functional(self.job, self.config, self.saved, self.parent, self.api)
        self.assertFalse(any('picrust2_pipeline.py' in args and '-o' in args for args in self.calls))
        self.config.pop('picrust_reference')
        self.raw_table = '#OTU ID\t001\t010\nASV_DIFFERENT\t10\t20\n'
        with self.assertRaisesRegex(ValueError, '标识不一致'):
            fp.run_functional(self.job, self.config, self.saved, self.parent, self.api)
        self.assertFalse(any('picrust2_pipeline.py' in args and '-o' in args for args in self.calls))


if __name__ == '__main__':
    unittest.main()
