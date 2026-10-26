"""Focused input-integrity and QC-report checks; QIIME processing is mocked."""
import csv
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import app
import test_quality_pipeline as fixtures
from qc_report import write_dada2_retention
from test_mapping import fastq


class AutomaticMappingIntegrityTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def test_two_column_mapping_rejects_reused_symlink_and_hardlink(self):
        for kind in ('symlink', 'hardlink'):
            with self.subTest(kind=kind):
                raw = self.root / kind
                raw.mkdir()
                for mate in (1, 2):
                    source = fastq(raw / f'A_R{mate}.fastq', mate=mate)
                    target = raw / f'B_R{mate}.fastq'
                    if kind == 'symlink': target.symlink_to(source.name)
                    else: os.link(source, target)
                sheet = self.root / f'{kind}.csv'
                sheet.write_text('SampleID,Group\nA,CON\nB,TREAT\n', encoding='utf8')
                with self.assertRaisesRegex(ValueError, '文件重复引用'):
                    app.discover({'input': str(raw), 'metadata': str(sheet), 'mode': 'paired'})

    def test_automatic_symlink_outside_input_is_rejected(self):
        raw = self.root / 'raw'
        raw.mkdir()
        source = fastq(self.root / 'outside.fastq')
        (raw / '001.fastq').symlink_to(source)
        sheet = self.root / 'groups.csv'
        sheet.write_text('SampleID,Group\n001,CON\n', encoding='utf8')
        with self.assertRaisesRegex(ValueError, '必须位于选择的原始目录内'):
            app.discover({'input': str(raw), 'metadata': str(sheet), 'mode': 'single'})

    def test_distinct_automatic_inputs_preserve_ids_and_raw_bytes(self):
        raw = self.root / 'raw'
        raw.mkdir()
        originals = {}
        for mate in (1, 2):
            filename = fastq(raw / f'001_R{mate}.fastq', mate=mate)
            originals[filename] = filename.read_bytes()
        sheet = self.root / 'groups.csv'
        sheet.write_text('SampleID,Group\n001,CON\n', encoding='utf8')
        samples = app.discover({'input': str(raw), 'metadata': str(sheet), 'mode': 'paired'})
        self.assertEqual(samples[0]['SampleID'], '001')
        self.assertEqual({filename: filename.read_bytes() for filename in originals}, originals)


class CurveDepthAuditTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.QualityPipelineTests('test_auto_parameters_explicit_chimera_rules_and_stage_outputs')
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.config = self.fixture.config | {'parameter_mode': 'manual', 'trunc_f': '300', 'trunc_r': '300',
                                            'trim_f': '0', 'trim_r': '0', 'max_depth': '1000'}

    def run_preprocess(self, conversion=None):
        calls = []
        def fake(job, config, *args):
            calls.append(args)
            fixtures.materialize_qiime_outputs(job, config, *args)
        job = {}
        with patch.object(app, 'q', side_effect=fake), patch.object(app, 'command', side_effect=conversion or self.fixture.conversion), patch('abundance_tables.run_abundance', return_value={}):
            app.preprocess(job, self.config)
        return job, calls

    def test_requested_depth_is_capped_with_explicit_saved_warning(self):
        job, calls = self.run_preprocess()
        run = Path(job['run'])
        saved = json.loads((run / 'parameters.json').read_text(encoding='utf8'))
        self.assertEqual(saved['max_depth'], '1000')
        self.assertEqual(saved['requested_curve_max_depth'], 1000)
        self.assertEqual(saved['effective_curve_max_depth'], 125)
        self.assertTrue(saved['curve_depth_warnings'])
        curve = next(args for args in calls if args[:2] == ('diversity', 'alpha-rarefaction'))
        self.assertEqual(curve[curve.index('--p-max-depth') + 1], 125)
        report = json.loads((run / '06_rarefaction/report.json').read_text(encoding='utf8'))
        self.assertEqual(report['results']['requested_max_depth'], 1000)
        self.assertEqual(report['results']['effective_max_depth'], 125)
        self.assertEqual(report['results']['warnings'], job['result']['curve_depth_warnings'])
        self.assertTrue((run / '03_dada2/dada2_retention.csv').is_file())

    def test_small_actual_maximum_rejects_before_curve_command(self):
        calls = []
        def fake(job, config, *args):
            calls.append(args)
            fixtures.materialize_qiime_outputs(job, config, *args)
        def conversion(job, args, cwd=None, process_env=None):
            self.fixture.conversion(job, args, cwd, process_env)
            if 'biom' in args:
                destination = Path(args[args.index('-o') + 1])
                destination.write_text('#OTU ID\t001\nASV1\t2\n', encoding='utf8')
        with patch.object(app, 'q', side_effect=fake), patch.object(app, 'command', side_effect=conversion):
            with self.assertRaisesRegex(ValueError, '最大样本计数为 2'):
                app.preprocess({}, self.config)
        self.assertFalse(any(args[:2] == ('diversity', 'alpha-rarefaction') for args in calls))


class Dada2RetentionReportTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def test_paired_adjacent_fractions_and_warning_do_not_remove_samples(self):
        source = {'per_sample': [
            {'SampleID': '001', 'input': 100, 'filtered': 80, 'denoised': 70, 'merged': 35, 'non-chimeric': 28},
            {'SampleID': 'low', 'input': 100, 'filtered': 20, 'denoised': 10, 'merged': 2, 'non-chimeric': 0},
        ]}
        result = write_dada2_retention(self.root, source, paired=True)
        report = json.loads(Path(result['json']).read_text(encoding='utf8'))
        good, low = report['per_sample']
        self.assertEqual([row['SampleID'] for row in report['per_sample']], ['001', 'low'])
        self.assertEqual((good['filter_retention_fraction'], good['merge_retention_fraction'], good['nonchimera_retention_fraction']), (0.8, 0.5, 0.8))
        self.assertFalse(good['warnings'])
        self.assertTrue(low['warnings'])
        self.assertFalse(report['criteria']['automatically_excludes_samples'])
        with Path(result['csv']).open(encoding='utf-8-sig', newline='') as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(rows[0]['SampleID'], '001')
        self.assertEqual(len(rows), 2)

    def test_single_end_and_zero_denominator_are_reported_honestly(self):
        source = {'per_sample': [{'SampleID': 'zero', 'input': 10, 'filtered': 0, 'denoised': 0, 'non-chimeric': 0}]}
        result = write_dada2_retention(self.root, source, paired=False)
        row = json.loads(Path(result['json']).read_text(encoding='utf8'))['per_sample'][0]
        self.assertNotIn('merge_retention_fraction', row)
        self.assertIsNone(row['nonchimera_retention_fraction'])
        self.assertEqual(row['filter_retention_fraction'], 0)
        self.assertTrue(row['warnings'])


if __name__ == '__main__':
    unittest.main()
