"""Synthetic QIIME outputs test orchestration, not real biological processing."""
import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import app
from run_layout import artifact_path, initialize_layout, quality_counts, table_summary
from test_mapping import fastq


def materialize_qiime_outputs(job, config, *args):
    """Shared mock writes the specific real export formats consumed by the app."""
    run = Path(job['run'])
    for index, value in enumerate(args[:-1]):
        if str(value).startswith('--o-'):
            target = artifact_path(run, args[index + 1])
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(b'placeholder: not a real QIIME artifact')
    if args[:2] != ('tools', 'export'):
        return
    source = str(args[args.index('--input-path') + 1])
    destination = artifact_path(run, args[args.index('--output-path') + 1])
    destination.mkdir(parents=True, exist_ok=True)
    saved = json.loads((run / 'parameters.json').read_text(encoding='utf8'))
    samples = saved.get('samples', saved.get('config', {}).get('samples', []))
    ids = [sample['SampleID'] for sample in samples]
    if source in ('demux.qzv', 'trimmed.qzv'):
        paired = saved['mode'] == 'paired'
        with (destination / 'per-sample-fastq-counts.tsv').open('w', newline='', encoding='utf8') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['sample ID', 'forward sequence count'] + (['reverse sequence count'] if paired else []))
            writer.writerows([sid, 1000] + ([1000] if paired else []) for sid in ids)
        for direction in ('forward', 'reverse') if paired else ('forward',):
            with (destination / (direction + '-seven-number-summaries.tsv')).open('w', newline='', encoding='utf8') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(['', *range(1, 301)])
                for key, number in [('count', 1000), ('2%', 30), ('9%', 32), ('25%', 33), ('50%', 35), ('75%', 38), ('91%', 40), ('98%', 40)]:
                    writer.writerow([key, *([number] * 300)])
    elif source == 'dada2-stats.qza':
        (destination / 'stats.tsv').write_text('sample-id\tinput\tfiltered\tdenoised\tmerged\tnon-chimeric\n#q2:types\tnumeric\tnumeric\tnumeric\tnumeric\tnumeric\n' + ''.join(f'{sid}\t1000\t950\t900\t800\t{125*(i+1)}\n' for i, sid in enumerate(ids)), encoding='utf8')
    elif source == 'alpha-rarefaction.qzv':
        for metric,value in [('observed_features',1),('shannon',0),('faith_pd',0)]:
            (destination/(metric+'.csv')).write_text('sample-id,depth-10_iter-1,depth-50_iter-1,depth-100_iter-1,Group\n'+''.join(f'{sample["SampleID"]},{value},{value},{value},{sample["Group"]}\n' for sample in samples),encoding='utf8')


class QualityPipelineTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.raw = self.root / 'raw'; self.raw.mkdir()
        for mate in (1, 2): fastq(self.raw / f'001_R{mate}.fastq.gz', mate=mate)
        self.metadata = self.root / 'groups.csv'
        self.metadata.write_text('SampleID,Group\n001,CON\n', encoding='utf8')
        classifier = self.root / 'classifier.qza'; classifier.write_bytes(b'placeholder')
        self.config = dict(input=str(self.raw), metadata=str(self.metadata), output=str(self.root / 'results'), mode='paired', primer_mode='removed', threads='1', env='mock', r_env='mock-r', conda='/fake/conda', classifier=str(classifier), font='DejaVu Serif', parameter_mode='auto', min_frequency='1', max_depth='', amplicon='')

    def conversion(self, job, args, cwd=None, process_env=None):
        if 'biom' in args:
            destination = Path(args[args.index('-o') + 1])
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text('#OTU ID\t001\nASV1\t125\n', encoding='utf8')

    def test_quality_exports_tables_and_recommends_without_dada2(self):
        job = {}; calls = []
        def fake(job, config, *args):
            calls.append(args); materialize_qiime_outputs(job, config, *args)
        with patch.object(app, 'q', side_effect=fake), patch.object(app, 'command'):
            app.preprocess(job, self.config, quality_only=True)
        run = Path(job['run'])
        for name in ('per-sample-fastq-counts.tsv', 'forward-seven-number-summaries.tsv', 'reverse-seven-number-summaries.tsv'):
            self.assertTrue((run / '01_import/quality_tables' / name).is_file())
        self.assertFalse(any(args[0] == 'dada2' for args in calls))
        self.assertEqual(job['result']['quality_advice']['parameters']['amplicon'], '')
        report = json.loads((run / '01_import/report.json').read_text())
        self.assertEqual(report['results']['totals']['forward sequence count'], 1000)
        self.assertEqual(report['results']['per_sample'][0]['SampleID'], '001')
        self.assertTrue((run / '处理说明.md').is_file())

    def test_auto_parameters_explicit_chimera_rules_and_stage_outputs(self):
        job = {}; calls = []
        def fake(job, config, *args):
            calls.append(args); materialize_qiime_outputs(job, config, *args)
        with patch.object(app, 'q', side_effect=fake), patch.object(app, 'command', side_effect=self.conversion):
            app.preprocess(job, self.config)
        run = Path(job['run'])
        saved = json.loads((run / 'parameters.json').read_text())
        self.assertEqual(saved['trunc_f'], 300)
        self.assertEqual(saved['amplicon'], '')
        dada = next(args for args in calls if args[:2] == ('dada2', 'denoise-paired'))
        for key, value in [('--p-chimera-method', 'consensus'), ('--p-min-overlap', '12'), ('--p-max-merge-mismatch', '0'), ('--p-max-ee-f', '2'), ('--p-max-ee-r', '2')]:
            self.assertEqual(dada[dada.index(key) + 1], value)
        for folder in ('00_inputs', '01_import', '02_primers', '03_dada2', '04_taxonomy_filter', '05_phylogeny', '06_rarefaction', '08_parameter_advice'):
            self.assertTrue((run / folder / 'report.json').is_file())
        self.assertTrue((run / '03_dada2/dada2-table.qza').is_file())
        self.assertFalse((run / 'dada2-table.qza').exists())
        self.assertEqual(app.totals(run), {'001': 125})

    def test_known_insufficient_overlap_does_not_run_dada2(self):
        job = {}; calls = []
        def fake(job, config, *args):
            calls.append(args); materialize_qiime_outputs(job, config, *args)
        with patch.object(app, 'q', side_effect=fake), patch.object(app, 'command'):
            with self.assertRaisesRegex(ValueError, '重叠不足'):
                app.preprocess(job, self.config | {'amplicon': '600'})
        self.assertFalse(any(args[0] == 'dada2' for args in calls))

    def test_existing_legacy_quality_can_export_without_reimport(self):
        run = self.root / 'old-run'; run.mkdir(); app.marker(run)
        (run / 'parameters.json').write_text(json.dumps(self.config | {'samples':[{'SampleID':'001','Group':'CON'}]}))
        (run / 'quality_complete.json').write_text('{}')
        calls = []; job = {}
        def fake(job, config, *args):
            calls.append(args); materialize_qiime_outputs(job, config, *args)
        with patch.object(app, 'q', side_effect=fake):
            app.refresh_quality_advice(job, {'quality_run': str(run)})
        self.assertEqual([args[:2] for args in calls], [('tools', 'export')])
        self.assertEqual(job['result']['quality_advice']['parameters']['trunc_f'], 300)
        self.assertFalse((run / 'layout.json').exists())

    def test_failure_report_preserves_completed_quality_stage(self):
        job = {'status':'running','action':'preprocess'}
        def fake(job, config, *args):
            if args[0] == 'dada2': raise RuntimeError('synthetic DADA2 failure')
            materialize_qiime_outputs(job, config, *args)
        with patch.object(app, 'q', side_effect=fake), patch.object(app, 'command'):
            app.worker(job, 'preprocess', self.config)
        run = Path(job['run'])
        self.assertEqual(job['status'], 'failed')
        self.assertEqual(json.loads((run / '03_dada2/report.json').read_text())['status'], 'failed')
        self.assertEqual(json.loads((run / '01_import/report.json').read_text())['status'], 'done')
        self.assertFalse((run / 'preprocess_complete.json').exists())

    def test_mapping_keeps_classifier_absolute_and_core_output_new(self):
        run = self.root / 'analysis'; run.mkdir(); initialize_layout(run, 'diversity')
        self.assertFalse(artifact_path(run, 'core-metrics-results').exists())
        with patch.object(app, 'command') as command:
            app.q({'run':str(run)}, self.config, 'tools', 'export', '--input-path', self.config['classifier'], '--output-path', 'exports/test')
        args = command.call_args.args[1]
        self.assertIn(self.config['classifier'], args)
        self.assertIn('02_exports/test', args)

    def test_tables_reject_invalid_counts(self):
        bad = self.root / 'bad.tsv'
        bad.write_text('#OTU ID\t001\nASV1\tNaN\n')
        with self.assertRaisesRegex(ValueError, '无效'): table_summary(bad)
        folder = self.root / 'counts'; folder.mkdir()
        (folder / 'per-sample-fastq-counts.tsv').write_text('sample ID\tforward sequence count\n001\t-2\n')
        with self.assertRaisesRegex(ValueError, '无效'): quality_counts(folder)

    def test_input_check_saves_a_report_and_keeps_raw_files(self):
        before = {p: p.read_bytes() for p in self.raw.iterdir()}
        job = {'status':'running', 'log':str(self.root / 'inspect.log')}
        app.worker(job, 'inspect', self.config)
        self.assertEqual(job['status'], 'done')
        report = json.loads((Path(job['run']) / '00_checks/report.json').read_text())
        self.assertEqual(report['results']['inspection'][0]['SampleID'], '001')
        self.assertEqual({p:p.read_bytes() for p in self.raw.iterdir()}, before)

    def test_single_end_runs_without_reverse_parameters(self):
        config=self.config | {'mode':'single','parameter_mode':'manual','trunc_f':'100','trim_f':'0'}
        self.metadata.write_text('SampleID,Group,forward,reverse\n001,CON,001_R1.fastq.gz,\n',encoding='utf8')
        job={}; calls=[]
        def fake(job, config, *args):
            calls.append(args); materialize_qiime_outputs(job,config,*args)
            if args[:2]==('tools','export') and args[args.index('--input-path')+1]=='dada2-stats.qza':
                stats=artifact_path(job['run'],'exports/denoising-stats')/'stats.tsv'
                stats.write_text('sample-id\tinput\tfiltered\tdenoised\tnon-chimeric\n001\t1000\t950\t900\t125\n',encoding='utf8')
        with patch.object(app,'q',side_effect=fake),patch.object(app,'command',side_effect=self.conversion):
            app.preprocess(job,config)
        dada=next(args for args in calls if args[:2]==('dada2','denoise-single'))
        self.assertIn('--p-max-ee',dada)
        self.assertNotIn('--p-min-overlap',dada)
        self.assertNotIn('--p-trunc-len-r',dada)
        report=json.loads((Path(job['run'])/'03_dada2/report.json').read_text())
        self.assertNotIn('merged',report['results']['denoising']['totals'])


if __name__ == '__main__': unittest.main()
