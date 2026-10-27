"""Synthetic command failures verify retry boundaries and preservation of inputs."""
import hashlib
import http.client
import json
import os
import tempfile
import threading
import unittest
from pathlib import Path
from unittest.mock import patch

import app
from recovery import persist_task, prepare_retry, recovery_plan, seal_checkpoint
from run_layout import artifact_path, initialize_layout, stage_report
from test_mapping import fastq
from test_quality_pipeline import materialize_qiime_outputs


def snapshot(directory):
    return {str(p.relative_to(directory)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in directory.rglob('*') if p.is_file()}


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.raw = self.root / 'raw'
        self.raw.mkdir()
        for mate in (1, 2):
            fastq(self.raw / f'001_R{mate}.fastq.gz', mate=mate)
        metadata = self.root / 'groups.csv'
        metadata.write_text('SampleID,Group\n001,CON\n')
        classifier = self.root / 'classifier.qza'
        classifier.write_bytes(b'synthetic classifier')
        self.config = dict(input=str(self.raw), metadata=str(metadata), output=str(self.root / 'results'),
                           mode='paired', primer_mode='removed', data_state='clean', threads='1',
                           env='mock', r_env='mock-r', conda='/fake/conda', classifier=str(classifier),
                           font='DejaVu Serif', parameter_mode='manual', min_frequency='1',
                           max_depth='', amplicon='', trunc_f='280', trunc_r='220', trim_f='0', trim_r='0',
                           batch_id='batch-a')

    def conversion(self, job, args, cwd=None, process_env=None):
        if 'biom' in args:
            destination = Path(args[args.index('-o') + 1])
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text('#OTU ID\t001\nASV1\t125\n')

    def failed_preprocess(self, failure='dada2'):
        job = {'id': 'original-' + failure, 'action': 'preprocess', 'status': 'running'}
        def fake(job, config, *args):
            if args[0] == failure:
                raise RuntimeError('synthetic ' + failure + ' failure')
            materialize_qiime_outputs(job, config, *args)
            if args[:2] == ('tools', 'import'):
                artifact_path(job['run'], 'demux.qza').write_bytes(b'complete synthetic import')
        with patch.object(app, 'q', side_effect=fake), patch.object(app, 'command', side_effect=self.conversion):
            app.worker(job, 'preprocess', self.config)
        self.assertEqual(job['status'], 'failed', job.get('error'))
        return Path(job['run'])

    def task(self, run, action, config, status='failed'):
        job = {'id': 'saved-task', 'action': action, 'run': str(run), 'status': status,
               'step': 'synthetic failure', '_task': {'action': action, 'config': config, 'instance_id': 'old-process'}}
        persist_task(job)

    def test_failed_dada2_forks_quality_accepts_numeric_correction_and_keeps_source_read_only(self):
        run = self.failed_preprocess()
        before, raw_before = snapshot(run), snapshot(self.raw)
        state = json.loads((run / 'task_state.json').read_text())
        self.assertEqual(state['action'], 'preprocess')
        self.assertEqual(state['status'], 'failed')
        plan = recovery_plan(run, app.module_api())
        self.assertTrue(plan['recoverable'])
        self.assertEqual(plan['strategy'], 'quality')
        self.assertEqual(plan['fields']['run'], '')
        self.assertEqual(plan['config']['trunc_f'], '280')
        self.assertEqual(before, snapshot(run))
        calls = []
        def fail_again(job, config, *args):
            calls.append(args)
            if args[0] == 'dada2':
                raise RuntimeError('synthetic corrected DADA2 failure')
            materialize_qiime_outputs(job, config, *args)
        retry = {'id': 'retry-1', 'action': 'recover', 'status': 'running'}
        with patch.object(app, 'q', side_effect=fail_again), patch.object(app, 'command', side_effect=self.conversion):
            app.worker(retry, 'recover', dict(plan['config'], recovery_run=str(run), trunc_f='250', r_env='fixed-r',
                                              forward='', reverse='', input_layout='unmerged'))
        self.assertEqual(retry['status'], 'failed')
        destination = Path(retry['run'])
        self.assertNotEqual(destination, run)
        self.assertEqual(retry['recovery_action'], 'preprocess')
        self.assertFalse(any(args[:2] == ('tools', 'import') or args[0] == 'cutadapt' for args in calls))
        dada2 = next(args for args in calls if args[0] == 'dada2')
        self.assertEqual(dada2[dada2.index('--p-trunc-len-f') + 1], '250')
        self.assertEqual(snapshot(run), before)
        self.assertEqual(snapshot(self.raw), raw_before)
        self.assertNotEqual(artifact_path(destination, 'demux.qza').stat().st_ino,
                            artifact_path(run, 'demux.qza').stat().st_ino)
        self.assertFalse((destination / 'preprocess_complete.json').exists())
        self.assertEqual(json.loads((destination / 'task_state.json').read_text())['action'], 'preprocess')

    def test_taxonomy_failure_restarts_at_dada2_not_raw_import(self):
        run = self.failed_preprocess('feature-classifier')
        self.assertTrue((run / '03_dada2/dada2-table.qza').is_file())
        plan = recovery_plan(run, app.module_api())
        job = {'_task': {}}
        action, config = prepare_retry(job, dict(plan['config'], recovery_run=str(run), trunc_r='210'), app.module_api())
        self.assertEqual(action, 'preprocess')
        destination = Path(job['run'])
        self.assertTrue(artifact_path(destination, 'demux.qza').is_file())
        self.assertFalse(artifact_path(destination, 'dada2-table.qza').exists())
        self.assertEqual(config['trunc_r'], '210')

    def test_changed_inputs_primers_mapping_and_checkpoint_are_rejected_before_copy(self):
        for changed in ('raw', 'raw_preserved_stat', 'primer', 'metadata', 'artifact'):
            with self.subTest(changed=changed):
                run = self.failed_preprocess()
                plan = recovery_plan(run, app.module_api())
                request = dict(plan['config'], recovery_run=str(run))
                if changed in ('raw', 'raw_preserved_stat'):
                    filename = self.raw / '001_R1.fastq.gz'
                    original_stat = filename.stat()
                    content = filename.read_bytes()
                    filename.write_bytes(content + b'changed' if changed == 'raw' else content[:-1] + bytes([content[-1] ^ 1]))
                    if changed == 'raw_preserved_stat':
                        os.utime(filename, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))
                elif changed == 'primer':
                    request.update(data_state='raw', primer_mode='341F')
                elif changed == 'metadata':
                    Path(self.config['metadata']).write_text('SampleID,Group\n001,OTHER\n')
                else:
                    artifact_path(run, 'demux.qza').write_bytes(b'changed artifact')
                directories = set(Path(self.config['output']).iterdir())
                with self.assertRaises((ValueError, OSError)):
                    prepare_retry({'_task': {}}, request, app.module_api())
                self.assertEqual(set(Path(self.config['output']).iterdir()), directories)
                # Restore real raw input for the next independent fixture.
                for mate in (1, 2):
                    fastq(self.raw / f'001_R{mate}.fastq.gz', mate=mate)
                Path(self.config['metadata']).write_text('SampleID,Group\n001,CON\n')

    def test_abundance_retry_forks_completed_outputs_and_never_copies_partial_export(self):
        run = self.failed_preprocess('feature-classifier')
        # Synthetic completed prerequisites; the retry itself still invokes the
        # real resume_abundance orchestration, with only its exporter replaced.
        for logical in ('filtered-table.qza', 'taxonomy.qza', 'rooted-tree.qza', 'taxa-barplot.qzv', 'alpha-rarefaction.qzv'):
            filename = artifact_path(run, logical)
            filename.parent.mkdir(parents=True, exist_ok=True)
            filename.write_bytes(('complete ' + logical).encode())
        table = artifact_path(run, 'exported-table/feature-table.tsv')
        table.parent.mkdir(parents=True, exist_ok=True)
        table.write_text('#OTU ID\t001\nASV1\t125\n')
        for stage in ('04_taxonomy_filter', '05_phylogeny', '06_rarefaction'):
            stage_report(run, stage, stage)
        seal_checkpoint(run, 'abundance')
        partial = run / '10_abundance/partial.qza'
        partial.parent.mkdir()
        partial.write_bytes(b'incomplete')
        before, raw_before = snapshot(run), snapshot(self.raw)
        plan = recovery_plan(run, app.module_api())
        self.assertEqual(plan['strategy'], 'abundance')
        retry = {'id': 'abundance-retry', 'action': 'recover', 'status': 'running'}
        with patch('abundance_tables.run_abundance', return_value={'summary': 'synthetic exporter'}), patch.object(app, 'q') as qiime:
            app.worker(retry, 'recover', dict(plan['config'], recovery_run=str(run)))
        self.assertEqual(retry['status'], 'done', retry.get('error'))
        self.assertEqual(retry['recovery_action'], 'resume-abundance')
        destination = Path(retry['run'])
        self.assertTrue((destination / 'preprocess_complete.json').is_file())
        self.assertFalse((destination / '10_abundance/partial.qza').exists())
        self.assertFalse(artifact_path(destination, 'dada2-table.qza').exists())
        qiime.assert_not_called()
        self.assertEqual(snapshot(run), before)
        self.assertEqual(snapshot(self.raw), raw_before)

    def child(self, action='diversity', cached=False):
        parent = self.root / 'completed-parent'
        parent.mkdir()
        app.marker(parent)
        initialize_layout(parent, 'preprocess')
        saved = dict(self.config, samples=app.discover(self.config))
        (parent / 'parameters.json').write_text(json.dumps(saved))
        (parent / 'preprocess_complete.json').write_text('{}')
        for logical in ('filtered-table.qza', 'taxonomy.qza', 'rooted-tree.qza'):
            artifact_path(parent, logical).write_bytes(b'completed synthetic parent')
        child = parent / '07_diversity' / 'failed-child'
        child.mkdir(parents=True)
        app.marker(child)
        initialize_layout(child, 'diversity')
        params = {'preprocess': str(parent), 'config': saved, 'depth': 100, 'permutations': 999}
        (child / 'parameters.json').write_text(json.dumps(params))
        config = dict(saved, run=str(parent), depth='100', permutations='999')
        if cached:
            (child / 'qiime_statistics_complete.json').write_text('{}')
            artifact_path(child, 'core-metrics-results').mkdir(exist_ok=True)
            artifact_path(child, 'core-metrics-results/rarefied_table.qza').write_bytes(b'saved random table')
        self.task(child, action, config)
        return parent, child

    def test_downstream_retries_only_analysis_from_saved_parent_and_rejects_changed_parent(self):
        parent, child = self.child(action='differential')
        plan = recovery_plan(child, app.module_api())
        self.assertEqual(plan['strategy'], 'analysis')
        action, config = prepare_retry({'_task': {}}, dict(plan['config'], recovery_run=str(child), run='/other', r_env='fixed-r'), app.module_api())
        self.assertEqual(action, 'differential')
        self.assertEqual(config['run'], str(parent))
        self.assertEqual(config['r_env'], 'fixed-r')
        artifact_path(parent, 'filtered-table.qza').write_bytes(b'changed parent')
        with self.assertRaisesRegex(ValueError, '检查点文件已改变'):
            prepare_retry({'_task': {}}, dict(plan['config'], recovery_run=str(child)), app.module_api())

    def test_diversity_figure_failure_uses_original_cache_and_changed_depth_restarts_analysis(self):
        parent, child = self.child(cached=True)
        before = snapshot(child)
        plan = recovery_plan(child, app.module_api())
        self.assertEqual(plan['strategy'], 'diversity-render')
        action, config = prepare_retry({'_task': {}}, dict(plan['config'], recovery_run=str(child), diversity_run='/other'), app.module_api())
        self.assertEqual(action, 'diversity-render')
        self.assertEqual(config['diversity_run'], str(child))
        action, config = prepare_retry({'_task': {}}, dict(plan['config'], recovery_run=str(child), depth='200'), app.module_api())
        self.assertEqual(action, 'diversity')
        self.assertEqual(config['run'], str(parent))
        self.assertEqual(snapshot(child), before)
        artifact_path(child, 'core-metrics-results/rarefied_table.qza').write_bytes(b'changed random table')
        with self.assertRaisesRegex(ValueError, '检查点文件已改变'):
            prepare_retry({'_task': {}}, dict(plan['config'], recovery_run=str(child)), app.module_api())

    def test_stale_running_record_becomes_interrupted_and_checks_have_no_recovery(self):
        run = self.failed_preprocess()
        self.task(run, 'preprocess', self.config, status='running')
        before = snapshot(run)
        plan = recovery_plan(run, app.module_api())
        self.assertEqual(plan['status'], 'interrupted')
        self.assertTrue(plan['recoverable'])
        active = [{'status': 'running', 'run': str(run)}]
        self.assertFalse(recovery_plan(run, app.module_api(), active)['recoverable'])
        self.assertEqual(snapshot(run), before)
        for action in ('inspect', 'preflight'):
            self.task(run, action, self.config)
            self.assertFalse(recovery_plan(run, app.module_api())['recoverable'])

    def test_external_import_retry_keeps_explicit_file_request_and_external_child_guard(self):
        run = self.root / 'failed-external'
        run.mkdir()
        app.marker(run)
        config = dict(self.config, external_module='differential', external_table='/saved/counts.tsv',
                      external_taxonomy='/saved/taxonomy.tsv', external_counts_confirmed=True)
        (run / 'parameters.json').write_text(json.dumps(config))
        self.task(run, 'external-inputs', config)
        before = snapshot(run)
        plan = recovery_plan(run, app.module_api())
        self.assertTrue(plan['recoverable'])
        action, retried = prepare_retry({'_task': {}}, dict(recovery_run=str(run), external_table='',
                                                           external_module='', metadata='', threads='3'), app.module_api())
        self.assertEqual(action, 'external-inputs')
        self.assertEqual(retried['external_table'], '/saved/counts.tsv')
        self.assertEqual(retried['external_module'], 'differential')
        self.assertEqual(retried['metadata'], self.config['metadata'])
        self.assertEqual(retried['threads'], '3')
        self.assertEqual(before, snapshot(run))
        parent, child = self.child(action='differential')
        (parent / 'preprocess_complete.json').unlink()
        (parent / 'external_inputs_complete.json').write_text('{}')
        params = json.loads((parent / 'parameters.json').read_text())
        params['supported_modules'] = ['differential']
        (parent / 'parameters.json').write_text(json.dumps(params))
        self.assertTrue(recovery_plan(child, app.module_api())['recoverable'])
        self.task(child, 'functional', dict(self.config, run=str(parent)))
        self.assertFalse(recovery_plan(child, app.module_api())['recoverable'])

    def test_task_api_is_small_and_omits_private_records(self):
        server = app.ThreadingHTTPServer(('127.0.0.1', 0), app.Handler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            jobs = {str(i): {'id': str(i), 'status': 'done', '_execution': {'sources': 'private'},
                             '_task': {'config': 'private'}, 'log': '/private'} for i in range(25)}
            jobs['active'] = {'id': 'active', 'status': 'running', 'action': 'recover', 'recovery_action': 'preprocess'}
            with patch.dict(app.JOBS, jobs, clear=True):
                conn = http.client.HTTPConnection('127.0.0.1', server.server_port)
                conn.request('GET', '/api/tasks', headers={'X-Workstation-Token': app.TOKEN})
                response = conn.getresponse()
                result = json.loads(response.read())
                self.assertEqual(response.status, 200)
                self.assertEqual(len(result['jobs']), 21)
                self.assertEqual(result['jobs'][0]['id'], 'active')
                self.assertNotIn('private', json.dumps(result))
                conn.close()
        finally:
            server.shutdown()
            server.server_close()
            thread.join(timeout=3)


if __name__ == '__main__':
    unittest.main()
