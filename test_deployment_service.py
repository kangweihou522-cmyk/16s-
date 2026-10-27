"""Exercise deployment defaults and the real HTTP lifecycle without analysis."""
import json
import re
import tempfile
import threading
import unittest
import urllib.error
import urllib.request
from pathlib import Path
from unittest.mock import patch

import app
from deployment_context import frontend_defaults, load_context


class DeploymentContextTests(unittest.TestCase):
    def test_legacy_launch_and_complete_environment_mapping(self):
        self.assertEqual(load_context(app.ROOT, ''), {})
        config = {'runs_wsl': '/mnt/d/16S/分析结果', 'conda_wsl': '/root/miniconda3/bin/conda',
                  'envs': {'qiime': 'q-env', 'r': 'r-env', 'ancom': 'a-env',
                           'lefse': 'l-env', 'picrust2': 'p-env'}}
        self.assertEqual(frontend_defaults(config), {
            'output': '/mnt/d/16S/分析结果', 'conda': '/root/miniconda3/bin/conda',
            'env': 'q-env', 'r_env': 'r-env', 'ancom_env': 'a-env',
            'lefse_env': 'l-env', 'picrust_env': 'p-env'})

    def test_reject_wrong_service_program_and_result_location(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / '配置').mkdir()
            config = {'schema': 1, 'service_id': '16s-workstation', 'root_wsl': str(root),
                      'program_wsl': str(app.ROOT), 'runs_wsl': str(root / '分析结果')}
            filename = root / '配置/workstation.json'
            filename.write_text(json.dumps(config))
            self.assertEqual(load_context(app.ROOT, root), config)
            for key, value in [('service_id', 'mouse-bulk-rnaseq'), ('root_wsl', '/other'),
                               ('program_wsl', '/old-version'), ('runs_wsl', '/other-results')]:
                filename.write_text(json.dumps(dict(config, **{key: value})))
                with self.subTest(key=key), self.assertRaises(ValueError):
                    load_context(app.ROOT, root)


class DeploymentServiceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name)
        self.config = {'schema': 1, 'service_id': '16s-workstation', 'root_wsl': str(root),
                       'root_windows': 'D:\\16S分析工作站', 'program_wsl': str(app.ROOT),
                       'runs_wsl': str(root / '分析结果'), 'runs_windows': 'D:\\16S分析工作站\\分析结果',
                       'envs': {'qiime': 'custom-q', 'r': 'custom-r'}, 'distro': 'Ubuntu-24.04'}
        for target, value in [('DEPLOYMENT', self.config), ('JOBS', {}), ('STOPPING', False)]:
            manager = patch.object(app, target, value)
            manager.start()
            self.addCleanup(manager.stop)
        self.server = app.ThreadingHTTPServer(('127.0.0.1', 0), app.Handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.addCleanup(self.close)
        self.base = 'http://127.0.0.1:' + str(self.server.server_port)
        self.opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))

    def close(self):
        self.server.shutdown()
        self.thread.join(timeout=3)
        self.server.server_close()

    def request(self, route, *, post=False, token=False, host=None):
        headers = {'X-Workstation-Token': app.TOKEN} if token else {}
        if host:
            headers['Host'] = host
        request = urllib.request.Request(self.base + route, headers=headers,
                                         data=b'{}' if post else None)
        try:
            with self.opener.open(request, timeout=3) as response:
                return response.status, response.read()
        except urllib.error.HTTPError as error:
            return error.code, error.read()

    def test_health_is_identifiable_readonly_and_host_restricted(self):
        status, body = self.request('/api/health')
        health = json.loads(body)
        self.assertEqual(status, 200)
        self.assertEqual(health['service_id'], '16s-workstation')
        self.assertEqual(health['install_root'], self.config['root_wsl'])
        self.assertFalse(health['busy'])
        self.assertNotIn('token', health)
        self.assertEqual(self.request('/api/health', host='untrusted.example')[0], 403)
        self.assertEqual(self.request('/api/deployment')[0], 403)

    def test_html_defaults_preserve_valid_json_and_configured_paths(self):
        self.config['envs']['qiime'] = 'q<&special'
        status, body = self.request('/')
        self.assertEqual(status, 200)
        html = body.decode()
        payload = re.search(r'id="deployment-defaults">(.*?)</script>', html, re.S)[1]
        defaults = json.loads(payload)
        self.assertEqual(defaults['env'], 'q<&special')
        self.assertEqual(defaults['output'], self.config['runs_wsl'])
        self.assertNotIn('<', payload)
        status, body = self.request('/api/deployment', token=True)
        self.assertEqual(status, 200)
        self.assertEqual(json.loads(body)['root_windows'], self.config['root_windows'])

    def test_stop_needs_token_and_rejects_busy_without_cancelling_task(self):
        self.assertEqual(self.request('/api/station-stop', post=True)[0], 403)
        app.JOBS['real-analysis'] = {'status': 'running'}
        status, body = self.request('/api/station-stop', post=True, token=True)
        self.assertEqual(status, 409)
        self.assertIn('仍在运行', json.loads(body)['error'])
        self.assertEqual(app.JOBS['real-analysis']['status'], 'running')
        self.assertFalse(app.STOPPING)
        self.assertTrue(json.loads(self.request('/api/health')[1])['busy'])

    def test_idle_stop_prevents_new_jobs_atomically(self):
        status, body = self.request('/api/station-stop', post=True, token=True)
        self.assertEqual(status, 200)
        self.assertTrue(json.loads(body)['stopping'])
        self.assertTrue(app.STOPPING)
        with self.assertRaisesRegex(ValueError, '正在停止或更新'):
            app.launch('inspect', {})
        self.thread.join(timeout=3)
        self.assertFalse(self.thread.is_alive())

    def test_finalizing_task_is_busy_until_archives_are_saved(self):
        app.JOBS['saving-code'] = {'status': 'done', '_execution': {}}
        self.assertTrue(json.loads(self.request('/api/health')[1])['busy'])
        self.assertEqual(self.request('/api/station-stop', post=True, token=True)[0], 409)
        with self.assertRaisesRegex(ValueError, '已有任务运行中'):
            app.launch('inspect', {})
        app.JOBS['saving-code'].pop('_execution')
        app.JOBS['saving-code']['_task'] = {}
        self.assertTrue(json.loads(self.request('/api/health')[1])['busy'])
        app.JOBS['saving-code'].pop('_task')
        self.assertFalse(json.loads(self.request('/api/health')[1])['busy'])

    def test_initial_job_logs_are_outside_the_immutable_program(self):
        with patch.object(app.threading, 'Thread'):
            job_id = app.launch('preflight', {})
        log = Path(app.JOBS[job_id]['log'])
        self.assertTrue(log.is_relative_to(Path(self.config['root_wsl']) / '日志/tasks'))
        self.assertFalse(log.is_relative_to(app.ROOT))
        self.assertTrue(log.parent.is_dir())

    def test_failed_inventory_reports_failure_and_releases_the_busy_guard(self):
        job = {'id': 'save-error', 'action': 'preflight', 'status': 'running',
               'run': self.config['root_wsl']}
        app.JOBS[job['id']] = job
        with patch.object(app, 'preflight'), patch.object(app, 'finalize_execution'), \
                patch.object(app, 'inventory', side_effect=OSError('disk write failed')):
            app.worker(job, 'preflight', {})
        self.assertEqual(job['status'], 'failed')
        self.assertIn('保存任务清单失败', job['error'])
        self.assertNotIn('_execution', job)
        self.assertFalse(app.jobs_busy())


if __name__ == '__main__':
    unittest.main()
