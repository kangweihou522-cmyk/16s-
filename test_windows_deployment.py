"""Deployment transactions, raw-data preservation and authenticated service ownership."""
from pathlib import Path
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import os
import subprocess
import sys
import threading

import contextlib
import io
import tempfile
import unittest
from unittest.mock import patch

from windows import deploy, station_runtime


def make_station(tmp_path, monkeypatch):
    source = tmp_path/'downloaded source'
    for name in deploy.REQUIRED:
        path = source/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('# fixture '+name, encoding='utf-8')
    root = tmp_path/"D 盘/16S O'Brien $(touch WRONG)"
    conda = tmp_path/'existing conda'
    conda.write_text('#!/bin/sh\nexit 0\n')
    conda.chmod(0o755)
    monkeypatch.setattr(deploy, 'find_conda', lambda preferred=None: str(conda))
    monkeypatch.setattr(station_runtime, 'probe', lambda *args: dict(state='offline', health={}))
    monkeypatch.setattr(station_runtime, 'stop', lambda *args: dict(state='offline', health={}))
    return source, root


def stage(source, root):
    return deploy.prepare(source, root, r"D:\实验 16S\工作站 O'Brien $(touch WRONG)",
                          'Ubuntu-24.04', 8767, 'researcher')


def install(source, root):
    plan = stage(source, root)
    config = deploy.activate(root, plan['plan_file'])
    deploy.complete(root, plan['plan_file'])
    return config


@contextlib.contextmanager
def make_health_server(tmp_path):
    root = tmp_path/'station'; program = root/'程序/release'; program.mkdir(parents=True)
    payload = dict(service_id='16s-workstation', install_root=str(root),
                   program_dir=str(program), pid=os.getpid(), instance_id='instance-one', busy=False)
    state = dict(payload=payload, posts=[], status=200, close=False)
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            if state.get('redirect'):
                self.send_response(302)
                self.send_header('Location', state['redirect'])
                self.end_headers()
                return
            self.send_response(200); self.end_headers()
            self.wfile.write(json.dumps(state['payload']).encode())
        def do_POST(self):
            state['posts'].append((self.path, self.headers.get('X-Workstation-Token')))
            self.send_response(state['status']); self.end_headers(); self.wfile.write(b'{}')
            if state['close']:
                def shutdown():
                    server.shutdown(); server.server_close()
                threading.Thread(target=shutdown, daemon=True).start()
        def log_message(self, *args):
            pass
    server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True); thread.start()
    config_dir = root/'配置'; config_dir.mkdir()
    record = dict(payload, token='root-owned-token', port=server.server_port)
    record.pop('service_id'); record.pop('busy')
    (config_dir/'service.json').write_text(json.dumps(record))
    yield root, server.server_port, state
    server.shutdown(); server.server_close(); thread.join()


class _Patches:
    def __init__(self, stack):
        self.stack = stack

    def setattr(self, target, name, value):
        self.stack.enter_context(patch.object(target, name, value))

    def setenv(self, name, value):
        self.stack.enter_context(patch.dict(os.environ, {name: value}))


class WindowsDeploymentTests(unittest.TestCase):
    def setUp(self):
        self.stack = contextlib.ExitStack()
        self.addCleanup(self.stack.close)
        self.tmp_path = Path(self.stack.enter_context(tempfile.TemporaryDirectory()))
        self.monkeypatch = _Patches(self.stack)
        # Runtime ownership tests use a real HTTP server, without deployment mocks.
        if self._testMethodName.startswith(('test_health_', 'test_stale_service_', 'test_authenticated_stop_')):
            self.health_server = self.stack.enter_context(make_health_server(self.tmp_path))
        elif not self._testMethodName.startswith('test_conda_detection_'):
            self.station = make_station(self.tmp_path, self.monkeypatch)

    def test_allowlist_and_layout_preserve_private_data(self):
        station = self.station
        source, root = station
        for name in ['sample.fastq.gz', 'customer.csv', '.env', 'runs/job/status.json',
                     'updates/old.bundle', 'assets/private.FASTQ', 'windows/runs/job.fastq',
                     '.git/objects/file', 'assets/font.woff2', 'fonts/Chinese.ttf',
                     'assets/view.js', 'future_analysis.py']:
            path = source/name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('# source')
        config = install(source, root)
        program = Path(config['program_wsl'])
        for name in ['sample.fastq.gz', 'customer.csv', '.env', 'runs', 'updates',
                     'assets/private.FASTQ', 'windows/runs', '.git']:
            assert not (program/name).exists()
        for name in ['assets/font.woff2', 'fonts/Chinese.ttf', 'assets/view.js', 'future_analysis.py']:
            assert (program/name).is_file()
        assert all((root/name).is_dir() for name in deploy.DIRECTORIES)
        assert config['service_id'] == '16s-workstation' and config['schema'] == 1
        assert config['wsl_user'] == 'researcher' and config['envs'] == deploy.DEFAULT_ENVS
        assert config['runs_wsl'] == str(root/'分析结果')
        assert config['version'].startswith('sha256:')
        assert (root/'启动工作站.bat').read_bytes().startswith(b'@echo off\r\n')
        assert b'Station.ps1' in (root/'停止工作站.bat').read_bytes()
        assert not (root.parent/'WRONG').exists()

    def test_optional_advanced_installer_targets_only_its_independent_environment(self):
        source, root = self.station
        config = install(source, root)
        before = (root/'配置/workstation.json').read_bytes()
        for exists in (False, True):
            envs = ['/existing/envs/'+name for name in config['envs'].values()]
            if exists:
                envs.append('/existing/envs/16s-workstation-advanced')
            def external(args, **kwargs):
                if args[1:4] == ['env', 'list', '--json']:
                    return subprocess.CompletedProcess(args, 0, json.dumps({'envs': envs}), '')
                return subprocess.CompletedProcess(args, 0)
            with self.subTest(exists=exists), patch.object(deploy.subprocess, 'run', side_effect=external) as command:
                result = deploy.install_advanced(root)
            assert result == {'environment': '16s-workstation-advanced', 'available': True}
            calls = [call.args[0] for call in command.call_args_list]
            assert len(calls) == 3
            assert calls[1][1] == ('install' if exists else 'create')
            assert calls[1][calls[1].index('-n')+1] == '16s-workstation-advanced'
            assert '--override-channels' in calls[1] and '--strict-channel-priority' in calls[1]
            assert 'r-vegan>=2.6,<3' in calls[1] and 'r-ppcor' in calls[1]
            assert calls[2][calls[2].index('-n')+1] == '16s-workstation-advanced'
            assert not any(name in calls[1] for name in config['envs'].values())
        assert (root/'配置/workstation.json').read_bytes() == before

    def test_optional_advanced_installer_refuses_busy_service_before_any_command(self):
        source, root = self.station
        install(source, root)
        with patch.object(deploy, 'ensure_idle', side_effect=ValueError('任务仍在运行')), \
                patch.object(deploy.subprocess, 'run') as command:
            with self.assertRaisesRegex(ValueError, '仍在运行'):
                deploy.install_advanced(root)
        command.assert_not_called()


    def test_updates_preserve_raw_results_and_changed_releases(self):
        station = self.station
        source, root = station
        first = install(source, root)
        program = Path(first['program_wsl'])
        raw = root/'原始数据/reads.fastq.gz'; raw.write_bytes(b'USER RAW')
        result = root/'分析结果/result.csv'; result.write_text('USER RESULT')
        # Offline status files are not deployment locks and raw directories are not scanned.
        job = root/'分析结果/old'; job.mkdir(); (job/'status.json').write_text('invalid JSON')
        assert install(source, root)['program_wsl'] == str(program)
        assert install(program, root)['program_wsl'] == str(program)
        (program/'app.py').write_text('MY LOCAL EDIT')
        second = install(source, root)
        assert second['program_wsl'] != str(program)
        assert (program/'app.py').read_text() == 'MY LOCAL EDIT'
        source.joinpath('app.py').write_text('# next source')
        third = install(source, root)
        assert third['program_wsl'] not in {first['program_wsl'], second['program_wsl']}
        assert raw.read_bytes() == b'USER RAW' and result.read_text() == 'USER RESULT'


    def test_edited_root_wrappers_refuse_activation(self):
        station = self.station
        name = 'Station.ps1'
        source, root = station
        first = install(source, root)
        wrapper = root/name; wrapper.write_text('USER WRAPPER')
        source.joinpath('app.py').write_text('# new source')
        plan = stage(source, root)
        with self.assertRaisesRegex(ValueError, '未覆盖'):
            deploy.activate(root, plan['plan_file'])
        assert wrapper.read_text() == 'USER WRAPPER'
        assert deploy.read_config(root) == first


    def test_malformed_source_and_config_fail_without_overwriting(self):
        source, root = self.station
        source.joinpath('app.py').write_text('def broken(:')
        with self.assertRaisesRegex(ValueError, '源码格式错误'):
            stage(source, root)
        assert not root.exists()
        source.joinpath('app.py').write_text('# valid source')
        install(source, root)
        config_file = root/'配置/workstation.json'
        config_file.write_text('{invalid JSON')
        with self.assertRaisesRegex(ValueError, '配置损坏'):
            stage(source, root)
        assert config_file.read_text() == '{invalid JSON'


    def test_tampered_plan_cannot_select_external_program(self):
        source, root = self.station
        first = install(source, root)
        plan = stage(source, root)
        path = Path(plan['plan_file'])
        value = json.loads(path.read_text())
        value['config']['program_wsl'] = str(root.parent/'external program')
        path.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError, '路径无效'):
            deploy.activate(root, path)
        assert deploy.read_config(root) == first


    def test_staged_hash_tampering_and_concurrent_activation_are_blocked(self):
        station = self.station
        source, root = station
        first = install(source, root)
        source.joinpath('app.py').write_text('# second')
        damaged = stage(source, root)
        Path(damaged['config']['program_wsl']).joinpath('app.py').write_text('TAMPERED')
        with self.assertRaisesRegex(ValueError, '改变'):
            deploy.activate(root, damaged['plan_file'])
        assert deploy.read_config(root) == first
        stale = stage(source, root)
        source.joinpath('app.py').write_text('# third')
        current = install(source, root)
        with self.assertRaisesRegex(ValueError, '另一部署'):
            deploy.activate(root, stale['plan_file'])
        assert deploy.read_config(root) == current


    def test_write_failure_rolls_back_wrappers_and_preserves_config(self):
        station = self.station
        monkeypatch = self.monkeypatch
        source, root = station
        first = install(source, root)
        originals = {name: (root/name).read_bytes() for name in first['managed_files']}
        source.joinpath('windows/Station.ps1').write_text('new PS helper')
        plan = stage(source, root)
        original_write = deploy.atomic_write
        def failing_write(path, content):
            if Path(path) == root/'配置/workstation.json':
                raise OSError('simulated storage failure')
            original_write(path, content)
        monkeypatch.setattr(deploy, 'atomic_write', failing_write)
        with self.assertRaisesRegex(OSError, 'storage failure'):
            deploy.activate(root, plan['plan_file'])
        assert deploy.read_config(root) == first
        assert all((root/name).read_bytes() == content for name, content in originals.items())
        assert not (root/'配置/.deployment-lock').exists()


    def test_start_failure_rollback_restores_previous_config_and_wrappers(self):
        station = self.station
        source, root = station
        first = install(source, root)
        original_ps = (root/'Station.ps1').read_bytes()
        raw = root/'原始数据/keep.fastq'; raw.write_text('KEEP')
        source.joinpath('windows/Station.ps1').write_text('new helper')
        source.joinpath('app.py').write_text('raise RuntimeError("broken startup")')
        plan = stage(source, root)
        second = deploy.activate(root, plan['plan_file'])
        assert second != first
        # Recovery also works if a failed candidate writes into its source directory.
        Path(second['program_wsl']).joinpath('unexpected-startup.log').write_text('failure')
        assert deploy.rollback(root, plan['plan_file']) == first
        assert deploy.read_config(root) == first
        assert (root/'Station.ps1').read_bytes() == original_ps
        assert raw.read_text() == 'KEEP'
        assert Path(second['program_wsl']).is_dir()


    def test_first_install_rollback_removes_only_managed_files(self):
        station = self.station
        source, root = station
        plan = stage(source, root)
        config = deploy.activate(root, plan['plan_file'])
        raw = root/'原始数据/keep.fastq'; raw.write_text('KEEP')
        assert deploy.rollback(root, plan['plan_file']) == {}
        assert deploy.read_config(root) == {}
        assert all(not (root/name).exists() for name in config['managed_files'])
        assert raw.read_text() == 'KEEP'


    def test_incomplete_source_and_invalid_paths_do_not_create_station(self):
        station = self.station
        source, root = station
        source.joinpath('deployment_context.py').unlink()
        with self.assertRaisesRegex(ValueError, '源码不完整'):
            stage(source, root)
        assert not root.exists()
        for path in ['D:\\', 'relative', 'D:\\bad\nname', 'D:\\bad*name', 'D:\\..\\station', r'D:\CON']:
            with self.assertRaisesRegex(ValueError, '路径|目录'):
                deploy.prepare(source, root, path, 'Ubuntu', 8767)
        for port in [80, 65536, True]:
            with self.assertRaisesRegex(ValueError, '端口'):
                deploy.prepare(source, root, r'D:\station', 'Ubuntu', port)
        with self.assertRaisesRegex(ValueError, '源码目录'):
            deploy.prepare(source, source/'nested install', r'D:\station', 'Ubuntu', 8767)


    def test_symlinked_source_and_destination_are_never_copied(self):
        station = self.station
        tmp_path = self.tmp_path
        source, root = station
        private = tmp_path/'private.py'; private.write_text('secret')
        (source/'assets').mkdir()
        (source/'assets/leak.py').symlink_to(private)
        with self.assertRaisesRegex(ValueError, '符号链接'):
            stage(source, root)
        assert not root.exists()
        root.parent.mkdir(parents=True)
        root.symlink_to(tmp_path/'other directory')
        with self.assertRaisesRegex(ValueError, '符号链接'):
            stage(source, root)


    def test_online_busy_or_foreign_port_blocks_update(self):
        station = self.station
        monkeypatch = self.monkeypatch
        source, root = station
        first = install(source, root)
        for state in ['busy', 'foreign']:
            monkeypatch.setattr(station_runtime, 'probe', lambda *args: dict(
                state='running' if state == 'busy' else 'occupied', health={'busy': True}))
            with self.assertRaisesRegex(ValueError, '未结束|其他程序'):
                stage(source, root)
        assert deploy.read_config(root) == first


    def test_candidate_check_reports_missing_envs_without_blocking_activation(self):
        station = self.station
        monkeypatch = self.monkeypatch
        source, root = station
        plan = stage(source, root)
        monkeypatch.setattr(deploy.subprocess, 'run', lambda *args, **kwargs:
                            subprocess.CompletedProcess(args, 0, json.dumps({'envs': [
                                '/existing/envs/qiime2-amplicon-2025.7', '/existing/envs/16s-workstation-r']})))
        report = deploy.check(root, plan['plan_file'])
        assert report['conda_available'] and not report['ready']
        assert report['envs']['qiime']['available'] and report['envs']['r']['available']
        assert report['missing_envs'] == ['16s-workstation-ancom', 'lefse', 'picrust2']
        assert deploy.activate(root, plan['plan_file'])['envs'] == deploy.DEFAULT_ENVS


    def test_conda_detection_reuses_hidden_miniconda_and_never_installs(self):
        tmp_path = self.tmp_path
        monkeypatch = self.monkeypatch
        home = tmp_path/'home'
        conda = home/'.miniconda3/bin/conda'
        conda.parent.mkdir(parents=True); conda.write_text('#!/bin/sh\nexit 0\n'); conda.chmod(0o755)
        monkeypatch.setattr(Path, 'home', lambda: home)
        monkeypatch.setattr(deploy.shutil, 'which', lambda name: None)
        assert deploy.find_conda() == str(conda)
        conda.unlink()
        with self.assertRaisesRegex(ValueError, '找不到已有的 Conda'):
            deploy.find_conda()


    def test_start_appends_console_log_and_sets_station_root(self):
        station = self.station
        source, root = station
        source.joinpath('app.py').write_text(
            'import os,sys\nprint("station="+os.environ["STATION_ROOT"])\n'
            'print("args="+repr(sys.argv[1:]))\nprint("bytecode="+str(sys.dont_write_bytecode))\n')
        install(source, root)
        log = root/'日志/server.log'; log.write_bytes(b'PREVIOUS LOG\n')
        console = io.StringIO()
        with contextlib.redirect_stdout(console):
            assert deploy.start(root) == 0
        text = log.read_text()
        assert text.startswith('PREVIOUS LOG\n') and 'station='+str(root) in text
        assert "['--port', '8767', '--no-browser']" in text and 'bytecode=True' in text
        assert 'station='+str(root) in console.getvalue()


    def test_health_rejects_foreign_root_and_legacy_service(self):
        health_server = self.health_server
        root, port, state = health_server
        assert station_runtime.probe(root, port)['state'] == 'running'
        state['payload']['install_root'] = str(root.parent/'another station')
        assert station_runtime.probe(root, port)['state'] == 'occupied'
        with self.assertRaisesRegex(ValueError, '其他程序'):
            station_runtime.stop(root, port)
        state['payload']['install_root'] = str(root)
        state['payload'].pop('instance_id')
        assert station_runtime.probe(root, port)['state'] == 'occupied'
        state['redirect'] = 'http://127.0.0.1:1/remote-health'
        assert station_runtime.probe(root, port)['state'] == 'occupied'
        assert not state['posts']


    def test_stale_service_record_and_busy_race_never_stop_service(self):
        health_server = self.health_server
        root, port, state = health_server
        state['payload']['instance_id'] = 'reused-pid-new-instance'
        with self.assertRaisesRegex(ValueError, '进程与记录不一致'):
            station_runtime.stop(root, port)
        assert not state['posts']
        state['payload']['instance_id'] = 'instance-one'
        state['payload']['busy'] = True
        with self.assertRaisesRegex(ValueError, '未结束'):
            station_runtime.stop(root, port)
        assert not state['posts']
        state['payload']['busy'] = False
        state['status'] = 409
        with self.assertRaisesRegex(ValueError, '未结束'):
            station_runtime.stop(root, port)
        assert state['posts'] == [('/api/station-stop', 'root-owned-token')]


    def test_authenticated_stop_waits_for_owned_service_to_exit(self):
        health_server = self.health_server
        monkeypatch = self.monkeypatch
        root, port, state = health_server
        monkeypatch.setenv('http_proxy', 'http://127.0.0.1:1')
        monkeypatch.setenv('HTTP_PROXY', 'http://127.0.0.1:1')
        monkeypatch.setenv('NO_PROXY', '')
        state['close'] = True
        assert station_runtime.stop(root, port)['state'] == 'offline'
        assert state['posts'] == [('/api/station-stop', 'root-owned-token')]


if __name__ == "__main__":
    unittest.main()
