"""Exercise text archives with real subprocesses; the fake CLI is not QIIME/R."""
import http.client
import json
import sys
import tempfile
import threading
import unittest
from pathlib import Path
from unittest.mock import patch

import app
from execution_record import prepare_execution, finalize_execution
from run_layout import initialize_layout


class ExecutionRecordTests(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup)
        self.root=Path(temporary.name);self.run=self.root/'run';self.run.mkdir()
        app.marker(self.run);initialize_layout(self.run,'diversity')
        self.cli=self.root/'fake-conda'
        self.cli.write_text('''#!/usr/bin/env python3
import sys
from pathlib import Path
if 'qiime' in sys.argv:
    option='--o-alpha-diversity'
    if option in sys.argv:
        Path(sys.argv[sys.argv.index(option)+1]).write_bytes(b'FAKE artifact')
    else:
        print('fake QIIME component available')
elif '--json' in sys.argv:
    print('{"envs":[]}')
''')
        self.cli.chmod(0o755)
        self.config={'conda':str(self.cli),'env':'fake-qiime','r_env':'fake-r','font':'Font with spaces'}

    def attach_run(self,job,parameters=None):
        job.update(run=str(self.run),log=str(self.run/'workflow.log'))
        (self.run/'parameters.json').write_text(json.dumps(parameters or {'trunc_f':275,'threads':1}),encoding='utf8')

    def report(self,job):
        return Path(job['execution_code_directory'])/'完整运行代码.txt'

    def test_worker_archives_actual_mapped_commands_effective_parameters_and_custom_r(self):
        assets=self.root/'custom plots';assets.mkdir()
        for name in ['qiime2_io.R','plot_alpha_original_style.R','plot_beta_original_style.R']:
            (assets/name).write_text('# custom executed asset '+name+'\n',encoding='utf8')
        script=self.root/'render.R';script.write_text('# rendering version BEFORE run\n',encoding='utf8')
        job={'id':'testactual','status':'running','action':'preprocess'}
        def workflow(job,config):
            self.attach_run(job,{'trunc_f':275,'threads':1})
            app.begin_stage(job,'02_exports','额外 α 指标',{'metric':'chao1'})
            app.q(job,self.config,'diversity','alpha','--i-table','core-metrics-results/rarefied_table.qza','--p-metric','chao1','--o-alpha-diversity','core-metrics-results/chao1_vector.qza')
            app.finish_stage(job,{'metric':'chao1'})
            app.begin_stage(job,'04_R_figures','R 绘图',{'font':self.config['font']})
            app.command(job,app.r_command(self.config,str(script),str(self.run),self.config['font'],str(self.root),str(assets)),self.run,process_env=app.r_execution_env())
            script.write_text('# rendering version AFTER run\n')
            (assets/'qiime2_io.R').write_text('# later modified asset\n')
            app.finish_stage(job)
        with patch.object(app,'preprocess',side_effect=workflow):app.worker(job,'preprocess',{'trunc_f':'auto'})
        self.assertEqual(job['status'],'done')
        report=self.report(job).read_text(encoding='utf8')
        self.assertIn('qiime diversity alpha --i-table 01_core_metrics/core-metrics-results/rarefied_table.qza',report)
        self.assertIn('--o-alpha-diversity 01_core_metrics/core-metrics-results/chao1_vector.qza',report)
        self.assertIn("'Font with spaces'",report)
        self.assertIn('状态：done；退出码：0',report)
        r_source=(self.report(job).parent/'R完整代码.txt').read_text(encoding='utf8')
        self.assertIn('BEFORE run',r_source);self.assertNotIn('AFTER run',r_source)
        self.assertIn('custom executed asset qiime2_io.R',r_source);self.assertNotIn('later modified asset',r_source)
        self.assertIn('def preprocess(',report);self.assertNotIn(app.TOKEN,report)
        parameters=json.loads((self.report(job).parent/'本次参数.txt').read_text())
        self.assertEqual(parameters['requested']['trunc_f'],'auto')
        self.assertEqual(parameters['saved_effective']['trunc_f'],275)
        status=json.loads((self.run/'job_status.json').read_text())
        self.assertNotIn('_execution',status)
        self.assertIn('09_execution_code', (self.run/'outputs.tsv').read_text())

    def test_failed_process_and_launch_failure_remain_distinguishable(self):
        for executable,expected in [([sys.executable,'-c','import sys;sys.exit(9)'],'failed'),
                                    ([str(self.root/'missing executable')],'launch_failed')]:
            with self.subTest(expected=expected):
                job={'status':'running','action':'preprocess'}
                def workflow(job,config):
                    self.attach_run(job)
                    app.begin_stage(job,'03_statistics','检验命令')
                    app.command(job,executable,self.run)
                with patch.object(app,'preprocess',side_effect=workflow):app.worker(job,'preprocess',{})
                self.assertEqual(job['status'],'failed')
                self.assertIn('状态：'+expected,self.report(job).read_text(encoding='utf8'))
                last=json.loads((self.run/'commands.jsonl').read_text().splitlines()[-1])
                self.assertEqual(last['status'],expected)
                self.assertEqual(last['exit_code'],9 if expected=='failed' else None)

    def test_stopped_task_saves_record_without_claiming_unrun_command(self):
        job={'status':'running','action':'preprocess'}
        def workflow(job,config):
            self.attach_run(job);job['cancel']=True
            app.command(job,[sys.executable,'-c','print("must not run")'],self.run)
        with patch.object(app,'preprocess',side_effect=workflow):app.worker(job,'preprocess',{})
        self.assertEqual(job['status'],'cancelled')
        self.assertIn('未通过命令执行器运行外部命令',self.report(job).read_text(encoding='utf8'))
        self.assertFalse((self.run/'commands.jsonl').exists())

    def test_reused_run_preserves_old_task_archive_and_does_not_relabel_commands(self):
        reports=[]
        for action,label in [('quality','original quality'),('preprocess','later preprocessing')]:
            job={'status':'done'};self.attach_run(job)
            prepare_execution(job,app.ROOT,action,{},app.RUNTIME)
            job['stage_title']=label
            app.command(job,[sys.executable,'-c','print("ok")'],self.run)
            finalize_execution(job);reports.append(self.report(job))
        self.assertNotEqual(*reports)
        original=reports[0].read_text(encoding='utf8');later=reports[1].read_text(encoding='utf8')
        self.assertIn('阶段：original quality',original);self.assertNotIn('阶段：later preprocessing',original)
        self.assertIn('阶段：later preprocessing',later);self.assertNotIn('阶段：original quality',later)

    def test_final_analysis_includes_parent_archive_and_honest_legacy_fallback(self):
        parent=self.root/'preprocess';parent.mkdir()
        old=(parent/'09_execution_code/preprocess_old');old.mkdir(parents=True)
        (old/'完整运行代码.txt').write_text('# actual historical preprocessing source\n',encoding='utf8')
        job={'status':'done'};self.attach_run(job,{'preprocess':str(parent),'depth':17100})
        prepare_execution(job,app.ROOT,'diversity',{'depth':17100},app.RUNTIME)
        finalize_execution(job)
        self.assertIn('actual historical preprocessing source',self.report(job).read_text(encoding='utf8'))
        self.assertTrue((self.report(job).parent/'关联预处理代码.txt').is_file())
        (old/'完整运行代码.txt').unlink()
        (parent/'commands.jsonl').write_text('{"args":["old-qiime"],"exit_code":0}\n')
        prepare_execution(job,app.ROOT,'diversity',{'depth':17100},app.RUNTIME)
        finalize_execution(job)
        text=self.report(job).read_text(encoding='utf8')
        self.assertIn('不能用当前源码冒充旧版实际代码',text);self.assertIn('old-qiime',text)

    def test_running_status_api_does_not_transmit_private_source_snapshots(self):
        job={'status':'running','log':str(self.run/'missing.log'),'_execution':{'sources':'private large snapshot'}}
        server=app.ThreadingHTTPServer(('127.0.0.1',0),app.Handler)
        thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
        try:
            with patch.dict(app.JOBS,{'source-test':job}):
                client=http.client.HTTPConnection('127.0.0.1',server.server_port,timeout=5)
                client.request('GET','/api/status?id=source-test',headers={'X-Workstation-Token':app.TOKEN})
                response=client.getresponse();payload=json.loads(response.read());client.close()
            self.assertEqual(response.status,200);self.assertNotIn('_execution',payload)
            self.assertEqual(payload['status'],'running')
        finally:
            server.shutdown();server.server_close();thread.join(timeout=2)

    def test_environment_probes_are_also_in_check_archive(self):
        job={'status':'running','action':'preflight','log':str(self.root/'probe.log')}
        app.worker(job,'preflight',self.config|{'output':str(self.root/'checks')})
        self.assertEqual(job['status'],'done')
        run=Path(job['run'])
        entries=[json.loads(line) for line in (run/'commands.jsonl').read_text().splitlines()]
        self.assertEqual(len(entries),10)
        report=self.report(job).read_text(encoding='utf8')
        self.assertIn('qiime info',report);self.assertIn('qiime dada2 --help',report)
        self.assertIn('check_r.R', (self.report(job).parent/'R完整代码.txt').read_text(encoding='utf8'))


if __name__=='__main__':unittest.main()
