"""History restores complete outputs from the same saved source, without mutation."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import app
from batch_context import list_batches, restore_batch


class BatchTests(unittest.TestCase):
    def setUp(self):
        t=tempfile.TemporaryDirectory(); self.addCleanup(t.cleanup); self.root=Path(t.name)
        self.run=self.root/'run_20260101_123456_id';self.run.mkdir();app.marker(self.run)
        self.config={'samples':[{'SampleID':'001','Group':'CON'}],'input':'/mnt/d/raw','metadata':'/mnt/d/group.csv',
                     'output':str(self.root),'primer_mode':'removed','batch_id':'batch-a','batch_name':'猪试验 A'}
        (self.run/'parameters.json').write_text(json.dumps(self.config))
        (self.run/'quality_complete.json').write_text('{}')

    def child(self,folder,marker,parent=None,batch_id='batch-a'):
        run=self.run/folder;run.mkdir(parents=True);app.marker(run)
        (run/'parameters.json').write_text(json.dumps({'preprocess':str(parent or self.run),'config':{'batch_id':batch_id}}))
        (run/marker).write_text('{}');return run

    def test_quality_resume_and_consumed_cache_detachment(self):
        fields=restore_batch(self.run,app.module_api())['fields']
        self.assertEqual(fields['quality_run'],str(self.run));self.assertEqual(fields['run'],'')
        (self.run/'preprocess_started.json').write_text('{}')
        (self.run/'preprocess_complete.json').write_text('{}')
        result=restore_batch(self.run,app.module_api())
        self.assertEqual(result['fields']['quality_run'],'')
        self.assertEqual(result['fields']['run'],str(self.run))
        self.assertEqual(result['config']['data_state'],'clean')
        self.assertEqual(result['config']['samples'][0]['SampleID'],'001')

    def test_child_outputs_parent_and_batch_identity_guard_no_mutation(self):
        good=self.child('07_diversity/diversity_100_good','qiime_statistics_complete.json')
        self.child('07_diversity/diversity_900_wrong','qiime_statistics_complete.json',batch_id='other')
        self.child('07_diversity/diversity_999_wrong','qiime_statistics_complete.json',parent=self.root/'other')
        self.child('07_diversity/diversity_9999_incomplete','metadata.tsv')
        lefse=self.child('11_extensions/lefse_01','lefse_complete.json')
        diff=self.child('11_extensions/differential_01','differential_complete.json')
        (diff/'differential_options.json').write_text('{}')
        before={p:hashlib.sha256(p.read_bytes()).hexdigest() for p in self.run.rglob('*') if p.is_file()}
        result=restore_batch(good,app.module_api())
        self.assertEqual(result['fields']['diversity_run'],str(good))
        self.assertEqual(result['fields']['lefse_run'],str(lefse))
        self.assertEqual(result['fields']['differential_run'],str(diff))
        self.assertEqual(result['batch_name'],'猪试验 A')
        self.assertEqual(before,{p:hashlib.sha256(p.read_bytes()).hexdigest() for p in self.run.rglob('*') if p.is_file()})

    def test_list_ignores_arbitrary_folders_and_existing_settings_retained(self):
        unrelated=self.root/'run_wrong';unrelated.mkdir();(unrelated/'parameters.json').write_text('{}')
        result=list_batches(self.root,app.module_api())
        self.assertEqual(len(result['runs']),1)
        self.assertEqual(result['runs'][0]['batch_name'],'猪试验 A')
        (self.run/'group_style.json').write_text(json.dumps({'group_order':['CON'],'group_colors':{'CON':'#123456'}}))
        (self.run/'plot_settings.json').write_text(json.dumps({'taxonomy':{'palette':'user_R'}}))
        restored=restore_batch(self.run,app.module_api())
        self.assertEqual(restored['config']['group_colors']['CON'],'#123456')
        self.assertEqual(restored['config']['plot_settings']['taxonomy']['palette'],'user_R')
        with self.assertRaises(ValueError):restore_batch(unrelated,app.module_api())
