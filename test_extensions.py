"""Verify integration preserves historical inputs and cached fit provenance."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import app
from test_quality_pipeline import materialize_qiime_outputs


class ExtensionTests(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup)
        self.root=Path(temporary.name);self.parent=self.root/'parent';self.parent.mkdir();app.marker(self.parent)
        self.samples=[{'SampleID':f'00{i+1}','Group':'A' if i<2 else 'B'} for i in range(4)]
        self.saved={'samples':self.samples,'env':'fake','conda':'fake-conda','font':'DejaVu Serif'}
        (self.parent/'parameters.json').write_text(json.dumps(self.saved))
        (self.parent/'preprocess_complete.json').write_text('{}')
        for name in ['filtered-table.qza','taxonomy.qza','taxa-barplot.qzv']:(self.parent/name).write_bytes(b'fake artifact')
        table=self.parent/'exported-table';table.mkdir()
        (table/'feature-table.tsv').write_text('#OTU ID\t001\t002\t003\t004\nASV1\t125\t250\t375\t500\n')
        (table/'feature-table.biom').write_bytes(b'fake biom')
        taxonomy=self.parent/'exported-taxonomy';taxonomy.mkdir()
        (taxonomy/'taxonomy.tsv').write_text('Feature ID\tTaxon\tConfidence\nASV1\td__Bacteria; p__P; c__C; o__O; f__F; g__G; s__S\t0.99\n')

    def convert(self,job,args,cwd=None,process_env=None):
        if 'biom' not in args:return
        output=Path(args[args.index('-o')+1]);output.parent.mkdir(parents=True,exist_ok=True)
        values=['1']*4 if 'relative_export' in output.parts else ['125','250','375','500']
        output.write_text('#OTU ID\t001\t002\t003\t004\nASV1\t'+'\t'.join(values)+'\n')

    def test_historical_abundance_is_isolated_and_does_not_run_dada2(self):
        original={file.relative_to(self.parent):file.read_bytes() for file in self.parent.rglob('*') if file.is_file()}
        calls=[];job={'status':'running','action':'abundance'}
        def q(job,saved,*args):calls.append(args);materialize_qiime_outputs(job,saved,*args)
        with patch.object(app,'q',side_effect=q),patch.object(app,'command',side_effect=self.convert):
            app.worker(job,'abundance',{'run':str(self.parent)})
        self.assertEqual(job['status'],'done',job.get('error'))
        run=Path(job['run']);self.assertTrue((run/'abundance_complete.json').is_file())
        self.assertEqual(len(job['result']['abundance']['levels']),6)
        self.assertTrue((run/'10_abundance/asv/relative_with_taxonomy.csv').is_file())
        self.assertFalse(any(args[0]=='dada2' for args in calls))
        for filename,content in original.items():self.assertEqual((self.parent/filename).read_bytes(),content)
        self.assertTrue(Path(job['execution_code_directory'],'完整运行代码.txt').is_file())

    def test_failed_extension_preserves_parent_and_has_no_false_completion(self):
        job={'status':'running','action':'functional'}
        with patch('functional_prediction.run_functional',side_effect=RuntimeError('missing local PICRUSt2')):
            app.worker(job,'functional',{'run':str(self.parent),'ko_mapping':'selected.tsv'})
        self.assertEqual(job['status'],'failed');self.assertIn('PICRUSt2',job['error'])
        self.assertFalse(Path(job['run'],'functional_complete.json').exists())
        self.assertTrue(Path(job['execution_code_directory'],'完整运行代码.txt').is_file())
        self.assertTrue((self.parent/'preprocess_complete.json').is_file())

    def test_lefse_fit_marker_hash_survives_wrapper_and_render_uses_new_directory(self):
        job={'status':'running','action':'lefse'}
        def fit(job,c,saved,parent,api):
            Path(job['run'],'lefse_complete.json').write_text(json.dumps({'result_sha256':'preserve-fit-hash'}))
            return {'lefse_run':job['run']}
        with patch('lefse_analysis.run_lefse',side_effect=fit):app.worker(job,'lefse',{'run':str(self.parent)})
        self.assertEqual(job['status'],'done',job.get('error'));analysis=Path(job['run'])
        self.assertEqual(json.loads((analysis/'lefse_complete.json').read_text())['result_sha256'],'preserve-fit-hash')
        original={file.relative_to(analysis):file.read_bytes() for file in analysis.rglob('*') if file.is_file()}
        render_job={'status':'running','action':'lefse-render'}
        def render(job,c,cached,api):
            self.assertEqual(cached,analysis);self.assertNotEqual(Path(job['run']),analysis)
            return {'lefse_run':str(cached),'statistics_refit':False,'figures':{}}
        with patch('lefse_analysis.render_lefse',side_effect=render),patch('lefse_analysis.run_lefse') as model:
            app.worker(render_job,'lefse-render',{'lefse_run':str(analysis),'lda_threshold':'4'})
        self.assertEqual(render_job['status'],'done',render_job.get('error'));model.assert_not_called()
        for filename,content in original.items():self.assertEqual((analysis/filename).read_bytes(),content)
        record=Path(render_job['execution_code_directory'],'完整运行代码.txt').read_text(encoding='utf8')
        self.assertIn('用于重绘的 LEfSe 原始检验',record)

    def test_whole_group_missing_from_filtered_table_cannot_be_silently_dropped(self):
        (self.parent/'exported-table/feature-table.tsv').write_text('#OTU ID\t001\t002\nASV1\t125\t250\n')
        job={'status':'running','action':'differential'}
        with patch('differential_analysis.run_differential') as model:
            app.worker(job,'differential',{'run':str(self.parent),'diff_method':'ancombc2'})
        self.assertEqual(job['status'],'failed');self.assertIn('某个组',job['error']);model.assert_not_called()


if __name__=='__main__':unittest.main()
