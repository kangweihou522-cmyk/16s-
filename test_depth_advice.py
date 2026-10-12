import csv
import tempfile
import unittest
import zipfile
import json
import http.client
import threading
from pathlib import Path

from depth_advice import read_curves, recommend_depth, THRESHOLDS
import app


class DepthAdviceTests(unittest.TestCase):
    def setUp(self):
        self.samples=[{'SampleID':sid,'Group':group} for sid,group in [('001','A'),('002','A'),('003','B'),('004','B')]]
        self.counts={sample['SampleID']:5000 for sample in self.samples}

    def plateau(self):
        values={'observed_features':[50,100,103,104],'shannon':[1,2,2.01,2.02],'faith_pd':[4,10,10.2,10.3]}
        return {metric:{sid:dict(zip([500,1000,2000,4000],series)) for sid in self.counts} for metric,series in values.items()}

    def test_earliest_supported_plateau_uses_all_groups_and_metrics(self):
        advice=recommend_depth(self.samples,self.counts,self.plateau())
        self.assertEqual(advice['recommended_depth'],1000)
        self.assertEqual(advice['status'],'plateau_candidate')
        self.assertEqual(advice['curve_evidence']['depths'],[1000,2000,4000])
        self.assertTrue(all(row['retained'] for row in advice['retention']))

    def test_one_unstable_group_prevents_claiming_platform(self):
        curves=self.plateau();curves['observed_features']['004'][4000]=200
        advice=recommend_depth(self.samples,self.counts,curves)
        self.assertEqual(advice['status'],'retention_only')
        self.assertIsNone(advice['curve_evidence'])
        self.assertEqual(advice['recommended_depth'],5000)

    def test_each_group_retention_is_not_hidden_by_large_other_group(self):
        samples=[{'SampleID':f'A{i}','Group':'A'} for i in range(20)]+[{'SampleID':f'B{i}','Group':'B'} for i in range(3)]
        counts={sample['SampleID']:10000 for sample in samples};counts.update(B0=700,B1=2000,B2=3000)
        advice=recommend_depth(samples,counts)
        self.assertEqual(advice['recommended_depth'],700)
        self.assertTrue(all(group['retained']>=group['required'] for group in advice['groups']))

    def test_one_low_sample_of_ten_can_be_dropped_at_ninety_percent(self):
        samples=[{'SampleID':f'{group}{i}','Group':group} for group in ['A','B'] for i in range(10)]
        counts={sample['SampleID']:6000 for sample in samples};counts['A0']=0
        advice=recommend_depth(samples,counts)
        self.assertEqual(advice['recommended_depth'],6000)
        self.assertEqual({g['Group']:g['retained'] for g in advice['groups']},{'A':9,'B':10})

    def test_missing_curve_or_metric_is_explicitly_preliminary(self):
        curves=self.plateau();del curves['faith_pd']
        advice=recommend_depth(self.samples,self.counts,curves)
        self.assertEqual(advice['status'],'retention_only')
        self.assertTrue(any('缺少' in note for note in advice['warnings']))

    def test_no_positive_feasible_depth_or_too_few_groups_returns_none(self):
        counts=self.counts|{'001':0}
        self.assertIsNone(recommend_depth(self.samples,counts)['recommended_depth'])
        self.assertIsNone(recommend_depth(self.samples[:2],self.counts)['recommended_depth'])

    def test_unknown_samples_and_invalid_counts_rejected(self):
        curves=self.plateau();curves['shannon']['other']={500:1}
        with self.assertRaisesRegex(ValueError,'分组不一致'): recommend_depth(self.samples,self.counts,curves)
        with self.assertRaisesRegex(ValueError,'非负整数'): recommend_depth(self.samples,self.counts|{'001':float('nan')})

    def test_actual_qiime_csv_shape_sparse_values_and_legacy_qzv(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder=Path(temporary);archive=folder/'curve.qzv'
            text='sample-id,depth-100_iter-1,depth-100_iter-2,depth-200_iter-1,Group\n001,100,120,,A\n'
            with zipfile.ZipFile(archive,'w') as zipped:
                for metric in THRESHOLDS:zipped.writestr('uuid/data/'+metric+'.csv',text)
            curves=read_curves(visualization=archive)
            self.assertEqual(curves['shannon']['001'],{100:110})
            (folder/'shannon.csv').write_text('sample-id,depth-100_iter-1\n001,inf\n')
            with self.assertRaisesRegex(ValueError,'无效数值'):read_curves(folder=folder)

    def test_blank_invalid_and_boundary_integer_errors_name_the_field(self):
        for raw in ['',None,' ','1.5',True,'1e3']:
            with self.subTest(raw=raw),self.assertRaisesRegex(ValueError,'本次抽平深度'):
                app.integer({'depth':raw},'depth',1,10000)
        self.assertEqual(app.integer({'depth':' 5000 '},'depth',1,10000),5000)
        with self.assertRaisesRegex(ValueError,'R1 截断位置'):app.integer({},'trunc_f',1,1000)

    def test_real_http_historical_recommendation_and_blank_depth_diagnostic(self):
        with tempfile.TemporaryDirectory() as temporary:
            run=Path(temporary);app.marker(run)
            (run/'parameters.json').write_text(json.dumps({'samples':self.samples}))
            table=run/'exported-table';table.mkdir()
            (table/'feature-table.tsv').write_text('#OTU ID\t001\t002\t003\t004\nASV1\t5000\t5000\t5000\t5000\n')
            curves=run/'exported-rarefaction';curves.mkdir()
            for metric,values in self.plateau().items():
                with (curves/(metric+'.csv')).open('w',newline='') as file:
                    writer=csv.writer(file);writer.writerow(['sample-id',*[f'depth-{depth}_iter-1' for depth in [500,1000,2000,4000]]])
                    for sid,series in values.items():writer.writerow([sid,*series.values()])
            server=app.ThreadingHTTPServer(('127.0.0.1',0),app.Handler)
            thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
            def get(path,authorized=True):
                connection=http.client.HTTPConnection('127.0.0.1',server.server_port,timeout=5)
                connection.request('GET',path,headers={'X-Workstation-Token':app.TOKEN} if authorized else {})
                response=connection.getresponse();status=response.status;payload=json.loads(response.read());connection.close();return status,payload
            try:
                status,advice=get('/api/depth-advice?run='+str(run))
                self.assertEqual(status,200);self.assertEqual(advice['recommended_depth'],1000)
                status,error=get('/api/retention?run='+str(run)+'&depth=')
                self.assertEqual(status,400);self.assertIn('本次抽平深度',error['error'])
                self.assertNotIn('invalid literal',error['error'])
                self.assertEqual(get('/api/depth-advice?run='+str(run),False)[0],403)
            finally:
                server.shutdown();server.server_close();thread.join(timeout=2)


if __name__=='__main__':unittest.main()
