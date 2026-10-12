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

    def test_supported_plateau_is_separate_from_minimum_read_recommendation(self):
        advice=recommend_depth(self.samples,self.counts,self.plateau())
        self.assertEqual(advice['recommended_depth'],4900)
        self.assertEqual(advice['status'],'minimum_read_rule')
        self.assertEqual(advice['curve_plateau_candidate_depth'],1000)
        self.assertEqual(advice['curve_evidence']['depths'],[1000,2000,4000])
        self.assertTrue(all(row['retained'] for row in advice['retention']))

    def test_one_unstable_group_prevents_platform_claim_but_not_numeric_recommendation(self):
        curves=self.plateau();curves['observed_features']['004'][4000]=200
        advice=recommend_depth(self.samples,self.counts,curves)
        self.assertEqual(advice['status'],'minimum_read_rule')
        self.assertIsNone(advice['curve_evidence'])
        self.assertIsNone(advice['curve_plateau_candidate_depth'])
        self.assertEqual(advice['recommended_depth'],4900)

    def test_recommendation_is_largest_hundred_strictly_below_smallest_original_sample(self):
        for minimum,expected in [(17105,17100),(17100,17000),(101,100)]:
            with self.subTest(minimum=minimum):
                counts={sid:30000 for sid in self.counts}|{'001':minimum}
                advice=recommend_depth(self.samples,counts)
                self.assertEqual(advice['recommended_depth'],expected)
                self.assertLess(expected,minimum)
                self.assertEqual(expected%100,0)
                self.assertTrue(all(row['retained'] for row in advice['retention']))

    def test_all_samples_in_small_and_large_groups_are_retained(self):
        samples=[{'SampleID':f'A{i}','Group':'A'} for i in range(20)]+[{'SampleID':f'B{i}','Group':'B'} for i in range(3)]
        counts={sample['SampleID']:10000 for sample in samples};counts.update(B0=700,B1=2000,B2=3000)
        advice=recommend_depth(samples,counts)
        self.assertEqual(advice['recommended_depth'],600)
        self.assertEqual(len(advice['retention']),len(samples))
        self.assertTrue(all(row['retained'] for row in advice['retention']))
        self.assertTrue(all(group['retained']==group['required']==group['original'] for group in advice['groups']))

    def test_zero_sample_of_ten_is_not_silently_dropped(self):
        samples=[{'SampleID':f'{group}{i}','Group':group} for group in ['A','B'] for i in range(10)]
        counts={sample['SampleID']:6000 for sample in samples};counts['A0']=0
        curves={metric:{sample['SampleID']:{1000:1,2000:1,4000:1} for sample in samples} for metric in THRESHOLDS}
        advice=recommend_depth(samples,counts,curves)
        self.assertIsNone(advice['recommended_depth'])
        self.assertEqual(advice['status'],'insufficient_depth')
        self.assertEqual(advice['curve_plateau_candidate_depth'],1000)
        self.assertEqual(len(advice['retention']),20)
        self.assertEqual(next(row['count'] for row in advice['retention'] if row['SampleID']=='A0'),0)
        self.assertTrue(any('不能自动丢弃' in note and '人工' in note for note in advice['warnings']))

    def test_missing_curve_or_metric_does_not_change_minimum_read_recommendation(self):
        curves=self.plateau();del curves['faith_pd']
        advice=recommend_depth(self.samples,self.counts,curves)
        self.assertEqual(advice['status'],'minimum_read_rule')
        self.assertEqual(advice['recommended_depth'],4900)
        self.assertIsNone(advice['curve_plateau_candidate_depth'])
        self.assertTrue(any('缺少' in note for note in advice['warnings']))

    def test_minimum_at_or_below_one_hundred_returns_none_and_manual_review_warning(self):
        for minimum in [0,1,99,100]:
            with self.subTest(minimum=minimum):
                advice=recommend_depth(self.samples,self.counts|{'001':minimum})
                self.assertIsNone(advice['recommended_depth'])
                self.assertEqual(advice['status'],'insufficient_depth')
                self.assertTrue(any('人工' in note for note in advice['warnings']))
        counts=self.counts.copy();del counts['001']
        self.assertIsNone(recommend_depth(self.samples,counts)['recommended_depth'])

    def test_unusable_groups_return_none_and_manual_review_warning(self):
        for samples in [[],self.samples[:2],self.samples[:3]]:
            with self.subTest(samples=samples):
                advice=recommend_depth(samples,self.counts)
                self.assertIsNone(advice['recommended_depth'])
                self.assertEqual(advice['status'],'insufficient_samples')
                self.assertTrue(any('人工' in note for note in advice['warnings']))

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
            (table/'feature-table.tsv').write_text('#OTU ID\t001\t002\t003\t004\nASV1\t17105\t20000\t20000\t20000\n')
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
                self.assertEqual(status,200);self.assertEqual(advice['recommended_depth'],17100)
                self.assertEqual(advice['curve_plateau_candidate_depth'],1000)
                self.assertTrue(all(row['retained'] for row in advice['retention']))
                (table/'feature-table.tsv').write_text('#OTU ID\t001\t002\t003\t004\nASV1\t17100\t20000\t20000\t20000\n')
                self.assertEqual(get('/api/depth-advice?run='+str(run))[1]['recommended_depth'],17000)
                (table/'feature-table.tsv').write_text('#OTU ID\t001\t002\t003\t004\nASV1\t0\t20000\t20000\t20000\n')
                status,advice=get('/api/depth-advice?run='+str(run))
                self.assertEqual(status,200);self.assertIsNone(advice['recommended_depth'])
                self.assertEqual(len(advice['retention']),4)
                status,error=get('/api/retention?run='+str(run)+'&depth=')
                self.assertEqual(status,400);self.assertIn('本次抽平深度',error['error'])
                self.assertNotIn('invalid literal',error['error'])
                self.assertEqual(get('/api/depth-advice?run='+str(run),False)[0],403)
            finally:
                server.shutdown();server.server_close();thread.join(timeout=2)


if __name__=='__main__':unittest.main()
