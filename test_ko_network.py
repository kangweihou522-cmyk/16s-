"""Real R module integration, preserved statistics and KGML-only recovery."""
import csv
import json
import os
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch

import app
from ko_network import import_differential
from recovery import recovery_plan

RSCRIPT='/workspace/r-conda-review/Rscript'
KGML='''<?xml version="1.0"?><pathway name="path:ko00010" org="ko" number="00010" title="Synthetic QA">
<entry id="1" name="ko:K00001 ko:K00002" type="ortholog"><graphics name="Enzyme group" type="line" coords="20,30,100,30,100,80" x="60" y="55" width="80" height="50"/></entry>
<entry id="2" name="ko:K00002" type="ortholog"><graphics name="Second" type="rectangle" x="130" y="80" width="40" height="20"/></entry>
<entry id="3" name="cpd:C00001" type="compound"><graphics name="Water" type="circle" x="20" y="30" width="8" height="8"/></entry></pathway>'''


class KOInputTests(unittest.TestCase):
    def test_import_preserves_direction_requires_adjusted_q_and_rejects_duplicates(self):
        with tempfile.TemporaryDirectory() as folder:
            source=Path(folder)/'results.csv';dest=Path(folder)/'canonical.tsv'
            source.write_text('KO,comparison,numerator,denominator,effect,q\nK00001,A–B,A,B,1.5,0.02\n')
            original=source.read_bytes();comparisons=import_differential(source,dest,'log2FC')
            self.assertEqual(comparisons[0]['numerator'],'A');self.assertEqual(source.read_bytes(),original)
            source.write_text(original.decode()+original.decode().splitlines()[1]+'\n')
            with self.assertRaisesRegex(ValueError,'重复'):import_differential(source,dest,'log2FC')
            source.write_text(original.decode().replace(',q\n',',pvalue\n'))
            with self.assertRaisesRegex(ValueError,'校正'):import_differential(source,dest,'log2FC')


@unittest.skipUnless(Path(RSCRIPT).exists(),'QA R runtime is unavailable')
class KONetworkRTests(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory();self.addCleanup(temp.cleanup);self.root=Path(temp.name)
        self.map=self.root/'map.xml';self.map.write_text(KGML)
        self.diff=self.root/'differential.tsv'
        self.diff.write_text('KO\tcomparison\tnumerator\tdenominator\teffect\tq\nK00001\tC001\tB\tA\t2\t0.01\nK00002\tC001\tB\tA\t-1\t0.02\nK99999\tC001\tB\tA\t3\t0.001\n')
        self.config={'ko_source':'differential','ko_differential_file':str(self.diff),'ko_effect_type':'log2FC',
                     'ko_kgml':str(self.map),'ko_map_id':'ko00010','output':str(self.root/'results'),
                     'r_env':'qa-r','font':'DejaVu Serif'}
        self.calls=[]
        def command(config,*args):
            self.calls.append(args)
            return [RSCRIPT,*args]
        self.addCleanup(patch.stopall)
        patch.object(app,'r_command',side_effect=command).start()
        patch.object(app,'r_execution_env',return_value=dict(os.environ)).start()

    def run_module(self,config=None):
        job={'status':'running','action':'ko-network'}
        app.worker(job,'ko-network',config or self.config)
        self.assertEqual(job['status'],'done',job.get('error'))
        return job,Path(job['run'])

    def test_import_map_mixed_unmapped_and_cached_single_figure_without_original_files(self):
        job,run=self.run_module()
        self.assertEqual(len(self.calls),1);self.assertIn('--render',self.calls[0])
        coverage=job['result']['ko_map_summary']['coverage']
        self.assertEqual(coverage['mixed_feature_groups'],1)
        self.assertEqual(coverage['unmapped_significant_kos'],1)
        self.assertIn('K99999',(run/'02_mapping/unmapped_ko_rows.tsv').read_text())
        before={p:p.read_bytes() for p in (run/'02_statistics').iterdir()}
        self.map.unlink();self.diff.unlink()
        preview={'status':'running'}
        app.worker(preview,'figure-render',{'image':str(run/'03_figures/metabolic_map.png'),
                     'plot_settings':{'metabolic_map':{'up_color':'#AA4400','width':9,'q_threshold':.01}}})
        self.assertEqual(preview['status'],'done',preview.get('error'))
        self.assertEqual(len(self.calls),2);self.assertTrue(all('--render' in call for call in self.calls))
        for file,raw in before.items():self.assertEqual(file.read_bytes(),raw)
        self.assertTrue((Path(preview['run'])/'03_figures/metabolic_map.svg').is_file())
        self.assertIn('K00001',(Path(preview['run'])/'02_mapping/ipath_selection.txt').read_text())

    def test_failed_map_download_keeps_completed_stats_and_recovers_only_map(self):
        request=dict(self.config,ko_kgml='')
        job={'status':'running','action':'ko-network'}
        with patch('ko_network.urllib.request.urlopen',side_effect=TimeoutError('synthetic timeout')):
            app.worker(job,'ko-network',request)
        self.assertEqual(job['status'],'failed')
        run=Path(job['run']);before={p:p.read_bytes() for p in (run/'02_statistics').iterdir()}
        plan=recovery_plan(run,app.module_api());self.assertEqual(plan['strategy'],'ko-map')
        retry={'status':'running','action':'recover'}
        app.worker(retry,'recover',dict(plan['config'],recovery_run=str(run),ko_kgml=str(self.map)))
        self.assertEqual(retry['status'],'done',retry.get('error'))
        self.assertNotEqual(Path(retry['run']),run)
        for file,raw in before.items():self.assertEqual(file.read_bytes(),raw)
        self.assertEqual(len(self.calls),1)

    def test_fractional_predicted_table_runs_model_once_then_reuses_stats(self):
        ids=[group+str(i) for group in ('A','B','C') for i in range(4)]
        metadata=self.root/'groups.csv';metadata.write_text('SampleID,Group\n'+''.join(sid+','+sid[0]+'\n' for sid in ids))
        table=self.root/'ko.tsv';rng=random.Random(77)
        with table.open('w',newline='') as stream:
            writer=csv.writer(stream,delimiter='\t');writer.writerow(['KO',*ids])
            for index,ko in enumerate(('K00001','K00002','K00003','K00004')):
                writer.writerow([ko,*[10.5+index*3+i*(index+1)+rng.random()*8 for i in range(len(ids))]])
        job,run=self.run_module(dict(self.config,ko_source='table',ko_table=str(table),ko_metadata=str(metadata)))
        self.assertEqual(len(self.calls),2);self.assertIn('--run',self.calls[0]);self.assertIn('--render',self.calls[1])
        comparisons=json.loads((run/'02_statistics/comparisons.json').read_text());self.assertEqual(len(comparisons),3)
        before={p:p.read_bytes() for p in (run/'02_statistics').iterdir()}
        retry={'status':'running'};app.worker(retry,'ko-network-render',{'ko_network_run':str(run),'ko_comparison':comparisons[1]['id']})
        self.assertEqual(retry['status'],'done',retry.get('error'))
        self.assertEqual(sum('--run' in call for call in self.calls),1)
        for file,raw in before.items():self.assertEqual(file.read_bytes(),raw)


if __name__=='__main__':unittest.main()
