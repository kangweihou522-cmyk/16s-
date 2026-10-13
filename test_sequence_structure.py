"""Real FASTQ inspection/export and opt-in official Cutadapt trimming checks."""
import csv
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import app
from sequence_structure import PRESETS, IUPAC, inspect_samples, primer_list, export_clean_fastq
from test_quality_pipeline import materialize_qiime_outputs


def fastq(filename, sequence, mate, n=40):
    filename.parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(filename,'wt',encoding='ascii') as stream:
        for i in range(n): stream.write(f'@read{i}/{mate}\n{sequence}\n+\n'+('I'*len(sequence))+'\n')


class StructureTests(unittest.TestCase):
    def setUp(self):
        t=tempfile.TemporaryDirectory(); self.addCleanup(t.cleanup); self.root=Path(t.name)
        self.raw=self.root/'raw'; self.raw.mkdir()
        self.f=PRESETS['341F'][0]; self.r=PRESETS['806R'][0]
        self.body='TTGAGCTAGCTACGAT'*8
        concrete=lambda p:''.join(IUPAC[b][0] for b in p)
        for mate,primer in [(1,self.f),(2,self.r)]: fastq(self.raw/f'001_R{mate}.fastq.gz','ACGTACGTACGT'+concrete(primer)+self.body,mate)
        self.metadata=self.root/'groups.csv'; self.metadata.write_text('SampleID,Group\n001,CON\n')
        self.config=dict(input=str(self.raw),metadata=str(self.metadata),output=str(self.root/'output'),mode='paired',data_state='unknown',primer_mode='341F',threads='1',env='mock',r_env='mock-r',conda='/fake/conda',font='DejaVu Sans',parameter_mode='auto',batch_id='batch1')
        self.samples=app.discover(self.config)
        self.before={p:hashlib.sha256(p.read_bytes()).hexdigest() for p in self.raw.iterdir()}

    def test_full_degenerate_primer_prefix_detected_without_barcode_identity_claim(self):
        report=inspect_samples(self.samples,self.config,app.reads)
        self.assertTrue(report['reliable'])
        for direction in ('R1','R2'):
            evidence=report['samples'][0]['structure'][direction]
            best=next(r for r in evidence['candidates'] if r['orientation']=='expected' and r['fraction']==1)
            self.assertEqual(best['modal_start'],12)
            self.assertEqual(best['prefix_reads'],40)
        self.assertFalse(report['barcode_identity_confirmed'])
        self.assertEqual(report['samples'][0]['SampleID'],'001')

    def test_no_match_and_reverse_orientation_do_not_become_clean(self):
        for mate in (1,2): fastq(self.raw/f'001_R{mate}.fastq.gz',self.body,mate)
        report=inspect_samples(app.discover(self.config),self.config,app.reads)
        self.assertFalse(report['reliable'])
        from sequence_structure import COMPLEMENT
        fastq(self.raw/'001_R1.fastq.gz',self.f.translate(COMPLEMENT)[::-1].replace('N','A').replace('W','A')+self.body,1)
        report=inspect_samples(app.discover(self.config),self.config,app.reads)
        self.assertFalse(report['reliable'])
        self.assertTrue(any(r['orientation']=='reverse_complement' and r['fraction']==1 for r in report['samples'][0]['structure']['R1']['candidates']))

    def test_unknown_failed_report_preserved_before_any_qiime(self):
        for mate in (1,2): fastq(self.raw/f'001_R{mate}.fastq.gz',self.body,mate)
        job={}
        with patch.object(app,'q') as qiime,self.assertRaisesRegex(ValueError,'未确认'):
            app.preprocess(job,self.config,quality_only=True)
        qiime.assert_not_called()
        self.assertTrue((Path(job['run'])/'00_checks/sequence_structure.json').is_file())

    def test_ambiguous_boundary_and_too_few_reads_require_review(self):
        config=dict(self.config,primer_mode='custom',forward=self.f+';'+self.f+'T',reverse=self.r)
        report=inspect_samples(self.samples,config,app.reads)
        self.assertFalse(report['reliable'])
        self.assertEqual(report['samples'][0]['structure']['R1']['status'],'ambiguous')
        fastq(self.raw/'001_R1.fastq.gz','CCTACGGGAGGCAGCAG'+self.body,1,n=2)
        fastq(self.raw/'001_R2.fastq.gz','GGACTACTAGGGTATCTAAT'+self.body,2,n=2)
        self.assertFalse(inspect_samples(app.discover(self.config),self.config,app.reads)['reliable'])
        for bad in ('NNNNNNNNNNNNNNNN','ACGT;touch /tmp/x','ACGT'):
            with self.assertRaises(ValueError):primer_list(bad)

    def test_clean_state_skips_cutadapt_and_merged_never_runs_dada2(self):
        config=dict(self.config,data_state='clean',primer_mode='auto')
        calls=[]
        def qiime(job,c,*args): calls.append(args); materialize_qiime_outputs(job,c,*args)
        with patch.object(app,'q',side_effect=qiime),patch.object(app,'command'):
            app.preprocess({},config,quality_only=True)
        self.assertFalse(any(args[0]=='cutadapt' for args in calls))
        with patch.object(app,'q') as q,self.assertRaisesRegex(ValueError,'已提前拼接'):
            app.preprocess({},dict(config,input_layout='merged'))
        q.assert_not_called()

    def test_clean_manifest_mapping_and_traversal_protection(self):
        export=self.root/'export'; export.mkdir()
        for mate in (1,2): fastq(export/f'seq{mate}.fastq.gz',self.body,mate)
        (export/'MANIFEST').write_text('sample-id,filename,direction\n001,seq1.fastq.gz,forward\n001,seq2.fastq.gz,reverse\n')
        result=export_clean_fastq(export,self.root/'clean',self.samples,self.config,app.reads)
        with Path(result['sample_sheet']).open(newline='') as stream: rows=list(csv.DictReader(stream))
        self.assertEqual(rows[0]['SampleID'],'001')
        self.assertEqual(next(app.reads(Path(rows[0]['forward'])))[1],self.body)
        self.assertFalse(result['residual_warning'])
        (export/'MANIFEST').write_text('sample-id,filename,direction\n001,../raw/001_R1.fastq.gz,forward\n')
        with self.assertRaises(ValueError): export_clean_fastq(export,self.root/'unsafe',self.samples,self.config,app.reads)
        self.assertEqual({p:hashlib.sha256(p.read_bytes()).hexdigest() for p in self.raw.iterdir()},self.before)

    @unittest.skipUnless(importlib.util.find_spec('cutadapt'),'Enable official Cutadapt via PYTHONPATH for real trimming')
    def test_real_cutadapt_export_and_quality_pipeline(self):
        calls=[]
        def qiime(job,c,*args):
            calls.append(args); materialize_qiime_outputs(job,c,*args)
            if args[:2]==('cutadapt','trim-paired'):
                clean=Path(job['run'])/'actual_trimmed';clean.mkdir()
                argv=[sys.executable,'-m','cutadapt','-e','.1','--discard-untrimmed','--pair-filter=any',
                      '-o',str(clean/'R1.fastq.gz'),'-p',str(clean/'R2.fastq.gz'),'--json',str(clean/'cutadapt_counts.json')]
                for flag,target in [('--p-front-f','-g'),('--p-front-r','-G')]:
                    i=args.index(flag)+1
                    while i<len(args) and not str(args[i]).startswith('--'):
                        argv.extend([target,str(args[i])]);i+=1
                argv.extend([str(self.raw/'001_R1.fastq.gz'),str(self.raw/'001_R2.fastq.gz')])
                result=subprocess.run(argv,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                report=json.loads((clean/'cutadapt_counts.json').read_text());self.assertEqual(report['read_counts']['output'],40)
            if args[:2]==('tools','export') and args[args.index('--input-path')+1]=='trimmed.qza':
                destination=Path(job['run'])/args[args.index('--output-path')+1];destination.mkdir(parents=True,exist_ok=True)
                import shutil
                for mate in (1,2):shutil.copyfile(Path(job['run'])/f'actual_trimmed/R{mate}.fastq.gz',destination/f'R{mate}.fastq.gz')
                (destination/'MANIFEST').write_text('sample-id,filename,direction\n001,R1.fastq.gz,forward\n001,R2.fastq.gz,reverse\n')
        job={}
        with patch.object(app,'q',side_effect=qiime),patch.object(app,'command'):
            app.preprocess(job,self.config,quality_only=True)
        run=Path(job['run']);saved=json.loads((run/'parameters.json').read_text())
        self.assertIn('1',saved['cleanup_plan'])
        self.assertTrue((run/'quality_complete.json').is_file())
        self.assertEqual(next(app.reads(run/'02_primers/clean_fastq/sample_0001_R1.fastq.gz'))[1],self.body)
        self.assertEqual(next(app.reads(run/'02_primers/clean_fastq/sample_0001_R2.fastq.gz'))[1],self.body)
        self.assertFalse(any(args[0]=='dada2' for args in calls))
        self.assertEqual({p:hashlib.sha256(p.read_bytes()).hexdigest() for p in self.raw.iterdir()},self.before)
