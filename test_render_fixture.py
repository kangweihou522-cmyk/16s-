#!/usr/bin/env python3
"""Reproducible real-R render validation. No QIIME calculations or app edits.
Artifacts: ./r-render-fixtures (synthetic exports).
Usage: python3 test_render_fixture.py; Rscript and plotting packages must be installed.
"""
from pathlib import Path
import csv, itertools, json, math, os, shutil, subprocess, sys
BASE=Path(__file__).resolve().parent
ROOT=BASE/'r-render-fixtures'
SCRIPT=BASE/'render.R'
RSCRIPT=Path(os.environ.get('WORKSTATION_TEST_RSCRIPT',shutil.which('Rscript') or 'Rscript'))
ASSETS=BASE/'plots'
IDS=[f'{i:03}' for i in range(1,9)]
GROUP={sid:'CON' if int(sid)<=4 else 'TRT' for sid in IDS}
ORDER=['008','003','001','006','002','007','004','005']
COORDS={'001':(-2.0,-1.0),'002':(-1.0,0.0),'003':(-1.5,1.0),'004':(-2.5,.5),'005':(1.0,-1.0),'006':(2.0,.2),'007':(1.8,1.2),'008':(2.7,.7)}
METRICS=['weighted_unifrac','unweighted_unifrac','bray_curtis','jaccard']
ALPHA=['observed_features','shannon','evenness','faith_pd','chao1','simpson','goods_coverage']
ALPHA_LABELS=['Observed features','Shannon diversity','Pielou evenness','Faith phylogenetic diversity','Chao1 richness','Simpson diversity',"Good's coverage"]

def write(path,text):
 path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')

def fixture(name,mutate=None):
 case=ROOT/name;run=case/'analysis';pre=case/'preprocess'
 if case.exists():shutil.rmtree(case)
 md='sample-id\tGroup\n#q2:types\tcategorical\n'+''.join(f'{sid}\t{GROUP[sid]}\n' for sid in IDS)
 write(run/'metadata.tsv',md);write(pre/'metadata.tsv',md)
 for j,metric in enumerate(ALPHA):
  write(run/'exports/alpha'/metric/'alpha-diversity.tsv','sample-id\t'+metric+'\n'+''.join(f'{sid}\t{(int(sid)+j)/10}\n' for sid in reversed(IDS)))
 for metric in METRICS:
  text='Eigvals\t2\n0.6\t0.4\n\nProportion explained\t2\n0.6\t0.4\n\nSpecies\t0\t0\n\nSite\t8\t2\n'
  text+=''.join(f'{sid}\t{COORDS[sid][0]}\t{COORDS[sid][1]}\n' for sid in ORDER)
  text+='\nBiplot\t0\t0\n\nSite constraints\t0\t0\n'
  write(run/'exports/pcoa'/metric/'ordination.txt',text)
  # Matrix intentionally ordered differently from saved ordination.
  text='\t'+'\t'.join(IDS)+'\n'
  for sid in IDS:
   text+=sid+'\t'+'\t'.join(f'{math.dist(COORDS[sid],COORDS[other]):.15g}' for other in IDS)+'\n'
  write(run/'exports/distance'/metric/'distance-matrix.tsv',text)
 a=[90,100,80,20,800,20,60,40];b=[10,900,20,80,200,80,40,60]
 write(pre/'exported-table/feature-table.tsv','# Constructed from biom file\n#OTU ID\t'+'\t'.join(IDS)+'\nASV1\t'+'\t'.join(map(str,a))+'\nASV2\t'+'\t'.join(map(str,b))+'\n')
 write(pre/'exported-taxonomy/taxonomy.tsv','Feature ID\tTaxon\tConfidence\nASV1\td__Bacteria; p__Firmicutes; c__Bacilli; o__OrdA; f__FamA; g__Alpha; s__A\t0.99\nASV2\td__Bacteria; p__Bacteroidota; c__Bacteroidia; o__OrdB; f__FamB; g__Beta; s__B\t0.98\n')
 if mutate:mutate(run,pre)
 return run,pre

def execute(name,assets,mutate=None,success=True):
 run,pre=fixture(name,mutate);env=os.environ.copy();env['XDG_CACHE_HOME']=str(ROOT/'cache')
 result=subprocess.run([str(RSCRIPT),str(SCRIPT),str(run),'DejaVu Serif',str(pre),str(assets)],capture_output=True,text=True,env=env)
 write(run.parent/'R_stdout.txt',result.stdout);write(run.parent/'R_stderr.txt',result.stderr)
 print(f'{name}: exit={result.returncode}',flush=True)
 if success and result.returncode:raise AssertionError(result.stderr)
 if not success and result.returncode==0:raise AssertionError(f'{name}: invalid fixture unexpectedly succeeded')
 return run,result

def read_csv(path):
 with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def validate_positive(run,original,out=None):
 out=out or run/'R_original_style';summary=json.loads((out/'render_summary.json').read_text())
 assert summary['original_assets']==original
 assert summary['source']=='saved_QIIME2_outputs'
 assert not (out/'alpha/alpha_KW_saved.csv').exists()
 assert not (out/'alpha/alpha_pairwise_KW_BH_saved.csv').exists()
 assert not (out/'beta/beta_QIIME_overall_saved.csv').exists()
 result_lines=(out/'results.tsv').read_text().splitlines()
 assert len(result_lines)-1>35, 'Artifact index must list generated results'
 for metric in METRICS:
  points=read_csv(out/'beta'/f'{metric}_PCoA_coordinates.csv')
  assert [p['SampleID'] for p in points]==ORDER
  assert all(float(p['PC1'])==COORDS[p['SampleID']][0] and float(p['PC2'])==COORDS[p['SampleID']][1] for p in points)
  axes=read_csv(out/'beta'/f'{metric}_PCoA_axes.csv')
  assert [float(p['ProportionExplained']) for p in axes]==[.6,.4]
  svg=(out/'beta'/f'PCoA_{metric}.svg').read_text()
  assert '60%' in svg or '60.0%' in svg
  assert '40%' in svg or '40.0%' in svg
  assert 'PERMANOVA:' not in svg,'No fabricated significance without saved tests'
 for rank,taxon in [('phylum','Firmicutes'),('genus','Alpha')]:
  rel=read_csv(out/'taxonomy'/f'{rank}_sample_relative_abundance_unrarefied.csv')
  assert set(rel[0])=={'Taxon',*IDS},'Leading-zero IDs must stay strings in header'
  for sid in IDS:assert abs(sum(float(r[sid]) for r in rel)-1)<1e-12
  means=read_csv(out/'taxonomy'/f'{rank}_group_mean_unrarefied.csv')
  selected=[m for m in means if m['Taxon']==taxon]
  assert len(selected)==2
  assert all(abs(float(m['Abundance'])-.5)<1e-12 for m in selected),'Means must not pool raw reads'
 annotations=read_csv(out/'alpha/alpha_letters_positions.csv')
 assert all(a['Letter']=='' for a in annotations),'No fabricated alpha annotations'
 alpha=read_csv(out/'alpha/alpha_plot_input.csv');assert set(a['Sample'] for a in alpha)==set(IDS);assert set(a['Metric'] for a in alpha)==set(ALPHA), 'All seven alpha metrics must render'
 for j,metric in enumerate(ALPHA):
  points=read_csv(out/'alpha'/f'{metric}_plot_input.csv')
  assert [p['Sample'] for p in points]==list(reversed(IDS)),f'{metric}: saved sample IDs/order changed'
  assert all(p['Group']==GROUP[p['Sample']] and p['Metric']==metric and abs(float(p['Value'])-(int(p['Sample'])+j)/10)<1e-12 for p in points),f'{metric}: saved values/group mapping changed'
  assert points==[p for p in alpha if p['Metric']==metric]
  letters=read_csv(out/'alpha'/f'{metric}_letters_positions.csv')
  assert {a['Group'] for a in letters}==set(GROUP.values()) and len(letters)==2
  assert all(a['Metric']==metric and a['Letter']=='' for a in letters),'No fabricated per-metric alpha annotations'
  for annotation in letters:
   expected=max(float(p['Value']) for p in points if p['Group']==annotation['Group'])*1.02
   assert abs(float(annotation['y'])-expected)<1e-12
  assert letters==[a for a in annotations if a['Metric']==metric]
  svg=(out/'alpha'/f'{metric}.svg').read_text()
  assert ALPHA_LABELS[j] in svg or ALPHA_LABELS[j].replace("'",'&apos;') in svg,f'{metric}: figure label missing'
 assert not any((out/'alpha'/f'alpha_diversity.{suffix}').exists() for suffix in ['svg','pdf','png']),'Combined alpha panel must not be generated'
 stems=[*[f'alpha/{metric}' for metric in ALPHA],'beta/PCoA_four_metrics',*[f'beta/PCoA_{metric}' for metric in METRICS],*[f'taxonomy/{rank}_{kind}_unrarefied' for rank in ['phylum','genus'] for kind in ['sample','group_mean']]]
 for stem in stems:
  for suffix in ['svg','pdf','png']:
   path=out/(stem+'.'+suffix);assert path.stat().st_size>500,str(path)
   prefix=path.read_bytes()[:8]
   if suffix=='pdf':assert prefix.startswith(b'%PDF-')
   if suffix=='png':assert prefix==b'\x89PNG\r\n\x1a\n'
 figures={p.relative_to(out).as_posix() for p in out.rglob('*') if p.suffix in {'.svg','.pdf','.png'}}
 assert figures=={f'{stem}.{suffix}' for stem in stems for suffix in ['svg','pdf','png']},'Expected exactly 48 named figure files'
 return len(stems)*3

def ord_replace(old,new):
 def mutate(run,pre):
  path=run/'exports/pcoa/weighted_unifrac/ordination.txt'
  write(path,path.read_text().replace(old,new,1))
 return mutate

def invalid_alpha(run,pre):
 path=run/'exports/alpha/observed_features/alpha-diversity.tsv';write(path,path.read_text().replace('008\t0.8','MISSING\t0.8'))

def invalid_distance(run,pre):
 path=run/'exports/distance/weighted_unifrac/distance-matrix.tsv';write(path,path.read_text().replace('001','999'))

def main():
 ROOT.mkdir(exist_ok=True)
 counts={}
 for original in [False,True]:
  name='original' if original else 'fallback';run,_=execute(name,ASSETS if original else ROOT/'missing-assets')
  counts[name]=validate_positive(run,original)
 for name,mutation in [('unknown_alpha_id',invalid_alpha),('unknown_pcoa_id',ord_replace('008\t2.7','999\t2.7')),('distance_id_mismatch',invalid_distance),('nonfinite_pcoa',ord_replace('008\t2.7','008\tNaN')),('negative_proportion',ord_replace('Proportion explained\t2\n0.6\t0.4','Proportion explained\t2\n-0.6\t0.4')),('nonfinite_proportion',ord_replace('Proportion explained\t2\n0.6\t0.4','Proportion explained\t2\nNaN\t0.4'))]:
  _,result=execute(name,ROOT/'missing-assets',mutation,success=False)
  assert 'Error' in result.stderr
 print(json.dumps({'positive_modes':counts,'negative_cases':6,'status':'PASS'},indent=2),flush=True)
if __name__=='__main__':main()
