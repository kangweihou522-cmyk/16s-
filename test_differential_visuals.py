"""Focused real-R cached plotting checks; synthetic saved results, no model fit."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent
RSCRIPT = os.environ.get('WORKSTATION_TEST_RSCRIPT') or shutil.which('Rscript')
if not RSCRIPT and Path('/workspace/r-conda-review/Rscript').is_file():
    RSCRIPT = '/workspace/r-conda-review/Rscript'


@unittest.skipUnless(RSCRIPT, 'A real R plotting runtime is not installed')
class DifferentialVisualTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='16s-diff-visual-')
        self.root = Path(self.temporary.name)
        self.source = self.root / 'source'
        self.folder = self.source / '02_differential' / 'L6_genus'
        self.folder.mkdir(parents=True)
        self.counts = self.source / 'counts.tsv'
        self.metadata = self.source / 'metadata.tsv'
        self.counts.write_text('FeatureID\t0001\t0002\t0003\t0004\t0005\t0006\n'
                               'taxon_A\t90\t100\t20\t0\t40\t60\n'
                               'taxon_B\t10\t900\t0\t20\t60\t40\n'
                               'taxon_SZ\t0\t0\t20\t20\t0\t0\n')
        self.metadata.write_text('sample-id\tGroup\n0006\tHDC2\n0001\tCON\n0003\tLDC_1\n'
                                 '0002\tCON\n0005\tHDC2\n0004\tLDC_1\n')
        self.options = {'method': 'ancombc2', 'reference_group': 'CON', 'groups': ['CON', 'LDC_1', 'HDC2'],
                        'metadata_tsv': str(self.metadata), 'alpha': .05, 'pseudo_sens': True,
                        'font': 'sans', 'tables': [{'level': 6, 'name': 'genus', 'counts_tsv': str(self.counts)}]}
        (self.source / 'differential_options.json').write_text(json.dumps(self.options))

    def tearDown(self):
        self.temporary.cleanup()

    def r(self, expression, *args):
        result = subprocess.run([RSCRIPT, '-e', expression, *map(str, args)], capture_output=True,
                                text=True, timeout=60)
        if 'there is no package called' in result.stderr:
            self.skipTest('The existing R runtime lacks plotting dependencies')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def write_fit(self):
        self.r("""a<-commandArgs(TRUE)
res<-data.frame(taxon=c('taxon_A','taxon_B'),lfc_GroupLDC_1=c(-.3,.3),se_GroupLDC_1=c(.1,.1),p_GroupLDC_1=c(.001,.02),q_GroupLDC_1=c(.006,.04),diff_GroupLDC_1=c(TRUE,TRUE),passed_ss_GroupLDC_1=c(TRUE,TRUE),lfc_GroupHDC2=c(.8,-.8),se_GroupHDC2=c(.2,.2),p_GroupHDC2=c(.002,.001),q_GroupHDC2=c(.009,.007),diff_GroupHDC2=c(TRUE,TRUE),passed_ss_GroupHDC2=c(TRUE,TRUE),check.names=FALSE)
pair<-res;pair$q_GroupLDC_1<-c(.03,.03);pair$passed_ss_GroupLDC_1<-c(TRUE,FALSE);pair$q_GroupHDC2<-c(.025,.005)
pair$lfc_GroupHDC2_GroupLDC_1<-c(1.1,-1.1);pair$se_GroupHDC2_GroupLDC_1<-c(.3,.3);pair$p_GroupHDC2_GroupLDC_1<-c(.003,.04);pair$q_GroupHDC2_GroupLDC_1<-c(.01,.2);pair$diff_GroupHDC2_GroupLDC_1<-c(TRUE,FALSE);pair$passed_ss_GroupHDC2_GroupLDC_1<-c(TRUE,TRUE)
zero<-data.frame(taxon=c('taxon_B','taxon_SZ','taxon_A'),CON=c(FALSE,TRUE,FALSE),LDC_1=c(FALSE,FALSE,FALSE),HDC2=c(FALSE,TRUE,FALSE),check.names=FALSE)
names(zero)[-1]<-paste0('structural_zero (Group = ',names(zero)[-1],')')
saveRDS(list(res=res,res_pair=pair,zero_ind=zero),a[1])""", self.folder / 'ancombc2_fit.rds')

    def redraw(self, **updates):
        target = self.root / ('redraw_' + str(len(list(self.root.glob('redraw_*')))))
        target.mkdir()
        request = {'analysis_run': str(self.source), 'diff_comparisons': 'all_pairs',
                   'diff_heatmap_scale': 'relative_abundance', 'bubble_color_mode': 'abundance',
                   'heatmap_low': '#000011', 'heatmap_mid': '#222233', 'heatmap_high': '#3333FF',
                   'bubble_low': '#112233', 'bubble_high': '#FFCC00', 'diff_plot_top_n': 30,
                   'group_style': {'group_order': ['HDC2', 'CON', 'LDC_1'],
                                   'group_colors': {'HDC2': '#F11010', 'CON': '#10F110', 'LDC_1': '#1010F1'}}}
        request.update(updates)
        parameters = target / 'options.json'
        parameters.write_text(json.dumps(request))
        result = subprocess.run([RSCRIPT, str(ROOT / 'differential_ancombc.R'), '--plot-ancom',
                                 str(parameters), str(target)], capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return target, json.loads((target / '02_differential' / 'summary.json').read_text())

    @staticmethod
    def table(path):
        with path.open(newline='') as stream:
            return list(csv.DictReader(stream, delimiter='\t'))

    def test_cached_all_pairs_abundance_and_structural_zero_are_exact(self):
        self.write_fit()
        before = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                  for path in self.source.rglob('*') if path.is_file()}
        target, summary = self.redraw()
        self.assertFalse(summary['statistics_refit']); self.assertFalse(summary['q_recomputed'])
        self.assertEqual(summary['original_analysis'], str(self.source))
        folder = target / '02_differential' / 'L6_genus'
        info = [json.loads(path.read_text()) for path in sorted((folder / 'comparisons').glob('*/comparison.json'))]
        self.assertEqual([item['contrast'] for item in info],
                         ['GroupHDC2', 'GroupHDC2_GroupLDC_1', 'GroupLDC_1'])
        self.assertTrue(all(item['model_result_source'] == 'res_pair' for item in info))
        pair = self.table(folder / 'comparisons' / 'pair_002' / 'comparison.tsv')
        a = next(row for row in pair if row['taxon'] == 'taxon_A')
        self.assertEqual(float(a['lfc']), 1.1); self.assertEqual(float(a['q']), .01)
        structural = next(row for row in pair if row['taxon'] == 'taxon_SZ')
        self.assertEqual(structural['source'], 'structural_zero_presence_absence')
        self.assertTrue(all(structural[field] == 'NA' for field in ('lfc', 'se', 'p', 'q')))
        control_pair = self.table(folder / 'comparisons' / 'pair_001' / 'comparison.tsv')
        self.assertNotIn('taxon_SZ', [row['taxon'] for row in control_pair], 'Both structural-zero groups do not define a differential finding')
        input_rows = self.table(folder / 'comparisons' / 'pair_003' / 'effect_plot_input.tsv')
        self.assertEqual([row['taxon'] for row in input_rows], ['taxon_A'], 'Failed sensitivity taxa must not be shown as robust')
        abundance = folder / 'abundance'
        means = self.table(abundance / 'group_mean_relative_abundance.tsv')
        self.assertEqual(float(next(row for row in means if row['taxon'] == 'taxon_A')['CON']), .5,
                         'Must average per-sample proportions rather than pool unequal count libraries')
        selected = self.table(abundance / 'taxon_selection.tsv')
        self.assertEqual(next(row for row in selected if row['taxon'] == 'taxon_SZ')['structural_zero'], 'TRUE')
        heat = (abundance / 'abundance_heatmap.svg').read_text()
        self.assertIn('[SZ]', heat)
        for color in ('#000011', '#3333FF', '#F11010', '#10F110', '#1010F1'):
            self.assertIn(color, heat)
        bubble = (abundance / 'abundance_bubble.svg').read_text()
        self.assertIn('#FFCC00', bubble)
        self.assertEqual(before, {path: hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in before},
                         'Cached plots must preserve original counts/metadata/options/model results')
        for filename in ('abundance_heatmap.png', 'abundance_heatmap.svg', 'abundance_bubble.png', 'abundance_bubble.svg'):
            self.assertGreater((abundance / filename).stat().st_size, 500)

    def test_two_group_primary_mapping_and_no_hit_selection(self):
        self.r("""a<-commandArgs(TRUE);source(a[1]);root<-a[2]
m<-read_count_input(file.path(root,'counts.tsv'));meta<-read_group_metadata(file.path(root,'metadata.tsv'),colnames(m),'CON',c('CON','LDC_1','HDC2'))
two<-droplevels(meta[meta$Group!='HDC2',,drop=FALSE]);plan<-official_contrast_map(two,'CON');stopifnot(nrow(plan)==1,plan$key=='GroupLDC_1',plan$positive_group=='LDC_1',plan$negative_group=='CON')
empty<-data.frame(taxon=c('taxon_A','taxon_B'),lfc_GroupLDC_1=c(NA_real_,1),se_GroupLDC_1=c(NA_real_,.2),p_GroupLDC_1=c(NA_real_,.3),q_GroupLDC_1=c(NA_real_,.7),diff_GroupLDC_1=c(FALSE,FALSE),passed_ss_GroupLDC_1=c(FALSE,TRUE),check.names=FALSE)
out<-file.path(root,'nohit');dir.create(out);opts<-list(reference_group='CON',alpha=.05,pseudo_sens=TRUE,font='sans')
report<-plot_ancom_level(list(res=empty),m[,rownames(two),drop=FALSE],two,out,opts,'synthetic two groups');stopifnot(report$result_source=='res',report$abundance$selected_total==0,report$abundance$shown==0)
map<-official_contrast_map(meta,'CON');stopifnot(map$key[3]=='GroupHDC2_GroupLDC_1')
zero<-data.frame(taxon=c('taxon_A','taxon_B'),CON=c(TRUE,TRUE),LDC_1=c(NA,TRUE),HDC2=c(FALSE,FALSE),check.names=FALSE)
names(zero)[-1]<-paste0('structural_zero (Group = ',names(zero)[-1],')')
sz<-structural_zero_data(list(zero_ind=zero),as.list(map[1,,drop=FALSE]),m,meta);stopifnot(is.null(sz));cat('Exact labels, two-group primary fallback and honest no-hit plots passed')""", ROOT / 'differential_ancombc.R', self.source)
        self.assertIn('No official differential findings selected',
                      (self.source / 'nohit' / 'abundance' / 'abundance_heatmap.svg').read_text())


if __name__ == '__main__':
    unittest.main()
