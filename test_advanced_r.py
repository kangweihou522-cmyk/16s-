"""Real vegan/ppcor analyses and cache-only rendering on independent synthetic samples."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parent
RSCRIPT = os.environ.get('WORKSTATION_TEST_RSCRIPT') or shutil.which('Rscript')
if not RSCRIPT and Path('/workspace/r-conda-review/Rscript').is_file():
    RSCRIPT = '/workspace/r-conda-review/Rscript'


def runtime_ready():
    if not RSCRIPT:
        return False
    probe = subprocess.run([RSCRIPT, '-e',
        "quit(status=if(all(vapply(c('vegan','permute','ppcor','ggplot2','jsonlite','svglite','systemfonts'),requireNamespace,logical(1),quietly=TRUE)))0 else 1)"],
        capture_output=True, timeout=30)
    return probe.returncode == 0


def expected_rank_correlation(values, traits, groups=None):
    """Independent rank/Pearson definition; Group adjustment centers ranks within group."""
    def ranks(items):
        result = [0.] * len(items)
        for value in sorted(set(items)):
            positions = [i for i, x in enumerate(items) if x == value]
            lower = sum(x < value for x in items)
            mean_rank = lower + (len(positions) + 1) / 2
            for i in positions:
                result[i] = mean_rank
        return result
    def centered(items):
        if groups is None:
            mean = sum(items) / len(items)
            return [x - mean for x in items]
        means = {g: sum(x for x, group in zip(items, groups) if group == g) / groups.count(g) for g in set(groups)}
        return [x - means[g] for x, g in zip(items, groups)]
    x, y = centered(ranks(values)), centered(ranks(traits))
    return sum(a * b for a, b in zip(x, y)) / math.sqrt(sum(a * a for a in x) * sum(b * b for b in y))


@unittest.skipUnless(runtime_ready(), 'Real R runtime with advanced packages is not installed')
class AdvancedRTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='16s-advanced-r-')
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.ids = [f'{group}{i:02}' for group in ('A', 'B', 'C') for i in range(1, 9)]
        self.features = ['k__Bacteria;f__One;g__Signal', 'k__Bacteria;f__One;g__shared',
                         'k__Bacteria;f__Two;g__shared', 'k__Bacteria;g__D',
                         'k__Bacteria;g__rare', 'k__Bacteria;g__excluded']
        self.counts = [[] for _ in self.features]
        rng = random.Random(53)
        self.metadata = []
        for i, sample in enumerate(self.ids):
            signal = rng.uniform(-2, 2)
            group = i // 8
            self.metadata.append([sample, ('A', 'B', 'C')[group], signal + rng.gauss(0, .07),
                                  30 + i % 7, ('A', 'B', 'C')[group], 0,
                                  signal + rng.gauss(0, .2)])
            values = [round(150 + 30 * signal + group * 30 + rng.gauss(0, 4)),
                      60 + rng.randrange(30), 50 + group * 20 + rng.randrange(15),
                      80 + rng.randrange(15), 2 if i < 3 else 0, 1000]
            for row, value in zip(self.counts, values):
                row.append(value)
        self.metadata_path = self.directory / 'metadata.csv'
        self.write_metadata()
        for filename, features in [('selected.tsv', 5), ('full.tsv', 6)]:
            with (self.directory / filename).open('w', newline='', encoding='utf8') as handle:
                writer = csv.writer(handle, delimiter='\t')
                writer.writerow(['FeatureID'] + self.ids)
                writer.writerows([[self.features[k]] + self.counts[k] for k in range(features)])
        self.options = {
            'counts_path': str(self.directory / 'selected.tsv'),
            'total_counts_path': str(self.directory / 'full.tsv'),
            'metadata_path': str(self.metadata_path),
            'modules': ['core', 'ordination', 'association'],
            'core_min_reads': 2, 'core_min_relative': .001, 'core_prevalence': .8,
            'feature_prevalence': .2, 'clr_pseudocount': .5,
            'predictors': ['Group'], 'covariates': ['Age'],
            'phenotypes': ['IL6', 'Constant'], 'association_covariates': ['Group'],
            'permutations': 49, 'seed': 42, 'font': 'DejaVu Serif',
            'group_order': ['C', 'A', 'B'],
            'group_colors': {'A': '#112233', 'B': '#445566', 'C': '#778899'},
        }

    def write_metadata(self):
        with self.metadata_path.open('w', newline='', encoding='utf8') as handle:
            writer = csv.writer(handle)
            writer.writerow(['SampleID', 'Group', 'IL6', 'Age', 'Alias', 'Constant', 'CRP (mg/L) + arbitrary'])
            writer.writerows(self.metadata)

    def execute(self, options=None, name='run', mode='--run', no_refit=False):
        options = self.options if options is None else options
        option_path = self.directory / (name + '.json')
        option_path.write_text(json.dumps(options), encoding='utf8')
        run = self.directory / name
        if no_refit:
            expression = "a<-commandArgs(TRUE);source(a[1]);adv_fit<-function(...)stop('FORBIDDEN REFIT');adv_main(a[2:4])"
            args = [RSCRIPT, '-e', expression, str(ROOT / 'advanced_analysis.R'), mode, str(option_path), str(run)]
        else:
            args = [RSCRIPT, str(ROOT / 'advanced_analysis.R'), mode, str(option_path), str(run)]
        result = subprocess.run(args, capture_output=True, text=True, timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return run

    @staticmethod
    def table(run, filename):
        with (run / '02_statistics' / filename).open(encoding='utf8', newline='') as handle:
            return list(csv.DictReader(handle, delimiter='\t'))

    def test_real_models_full_denominators_clr_and_known_adjusted_association(self):
        run = self.execute()
        diagnostics = {x['Module']: x for x in self.table(run, 'diagnostics.tsv')}
        self.assertTrue(all(x['Status'] == 'complete' for x in diagnostics.values()), diagnostics)
        audit = self.table(run, 'sample_audit.tsv')
        self.assertEqual([x['SampleID'] for x in audit], self.ids)
        self.assertEqual(float(audit[0]['FullLibrarySize']), sum(row[0] for row in self.counts))
        relative = self.table(run, 'sample_relative_abundance.tsv')
        denominator = sum(row[0] for row in self.counts)
        self.assertAlmostEqual(float(relative[0][self.ids[0]]), self.counts[0][0] / denominator, places=12)
        self.assertLess(sum(float(row[self.ids[0]]) for row in relative), .5)
        filtered = self.table(run, 'feature_filter.tsv')
        self.assertEqual(filtered[-1]['RetainedForMultivariateAndAssociation'], 'FALSE')
        self.assertNotEqual(filtered[1]['Label'], filtered[2]['Label'])
        core = self.table(run, 'core_by_group.tsv')
        rare_a = next(x for x in core if x['FeatureID'] == self.features[4] and x['Group'] == 'A')
        self.assertEqual(int(rare_a['DetectedSamples']), 3)
        self.assertAlmostEqual(float(rare_a['Prevalence']), 3 / 8)
        clr = self.table(run, 'clr_transformed.tsv')
        self.assertEqual(len(clr), 4)
        for sample in self.ids:
            self.assertAlmostEqual(sum(float(x[sample]) for x in clr), 0, places=10)
        raw_log = [math.log(row[0] + .5) for row in self.counts[:4]]
        self.assertAlmostEqual(float(clr[0][self.ids[0]]), raw_log[0] - sum(raw_log) / 4, places=11)
        assoc = self.table(run, 'phenotype_associations.tsv')
        signal = next(x for x in assoc if x['FeatureID'] == self.features[0] and x['Phenotype'] == 'IL6')
        transformed = self.table(run, 'association_transformed.tsv')
        self.assertEqual(len(transformed), 4)
        for sample_index, sample in enumerate(self.ids):
            logs = [math.log(row[sample_index] + .5) for row in self.counts[:4]]
            for row, logged in zip(transformed, logs):
                self.assertAlmostEqual(float(row[sample]), logged - sum(logs) / len(logs), places=11)
        expected = expected_rank_correlation([float(transformed[0][sample]) for sample in self.ids],
                                             [row[2] for row in self.metadata], [row[1] for row in self.metadata])
        self.assertAlmostEqual(float(signal['Rho']), expected, places=11)
        self.assertEqual(signal['Transform'], 'clr')
        self.assertGreater(float(signal['Rho']), .9)
        self.assertLess(float(signal['Q']), .001)
        self.assertEqual(signal['N'], '24')
        self.assertEqual(signal['CovariateDF'], '2')
        self.assertEqual(signal['SampleIDs'].split(';'), self.ids)
        self.assertTrue(all(x['Status'] == 'skipped' and 'Constant' in x['Reason'] for x in assoc if x['Phenotype'] == 'Constant'))
        tested = sorted((float(x['P']), float(x['Q'])) for x in assoc if x['Status'] == 'tested')
        expected = [min(1, min(tested[j][0] * len(tested) / (j + 1) for j in range(i, len(tested)))) for i in range(len(tested))]
        for (_, actual), correct in zip(tested, expected):
            self.assertAlmostEqual(actual, correct, places=12)
        permanova = self.table(run, 'permanova.tsv')
        self.assertEqual({x['Term'] for x in permanova if x['Status'] == 'tested'}, {'Group', 'Age'})
        distance = self.table(run, 'bray_curtis_distance.tsv')
        self.assertEqual([x['SampleID'] for x in distance], self.ids)
        self.assertTrue(all(float(row[sample]) == 0 for row, sample in zip(distance, self.ids)))
        for filename in ('permanova.tsv', 'permdisp.tsv', 'constrained_overall.tsv', 'constrained_marginal.tsv'):
            for row in self.table(run, filename):
                if row['Status'] == 'tested':
                    self.assertTrue(0 <= float(row['P']) <= 1)
                    self.assertTrue(math.isfinite(float(row['F'])))
        summary = json.loads((run / '03_figures/render_summary.json').read_text())
        self.assertEqual(summary['source'], 'advanced_analysis')
        self.assertEqual(summary['group_style']['group_order'], ['C', 'A', 'B'])
        self.assertEqual(len(summary['figures']), 5)
        stats_summary = json.loads((run / '02_statistics/analysis_summary.json').read_text())
        self.assertEqual(stats_summary['association_transform'], 'clr')
        self.assertIn('not absolute abundance', stats_summary['association_definition'])
        for figure in summary['figures']:
            for suffix in ('png', 'svg', 'pdf'):
                self.assertGreater((run / (figure['stem'] + '.' + suffix)).stat().st_size, 1000)

    def test_single_cached_render_uses_no_inputs_or_fit_and_preserves_every_statistic(self):
        run = self.execute()
        statistics = run / '02_statistics'
        before = {p.name: (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns) for p in statistics.iterdir()}
        Path(self.options['counts_path']).unlink()
        self.metadata_path.unlink()
        render = self.execute({
            'cache_run': str(run), 'plot_target': {'stem': '03_figures/CLR_PCA', 'profile': 'beta'},
            'group_order': ['B', 'A', 'C'], 'group_colors': {'A': '#FF0000', 'B': '#00FF00', 'C': '#0000FF'},
            'plot_settings': {'general': {'formats': 'png_svg_pdf', 'font_family': 'DejaVu Sans',
                                         'base_size': 14, 'theme': 'minimal'},
                              'beta': {'width': 7, 'height': 4, 'title': 'Saved CLR view', 'point_size': 5,
                                       'show_ellipse': False, 'show_sample_labels': True}},
        }, name='render', mode='--render', no_refit=True)
        after = {p.name: (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns) for p in statistics.iterdir()}
        self.assertEqual(before, after)
        self.assertFalse((render / '02_statistics').exists())
        self.assertEqual({p.name for p in (render / '03_figures').glob('*.png')}, {'CLR_PCA.png'})
        svg = (render / '03_figures/CLR_PCA.svg').read_text()
        self.assertIn('Saved CLR view', svg)
        self.assertIn('A01', svg)
        self.assertIn('#FF0000', svg)
        effective = json.loads((render / '03_figures/CLR_PCA_plot_settings.json').read_text())['effective']
        self.assertEqual((effective['width'], effective['height'], effective['font']), (7, 4, 'DejaVu Sans'))
        summary = json.loads((render / '03_figures/render_summary.json').read_text())
        self.assertFalse(summary['statistics_refit'])
        self.assertEqual(summary['advanced_run'], str(run))
        self.assertEqual(len(summary['figures']), 1)
        heatmap = self.execute({
            'cache_run': str(run), 'plot_target': {'stem': '03_figures/association_heatmap', 'profile': 'heatmap'},
            'association_transform': 'relative',  # Inferential settings cannot alter a cached view.
            'heatmap_low': '#006600', 'heatmap_mid': '#FFFFCC', 'heatmap_high': '#990099',
            'plot_settings': {'heatmap': {'top_n': 1, 'show_values': True}},
        }, name='heatmap', mode='--render', no_refit=True)
        heat_settings = json.loads((heatmap / '03_figures/association_heatmap_plot_settings.json').read_text())['effective']
        self.assertEqual(heat_settings['profile_settings']['top_n'], 1)
        self.assertEqual(heat_settings['colors']['heatmap_high'], '#990099')
        heat_svg = (heatmap / '03_figures/association_heatmap.svg').read_text()
        self.assertIn('g__Signal', heat_svg)
        self.assertNotIn('g__shared', heat_svg)
        self.assertIn('CLR:', heat_svg)
        self.assertEqual(before, {p.name: (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns) for p in statistics.iterdir()})

    def test_rank_deficiency_skips_inference_and_retains_descriptive_outputs(self):
        run = self.execute(self.options | {'covariates': ['Alias'], 'association_covariates': ['Group', 'Alias']})
        diagnostics = {x['Module']: x for x in self.table(run, 'diagnostics.tsv')}
        for module in ('permanova', 'constrained'):
            self.assertEqual(diagnostics[module]['Status'], 'skipped')
            self.assertIn('Rank deficient', diagnostics[module]['Reason'])
        self.assertEqual(diagnostics['core']['Status'], 'complete')
        self.assertEqual(diagnostics['clr_pca']['Status'], 'complete')
        rows = self.table(run, 'phenotype_associations.tsv')
        self.assertTrue(all(x['Status'] == 'skipped' and x['P'] == 'NA' and x['Q'] == 'NA' for x in rows))
        self.assertTrue(all('Rank deficient' in x['Reason'] for x in rows))
        self.assertTrue((run / '03_figures/CLR_PCA.svg').is_file())
        self.assertFalse((run / '03_figures/association_heatmap.svg').exists())

    def test_repeated_subjects_allow_core_and_pca_but_no_inferential_tests(self):
        run = self.execute(self.options | {'repeated_subjects': True})
        diagnostics = {x['Module']: x for x in self.table(run, 'diagnostics.tsv')}
        for module in ('core', 'clr_pca'):
            self.assertEqual(diagnostics[module]['Status'], 'complete')
        for module in ('permanova', 'permdisp', 'constrained', 'association'):
            self.assertEqual(diagnostics[module]['Status'], 'skipped')
        for filename in ('permanova.tsv', 'permdisp.tsv', 'constrained_overall.tsv'):
            self.assertFalse((run / '02_statistics' / filename).exists())
        rows = self.table(run, 'phenotype_associations.tsv')
        self.assertTrue(all(x['Status'] == 'skipped' and 'Independent samples only' in x['Reason'] and x['P'] == 'NA' for x in rows))

    def test_bray_dbrda_safe_variable_alias_and_trait_complete_cases(self):
        self.metadata[0][2] = ''
        self.write_metadata()
        run = self.execute(self.options | {
            'constrained_method': 'bray_dbrda', 'predictors': ['CRP (mg/L) + arbitrary'],
            'covariates': ['Group'], 'phenotypes': ['IL6'], 'association_method': 'spearman',
            'association_transform': 'relative',
        })
        diagnostics = {x['Module']: x for x in self.table(run, 'diagnostics.tsv')}
        self.assertEqual(diagnostics['constrained']['Status'], 'complete', diagnostics['constrained'])
        coordinates = self.table(run, 'constrained_coordinates.tsv')
        self.assertEqual(len(coordinates), 24)
        self.assertTrue(all(float(x['Axis2']) == 0 for x in coordinates))
        terms = self.table(run, 'constrained_marginal.tsv')
        self.assertIn('CRP (mg/L) + arbitrary', {x['Term'] for x in terms})
        assoc = self.table(run, 'phenotype_associations.tsv')
        self.assertTrue(all(x['N'] == '23' and x['CovariateDF'] == '0' and x['Covariates'] == '' for x in assoc))
        self.assertTrue(all(x['Transform'] == 'relative' for x in assoc))
        transformed = self.table(run, 'association_transformed.tsv')
        expected_values = [self.counts[0][i] / sum(row[i] for row in self.counts) for i in range(24)]
        for sample, expected_value in zip(self.ids, expected_values):
            self.assertAlmostEqual(float(transformed[0][sample]), expected_value, places=12)
        expected_rho = expected_rank_correlation(expected_values[1:], [row[2] for row in self.metadata[1:]])
        self.assertAlmostEqual(float(assoc[0]['Rho']), expected_rho, places=11)
        self.assertTrue(all('A01' not in x['SampleIDs'].split(';') for x in assoc))
        inclusion = self.table(run, 'association_sample_inclusion.tsv')
        self.assertEqual(next(x for x in inclusion if x['SampleID'] == 'A01')['Included'], 'FALSE')

    def test_clr_association_only_caches_inputs_and_requires_two_features(self):
        run = self.execute(self.options | {'modules': ['association'], 'phenotypes': ['IL6']})
        self.assertFalse((run / '02_statistics/clr_transformed.tsv').exists())
        transformed = self.table(run, 'association_transformed.tsv')
        self.assertEqual(len(transformed), 4)
        for i, sample in enumerate(self.ids):
            logs = [math.log(row[i] + .5) for row in self.counts[:4]]
            for row, logged in zip(transformed, logs):
                self.assertAlmostEqual(float(row[sample]), logged - sum(logs) / 4, places=11)
        with (self.directory / 'selected.tsv').open('w', newline='', encoding='utf8') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['FeatureID'] + self.ids)
            writer.writerow([self.features[0]] + self.counts[0])
        single = self.execute(self.options | {'modules': ['association'], 'phenotypes': ['IL6']}, name='single')
        diagnostic = self.table(single, 'diagnostics.tsv')[0]
        self.assertEqual(diagnostic['Status'], 'skipped')
        self.assertIn('at least two prevalent features', diagnostic['Reason'])

    def test_literal_identifiers_alias_collisions_and_nonfinite_test_guard(self):
        metadata = self.directory / 'literal.csv'
        metadata.write_text('SampleID,Group,Age\nNA,NA,32\n001,A,NA\n', encoding='utf8')
        expression = """a<-commandArgs(TRUE);source(a[1])
o<-adv_defaults(list(predictors=character(),association_covariates=character()))
stopifnot(o$association_transform=='clr',inherits(try(adv_defaults(list(association_transform='counts')),silent=TRUE),'try-error'))
m<-adv_metadata(a[2],c('NA','001'),o)
stopifnot(identical(rownames(m),c('NA','001')),as.character(m$Group[1])=='NA',is.na(m$Age[2]))
x<-data.frame(Df=c(1,1,1),F=c(3,4,Inf),row.names=c('v1','v2','Residual'))
x[['Pr(>F)']]<-c(.02,.04,.001)
d<-adv_anova_table(x,c(v2='v1',Age='v2'))
stopifnot(identical(d$Term,c('v2','Age','Residual')),is.na(d$P[3]),d$Status[3]=='not_tested',all(is.finite(d$Q[1:2])))
"""
        process = subprocess.run([RSCRIPT, '-e', expression, str(ROOT / 'advanced_analysis.R'), str(metadata)],
                                 capture_output=True, text=True, timeout=30)
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr)

    def test_degenerate_distances_cannot_produce_finite_inference(self):
        self.counts = [[10] * 24, [5] * 24, [0] * 24, [0] * 24, [0] * 24, [1000] * 24]
        for filename, features in [('selected.tsv', 5), ('full.tsv', 6)]:
            with (self.directory / filename).open('w', newline='', encoding='utf8') as handle:
                writer = csv.writer(handle, delimiter='\t')
                writer.writerow(['FeatureID'] + self.ids)
                writer.writerows([[self.features[k]] + self.counts[k] for k in range(features)])
        run = self.execute()
        diagnostics = {x['Module']: x for x in self.table(run, 'diagnostics.tsv')}
        self.assertEqual(diagnostics['core']['Status'], 'complete')
        for module in ('clr_pca', 'permanova', 'permdisp', 'constrained', 'association'):
            self.assertEqual(diagnostics[module]['Status'], 'skipped', diagnostics[module])
        for row in self.table(run, 'phenotype_associations.tsv'):
            self.assertEqual((row['Status'], row['P'], row['Q']), ('skipped', 'NA', 'NA'))
        self.assertTrue((run / '03_figures/core_heatmap.svg').is_file())


if __name__ == '__main__':
    unittest.main()
