"""Input/orchestration checks and real R plotting, not synthetic ANCOM validation."""
import csv
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import differential_analysis as differential

ROOT = Path(__file__).resolve().parent
SAMPLES = [{'SampleID': '0002', 'Group': 'CON'}, {'SampleID': '0001', 'Group': 'CON'},
           {'SampleID': 'CAP01', 'Group': 'A10'}, {'SampleID': 'CAP02', 'Group': 'A10'}]
RSCRIPT = os.environ.get('WORKSTATION_TEST_RSCRIPT') or shutil.which('Rscript')
if not RSCRIPT and Path('/workspace/r-conda-review/Rscript').is_file():
    RSCRIPT = '/workspace/r-conda-review/Rscript'


class DifferentialTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / 'source'; self.source.mkdir()
        self.run = self.root / 'differential'; self.run.mkdir()
        self.rows = [['d__Bacteria;p__Firmicutes', 10, 11, 3, 4],
                     ['d__Bacteria;p__Other', 20, 21, 23, 24]]
        self.ids = ['0001', '0002', 'CAP02', 'CAP01']

    def counts(self, filename, rows=None):
        filename.parent.mkdir(parents=True, exist_ok=True)
        with filename.open('w', encoding='utf8', newline='') as f:
            writer = csv.writer(f, delimiter='\t'); writer.writerow(['FeatureID', *self.ids])
            writer.writerows(self.rows if rows is None else rows)

    def manifest(self):
        taxonomy = ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Order;f__Family;g__Genus',
                    'd__Bacteria;p__Other;c__Other;o__Other;f__Other;g__Other']
        raw = self.source / 'exported-table' / 'feature-table.tsv'
        self.counts(raw, [['ASV' + str(index + 1), *row[1:]] for index, row in enumerate(self.rows)])
        taxfile = self.source / 'exported-taxonomy' / 'taxonomy.tsv'; taxfile.parent.mkdir()
        taxfile.write_text('FeatureID\tTaxon\tConfidence\n' + '\n'.join(
            'ASV' + str(index + 1) + '\t' + taxon + '\t0.99' for index, taxon in enumerate(taxonomy)) + '\n')
        entries = []
        for level in (2, 6):
            directory = self.source / '10_abundance' / ('L' + str(level))
            self.counts(directory / 'counts.tsv', [[(';'.join(taxonomy[index].split(';')[:level])), *row[1:]]
                                                 for index, row in enumerate(self.rows)])
            (directory / 'counts.qza').write_text('fixture artifact, never passed to real QIIME')
            entries.append({'level': level, 'counts': {'tsv': str((directory / 'counts.tsv').relative_to(self.source))},
                            'counts_qza': str((directory / 'counts.qza').relative_to(self.source))})
        (self.source / '10_abundance' / 'manifest.json').write_text(json.dumps({'levels': entries,
            'source': 'filtered_unrarefied', 'source_run': str(self.source),
            'source_sha256': {'filtered-table.tsv': hashlib.sha256(raw.read_bytes()).hexdigest(),
                              'taxonomy.tsv': hashlib.sha256(taxfile.read_bytes()).hexdigest()}}))

    def api(self, r_result=False):
        calls = []
        def q(job, config, *args):
            calls.append(('q', args))
        def export(job, config, input_path, output_path):
            calls.append(('export', str(input_path), str(output_path)))
            output = Path(output_path); output.mkdir(parents=True)
            (output / 'ancom.tsv').write_text('\tW\tReject null hypothesis\nd__Bacteria;p__Firmicutes\t1\tFalse\n')
        def command(job, args, cwd=None, process_env=None):
            calls.append(('command', args, process_env))
            if r_result and '--run' in args:
                directory = Path(job['run']) / '02_differential'; directory.mkdir(exist_ok=True)
                (directory / 'summary.json').write_text(json.dumps({'method': 'ANCOM-BC2', 'levels': []}))
        api = SimpleNamespace(ROOT=ROOT, q=q, export_artifact=export, command=command,
              conda_cmd=lambda c: ['conda', 'run', '-n', c['env']],
              r_command=lambda c, *args: ['conda', 'run', '-n', c['r_env'], 'Rscript', *map(str, args)],
              r_execution_env=lambda c: {'R_LIBS_USER': '/dev/null'},
              artifact_path=lambda p, n: Path(p) / n,
              begin_stage=lambda *a: calls.append(('begin', a[1])),
              finish_stage=lambda *a: calls.append(('finish', a[1])))
        return api, calls

    def test_literal_ids_and_groups_align_to_counts_without_numeric_order_inference(self):
        rows, sizes, excluded = differential.aligned_metadata(self.ids, self.rows, SAMPLES)
        self.assertEqual(rows, [('0001', 'CON'), ('0002', 'CON'), ('CAP02', 'A10'), ('CAP01', 'A10')])
        self.assertEqual(sizes, {'CON': 2, 'A10': 2}); self.assertEqual(excluded, [])
        opts = differential.options({'diff_method': 'ancombc2', 'diff_alpha': .02, 'diff_prevalence': .2}, SAMPLES)
        self.assertEqual(opts['reference_group'], 'CON'); self.assertEqual(opts['groups'], ['CON', 'A10'])
        self.assertEqual((opts['alpha'], opts['prv_cut']), (.02, .2))
        self.assertFalse(opts['global']); self.assertFalse(opts['pairwise'])

    def test_reject_relative_counts_unknown_samples_zero_libraries_and_unsupported_methods(self):
        counts = self.root / 'invalid.tsv'; self.counts(counts, [['taxon', .5, .2, 3, 4]])
        with self.assertRaisesRegex(ValueError, '整数'):
            differential.read_counts(counts)
        with self.assertRaisesRegex(ValueError, '分组'):
            differential.aligned_metadata(['unknown', *self.ids[1:]], self.rows, SAMPLES)
        with self.assertRaisesRegex(ValueError, '零文库'):
            differential.aligned_metadata(self.ids, [['taxon', 0, 1, 1, 1]], SAMPLES)
        for method in ('ancom2', 'original_ancom_ii', '', 'ANCOMBC2'):
            with self.assertRaises(ValueError):
                differential.options({'diff_method': method}, SAMPLES)
        with self.assertRaisesRegex(ValueError, '参照'):
            differential.options({'diff_method': 'ancombc2', 'reference_group': 'fake'}, SAMPLES)

    def test_omitted_saved_sample_is_reported_without_dropping_unknown_table_ids(self):
        _, _, excluded = differential.aligned_metadata(self.ids, self.rows,
            SAMPLES + [{'SampleID': 'dropped', 'Group': 'CON'}])
        self.assertEqual(excluded[0]['SampleID'], 'dropped')

    def test_three_nominal_groups_request_global_and_pairwise(self):
        opts = differential.options({'diff_method': 'ancombc2', 'reference_group': 'A10'},
            SAMPLES + [{'SampleID': 'E1', 'Group': 'A2'}, {'SampleID': 'E2', 'Group': 'A2'}])
        self.assertEqual(opts['groups'], ['CON', 'A10', 'A2'])
        self.assertEqual(opts['reference_group'], 'A10')
        self.assertTrue(opts['global']); self.assertTrue(opts['pairwise'])

    def test_qiime_ancom_uses_unrarefied_saved_counts_and_official_flags(self):
        self.manifest(); api, calls = self.api()
        config = {'diff_method': 'qiime_ancom', 'diff_levels': [2, 6], 'abundance_run': str(self.source)}
        saved = {'samples': SAMPLES, 'env': 'qiime-2025.7', 'r_env': '16s-workstation-r'}
        result = differential.run_differential({'run': str(self.run)}, config, saved, self.source, api)
        qcalls = [x[1] for x in calls if x[0] == 'q']
        self.assertEqual(sum(x[:2] == ('composition', 'add-pseudocount') for x in qcalls), 2)
        ancom = [x for x in qcalls if x[:2] == ('composition', 'ancom')]
        self.assertEqual(len(ancom), 2)
        self.assertTrue(all(x[x.index('--m-metadata-column') + 1] == 'Group' for x in ancom))
        self.assertTrue(all('raref' not in str(x).lower() for x in qcalls))
        opts = json.loads((self.run / 'differential_options.json').read_text())
        self.assertEqual(opts['group_sizes'], {'CON': 2, 'A10': 2})
        self.assertEqual(opts['p_adj_method'], 'official ANCOM default')
        self.assertIsNone(opts['prv_cut']); self.assertFalse(opts['pairwise'])
        ids, counts = differential.read_counts(self.run / '01_counts' / 'L2_phylum' / 'counts.tsv')
        self.assertEqual(ids, self.ids); self.assertEqual(counts, self.rows)
        self.assertEqual(result['method'], 'QIIME 2 ANCOM')
        self.assertTrue(any('--plot-qiime' in x[1] for x in calls if x[0] == 'command'))

    def test_ancombc2_uses_dedicated_ancom_env_and_preserves_raw_counts(self):
        self.manifest(); api, calls = self.api(r_result=True)
        config = {'diff_method': 'ancombc2', 'diff_levels': [2], 'abundance_run': str(self.source), 'ancom_env': 'my-ancom'}
        saved = {'samples': SAMPLES, 'env': 'qiime-2025.7', 'r_env': '16s-workstation-r'}
        differential.run_differential({'run': str(self.run)}, config, saved, self.source, api)
        commands = [x for x in calls if x[0] == 'command']
        self.assertEqual(len(commands), 2)
        self.assertEqual([x[1][3] for x in commands], ['my-ancom', 'my-ancom'])
        self.assertEqual([x[1][6] for x in commands], ['--check', '--run'])
        self.assertTrue(all(x[2]['R_LIBS_USER'] == '/dev/null' for x in commands))
        self.assertEqual(differential.read_counts(self.run / '01_counts' / 'L2_phylum' / 'counts.tsv')[1], self.rows)

    def test_manifest_paths_cannot_escape_saved_output(self):
        self.manifest(); manifest = self.source / '10_abundance' / 'manifest.json'
        content = json.loads(manifest.read_text()); content['levels'][0]['counts']['tsv'] = '../../outside.tsv'
        manifest.write_text(json.dumps(content))
        with self.assertRaisesRegex(ValueError, '超出'):
            differential._manifest_tables({'abundance_run': str(self.source)}, self.source, [2])

    def test_cache_rejects_wrong_batch_and_count_redistribution_with_matching_totals(self):
        self.manifest()
        wrong_parent = self.root / 'other_parent'; wrong_parent.mkdir()
        with self.assertRaisesRegex(ValueError, '批次'):
            differential._manifest_tables({'abundance_run': str(self.source)}, wrong_parent, [2])
        filename = self.source / '10_abundance' / 'L2' / 'counts.tsv'
        self.counts(filename, [[self.rows[0][0], 11, 11, 3, 4], [self.rows[1][0], 19, 21, 23, 24]])
        with self.assertRaisesRegex(ValueError, '内容不一致'):
            differential._manifest_tables({'abundance_run': str(self.source)}, self.source, [2])

    def test_cache_rejects_modified_parent_count_export(self):
        self.manifest()
        raw = self.source / 'exported-table' / 'feature-table.tsv'
        raw.write_text(raw.read_text().replace('10\t11', '12\t11'))
        with self.assertRaisesRegex(ValueError, '来源'):
            differential._manifest_tables({'abundance_run': str(self.source)}, self.source, [2])

    @unittest.skipUnless(RSCRIPT, 'A real R runtime is not installed')
    def test_real_r_parsing_count_and_metadata_validation(self):
        counts = self.root / 'counts.tsv'; self.counts(counts)
        metadata = self.root / 'metadata.tsv'
        metadata.write_text('sample-id\tGroup\n#q2:types\tcategorical\n0002\tCON\n0001\tCON\nCAP01\tA10\nCAP02\tA10\n')
        expr = "source(commandArgs(TRUE)[1]); m<-read_count_input(commandArgs(TRUE)[2]); d<-read_group_metadata(commandArgs(TRUE)[3],colnames(m),'CON',c('CON','A10')); stopifnot(identical(rownames(d),c('0001','0002','CAP02','CAP01')),identical(levels(d$Group),c('CON','A10')),!is.ordered(d$Group),m[1,1]==10); cat('real R input checks passed')"
        result = subprocess.run([RSCRIPT, '-e', expr, str(ROOT / 'differential_ancombc.R'), str(counts), str(metadata)],
                                capture_output=True, text=True, timeout=45)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('real R input checks passed', result.stdout)

    @unittest.skipUnless(RSCRIPT, 'A real R runtime is not installed')
    def test_real_r_model_call_uses_matrix_with_explicit_aligned_metadata(self):
        counts=self.root/'counts.tsv';self.counts(counts)
        metadata=self.root/'metadata.tsv'
        metadata.write_text('sample-id\tGroup\n0002\tCON\n0001\tCON\nCAP01\tA10\nCAP02\tA10\n')
        opts=differential.options({'diff_method':'ancombc2','threads':'1'},SAMPLES)
        options_file=self.root/'model_options.json';options_file.write_text(json.dumps(opts))
        expr="""a<-commandArgs(TRUE);source(a[1]);m<-read_count_input(a[2]);d<-read_group_metadata(a[3],colnames(m),'CON',c('CON','A10'));o<-jsonlite::fromJSON(a[4],simplifyVector=FALSE)
seen<-NULL;spy<-function(...){seen<<-list(...);'interface-only, not fitted statistics'}
x<-fit_ancom_model(m,d,o,fit_function=spy)
stopifnot(is.matrix(seen$data),identical(seen$data,m),identical(seen$meta_data,d),identical(seen$taxa_are_rows,TRUE),identical(seen$fix_formula,'Group'),is.null(seen$tax_level),identical(seen$p_adj_method,o$p_adj_method),identical(seen$pseudo_sens,o$pseudo_sens),identical(seen$prv_cut,o$prv_cut),identical(seen$lib_cut,o$lib_cut),identical(seen$struc_zero,o$struc_zero),identical(seen$neg_lb,o$neg_lb),identical(seen$alpha,o$alpha),identical(seen$n_cl,o$threads),identical(seen$global,o$global),identical(seen$pairwise,o$pairwise),identical(seen$iter_control,o$iter_control),identical(seen$em_control,o$em_control),identical(seen$mdfdr_control,o$mdfdr_control))
bad<-d[rev(rownames(d)),,drop=FALSE];failure<-tryCatch(fit_ancom_model(m,bad,o,fit_function=spy),error=function(e)conditionMessage(e));stopifnot(grepl('match exactly',failure));cat('matrix metadata interface checked; no model fitted')"""
        result=subprocess.run([RSCRIPT,'-e',expr,str(ROOT/'differential_ancombc.R'),str(counts),str(metadata),str(options_file)],capture_output=True,text=True,timeout=45)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertIn('no model fitted',result.stdout)

    @unittest.skipUnless(RSCRIPT, 'A real R runtime is not installed')
    def test_real_r_dependency_failure_is_explicit_without_fake_statistics(self):
        probe = subprocess.run([RSCRIPT, '-e', 'quit(status=if(requireNamespace("ANCOMBC",quietly=TRUE)) 0 else 1)'],
                               capture_output=True, text=True, timeout=30)
        if probe.returncode == 0:
            self.skipTest('ANCOMBC exists in selected runtime; dependency failure fixture is not applicable')
        result = subprocess.run([RSCRIPT, str(ROOT / 'differential_ancombc.R'), '--check', 'unused.json', str(self.run)],
                                capture_output=True, text=True, timeout=45)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Missing isolated ANCOM-BC2 dependencies', result.stderr)
        self.assertIn('install_differential.sh', result.stderr)
        self.assertFalse((self.run / '02_differential').exists())

    @unittest.skipUnless(RSCRIPT, 'A real R runtime is not installed')
    def test_real_r_plot_uses_official_w_fixture_and_outputs_png_svg(self):
        self.manifest(); api, _ = self.api()
        differential.run_differential({'run': str(self.run)},
            {'diff_method': 'qiime_ancom', 'diff_levels': [2], 'abundance_run': str(self.source)},
            {'samples': SAMPLES, 'env': 'fake-qiime', 'r_env': 'fake-r', 'font': 'sans'}, self.source, api)
        result = subprocess.run([RSCRIPT, str(ROOT / 'differential_ancombc.R'), '--plot-qiime',
                                str(self.run / 'differential_options.json'), str(self.run)],
                                capture_output=True, text=True, timeout=45)
        if 'Missing isolated R plotting dependencies' in result.stderr:
            self.skipTest('Real R plotting dependencies are not installed')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        folder = self.run / '02_differential' / 'L2_phylum'
        self.assertGreater((folder / 'ancom_W.png').stat().st_size, 1000)
        self.assertGreater((folder / 'ancom_W.svg').stat().st_size, 1000)

    @unittest.skipUnless(RSCRIPT, 'A real R runtime is not installed')
    def test_real_r_na_result_preservation_and_empty_effect_plot(self):
        expr = "source(commandArgs(TRUE)[1]); out<-commandArgs(TRUE)[2]; dir.create(out,showWarnings=FALSE); r<-data.frame(taxon=c('lineage;g__Unknown','lineage;g__Known'),lfc_GroupA10=c(NA,NA),se_GroupA10=c(NA,NA),q_GroupA10=c(NA,NA),check.names=FALSE); write_result(r,file.path(out,'res.tsv')); plot_effects(r,out,'CON','sans'); stopifnot(any(grepl('NA',readLines(file.path(out,'res.tsv'))))); cat('NA preserved and empty plot rendered')"
        result = subprocess.run([RSCRIPT, '-e', expr, str(ROOT / 'differential_ancombc.R'), str(self.run)],
                                capture_output=True, text=True, timeout=45)
        if 'there is no package called' in result.stderr:
            self.skipTest('Real R plotting dependencies are not installed')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('NA preserved', result.stdout)
        self.assertGreater((self.run / 'no_estimable_effects.png').stat().st_size, 1000)
        self.assertGreater((self.run / 'no_estimable_effects.svg').stat().st_size, 1000)


if __name__ == '__main__':
    unittest.main()
