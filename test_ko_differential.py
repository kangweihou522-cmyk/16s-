"""Real base-R KO inference checked against independent standard-library algebra."""
import csv
import itertools
import json
import math
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


def inverse(matrix):
    """Gauss-Jordan inverse, independent of R's QR-based lm implementation."""
    size = len(matrix)
    work = [list(row) + [float(i == j) for j in range(size)] for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(work[row][column]))
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [value / scale for value in work[column]]
        for row in range(size):
            if row != column:
                scale = work[row][column]
                work[row] = [a - scale * b for a, b in zip(work[row], work[column])]
    return [row[size:] for row in work]


def manual_contrast(values, design, contrast):
    """Normal equations and exact even-df Student-t CDF polynomial."""
    columns = list(zip(*design))
    covariance = inverse([[sum(a * b for a, b in zip(left, right)) for right in columns] for left in columns])
    xy = [sum(a * b for a, b in zip(column, values)) for column in columns]
    beta = [sum(a * b for a, b in zip(row, xy)) for row in covariance]
    residuals = [value - sum(a * b for a, b in zip(row, beta)) for row, value in zip(design, values)]
    df = len(values) - len(columns)
    variance = sum(value ** 2 for value in residuals) / df
    effect = sum(a * b for a, b in zip(contrast, beta))
    se = math.sqrt(variance * sum(contrast[i] * covariance[i][j] * contrast[j]
                                 for i in range(len(contrast)) for j in range(len(contrast))))
    assert df % 2 == 0
    r = abs(effect / se) / math.sqrt((effect / se) ** 2 + df)
    terms = [(math.comb((df - 2) // 2, k) * (-1) ** k) / (2 * k + 1) for k in range(df // 2)]
    p = 1 - sum(value * r ** (2 * k + 1) for k, value in enumerate(terms)) / sum(terms)
    return effect, se, p


def bh(values):
    ordered = sorted(range(len(values)), key=values.__getitem__)
    adjusted = [None] * len(values)
    running = 1.
    for rank in range(len(values), 0, -1):
        index = ordered[rank - 1]
        running = min(running, values[index] * len(values) / rank)
        adjusted[index] = running
    return adjusted


@unittest.skipUnless(RSCRIPT, 'Real Rscript is required')
class KoDifferentialTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.ids = ['001', 'NA', 'a3', 'b1', 'b2', 'b3', 'c1', 'c2', 'c3']
        self.groups = ['control (ref)'] * 3 + ['case + treatment'] * 3 + ['third / arm'] * 3
        self.rows = [
            ['K00001', .1, .2, .31, 1.2, .9, 1.7, .55, .8, .4],
            ['K00002', 1.8, 1.4, 1.7, .5, .7, .8, 1.1, .7, 1.5],
            ['K00003', .7, 0, .8, 1.4, .3, 0, .6, 1.1, .9],
            ['K00004', 0, 0, 0, 0, 0, 0, 0, 0, 0],
            ['K00005', .02, 0, 0, 0, 0, 0, 0, 0, 0],
        ]
        self.meta = [{'SampleID': sid, 'Group': group} for sid, group in zip(self.ids, self.groups)]

    def run_model(self, *, options=None, expect_success=True, run_name='run', rows=None, meta=None, ids=None):
        counts = self.root / 'counts.tsv'
        metadata = self.root / 'metadata.csv'
        samples = self.ids if ids is None else ids
        with counts.open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['FeatureID', *samples])
            writer.writerows(self.rows if rows is None else rows)
        records = self.meta if meta is None else meta
        fields = list(dict.fromkeys(field for record in records for field in record))
        with metadata.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader(); writer.writerows(records[::-1])
        opts = {'counts_path': str(counts), 'metadata_path': str(metadata), 'feature_prevalence': .2,
                'group_order': list(dict.fromkeys(self.groups))}
        opts.update(options or {})
        config = self.root / (run_name + '.json')
        config.write_text(json.dumps(opts))
        result = subprocess.run([RSCRIPT, str(ROOT / 'ko_differential.R'), '--run', str(config), str(self.root / run_name)],
                                cwd=ROOT, capture_output=True, text=True, timeout=60)
        if expect_success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
        return self.root / run_name / '02_statistics', result

    def table(self, path):
        with path.open(newline='') as handle:
            return list(csv.DictReader(handle, delimiter='\t'))

    def test_fractional_zero_inputs_all_pairs_bh_literal_ids_and_independent_ols(self):
        folder, _ = self.run_model()
        rows = self.table(folder / 'ko_differential.tsv')
        comparisons = json.loads((folder / 'comparisons.json').read_text())
        self.assertEqual([item['id'] for item in comparisons], ['C001', 'C002', 'C003'])
        self.assertEqual([(item['numerator'], item['denominator']) for item in comparisons],
                         [(self.groups[b * 3], self.groups[a * 3]) for a, b in itertools.combinations(range(3), 2)])
        self.assertEqual([row['SampleID'] for row in self.table(folder / 'sample_inclusion.tsv')], self.ids)
        raw = self.table(folder / 'raw_abundances.tsv')
        self.assertEqual(float(raw[0]['001']), .1)
        self.assertEqual(float(raw[0]['NA']), .2)
        totals = [sum(row[index + 1] for row in self.rows[:3]) for index in range(9)]
        logs = [[math.log2(row[i + 1] / totals[i] + 1e-6) for i in range(9)] for row in self.rows[:3]]
        expected_clr = [[value - sum(logs[k][i] for k in range(3)) / 3 for i, value in enumerate(row)] for row in logs]
        design = [[1., float(i >= 3 and i < 6), float(i >= 6)] for i in range(9)]
        contrasts = [[0., 1., 0.], [0., 0., 1.], [0., -1., 1.]]
        valid = [row for row in rows if row['Status'] == 'ok']
        self.assertEqual(len(valid), 9)
        for i, values in enumerate(expected_clr):
            for j, contrast in enumerate(contrasts):
                actual = next(row for row in valid if row['KO'] == self.rows[i][0] and row['comparison'] == f'C{j + 1:03d}')
                effect, se, p = manual_contrast(values, design, contrast)
                self.assertAlmostEqual(float(actual['effect']), effect, places=11)
                self.assertAlmostEqual(float(actual['SE']), se, places=11)
                self.assertAlmostEqual(float(actual['P']), p, places=11)
        adjusted = bh([float(row['P']) for row in valid])
        for row, q in zip(valid, adjusted):
            self.assertAlmostEqual(float(row['q']), q, places=12)
        self.assertEqual({row['Reason'] for row in rows if row['KO'] == 'K00004'}, {'all_zero'})
        self.assertEqual({row['Reason'] for row in rows if row['KO'] == 'K00005'}, {'below_prevalence'})
        report = json.loads((folder / 'transformation_report.json').read_text())
        self.assertEqual(report['multiplicity']['tested_hypotheses'], 9)
        self.assertEqual(report['clr_feature_universe'], ['K00001', 'K00002', 'K00003'])
        self.assertEqual(report['sample_ids'][:2], ['001', 'NA'])
        self.assertEqual(report['pseudocount'], 1e-6)
        self.assertIn('retained KOs', report['tss_denominator'])
        self.assertTrue((folder / 'sessionInfo.txt').is_file())

    def test_adjusted_numeric_and_categorical_covariates_use_safe_aliases(self):
        ages = [12, 15, 13, 17, 11, 16, 14, 19, 10]
        batches = ['lot + one', 'lot/two', 'lot + one', 'lot/two', 'lot + one', 'lot/two', 'lot + one', 'lot + one', 'lot/two']
        meta = [dict(row, **{'age (days)': age, 'batch + term': batch}) for row, age, batch in zip(self.meta, ages, batches)]
        folder, _ = self.run_model(meta=meta, options={'covariates': ['age (days)', {'name': 'batch + term', 'type': 'categorical'}]})
        report = json.loads((folder / 'transformation_report.json').read_text())
        self.assertEqual(report['effective_formula'], 'response ~ Group + V001 + V002')
        self.assertEqual(report['residual_df'], 4)
        clr = self.table(folder / 'clr_log2.tsv')
        design = [[1., float(3 <= i < 6), float(i >= 6), ages[i], float(batches[i] == 'lot/two')] for i in range(9)]
        expected = manual_contrast([float(clr[0][sid]) for sid in self.ids], design, [0, -1, 1, 0, 0])
        actual = next(row for row in self.table(folder / 'ko_differential.tsv') if row['KO'] == 'K00001' and row['comparison'] == 'C003')
        for key, value in zip(('effect', 'SE', 'P'), expected):
            self.assertAlmostEqual(float(actual[key]), value, places=10)

    def test_rank_deficient_and_saturated_designs_do_not_fabricate_pvalues(self):
        meta = [dict(row, confounded=int(i >= 3)) for i, row in enumerate(self.meta)]
        folder, _ = self.run_model(meta=meta, options={'covariates': ['confounded']})
        kept = [row for row in self.table(folder / 'ko_differential.tsv') if row['Status'] != 'excluded']
        self.assertEqual({row['Reason'] for row in kept}, {'rank_deficient_design'})
        self.assertTrue(all(row['P'] == row['q'] == row['effect'] == 'NA' for row in kept))
        meta = [dict(row, identity='person' + str(i)) for i, row in enumerate(self.meta)]
        folder, _ = self.run_model(meta=meta, options={'covariates': ['identity']}, run_name='saturated')
        self.assertTrue(all(row['P'] == 'NA' for row in self.table(folder / 'ko_differential.tsv')))
        # Two groups with one independent numeric covariate and zero residual df.
        ids = self.ids[:4]; groups = ['A', 'A', 'B', 'B']
        meta = [{'SampleID': sid, 'Group': group, 'x': i, 'z': i * i} for i, (sid, group) in enumerate(zip(ids, groups))]
        folder, _ = self.run_model(ids=ids, meta=meta, rows=[row[:5] for row in self.rows],
                                   options={'group_order': ['A', 'B'], 'covariates': ['x', 'z']}, run_name='no_df')
        self.assertEqual({row['Reason'] for row in self.table(folder / 'ko_differential.tsv') if row['Status'] != 'excluded'},
                         {'insufficient_residual_df'})

    def test_constant_and_zero_residual_features_have_na_tests(self):
        ids = self.ids[:6]
        meta = self.meta[:6]
        rows = [['K00001', 1, 1, 1, 1, 1, 1], ['K00002', 1, 1, 1, 2, 2, 2]]
        folder, _ = self.run_model(rows=rows, ids=ids, meta=meta, options={'group_order': list(dict.fromkeys(self.groups[:6]))})
        self.assertEqual({row['Reason'] for row in self.table(folder / 'ko_differential.tsv')}, {'zero_residual_variance'})
        folder, _ = self.run_model(rows=[['K00001', 1, 1, 1, 1, 1, 1], ['K00002', 2, 2, 2, 2, 2, 2]],
                                   ids=ids, meta=meta, options={'group_order': list(dict.fromkeys(self.groups[:6]))}, run_name='constant')
        self.assertEqual({row['Reason'] for row in self.table(folder / 'ko_differential.tsv')}, {'constant_response'})

    def test_repeated_subjects_sample_mismatch_and_invalid_counts_are_rejected(self):
        for field in ['AnimalID', 'SubjectID']:
            meta = [dict(row, **{field: 'same' if i < 2 else 'unique' + str(i)}) for i, row in enumerate(self.meta)]
            folder, result = self.run_model(meta=meta, expect_success=False, run_name=field)
            self.assertIn('repeated ' + field, result.stderr)
            self.assertFalse(folder.exists())
        _, result = self.run_model(meta=self.meta[:-1], expect_success=False, run_name='mismatch')
        self.assertIn('exact count sample set', result.stderr)
        for value in [-.1, float('nan'), float('inf')]:
            rows = [row[:] for row in self.rows]; rows[0][1] = value
            _, result = self.run_model(rows=rows, expect_success=False, run_name='invalid')
            self.assertIn('finite nonnegative', result.stderr)
        rows = [row[:] for row in self.rows]
        for row in rows: row[1] = 0
        _, result = self.run_model(rows=rows, expect_success=False, run_name='zero_total')
        self.assertIn('positive finite total', result.stderr)

    def test_two_per_group_warns_missing_covariates_reject_and_saved_statistics_are_immutable(self):
        indices = [0, 1, 3, 4]
        folder, result = self.run_model(ids=[self.ids[i] for i in indices], meta=[self.meta[i] for i in indices],
                                        rows=[[row[0], *(row[i + 1] for i in indices)] for row in self.rows],
                                        options={'group_order': list(dict.fromkeys(self.groups[:6]))})
        self.assertIn('only two samples', result.stderr)
        before = (folder / 'ko_differential.tsv').read_bytes()
        _, result = self.run_model(expect_success=False)
        self.assertIn('immutable', result.stderr)
        self.assertEqual((folder / 'ko_differential.tsv').read_bytes(), before)
        meta = [dict(row, age='' if i == 0 else i) for i, row in enumerate(self.meta)]
        _, result = self.run_model(meta=meta, options={'covariates': ['age']}, expect_success=False, run_name='missing')
        self.assertIn('Missing selected covariate', result.stderr)


if __name__ == '__main__':
    unittest.main()
