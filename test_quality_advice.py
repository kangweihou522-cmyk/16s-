"""Meaningful fixtures for the cautious quality-summary recommendation."""
import csv
import json
import tempfile
import unittest
from pathlib import Path

from quality_advice import recommend_parameters


class QualityAdviceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def summary(self, name="forward.tsv", *, q25=None, counts=None, row_oriented=False):
        q25 = q25 if q25 is not None else [35] * 120
        counts = counts if counts is not None else [1000] * len(q25)
        rows = [[""] + list(range(1, len(q25) + 1)), ["count"] + counts]
        for label, offset in (("2%", -5), ("9%", -2), ("25%", 0), ("50%", 1),
                              ("75%", 2), ("91%", 3), ("98%", 4)):
            rows.append([label] + [value + offset for value in q25])
        if row_oriented:
            rows = list(zip(*rows))
        path = self.root / name
        with path.open("w", newline="", encoding="utf-8") as stream:
            csv.writer(stream, delimiter="\t").writerows(rows)
        return path

    def test_read_length_coverage_uses_sampled_counts(self):
        result = recommend_parameters(self.summary(counts=[1000] * 95 + [990] * 5 + [989] * 20))
        self.assertEqual(result["parameters"]["trunc_f"], 100)
        evidence = result["evidence"]["forward"]
        self.assertEqual(evidence["coverage_limit"], 100)
        self.assertEqual(evidence["retained_read_fraction"], 0.99)

    def test_sustained_deterioration_excludes_first_low_median_window(self):
        result = recommend_parameters(self.summary(q25=[35] * 90 + [20] * 30))
        self.assertEqual(result["parameters"]["trunc_f"], 85)
        evidence = result["evidence"]["forward"]
        self.assertEqual(evidence["quality_window_start"], 86)
        self.assertEqual(evidence["quality_window_median_q25"], 27.5)

    def test_isolated_dips_and_windows_at_threshold_preserve_length(self):
        values = [35] * 120
        values[20] = values[60] = values[90] = 10
        self.assertEqual(recommend_parameters(self.summary(q25=values))["parameters"]["trunc_f"], 120)
        self.assertEqual(recommend_parameters(self.summary(q25=[30] * 120))["parameters"]["trunc_f"], 120)

    def test_known_overlap_failure_does_not_extend_candidates(self):
        forward = self.summary()
        reverse = self.summary("reverse.tsv", q25=[35] * 90 + [20] * 30)
        result = recommend_parameters(forward, reverse, amplicon="200", cpu_count=16)
        self.assertEqual(result["parameters"]["trunc_r"], 85)
        self.assertEqual(result["evidence"]["overlap"]["predicted_bases"], 5)
        self.assertTrue(result["requires_review"])
        self.assertIn("5 bp", result["warnings"][0])
        self.assertEqual(result["parameters"]["amplicon"], "200")
        self.assertEqual(result["parameters"]["threads"], 4)

    def test_known_overlap_at_minimum_is_accepted(self):
        result = recommend_parameters(self.summary(), self.summary("reverse.tsv"), amplicon=228)
        self.assertEqual(result["evidence"]["overlap"]["predicted_bases"], 12)
        self.assertTrue(result["evidence"]["overlap"]["sufficient"])
        self.assertFalse(result["requires_review"])

    def test_unknown_amplicon_never_fabricates_length_or_depth(self):
        result = recommend_parameters(self.summary(), self.summary("reverse.tsv"), cpu_count=2)
        parameters = result["parameters"]
        self.assertEqual(parameters["amplicon"], "")
        self.assertEqual(parameters["max_depth"], "")
        self.assertEqual(parameters["min_frequency"], 1)
        self.assertEqual(parameters["threads"], 2)
        self.assertEqual((parameters["trim_f"], parameters["trim_r"]), (0, 0))
        self.assertIsNone(result["evidence"]["overlap"]["predicted_bases"])
        self.assertTrue(result["requires_review"])
        self.assertIn("实际拼接保留率", result["warnings"][0])
        json.dumps(result, allow_nan=False)

    def test_single_end_does_not_invent_reverse_or_merging(self):
        result = recommend_parameters(self.summary(row_oriented=True), cpu_count=0)
        self.assertEqual(result["parameters"]["trunc_r"], "")
        self.assertEqual(result["parameters"]["threads"], 1)
        self.assertIsNone(result["evidence"]["reverse"])
        self.assertFalse(result["requires_review"])

    def test_untrustworthy_summaries_are_rejected_with_diagnostics(self):
        variants = {
            "zero": dict(counts=[0] * 120),
            "increase": dict(counts=[1000] * 100 + [1001] * 20),
            "nan": dict(q25=[35] * 119 + [float("nan")]),
            "short": dict(q25=[35] * 40),
            "early_deterioration": dict(q25=[20] * 120),
        }
        for name, options in variants.items():
            with self.subTest(name=name), self.assertRaises(ValueError):
                recommend_parameters(self.summary(name + ".tsv", **options))
        bad = self.root / "bad.tsv"
        for contents in ("", "\t1\t2\ncount\t1000\n", "\t1\ncount\t1000\n25%\t30\n"):
            bad.write_text(contents)
            with self.subTest(contents=contents), self.assertRaises(ValueError):
                recommend_parameters(bad)
        with self.assertRaisesRegex(ValueError, "cannot read"):
            recommend_parameters(self.root / "missing.tsv")

    def test_invalid_positions_and_unordered_quantiles_reject(self):
        path = self.summary()
        path.write_text(path.read_text().replace("\t1\t2\t3", "\t1\t3\t2", 1))
        with self.assertRaisesRegex(ValueError, "positions"):
            recommend_parameters(path)
        path = self.summary()
        path.write_text(path.read_text().replace("50%\t36", "50%\t20", 1))
        with self.assertRaisesRegex(ValueError, "unordered"):
            recommend_parameters(path)

    def test_invalid_amplicon_rejects_instead_of_guessing(self):
        forward = self.summary()
        for length in (0, -1, "unknown", "nan", 150.5, True):
            with self.subTest(length=length), self.assertRaises(ValueError):
                recommend_parameters(forward, amplicon=length)


if __name__ == "__main__":
    unittest.main()
