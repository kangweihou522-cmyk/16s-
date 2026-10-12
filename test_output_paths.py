"""Output isolation checks using real files and mocked QIIME commands."""
import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from test_mapping import fastq
from test_quality_pipeline import materialize_qiime_outputs

SPEC = importlib.util.spec_from_file_location("output_workstation", Path(__file__).with_name("app.py"))
app = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(app)


class OutputDirectoryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.raw = self.root / "原始数据"
        self.raw.mkdir()
        self.meta = self.raw / "metadata.csv"

    def sheet(self, header, rows):
        with self.meta.open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(header)
            writer.writerows(rows)

    def config(self, output):
        return dict(input=str(self.raw), metadata=str(self.meta), mode="paired",
                    output=str(output), threads="1", primer_mode="removed")

    def original_files(self):
        return {file: file.read_bytes() for file in self.raw.rglob("*") if file.is_file()}

    def assert_originals_preserved(self, original):
        self.assertEqual({file: file.read_bytes() for file in original}, original)

    def test_two_column_sheet_allows_results_subdirectory_and_independent_runs(self):
        forward = fastq(self.raw / "A1_R1.fastq.gz", mate=1)
        reverse = fastq(self.raw / "A1_R2.fastq.gz", mate=2)
        self.sheet(["SampleID", "Group"], [["A1", "CON"]])
        output = self.raw / "results"
        output.mkdir()
        (output / "user-note.txt").write_text("Preserve this user's existing output.", encoding="utf8")
        previous = output / "previous-analysis"
        previous.mkdir()
        (previous / "unused_R1.fastq").write_text("Invalid FASTQ must never enter new input discovery.")
        originals = self.original_files()
        config = self.config(output)
        jobs = [{}, {}]
        with patch.object(app, "q", side_effect=materialize_qiime_outputs), patch.object(app, "command"):
            for job in jobs:
                app.preprocess(job, config, quality_only=True)
        runs = [Path(job["run"]) for job in jobs]
        self.assertNotEqual(runs[0], runs[1])
        for run in runs:
            self.assertEqual(run.parent, output)
            self.assertTrue(run.name.startswith("run_"))
            self.assertTrue((run / "quality_complete.json").is_file())
            with app.artifact_path(run, "manifest.tsv").open(encoding="utf8") as stream:
                rows = list(csv.reader(stream, delimiter="\t"))
            self.assertEqual(rows[1], ["A1", str(forward), str(reverse)])
        # Automatic discovery still inspects only raw-folder files even after
        # generating nested result directories and manifests there.
        found = app.discover(config)
        self.assertEqual([row["SampleID"] for row in found], ["A1"])
        self.assertEqual(found[0]["files"], {"1": str(forward), "2": str(reverse)})
        self.assert_originals_preserved(originals)

    def test_four_column_sheet_allows_nested_results_and_explicit_input_selection(self):
        forward = fastq(self.raw / "reads" / "arbitrary-left.fq.gz", mate=1)
        reverse = fastq(self.raw / "reads" / "arbitrary-right.fastq", mate=2)
        (self.raw / "unlisted.fastq").write_text("Unlisted input should not be parsed.")
        self.sheet(["SampleID", "Group", "forward", "reverse"],
                   [["CAP1", "处理", "reads/arbitrary-left.fq.gz", "reads/arbitrary-right.fastq"]])
        output = self.raw / "analysis" / "results"
        originals = self.original_files()
        job = {}
        with patch.object(app, "q", side_effect=materialize_qiime_outputs), patch.object(app, "command"):
            app.preprocess(job, self.config(output), quality_only=True)
        run = Path(job["run"])
        self.assertEqual(run.parent, output)
        self.assertTrue((run / "quality_complete.json").is_file())
        with app.artifact_path(run, "manifest.tsv").open(encoding="utf8") as stream:
            rows = list(csv.reader(stream, delimiter="\t"))
        self.assertEqual(rows[1], ["CAP1", str(forward), str(reverse)])
        with app.artifact_path(run, "metadata.tsv").open(encoding="utf8") as stream:
            groups = list(csv.reader(stream, delimiter="\t"))
        self.assertEqual(groups[-1], ["CAP1", "处理"])
        self.assert_originals_preserved(originals)

    def test_same_ancestor_and_symlink_to_raw_outputs_are_rejected_without_changes(self):
        fastq(self.raw / "A1_R1.fastq", mate=1)
        fastq(self.raw / "A1_R2.fastq", mate=2)
        self.sheet(["SampleID", "Group"], [["A1", "CON"]])
        alias = self.root / "raw-alias"
        alias.symlink_to(self.raw, target_is_directory=True)
        originals = self.original_files()
        for output in [self.raw, self.root, alias]:
            with self.subTest(output=str(output)), \
                 patch.object(app, "q") as qiime, patch.object(app, "command") as command:
                with self.assertRaises(ValueError):
                    app.preprocess({}, self.config(output), quality_only=True)
                qiime.assert_not_called()
                command.assert_not_called()
                self.assert_originals_preserved(originals)
                self.assertEqual(set(self.original_files()), set(originals))


if __name__ == "__main__":
    unittest.main()
