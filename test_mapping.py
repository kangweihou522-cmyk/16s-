"""Checks for explicit SampleID/Group/FASTQ sample sheets.

Real small FASTQ files verify discovery and manifest routing; mocked external
commands do not establish that QIIME 2 or R is installed or operational.
"""
import csv
import gzip
import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location("mapping_workstation", Path(__file__).with_name("app.py"))
app = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(app)


def fastq(filename, mate=1, ids=("read1", "read2")):
    filename.parent.mkdir(parents=True, exist_ok=True)
    sequence = "ACGT" * 30
    text = "".join(f"@{rid}/{mate}\n{sequence}\n+\n{'I' * len(sequence)}\n" for rid in ids)
    if filename.name.lower().endswith(".gz"):
        with gzip.open(filename, "wt", encoding="ascii") as stream:
            stream.write(text)
    else:
        filename.write_text(text, encoding="ascii")
    return filename


class ExplicitMappingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.raw = self.root / "原始数据 有空格"
        self.raw.mkdir()
        self.meta = self.root / "分组.csv"
        self.forward = fastq(self.raw / "测序批次" / "unrelated-left.fq.gz", mate=1)
        self.reverse = fastq(self.raw / "测序批次" / "anything-right.fastq", mate=2)

    def sheet(self, rows, header=("SampleID", "Group", "forward", "reverse")):
        with self.meta.open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(header)
            writer.writerows(rows)

    def default_row(self, sid="0001", group="对照"):
        return (sid, group, "测序批次/unrelated-left.fq.gz", "测序批次/anything-right.fastq")

    def config(self, mode="paired"):
        return {"input": str(self.raw), "metadata": str(self.meta), "mode": mode}

    def processing_config(self):
        classifier = self.root / "classifier.qza"
        classifier.write_bytes(b"placeholder")
        return self.config() | dict(output=str(self.root / "results"), threads="1",
            min_frequency="1", max_depth="100", trunc_f="250", trunc_r="230",
            trim_f="0", trim_r="0", amplicon="460", primer_mode="removed",
            classifier=str(classifier), font="DejaVu Sans", env="fake-qiime",
            r_env="fake-r", conda="/fake/conda")

    def test_sheet_preserves_ids_and_columns_and_legacy_metadata_contract(self):
        self.sheet([self.default_row()])
        rows = app.read_sample_sheet(self.meta)
        self.assertEqual(rows[0]["SampleID"], "0001")
        self.assertEqual(rows[0]["Group"], "对照")
        self.assertEqual(rows[0]["forward"], self.default_row()[2])
        self.assertEqual(rows[0]["reverse"], self.default_row()[3])
        self.assertEqual(app.metadata(self.meta), [("0001", "对照")])

    def test_arbitrary_names_subdirectories_and_unlisted_files_are_supported(self):
        self.sheet([self.default_row()])
        # This invalid, unlisted FASTQ would fail both automatic naming and
        # sequence parsing if explicit selection were ignored.
        (self.raw / "unlisted-no-direction.fastq").write_text("not FASTQ")
        originals = {p: p.read_bytes() for p in self.raw.rglob("*") if p.is_file()}
        samples = app.discover(self.config())
        self.assertEqual(len(samples), 1)
        self.assertEqual(samples[0]["SampleID"], "0001")
        self.assertEqual(samples[0]["Group"], "对照")
        self.assertEqual(samples[0]["files"], {"1": str(self.forward), "2": str(self.reverse)})
        self.assertEqual(samples[0]["evidence"]["R1"]["sampled_reads"], 2)
        self.assertEqual({p: p.read_bytes() for p in originals}, originals)

    def test_absolute_files_inside_selected_directory_are_supported(self):
        self.sheet([("A1", "CON", str(self.forward), str(self.reverse))])
        self.assertEqual(app.discover(self.config())[0]["files"]["1"], str(self.forward))

    def test_windows_relative_separators_resolve_inside_directory(self):
        self.sheet([("A1", "CON", r"测序批次\unrelated-left.fq.gz", r"测序批次\anything-right.fastq")])
        self.assertEqual(app.discover(self.config())[0]["files"], {"1": str(self.forward), "2": str(self.reverse)})

    def test_single_sheet_uses_forward_only(self):
        self.sheet([("A1", "CON", self.default_row()[2], "")])
        sample = app.discover(self.config("single"))[0]
        self.assertEqual(sample["files"], {"1": str(self.forward)})
        self.assertEqual(set(sample["evidence"]), {"R1"})

    def test_invalid_direction_columns_do_not_silently_drop_input(self):
        cases = [
            ("paired", ("A1", "CON", self.default_row()[2], "")),
            ("paired", ("A1", "CON", "", self.default_row()[3])),
            ("paired", ("A1", "CON", "", "")),
            ("single", self.default_row()),
            ("single", ("A1", "CON", "", "")),
        ]
        for mode, row in cases:
            with self.subTest(mode=mode, row=row):
                self.sheet([row])
                with self.assertRaises(ValueError):
                    app.discover(self.config(mode))

    def test_invalid_sheet_shapes_and_duplicate_ids_are_rejected(self):
        cases = [
            (("SampleID", "Group", "forward"), [("A1", "CON", "x.fastq")]),
            (("SampleID", "Group", "reverse"), [("A1", "CON", "x.fastq")]),
            (("SampleID", "Group", "forward", "reverse", "other"), [self.default_row() + ("x",)]),
            (("SampleID", "Group", "forward", "reverse"), [("A1", "CON", "x.fastq")]),
            (("SampleID", "Group", "forward", "reverse"), [self.default_row() + ("extra",)]),
            (("SampleID", "Group", "forward", "reverse"), [self.default_row(), self.default_row(group="other")]),
            (("SampleID", "Group", "forward", "reverse"), [self.default_row(sid=" A1")]),
            (("SampleID", "Group", "forward", "reverse"), [self.default_row(group="")]),
            (("SampleID", "Group", "forward", "reverse"), []),
        ]
        for header, rows in cases:
            with self.subTest(header=header, rows=rows):
                self.sheet(rows, header)
                with self.assertRaises(ValueError):
                    app.read_sample_sheet(self.meta)

    def test_missing_or_non_fastq_file_is_rejected(self):
        wrong = self.raw / "not-sequencing.csv"
        wrong.write_text("placeholder")
        for name in ("missing.fastq", "not-sequencing.csv", "测序批次"):
            with self.subTest(name=name):
                self.sheet([("A1", "CON", name, self.default_row()[3])])
                with self.assertRaises(ValueError):
                    app.discover(self.config())

    def test_relative_and_absolute_files_outside_directory_are_rejected(self):
        outside = fastq(self.root / "outside.fastq")
        for name in ("../outside.fastq", str(outside)):
            with self.subTest(name=name):
                self.sheet([("A1", "CON", name, self.default_row()[3])])
                with self.assertRaises(ValueError):
                    app.discover(self.config())

    def test_symlink_pointing_outside_directory_is_rejected(self):
        outside = fastq(self.root / "outside.fastq")
        link = self.raw / "linked.fastq"
        link.symlink_to(outside)
        self.sheet([("A1", "CON", link.name, self.default_row()[3])])
        with self.assertRaises(ValueError):
            app.discover(self.config())

    def test_same_file_for_two_directions_is_rejected(self):
        self.sheet([("A1", "CON", self.default_row()[2], self.default_row()[2])])
        with self.assertRaises(ValueError):
            app.discover(self.config())

    def test_same_file_for_different_samples_is_rejected(self):
        self.sheet([self.default_row("A1"), self.default_row("B1")])
        with self.assertRaises(ValueError):
            app.discover(self.config())

    def test_symlink_alias_to_same_file_is_rejected(self):
        alias = self.raw / "alias.fastq.gz"
        alias.symlink_to(self.forward)
        self.sheet([("A1", "CON", self.default_row()[2], alias.name)])
        with self.assertRaises(ValueError):
            app.discover(self.config())

    def test_hardlink_alias_to_same_file_is_rejected(self):
        alias = self.raw / "alias.fastq.gz"
        os.link(self.forward, alias)
        self.sheet([("A1", "CON", self.default_row()[2], alias.name)])
        with self.assertRaises(ValueError):
            app.discover(self.config())

    def test_explicit_mapping_still_validates_paired_read_ids(self):
        fastq(self.reverse, mate=2, ids=("read1", "incorrect"))
        self.sheet([self.default_row()])
        with self.assertRaisesRegex(ValueError, "read ID"):
            app.discover(self.config())

    def test_explicit_sheet_routes_manifest_groups_and_saved_parameters(self):
        self.sheet([self.default_row()])
        config, job, calls = self.processing_config(), {}, []
        originals = {p: p.read_bytes() for p in (self.forward, self.reverse, self.meta)}
        with patch.object(app, "q", side_effect=lambda job, c, *args: calls.append(args)), patch.object(app, "command"):
            app.preprocess(job, config, quality_only=True)
        run = Path(job["run"])
        with (run / "manifest.tsv").open(encoding="utf8") as stream:
            manifest = list(csv.reader(stream, delimiter="\t"))
        self.assertEqual(manifest, [["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"],
                                    ["0001", str(self.forward), str(self.reverse)]])
        with (run / "metadata.tsv").open(encoding="utf8") as stream:
            groups = list(csv.reader(stream, delimiter="\t"))
        self.assertEqual(groups, [["sample-id", "Group"], ["#q2:types", "categorical"], ["0001", "对照"]])
        saved = json.loads((run / "parameters.json").read_text(encoding="utf8"))
        self.assertEqual(saved["samples"][0]["files"], {"1": str(self.forward), "2": str(self.reverse)})
        imported = next(c for c in calls if c[:2] == ("tools", "import"))
        self.assertIn("PairedEndFastqManifestPhred33V2", imported)
        self.assertEqual({p: p.read_bytes() for p in originals}, originals)

    def test_single_explicit_sheet_produces_single_manifest(self):
        self.sheet([("A1", "CON", self.default_row()[2], "")])
        job, calls = {}, []
        config = self.processing_config() | {"mode": "single"}
        with patch.object(app, "q", side_effect=lambda job, c, *args: calls.append(args)), patch.object(app, "command"):
            app.preprocess(job, config, quality_only=True)
        with (Path(job["run"]) / "manifest.tsv").open(encoding="utf8") as stream:
            self.assertEqual(list(csv.reader(stream, delimiter="\t")),
                             [["sample-id", "absolute-filepath"], ["A1", str(self.forward)]])
        imported = next(c for c in calls if c[:2] == ("tools", "import"))
        self.assertIn("SingleEndFastqManifestPhred33V2", imported)

    def test_quality_checkpoint_rejects_mapping_swap_with_unchanged_inventory(self):
        other_f = fastq(self.raw / "second-forward.fastq", mate=1)
        other_r = fastq(self.raw / "second-reverse.fastq", mate=2)
        self.sheet([self.default_row("A1"), ("B1", "处理", other_f.name, other_r.name)])
        config, job = self.processing_config(), {}
        with patch.object(app, "q"), patch.object(app, "command"):
            app.preprocess(job, config, quality_only=True)
        run = Path(job["run"])
        inventory_before = (run / "input_inventory.json").read_bytes()
        config["quality_run"] = job["run"]
        self.sheet([("A1", "对照", other_f.name, other_r.name),
                    ("B1", "处理", self.default_row()[2], self.default_row()[3])])
        with patch.object(app, "q") as qiime, patch.object(app, "command") as command:
            with self.assertRaises(ValueError):
                app.preprocess({}, config)
            qiime.assert_not_called()
            command.assert_not_called()
        self.assertEqual((run / "input_inventory.json").read_bytes(), inventory_before)
        self.assertFalse((run / "preprocess_started.json").exists())


if __name__ == "__main__":
    unittest.main()
