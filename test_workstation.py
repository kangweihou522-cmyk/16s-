"""Meaningful backend checks; no real QIIME 2 execution is simulated as validation."""
import csv
import gzip
import http.client
import importlib.util
import json
import os
import tempfile
import threading
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.parse import urlencode

SPEC = importlib.util.spec_from_file_location("workstation", Path(__file__).with_name("app.py"))
app = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(app)


def fastq(p, ids=("read01", "read02"), sequence="ACGT" * 30, mate=None):
    suffix = f"/{mate}" if mate else ""
    content = "".join(f"@{i}{suffix}\n{sequence}\n+\n{'I' * len(sequence)}\n" for i in ids)
    if str(p).endswith(".gz"):
        with gzip.open(p, "wt", encoding="ascii") as f:
            f.write(content)
    else:
        p.write_text(content, encoding="ascii")


from test_quality_pipeline import materialize_qiime_outputs


class BackendTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.raw = self.root / "raw"
        self.raw.mkdir()
        self.meta = self.root / "groups.csv"

    def metadata(self, rows):
        with self.meta.open("w", newline="", encoding="utf-8-sig") as f:
            w = csv.writer(f)
            w.writerow(["SampleID", "Group"])
            w.writerows(rows)

    def pair(self, sid, ids=("read01", "read02")):
        for mate in (1, 2):
            fastq(self.raw / f"{sid}_R{mate}.fastq.gz", ids=ids, mate=mate)

    def config(self):
        return {"input": str(self.raw), "metadata": str(self.meta), "mode": "paired"}

    def processing_config(self):
        classifier = self.root / "classifier.qza"
        classifier.write_bytes(b"mock-only")
        return self.config() | dict(output=str(self.root / "results"), threads="1",
            min_frequency="1", max_depth="100", trunc_f="250", trunc_r="230",
            trim_f="0", trim_r="0", amplicon="460", primer_mode="338F",
            classifier=str(classifier), font="DejaVu Sans", env="fake-qiime",
            r_env="custom-r-plotting", conda="/fake/conda")

    def record_q(self, calls, job, config, *args):
        calls.append(args)
        materialize_qiime_outputs(job, config, *args)

    def mock_biom_conversion(self, job, args, cwd=None, process_env=None):
        """Materialize only synthetic TSV counts needed by downstream checks."""
        if "biom" not in args or "convert" not in args:
            return
        run = Path(cwd or job["run"])
        output = run / args[args.index("-o") + 1]
        output.parent.mkdir(parents=True, exist_ok=True)
        config = json.loads((run / "parameters.json").read_text(encoding="utf8"))
        samples = config.get("samples", config.get("config", {}).get("samples", []))
        ids = [x["SampleID"] for x in samples]
        output.write_text("#OTU ID\t" + "\t".join(ids) + "\nASV1\t" +
                          "\t".join(str(125 * (i + 1)) for i in range(len(ids))) + "\n", encoding="utf8")

    def test_paired_ids_and_leading_zero_sample_names(self):
        self.metadata([("0001", "A"), ("0002", "B")])
        self.pair("0001")
        self.pair("0002")
        found = app.discover(self.config())
        self.assertEqual([x["SampleID"] for x in found], ["0001", "0002"])
        self.assertEqual(found[0]["evidence"]["R1"]["sampled_reads"], 2)
        self.assertEqual(Path(found[0]["files"]["2"]).name, "0001_R2.fastq.gz")

    def test_preserve_real_sample_id_ending_in_illumina_style_suffix(self):
        self.metadata([("mouse_S1", "A")])
        self.pair("mouse_S1")
        self.assertEqual(app.discover(self.config())[0]["SampleID"], "mouse_S1")

    def test_illumina_lane_names_match_metadata(self):
        self.metadata([("0001", "A")])
        for mate in (1, 2):
            fastq(self.raw / f"0001_S7_L001_R{mate}_001.fastq.gz", mate=mate)
        self.assertEqual(app.discover(self.config())[0]["SampleID"], "0001")

    def test_raw_dot_filename_prefix_maps_to_metadata_sample(self):
        self.metadata([("A1", "CON")])
        for mode in ("paired", "single"):
            with self.subTest(mode=mode):
                raw = self.root / mode
                raw.mkdir()
                if mode == "paired":
                    for mate in (1, 2):
                        fastq(raw / f"raw.A1.{mate}.fq.gz", mate=mate)
                else:
                    fastq(raw / "raw.A1.fq.gz")
                config = self.config() | {"input": str(raw), "mode": mode}
                sample = app.discover(config)[0]
                self.assertEqual(sample["SampleID"], "A1")
                self.assertTrue(Path(sample["files"]["1"]).name.startswith("raw.A1"))

    def test_raw_prefix_is_rejected_when_metadata_mapping_is_ambiguous(self):
        self.metadata([("A1", "CON"), ("raw.A1", "LDC")])
        for mode in ("paired", "single"):
            with self.subTest(mode=mode):
                raw = self.root / mode
                raw.mkdir()
                if mode == "paired":
                    for mate in (1, 2):
                        fastq(raw / f"raw.A1.{mate}.fq.gz", mate=mate)
                else:
                    fastq(raw / "raw.A1.fq.gz")
                with self.assertRaisesRegex(ValueError, "多个 SampleID"):
                    app.discover(self.config() | {"input": str(raw), "mode": mode})

    def test_sample_raw_suffix_preserves_sample_digit_and_pairs_mates(self):
        self.metadata([("CAP1", "CON")])
        for mate in (1, 2):
            fastq(self.raw / f"CAP1.raw_{mate}.fastq", mate=mate)
        sample = app.discover(self.config())[0]
        self.assertEqual(sample["SampleID"], "CAP1")
        self.assertEqual(Path(sample["files"]["1"]).name, "CAP1.raw_1.fastq")
        self.assertEqual(Path(sample["files"]["2"]).name, "CAP1.raw_2.fastq")

    def test_sample_raw_suffix_rejects_ambiguous_metadata(self):
        self.metadata([("CAP1", "CON"), ("CAP1.raw", "LDC")])
        for mate in (1, 2):
            fastq(self.raw / f"CAP1.raw_{mate}.fastq", mate=mate)
        with self.assertRaisesRegex(ValueError, "多个 SampleID"):
            app.discover(self.config())

    def test_two_forward_files_or_missing_reverse_never_form_fake_pair(self):
        self.metadata([("CAP1", "CON")])
        fastq(self.raw / "CAP1.raw_1.fastq", mate=1)
        with self.assertRaisesRegex(ValueError, "缺少 R1 或 R2"):
            app.discover(self.config())
        # A second real filename containing direction 1 is a duplicate R1,
        # even when it has another compression or FASTQ extension.
        fastq(self.raw / "CAP1.raw_1.fq.gz", mate=1)
        with self.assertRaisesRegex(ValueError, "多个同方向"):
            app.discover(self.config())

    def test_reject_unmatched_metadata(self):
        self.metadata([("0002", "A")])
        self.pair("0001")
        with self.assertRaisesRegex(ValueError, "样本不匹配"):
            app.discover(self.config())

    def test_reject_duplicate_metadata(self):
        self.metadata([("0001", "A"), ("0001", "B")])
        with self.assertRaises(ValueError):
            app.metadata(self.meta)

    def test_reject_extra_csv_cells(self):
        self.metadata([("0001", "A", "unexpected")])
        with self.assertRaises(ValueError):
            app.metadata(self.meta)

    def test_reject_pair_id_mismatch(self):
        self.metadata([("0001", "A")])
        fastq(self.raw / "0001_R1.fastq", ids=("read01", "read02"), mate=1)
        fastq(self.raw / "0001_R2.fastq", ids=("read01", "wrong"), mate=2)
        with self.assertRaisesRegex(ValueError, "read ID"):
            app.discover(self.config())

    def test_missing_mate_is_rejected(self):
        self.metadata([("0001", "A")])
        fastq(self.raw / "0001_R1.fastq", mate=1)
        with self.assertRaisesRegex(ValueError, "缺少"):
            app.discover(self.config())

    def test_multilane_is_rejected_without_silent_partial_use(self):
        self.metadata([("0001", "A")])
        for lane in (1, 2):
            for mate in (1, 2):
                fastq(self.raw / f"0001_S7_L{lane:03}_R{mate}_001.fastq", mate=mate)
        with self.assertRaisesRegex(ValueError, "多个同方向"):
            app.discover(self.config())

    def test_truncated_gzip_is_rejected(self):
        p = self.raw / "sample.fastq.gz"
        fastq(p)
        p.write_bytes(p.read_bytes()[:-8])
        with self.assertRaises((EOFError, gzip.BadGzipFile, OSError, ValueError)):
            list(app.reads(p))

    def test_fastq_quality_length_is_validated(self):
        p = self.raw / "sample.fastq"
        p.write_text("@r1\nACGT\n+\nIII\n", encoding="ascii")
        with self.assertRaisesRegex(ValueError, "FASTQ"):
            list(app.reads(p))

    def test_primer_detection_does_not_distinguish_shared_core(self):
        # 338F contains a 341F-compatible core. Both positive matches are
        # evidence of ambiguity, never proof that a read used 341F.
        seq = app.PRIMERS["338F"] + "A" * 50
        self.assertTrue(app.match_primer(seq, app.PRIMERS["338F"]))
        self.assertTrue(app.match_primer(seq, app.PRIMERS["341F"]))
        self.assertFalse(app.match_primer("T" * 120, app.PRIMERS["341F"]))
        self.assertFalse(app.match_primer("ACG", app.PRIMERS["341F"]))

    def test_overlap_rejects_insufficient_shared_bases(self):
        c = dict(trunc_f=240, trunc_r=220, trim_f=0, trim_r=0, amplicon=450)
        with self.assertRaises(ValueError):
            app.check_overlap(c)

    def test_overlap_keeps_normal_trim_geometry(self):
        # Removing outer bases from both reads does not reduce their middle
        # overlap if neither left trim reaches that overlap.
        c = dict(trunc_f=250, trunc_r=230, trim_f=17, trim_r=20, amplicon=460)
        self.assertEqual(app.check_overlap(c), (250, 230, 17, 20))

    def test_overlap_rejects_left_trim_entering_shared_region(self):
        c = dict(trunc_f=100, trunc_r=100, trim_f=90, trim_r=0, amplicon=180)
        with self.assertRaises(ValueError):
            app.check_overlap(c)

    def test_retention_keeps_boundary_depth_and_reports_missing_samples(self):
        run = self.root / "run"
        (run / "exported-table").mkdir(parents=True)
        (run / "exported-table" / "feature-table.tsv").write_text(
            "# Constructed from biom file\n#OTU ID\t0001\t0002\nASV1\t70\t20\nASV2\t30\t30\n",
            encoding="utf8",
        )
        (run / "parameters.json").write_text(json.dumps({"samples": [
            {"SampleID": "0001", "Group": "A"}, {"SampleID": "0002", "Group": "B"},
            {"SampleID": "0003", "Group": "B"}]}), encoding="utf8")
        result = app.retention(run, 100)
        self.assertEqual([x["retained"] for x in result], [True, False, False])
        self.assertEqual([x["count"] for x in result], [100, 50, 0])

    def test_windows_drive_path_mapping(self):
        self.assertEqual(str(app.path(r"D:\测序数据\raw")), "/mnt/d/测序数据/raw")

    def test_fake_orchestration_paired_manifest_and_saved_inputs(self):
        # This validates argument assembly and output isolation only. Mocked
        # subprocesses provide no evidence that QIIME 2 or R actually works.
        self.metadata([("0001", "A")])
        self.pair("0001")
        original = {p: p.read_bytes() for p in self.raw.iterdir()}
        config = self.processing_config()
        qiime_args, command_args, command_envs = [], [], []
        job = {}
        def fake_command(job, args, cwd=None, process_env=None):
            command_args.append(args)
            command_envs.append(process_env)
            self.mock_biom_conversion(job, args, cwd)
        with patch.object(app, "q", side_effect=lambda job, c, *args: self.record_q(qiime_args, job, c, *args)), \
             patch.object(app, "command", side_effect=fake_command):
            app.preprocess(job, config)
        run = Path(job["run"])
        with app.artifact_path(run, "manifest.tsv").open(encoding="utf8") as f:
            manifest = list(csv.reader(f, delimiter="\t"))
        self.assertEqual(manifest[0], ["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"])
        self.assertEqual(manifest[1][0], "0001")
        import_args = next(a for a in qiime_args if a[:2] == ("tools", "import"))
        self.assertIn("PairedEndFastqManifestPhred33V2", import_args)
        self.assertTrue(any(a[:2] == ("cutadapt", "trim-paired") for a in qiime_args))
        self.assertTrue(any(a[:2] == ("dada2", "denoise-paired") for a in qiime_args))
        r_calls = [(a, e) for a, e in zip(command_args, command_envs) if "Rscript" in a]
        self.assertEqual(len(r_calls), 2)
        for args, process_env in r_calls:
            self.assertEqual(args[args.index("-n") + 1], "custom-r-plotting")
            if args[-1] == "--version":
                self.assertNotIn("--vanilla", args)
            else:
                self.assertIn("--vanilla", args)
            self.assertEqual(process_env["R_LIBS_USER"], os.devnull)
        biom = next(a for a in command_args if "biom" in a)
        self.assertEqual(biom[biom.index("-n") + 1], "fake-qiime")
        self.assertTrue((run / "preprocess_complete.json").exists())
        self.assertTrue((run / "input_inventory.json").exists())
        self.assertEqual({p: p.read_bytes() for p in self.raw.iterdir()}, original)

    def test_blank_curve_max_depth_uses_filtered_sample_totals(self):
        self.metadata([("0001", "CON"), ("0002", "LDC")])
        self.pair("0001")
        self.pair("0002")
        config = self.processing_config() | {"max_depth": ""}
        job, calls = {}, []
        with patch.object(app, "q", side_effect=lambda job, c, *args: self.record_q(calls, job, c, *args)), \
             patch.object(app, "command", side_effect=self.mock_biom_conversion):
            app.preprocess(job, config)
        run = Path(job["run"])
        saved = json.loads((run / "parameters.json").read_text(encoding="utf8"))
        self.assertEqual(saved["effective_curve_max_depth"], 250)
        rarefaction = next(a for a in calls if a[:2] == ("diversity", "alpha-rarefaction"))
        self.assertEqual(rarefaction[rarefaction.index("--p-max-depth") + 1], 250)

    def test_mock_quality_phase_needs_no_truncation_or_classifier(self):
        self.metadata([("0001", "A")])
        self.pair("0001")
        config = self.processing_config()
        for k in ("trunc_f", "trunc_r", "trim_f", "trim_r", "amplicon", "classifier"):
            del config[k]
        job, calls = {}, []
        with patch.object(app, "q", side_effect=lambda job, c, *args: self.record_q(calls, job, c, *args)), \
             patch.object(app, "command"):
            app.preprocess(job, config, quality_only=True)
        run = Path(job["run"])
        self.assertTrue((run / "quality_complete.json").exists())
        self.assertFalse((run / "preprocess_complete.json").exists())
        self.assertTrue(any(a[:2] == ("tools", "import") for a in calls))
        self.assertTrue(any(a[:2] == ("cutadapt", "trim-paired") for a in calls))
        self.assertFalse(any(a[0] in ("dada2", "feature-classifier") for a in calls))

    def test_mock_quality_reuse_skips_import_and_cutadapt(self):
        self.metadata([("0001", "A")])
        self.pair("0001")
        config, job, calls = self.processing_config(), {}, []
        with patch.object(app, "q", side_effect=lambda job, c, *args: self.record_q(calls, job, c, *args)), \
             patch.object(app, "command", side_effect=self.mock_biom_conversion):
            app.preprocess(job, config, quality_only=True)
            config["quality_run"] = job["run"]
            calls.clear()
            app.preprocess(job, config)
        self.assertFalse(any(a[:2] == ("tools", "import") for a in calls))
        self.assertFalse(any(a[0] == "cutadapt" for a in calls))
        dada2 = next(a for a in calls if a[:2] == ("dada2", "denoise-paired"))
        self.assertEqual(dada2[dada2.index("--i-demultiplexed-seqs") + 1], "trimmed.qza")
        self.assertTrue((Path(job["run"]) / "preprocess_complete.json").exists())

    def test_quality_reuse_rejects_changed_input_before_running_commands(self):
        self.metadata([("0001", "A")])
        self.pair("0001")
        config, job = self.processing_config(), {}
        with patch.object(app, "q", side_effect=lambda job,c,*args: self.record_q([],job,c,*args)), patch.object(app, "command"):
            app.preprocess(job, config, quality_only=True)
        config["quality_run"] = job["run"]
        p = self.raw / "0001_R1.fastq.gz"
        st = p.stat()
        os.utime(p, ns=(st.st_atime_ns, st.st_mtime_ns + 10_000_000))
        with patch.object(app, "q") as qiime, patch.object(app, "command") as command:
            with self.assertRaisesRegex(ValueError, "原始文件已改变"):
                app.preprocess(job, config)
            qiime.assert_not_called()
            command.assert_not_called()

    def test_cancel_prevents_subprocess_start(self):
        with patch.object(app.subprocess, "Popen") as popen:
            with self.assertRaisesRegex(RuntimeError, "已停止"):
                app.command({"cancel": True}, ["echo", "unused"])
            popen.assert_not_called()

    def test_constant_alpha_is_skipped_but_variable_data_still_tested(self):
        filename = self.root / "alpha-diversity.tsv"
        rows = [("A1", "CON"), ("A2", "CON"), ("B1", "LDC"), ("B2", "LDC")]
        filename.write_text("sample-id\tvalue\nA1\t1\nA2\t1\nB1\t1\nB2\t1\n", encoding="utf8")
        self.assertIn("所有样本", app.alpha_skip_reason(filename, rows))
        filename.write_text("sample-id\tvalue\nA1\t1\nA2\t2\nB1\t3\nB2\t4\n", encoding="utf8")
        self.assertEqual(app.alpha_skip_reason(filename, rows), "")
        filename.write_text("sample-id\tvalue\nA1\t1\nA2\t1\nB1\t1\nB2\t1\nC1\t2\nC2\t3\n", encoding="utf8")
        self.assertIn("两两", app.alpha_skip_reason(filename, rows + [("C1", "HDC"), ("C2", "HDC")]))

    def test_zero_distances_are_skipped_but_nonzero_distances_not_skipped(self):
        filename = self.root / "distance-matrix.tsv"
        filename.write_text("\tA1\tA2\nA1\t0\t0\nA2\t0\t0\n", encoding="utf8")
        self.assertIn("所有距离为0", app.distance_skip_reason(filename))
        filename.write_text("\tA1\tA2\nA1\t0\t1\nA2\t1\t0\n", encoding="utf8")
        self.assertEqual(app.distance_skip_reason(filename), "")

    def test_mock_diversity_uses_one_rarefaction_and_saved_qiime_statistics(self):
        # Synthetic counts and placeholder artifacts test the orchestration.
        # Their presence is not a substitute for a real QIIME 2/R run.
        parent = self.root / "preprocessed"
        (parent / "exported-table").mkdir(parents=True)
        (parent / "exported-table" / "feature-table.tsv").write_text(
            "#OTU ID\t0001\t0002\t0003\t0004\nASV1\t100\t100\t100\t100\n",
            encoding="utf8",
        )
        samples = [{"SampleID": f"{i:04}", "Group": "A" if i < 3 else "B"}
                   for i in range(1, 5)]
        saved = self.processing_config() | {"samples": samples}
        (parent / "parameters.json").write_text(json.dumps(saved), encoding="utf8")
        (parent / "preprocess_complete.json").write_text("{}", encoding="utf8")
        for name in ("rooted-tree.qza", "filtered-table.qza"):
            (parent / name).write_bytes(b"placeholder-for-command-tests")
        qiime_calls, commands, command_envs, events, job = [], [], [], [], {}

        def fake_q(job, config, *args):
            qiime_calls.append(tuple(args))
            events.append(("qiime", tuple(args)))
            run = Path(job["run"])
            if args[:2] == ("diversity", "core-metrics-phylogenetic"):
                core = app.artifact_path(run, "core-metrics-results")
                core.mkdir()
                names = ["rarefied_table.qza"]
                names += [f"{m}_vector.qza" for m in
                          ("observed_features", "shannon", "evenness", "faith_pd")]
                for metric in ("weighted_unifrac", "unweighted_unifrac", "bray_curtis", "jaccard"):
                    names += [f"{metric}_distance_matrix.qza", f"{metric}_pcoa_results.qza",
                              f"{metric}_emperor.qzv"]
                for name in names:
                    (core / name).write_bytes(b"placeholder-for-command-tests")
            for idx, value in enumerate(args[:-1]):
                if str(value).startswith("--o-"):
                    output = app.artifact_path(run, str(args[idx + 1]))
                    output.parent.mkdir(parents=True, exist_ok=True)
                    output.write_bytes(b"placeholder-for-command-tests")
            if args[:2] == ("tools", "export"):
                destination = app.artifact_path(run, str(args[args.index("--output-path") + 1]))
                destination.mkdir(parents=True, exist_ok=True)
                source = str(args[args.index("--input-path") + 1])
                if "_vector.qza" in source:
                    (destination / "alpha-diversity.tsv").write_text(
                        "sample-id\tvalue\n0001\t1\n0002\t2\n0003\t3\n0004\t4\n", encoding="utf8")
                elif "_distance_matrix.qza" in source:
                    (destination / "distance-matrix.tsv").write_text(
                        "\t0001\t0002\t0003\t0004\n0001\t0\t1\t2\t3\n"
                        "0002\t1\t0\t3\t2\n0003\t2\t3\t0\t1\n0004\t3\t2\t1\t0\n", encoding="utf8")

        def fake_command(job, args, cwd=None, process_env=None):
            commands.append(args)
            command_envs.append(process_env)
            events.append(("command", tuple(args)))
            self.mock_biom_conversion(job,args,cwd,process_env)
            if "Rscript" in args:
                output=app.artifact_path(Path(job["run"]),"R_original_style")
                output.mkdir(parents=True,exist_ok=True)
                (output/"render_summary.json").write_text(json.dumps({"source":"saved_QIIME2_outputs","notes":[]}))

        with patch.object(app, "q", side_effect=fake_q), \
             patch.object(app, "command", side_effect=fake_command):
            app.diversity(job, {"run": str(parent), "depth": "100", "permutations": "999"})
        self.assertEqual(sum(a[:2] == ("diversity", "core-metrics-phylogenetic") for a in qiime_calls), 1)
        self.assertFalse(any(a[:2] == ("feature-table", "rarefy") for a in qiime_calls))
        alpha = [a for a in qiime_calls if a[:2] == ("diversity", "alpha-group-significance")]
        beta = [a for a in qiime_calls if a[:2] == ("diversity", "beta-group-significance")]
        self.assertEqual(len(alpha), 7)
        self.assertEqual(len(beta), 8)
        self.assertEqual({a[a.index("--p-method") + 1] for a in beta}, {"permanova", "permdisp"})
        self.assertTrue(all("--p-pairwise" in a for a in beta))
        self.assertTrue(any("Rscript" in a for a in commands))
        render_args, render_env = next((a, e) for a, e in zip(commands, command_envs) if "Rscript" in a)
        self.assertEqual(render_args[render_args.index("-n") + 1], "custom-r-plotting")
        self.assertEqual(render_env["R_LIBS_USER"], os.devnull)
        self.assertIn("--vanilla", render_args)
        biom = next(a for a in commands if "biom" in a)
        self.assertEqual(biom[biom.index("-n") + 1], "fake-qiime")
        self.assertTrue((Path(job["run"]) / "complete.json").exists())
        self.assertTrue((Path(job["run"]) / "qiime_outputs_complete.json").exists())
        render_index = next(i for i, (kind, args) in enumerate(events)
                            if kind == "command" and "Rscript" in args)
        export_events = [i for i, (kind, args) in enumerate(events)
                         if kind == "qiime" and args[:2] == ("tools", "export")]
        self.assertTrue(export_events)
        self.assertLess(max(export_events), render_index)
        exported = {str(args[args.index("--input-path") + 1]) for args in qiime_calls
                    if args[:2] == ("tools", "export")}
        self.assertIn("core-metrics-results/rarefied_table.qza", exported)
        for metric in ("observed_features", "shannon", "evenness", "faith_pd", "chao1", "simpson", "goods_coverage"):
            self.assertIn(f"core-metrics-results/{metric}_vector.qza", exported)
        for metric in ("weighted_unifrac", "unweighted_unifrac", "bray_curtis", "jaccard"):
            self.assertIn(f"core-metrics-results/{metric}_pcoa_results.qza", exported)
            self.assertIn(f"core-metrics-results/{metric}_distance_matrix.qza", exported)

        # A plotting failure must leave the underlying QIIME products usable.
        # Run again in a new result directory with an intentionally failed R.
        first_sub = Path(job["run"])
        first_parameters = (first_sub / "parameters.json").read_bytes()

        def failed_render(job, args, cwd=None, process_env=None):
            if "Rscript" in args:
                raise RuntimeError("intentional mocked R rendering failure")
            fake_command(job, args, cwd, process_env)

        with patch.object(app, "q", side_effect=fake_q), \
             patch.object(app, "command", side_effect=failed_render):
            with self.assertRaisesRegex(RuntimeError, "mocked R rendering failure"):
                app.diversity(job, {"run": str(parent), "depth": "100", "permutations": "999"})
        failed_sub = Path(job["run"])
        self.assertNotEqual(first_sub, failed_sub)
        self.assertTrue((failed_sub / "qiime_outputs_complete.json").exists())
        self.assertTrue((failed_sub / "qiime_statistics_complete.json").exists())
        self.assertTrue(app.artifact_path(failed_sub, "core-metrics-results/rarefied_table.qza").exists())
        self.assertFalse((failed_sub / "complete.json").exists())
        self.assertEqual((first_sub / "parameters.json").read_bytes(), first_parameters)

    def test_preflight_discovers_environments_without_requiring_one_selected(self):
        # This is a diagnostic mock, not evidence of installed tools.
        job = {"log": str(self.root / "preflight.log")}

        def fake_run(args, **kwargs):
            if args[:3] == ["/mock/conda", "env", "list"]:
                return app.subprocess.CompletedProcess(args, 0, '{"envs":["/mock/envs/qiime-installed"]}', "")
            if "Rscript" in args:
                return app.subprocess.CompletedProcess(args, 1, "", "EnvironmentNameNotFound: 16s-workstation-r")
            self.fail(f"Unexpected process with no environment selected: {args}")

        with patch.object(app, "find_conda", return_value="/mock/conda"), \
             patch.object(app.subprocess, "run", side_effect=fake_run):
            app.preflight(job, {"env": "", "font": "DejaVu Sans"})
        result = job["result"]
        self.assertEqual(result["environments"], ["qiime-installed"])
        checks = {x["component"]: x for x in result["checks"]}
        self.assertTrue(checks["conda"]["ok"])
        self.assertFalse(checks["QIIME 2"]["ok"])
        self.assertFalse(checks["专用 R 环境与绘图依赖"]["ok"])
        self.assertFalse(checks["用户选择的分类器"]["ok"])

    def test_preflight_stops_before_start_and_between_checks(self):
        job = {"log": str(self.root / "preflight.log"), "cancel": True}
        with patch.object(app.subprocess, "run") as run:
            with self.assertRaisesRegex(RuntimeError, "已停止环境检查"):
                app.preflight(job, {"env": "fake-qiime"})
            run.assert_not_called()
        job["cancel"] = False

        def cancel_after_first(args, **kwargs):
            job["cancel"] = True
            return app.subprocess.CompletedProcess(args, 0, '{"envs":[]}', "")

        with patch.object(app, "find_conda", return_value="/mock/conda"), \
             patch.object(app.subprocess, "run", side_effect=cancel_after_first) as run:
            with self.assertRaisesRegex(RuntimeError, "已停止环境检查"):
                app.preflight(job, {"env": "fake-qiime"})
            self.assertEqual(run.call_count, 1)


class HTTPBoundaryTests(unittest.TestCase):
    """Real local HTTP requests with fake files, never actual analysis."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.server = app.ThreadingHTTPServer(("127.0.0.1", 0), app.Handler)
        thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        thread.start()

        def close():
            self.server.shutdown()
            self.server.server_close()
            thread.join(timeout=2)

        self.addCleanup(close)

    def get(self, route, authorized=True):
        conn = http.client.HTTPConnection("127.0.0.1", self.server.server_port, timeout=5)
        try:
            headers = {"X-Workstation-Token": app.TOKEN} if authorized else {}
            conn.request("GET", route, headers=headers)
            reply = conn.getresponse()
            return reply.status, reply.read(), dict(reply.getheaders())
        finally:
            conn.close()

    def test_browse_skips_windows_protected_file(self):
        protected = self.root / "DumpStack.log.tmp"
        protected.touch()
        (self.root / "reads").mkdir()
        (self.root / "groups.csv").write_text("SampleID,Group\n", encoding="utf8")
        real_is_dir = Path.is_dir

        def guarded_is_dir(candidate):
            if candidate == protected:
                raise PermissionError(13, "Permission denied", str(protected))
            return real_is_dir(candidate)

        with patch.object(Path, "is_dir", guarded_is_dir):
            status, body, _ = self.get("/api/browse?" + urlencode({"path": str(self.root)}))
        self.assertEqual(status, 200, body)
        entries = json.loads(body)["entries"]
        self.assertEqual({x["name"] for x in entries}, {"reads", "groups.csv"})
        self.assertTrue(next(x for x in entries if x["name"] == "reads")["directory"])

    def test_browse_requires_local_session_token(self):
        status, _, _ = self.get("/api/browse?" + urlencode({"path": str(self.root)}), authorized=False)
        self.assertEqual(status, 403)

    def test_download_and_listing_are_limited_to_marked_run_outputs(self):
        run = self.root / "run"
        run.mkdir()
        app.marker(run)
        output = run / "图像.svg"
        output.write_bytes(b"synthetic-test-output")
        secret = self.root / "private.txt"
        secret.write_bytes(b"must-not-be-served")
        (run / "escaped.txt").symlink_to(secret)
        status, body, _ = self.get("/api/results?" + urlencode({"run": str(run)}))
        self.assertEqual(status, 200, body)
        names = {x["name"] for x in json.loads(body)["files"]}
        self.assertIn("图像.svg", names)
        self.assertNotIn("escaped.txt", names)
        status, body, headers = self.get("/api/file?" + urlencode({"path": str(output)}))
        self.assertEqual(status, 200, body)
        self.assertEqual(body, output.read_bytes())
        self.assertEqual(headers["X-Content-Type-Options"], "nosniff")
        self.assertTrue(headers["Content-Disposition"].startswith("attachment;"))
        for filename in (secret, run / "escaped.txt", run / ".." / "private.txt"):
            with self.subTest(filename=filename):
                status, body, _ = self.get("/api/file?" + urlencode({"path": str(filename)}))
                self.assertEqual(status, 400, body)
                self.assertNotIn(secret.read_bytes(), body)
        status, _, _ = self.get("/api/results?" + urlencode({"run": str(self.root)}))
        self.assertEqual(status, 400)
        status, _, _ = self.get("/api/file?" + urlencode({"path": str(output)}), authorized=False)
        self.assertEqual(status, 403)


if __name__ == "__main__":
    unittest.main()
