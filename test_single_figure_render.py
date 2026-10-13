"""Opt-in real R checks for cached per-figure editing; no model fitting.

Run with WORKSTATION_RSCRIPT=/path/to/Rscript python -m unittest test_single_figure_render.
Only synthetic saved values and official-shaped result objects are used.
"""
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parent
RSCRIPT = os.environ.get("WORKSTATION_RSCRIPT")


@unittest.skipUnless(RSCRIPT, "Set WORKSTATION_RSCRIPT to enable focused real R rendering")
class SingleFigureRendering(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="16s-single-figure-")
        cls.base = Path(cls.temporary.name)
        cls.source = cls.base / "source"
        cls.analysis = cls.source / "analysis"
        cls.preprocess = cls.source / "preprocess"
        cls.ids = [f"S{i}" for i in range(1, 7)]
        cls.groups = ["CON", "CON", "LDC", "LDC", "HDC", "HDC"]
        metadata = "sample-id\tGroup\n" + "".join(
            f"{sample}\t{group}\n" for sample, group in zip(cls.ids, cls.groups)
        )
        for directory in (cls.analysis, cls.preprocess):
            cls.write(directory / "metadata.tsv", metadata)
        cls.write(cls.analysis / "exports/alpha/shannon/alpha-diversity.tsv",
                  "sample-id\tshannon\n" + "".join(f"{sample}\t{index / 10}\n" for index, sample in enumerate(cls.ids, 1)))
        coords = [(-2, -1), (-1, 1), (0, -.5), (.5, .5), (2, -1), (3, 1)]
        ordination = ("Eigvals\t2\n0.6\t0.4\n\nProportion explained\t2\n0.6\t0.4\n\n"
                      "Species\t0\t0\n\nSite\t6\t2\n" +
                      "".join(f"{sample}\t{x}\t{y}\n" for sample, (x, y) in zip(cls.ids, coords)))
        distance = "\t" + "\t".join(cls.ids) + "\n" + "".join(
            sample + "\t" + "\t".join(str(math.dist(xy, other)) for other in coords) + "\n"
            for sample, xy in zip(cls.ids, coords)
        )
        for metric in ("weighted_unifrac", "unweighted_unifrac", "bray_curtis", "jaccard"):
            cls.write(cls.analysis / f"exports/pcoa/{metric}/ordination.txt", ordination)
            cls.write(cls.analysis / f"exports/distance/{metric}/distance-matrix.tsv", distance)
        counts = "\t".join(cls.ids) + "\ntaxon_A\t10\t20\t30\t40\t50\t60\ntaxon_B\t60\t50\t40\t30\t20\t10\n"
        cls.write(cls.preprocess / "exported-table/feature-table.tsv", "#OTU ID\t" + counts)
        cls.write(cls.preprocess / "exported-taxonomy/taxonomy.tsv",
                  "Feature ID\tTaxon\ntaxon_A\tk__Bacteria; p__P1; c__C1; o__O1; f__F1; g__G1\n"
                  "taxon_B\tk__Bacteria; p__P2; c__C2; o__O2; f__F2; g__G2\n")
        cls.diff = cls.source / "differential"
        cls.write(cls.diff / "metadata.tsv", metadata)
        cls.write(cls.diff / "counts.tsv", "FeatureID\t" + counts)
        tables = [{"level": 2, "name": "phylum", "counts_tsv": str(cls.diff / "unavailable_counts.tsv")},
                  {"level": 6, "name": "genus", "counts_tsv": str(cls.diff / "counts.tsv")}]
        options = {"method": "ancombc2", "reference_group": "CON", "groups": ["CON", "LDC", "HDC"],
                   "metadata_tsv": str(cls.diff / "metadata.tsv"), "alpha": .05, "pseudo_sens": True,
                   "font": "sans", "tables": tables}
        cls.write(cls.diff / "differential_options.json", json.dumps(options))
        fit_dir = cls.diff / "02_differential/L6_genus"
        fit_dir.mkdir(parents=True)
        # This is a cached fixture, not an ANCOM-BC2 fit or an inference test.
        script = f"""
source({json.dumps(str(ROOT / 'differential_ancombc.R'))})
counts <- read_count_input({json.dumps(str(cls.diff / 'counts.tsv'))})
meta <- read_group_metadata({json.dumps(str(cls.diff / 'metadata.tsv'))},colnames(counts),'CON',c('CON','LDC','HDC'))
plan <- official_contrast_map(meta,'CON')
res <- data.frame(taxon=rownames(counts),stringsAsFactors=FALSE)
for (key in plan$key) {{
  res[[paste0('lfc_',key)]] <- c(1,-2)
  res[[paste0('se_',key)]] <- c(.2,.3)
  res[[paste0('p_',key)]] <- c(.005,.008)
  res[[paste0('q_',key)]] <- c(.01,.02)
  res[[paste0('diff_',key)]] <- TRUE
  res[[paste0('passed_ss_',key)]] <- TRUE
}}
saveRDS(list(res=res,res_pair=res),{json.dumps(str(fit_dir / 'ancombc2_fit.rds'))})
"""
        cls.execute(["-e", script])
        cls.qiime = cls.source / "qiime"
        cls.qiime.mkdir()
        options["method"] = "qiime_ancom"
        cls.write(cls.qiime / "differential_options.json", json.dumps(options))
        cls.write(cls.qiime / "02_differential/L6_genus/ancom_results.tsv",
                  "taxon\tW\tReject null hypothesis\ntaxon_A\t20\tTrue\ntaxon_B\t4\tFalse\n")
        cls.source_hashes = cls.hashes(cls.source)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    @staticmethod
    def write(path, text):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")

    @staticmethod
    def hashes(path):
        return {str(item.relative_to(path)): hashlib.sha256(item.read_bytes()).hexdigest()
                for item in path.rglob("*") if item.is_file()}

    @classmethod
    def execute(cls, arguments):
        result = subprocess.run([RSCRIPT, "--vanilla", *arguments], cwd=ROOT,
                                capture_output=True, text=True, timeout=60)
        if result.returncode:
            raise AssertionError(result.stdout + result.stderr)

    def settings(self):
        return {"general": {"dpi": 72, "formats": "png_svg", "show_caption": False},
                "alpha": {"width": 4, "height": 3, "show_letters": False},
                "beta": {"width": 4, "height": 3, "show_ellipse": False},
                "taxonomy": {"width": 4, "height": 3},
                "heatmap": {"width": 4, "height": 3},
                "bubble": {"width": 4, "height": 3},
                "effect": {"width": 4, "height": 3}}

    def assert_single(self, output, profile, stem):
        manifest = json.loads((output / "figure_manifest.json").read_text())
        self.assertEqual([(row["profile"], row["stem"]) for row in manifest["figures"]], [(profile, stem)])
        images = sorted(str(path.relative_to(output)) for path in output.rglob("*")
                        if path.suffix in (".png", ".svg", ".pdf"))
        self.assertEqual(images, sorted([stem + ".png", stem + ".svg"]))
        self.assertEqual(self.hashes(self.source), self.source_hashes)

    def render_diversity(self, profile, stem):
        output = self.base / self.id().split(".")[-1]
        settings = self.base / "settings.json"
        target = self.base / "target.json"
        self.write(settings, json.dumps(self.settings()))
        self.write(target, json.dumps({"profile": profile, "stem": stem}))
        self.execute([str(ROOT / "render.R"), str(self.analysis), "sans", str(self.preprocess),
                      str(ROOT / "plots"), "", str(output), str(settings), str(target)])
        self.assert_single(output, profile, stem)
        return output

    def test_alpha_only_does_not_require_other_alpha_exports(self):
        self.render_diversity("alpha", "alpha/shannon")

    def test_single_pcoa_does_not_save_other_metrics_or_combined(self):
        self.render_diversity("beta", "beta/PCoA_bray_curtis")

    def test_combined_pcoa_does_not_save_individual_figures(self):
        output = self.render_diversity("beta", "beta/PCoA_four_metrics")
        # PNG IHDR stores the exact selected final canvas dimensions at 72 DPI.
        import struct
        dimensions = struct.unpack(">II", (output / "beta/PCoA_four_metrics.png").read_bytes()[16:24])
        self.assertEqual(dimensions, (288, 216))

    def test_taxonomy_target_overrides_rank_view_defaults(self):
        self.render_diversity("taxonomy", "taxonomy/genus_group_mean_unrarefied")

    def render_differential(self, profile, suffix, source=None, comparison_key=None):
        stem = "02_differential/L6_genus/" + suffix
        output = self.base / self.id().split(".")[-1]
        output.mkdir()
        target = {"profile": profile, "stem": stem}
        if comparison_key:
            target["comparison_key"] = comparison_key
        request = {"analysis_run": str(source or self.diff), "plot_target": target,
                   "plot_settings": self.settings(),
                   "group_style": {"group_order": ["HDC", "LDC", "CON"]}}
        self.write(output / "options.json", json.dumps(request))
        self.execute([str(ROOT / "differential_ancombc.R"), "--plot-ancom",
                      str(output / "options.json"), str(output)])
        self.assert_single(output, profile, stem)
        self.assertFalse((output / "02_differential/L2_phylum").exists())
        return output

    def test_heatmap_only_preserves_joint_all_group_data(self):
        output = self.render_differential("heatmap", "abundance/abundance_heatmap")
        self.assertFalse((output / "02_differential/L6_genus/comparisons").exists())
        self.assertFalse((output / "02_differential/L6_genus/abundance/bubble_plot_input.tsv").exists())
        rows = (output / "02_differential/L6_genus/abundance/heatmap_plot_input.tsv").read_text()
        for group in ("CON", "LDC", "HDC"):
            self.assertIn(group, rows)

    def test_bubble_only_does_not_construct_heatmap_outputs(self):
        output = self.render_differential("bubble", "abundance/abundance_bubble")
        self.assertFalse((output / "02_differential/L6_genus/abundance/heatmap_plot_input.tsv").exists())

    def test_effect_identity_stays_selected_after_group_order_changes(self):
        output = self.render_differential("effect", "comparisons/pair_001/effect", comparison_key="GroupLDC")
        info = json.loads((output / "02_differential/L6_genus/comparisons/pair_001/comparison.json").read_text())
        self.assertEqual(info["contrast"], "GroupLDC")
        self.assertEqual(info["positive_group"], "LDC")
        self.assertEqual(info["negative_group"], "CON")
        self.assertFalse((output / "02_differential/L6_genus/abundance").exists())

    def test_qiime_ancom_heatmap_uses_only_saved_reject_decision(self):
        output = self.render_differential("heatmap", "abundance/abundance_heatmap", source=self.qiime)
        self.assertFalse((output / "02_differential/L6_genus/ancom_W_plot_settings.json").exists())
        rows = (output / "02_differential/L6_genus/abundance/heatmap_plot_input.tsv").read_text()
        self.assertIn("taxon_A", rows)
        self.assertNotIn("taxon_B", rows)

    def test_qiime_ancom_W_only_skips_abundance(self):
        output = self.render_differential("effect", "ancom_W", source=self.qiime)
        self.assertFalse((output / "02_differential/L6_genus/abundance").exists())


if __name__ == "__main__":
    unittest.main()
