"""Focused statistics/reader checks; synthetic QZV files do not run QIIME."""
import csv
import importlib.util
import io
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import zipfile

BASE = Path(__file__).resolve().parent
RSCRIPT = os.environ.get("WORKSTATION_TEST_RSCRIPT") or shutil.which("Rscript")
SPEC = importlib.util.spec_from_file_location("statistics_app", BASE / "app.py")
app = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(app)


class StatisticsAuditTests(unittest.TestCase):
    def test_constant_pair_keeps_estimable_overall_alpha(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "alpha.tsv"
            path.write_text("id\tvalue\na1\t1\na2\t1\nb1\t1\nb2\t1\nc1\t2\nc2\t3\n")
            rows = [(f"{g}{i}", g) for g in "abc" for i in (1, 2)]
            self.assertEqual(app.alpha_skip_reason(path, rows), "")

    @unittest.skipUnless(RSCRIPT, "Set WORKSTATION_TEST_RSCRIPT for real R reader checks")
    def test_real_r_saved_statistics_keep_valid_results_and_mark_invalid_pairs(self):
        with tempfile.TemporaryDirectory(prefix="statistics-qzv-") as directory:
            root = Path(directory)
            ids = [f"{g}{i}" for g in "ABC" for i in (1, 2, 3)]
            groups = {sid: sid[0] for sid in ids}
            raw = "SubjectID1\tSubjectID2\tGroup1\tGroup2\tDistance\n"
            for i, first in enumerate(ids):
                for j in range(i + 1, len(ids)):
                    second = ids[j]
                    raw += f"{first}\t{second}\t{groups[first]}\t{groups[second]}\t{j-i}\n"
            overall = {"method name": "PERMANOVA", "sample size": 9,
                       "number of groups": 3, "test statistic": 1.4,
                       "p-value": .04, "number of permutations": 999}

            def beta_archive(name, invalid=False, stale=False, invalid_global=False):
                overview = overall | ({"test statistic": "nan"} if invalid_global else {})
                html = "<table>" + "".join(f"<tr><th>{k}</th><td>{v}</td></tr>" for k, v in overview.items()) + "</table>"
                table = io.StringIO()
                writer = csv.writer(table)
                writer.writerow(["Group 1", "Group 2", "Sample size", "Permutations", "pseudo-F", "p-value", "q-value"])
                writer.writerow(["A", "B", 6, 999, 1.3, .03, .06])
                writer.writerow(["A", "C", 6, 999, "nan" if invalid else 1.2, .001 if invalid else .2, .003 if invalid else .2])
                writer.writerow(["B", "C", 6, 999, 2, .04, .06])
                with zipfile.ZipFile(root / name, "w") as archive:
                    archive.writestr("fixture/metadata.yaml", "uuid: fixture\ntype: Visualization\n")
                    archive.writestr("fixture/data/index.html", html)
                    archive.writestr("fixture/data/raw_data.tsv", raw.replace("A1\tA2\tA\tA", "A1\tA2\tB\tA") if stale else raw)
                    archive.writestr("fixture/data/permanova-pairwise.csv", table.getvalue())

            beta_archive("valid.qzv")
            beta_archive("invalid_pair.qzv", invalid=True)
            beta_archive("stale.qzv", stale=True)
            beta_archive("invalid_global.qzv", invalid_global=True)
            ordination = "Eigvals\t2\n0.6\t0.4\n\nProportion explained\t2\n0.6\t0.4\n\nSite\t9\t2\n" + "".join(
                f"{sid}\t{i}\t{i*i}\n" for i, sid in enumerate(ids))
            with zipfile.ZipFile(root / "pcoa.qza", "w") as archive:
                archive.writestr("fixture/metadata.yaml", "uuid: fixture\ntype: PCoAResults\n")
                archive.writestr("fixture/data/ordination.txt", ordination)
            alpha_metadata = "id\tGroup\tshannon\n#q2:types\tcategorical\tnumeric\n" + "".join(
                f"{sid}\t{groups[sid]}\t{1 if sid[0] in 'AB' else int(sid[1])+1}\n" for sid in ids)
            with zipfile.ZipFile(root / "alpha.qzv", "w") as archive:
                archive.writestr("fixture/metadata.yaml", "uuid: fixture\ntype: Visualization\n")
                archive.writestr("fixture/data/column-Group.jsonp", 'load_data("Group", {"H": 6, "p": 0.05});')
                archive.writestr("fixture/data/metadata.tsv", alpha_metadata)
                archive.writestr("fixture/data/kruskal-wallis-pairwise-Group.csv",
                    "Group 1,Group 2,H,p-value,q-value\nA (n=3),C (n=3),4,0.04,0.04\nB (n=3),C (n=3),4,0.04,0.04\n")
            script = root / "check.R"
            script.write_text(r'''
args <- commandArgs(TRUE)
source(args[1]); root <- args[2]
source(file.path(dirname(args[1]),'plot_beta_original_style.R'))
ids <- unlist(lapply(c('A','B','C'), function(g) paste0(g,1:3)))
meta <- data.frame(SampleID=ids,Group=factor(substr(ids,1,1)))
dm <- abs(outer(seq_along(ids),seq_along(ids),'-')); dimnames(dm) <- list(ids,ids)
ord <- q2_read_pcoa(file.path(root,'pcoa.qza'))
stopifnot(identical(rownames(ord$coords),ids),identical(ord$proportion,c(.6,.4)),
          identical(unname(ord$coords[,1]),as.numeric(0:8)))
plot <- original_pcoa_plot(data.frame(ord$coords,Group=meta$Group),ord$proportion,
                          group_colors=c(A='red',B='blue',C='black'),ellipse=FALSE)
stopifnot(plot$coordinates$ratio == 1)
read <- function(name) q2_read_qiime_stats(file.path(root,name),'permanova',dm,meta,999)
ok <- read('valid.qzv')
stopifnot(ok$overall$P == .04, all(ok$pairwise$Status == 'ok'),
          identical(ok$pairwise$q_within_metric,c(.06,.2,.06)))
bad <- read('invalid_pair.qzv')
stopifnot(bad$overall$P == .04, bad$pairwise$P[1] == .03,
          is.na(bad$pairwise$P[2]),bad$pairwise$P_saved[2] == .001,
          all(is.na(bad$pairwise$q_within_metric)),
          bad$pairwise$Status[2] == 'invalid_saved_test')
stopifnot(inherits(try(read('stale.qzv'),silent=TRUE),'try-error'),
          inherits(try(read('invalid_global.qzv'),silent=TRUE),'try-error'))
values <- c(rep(1,6),2:4)
alpha <- q2_read_alpha_stats(file.path(root,'alpha.qzv'),'shannon',
             data.frame(Sample=ids,Group=meta$Group,Value=values))
stopifnot(alpha$overall$KW_p == .05,
          alpha$pairwise$Status[1] == 'not_in_qiime_output',
          is.na(alpha$pairwise$P[1]),all(alpha$pairwise$q_within_metric[-1] == .04))
cat('PASS: real R readers, missing/invalid pairs, original q values, stale and invalid overall rejected\n')
''')
            result = subprocess.run([RSCRIPT, str(script), str(BASE / "plots/qiime2_io.R"), str(root)], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("PASS:", result.stdout)

    @unittest.skipUnless(RSCRIPT, "Set WORKSTATION_TEST_RSCRIPT for real R rendering")
    def test_real_r_zero_distance_placeholder_preserves_other_metrics(self):
        import test_render_fixture as fixture
        with tempfile.TemporaryDirectory(prefix="statistics-render-") as directory:
            with patch.object(fixture, "ROOT", Path(directory)):
                def mutate(run, parent):
                    metric = "weighted_unifrac"
                    fixture.write(run / f"exports/distance/{metric}/distance-matrix.tsv",
                        "\t" + "\t".join(fixture.IDS) + "\n" + "".join(
                            sid + "\t" + "\t".join(["0"] * len(fixture.IDS)) + "\n" for sid in fixture.IDS))
                    path = run / f"exports/pcoa/{metric}/ordination.txt"
                    path.write_text(path.read_text().replace("Proportion explained\t2\n0.6\t0.4", "Proportion explained\t2\nNaN\tNaN"))
                run, parent = fixture.fixture("zero", mutate)
                result = subprocess.run([RSCRIPT, str(BASE / "render.R"), str(run), "DejaVu Serif", str(parent), str(BASE / "plots")], text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                out = run / "R_original_style/beta"
                self.assertIn("all pairwise distances are zero", (out / "PCoA_weighted_unifrac.svg").read_text())
                self.assertFalse((out / "weighted_unifrac_PCoA_coordinates.csv").exists())
                self.assertIn("undefined_zero_distances", (out / "weighted_unifrac_PCoA_status.csv").read_text())
                for metric in fixture.METRICS[1:]:
                    rows = fixture.read_csv(out / f"{metric}_PCoA_coordinates.csv")
                    self.assertEqual([row["SampleID"] for row in rows], fixture.ORDER)
                    self.assertTrue(all(float(row["PC1"]) == fixture.COORDS[row["SampleID"]][0] for row in rows))
                self.assertTrue((run / "R_original_style/alpha/shannon.svg").exists())


if __name__ == "__main__":
    unittest.main()
