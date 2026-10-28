"""KGML geometry and cached-statistic rendering checks using synthetic KO values.

Run real R checks with WORKSTATION_RSCRIPT=/path/to/Rscript python -m unittest
test_ko_map_render. The fixture models graphics only, not biological inference.
"""
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess
import tempfile
import unittest

from ko_map_geometry import normalize_map_id, parse_kgml, parse_kgml_text, write_geometry


ROOT = Path(__file__).resolve().parent
RSCRIPT = os.environ.get("WORKSTATION_RSCRIPT")
KGML = '''<?xml version="1.0"?>
<pathway name="path:ko00010" org="ko" number="00010" title="Synthetic positioning fixture">
  <entry id="1" name="ko:K00001 ko:k00002 ko:K00001" type="ortholog" reaction="rn:R00001">
    <graphics name="" type="line" coords="10,20,40,20,40,55"/>
  </entry>
  <entry id="2" name="KO:k00001" type="ortholog">
    <graphics name="Synthetic enzyme A" type="rectangle" x="65" y="20" width="40" height="12"/>
  </entry>
  <entry id="3" name="ko:K00003 ko:K00004" type="ortholog">
    <graphics name="Synthetic group B" type="roundrectangle" x="65" y="55" width="40" height="12"/>
  </entry>
  <entry id="4" name="cpd:C00001" type="compound">
    <graphics name="Synthetic compound" type="circle" x="10" y="55" width="8" height="8"/>
  </entry>
  <entry id="5" name="path:ko00020" type="map">
    <graphics name="Source pathway label" type="roundrectangle" x="65" y="85" width="40" height="12"/>
  </entry>
  <relation entry1="1" entry2="2" type="ECrel"><subtype name="compound" value="4"/></relation>
  <reaction id="1" name="rn:R00001" type="reversible"><substrate id="4" name="cpd:C00001"/><product id="4" name="cpd:C00001"/></reaction>
</pathway>'''
STATS = ("KO\teffect\tq\tcomparison\n"
         "ko:k00001\t2\t0.01\tC001\n"
         "K00001\t2\t0.01\tC001\n"  # identical copy retains row provenance
         "K00002\t-3\t0.02\tC001\n"
         "K00003\t0.4\t0.001\tC001\n"  # filtered by magnitude
         "K00004\t-5\t0.9\tC001\n"  # filtered by q
         "K99999\t4\t0.01\tC001\n"  # valid differential KO absent from map
         "K00004\t-2\t0.01\tC002\n")


class KGMLGeometry(unittest.TestCase):
    def test_original_coordinates_members_and_metadata(self):
        geometry = parse_kgml_text(KGML, "map00010")
        self.assertEqual(geometry["map_id"], "ko00010")
        self.assertEqual(geometry["entries"][0]["ko_ids"], ["K00001", "K00002"])
        self.assertEqual(geometry["entries"][1]["ko_ids"], ["K00001"])
        self.assertEqual(geometry["entries"][0]["graphics"][0]["coords"], [[10, 20], [40, 20], [40, 55]])
        self.assertEqual(geometry["entries"][1]["graphics"][0]["x"], 65)
        self.assertEqual(geometry["counts"], {"entries": 5, "ko_feature_groups": 3, "ko_graphics": 3, "unique_kos": 4})
        self.assertEqual(geometry["relations"][0]["entry2"], "2")
        self.assertEqual(geometry["reactions"][0]["substrates"][0]["id"], "4")
        self.assertEqual(geometry["bounds"], {"xmin": 6, "xmax": 85, "ymin": 14, "ymax": 91})

    def test_display_id_canonicalization_and_source_identity(self):
        self.assertEqual(normalize_map_id("MAP01100"), "ko01100")
        for invalid in ("hsa00010", "00010", "ko0001", "ko000100", "ko00010/../../a"):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                normalize_map_id(invalid)
        with self.assertRaisesRegex(ValueError, "does not match"):
            parse_kgml_text(KGML, "ko01100")
        for content in (KGML.replace('name="path:ko00010"', 'name="path:map00010"'),
                        KGML.replace('org="ko"', 'org="hsa"'), KGML.replace('number="00010"', 'number="01100"')):
            with self.subTest(content=content[:40]), self.assertRaises(ValueError):
                parse_kgml_text(content)

    def test_official_external_doctype_stripped_without_resolving(self):
        for host in ("www.genome.jp", "www.kegg.jp"):
            text = KGML.replace('<pathway ', f'<!DOCTYPE pathway SYSTEM "http://{host}/kegg/xml/KGML_v0.7.1_.dtd">\n<pathway ', 1)
            self.assertEqual(parse_kgml_text(text)["map_id"], "ko00010")

    def test_reject_entity_dtd_and_xml_include(self):
        prefixes = ('<!DOCTYPE pathway SYSTEM "file:///etc/passwd">',
                    '<!DOCTYPE pathway [<!ENTITY a "expanded">]>',
                    '<!DOCTYPE pathway SYSTEM "https://example.com/map.dtd">')
        for prefix in prefixes:
            with self.subTest(prefix=prefix), self.assertRaises(ValueError):
                parse_kgml_text(KGML.replace('<pathway ', prefix+'\n<pathway ', 1))
        with self.assertRaises(ValueError):
            parse_kgml_text(KGML.replace('<pathway ', '<pathway xmlns:xi="http://www.w3.org/2001/XInclude" ', 1))

    def test_invalid_or_absent_ko_graphics_fail_clearly(self):
        with self.assertRaisesRegex(ValueError, "no valid KO"):
            parse_kgml_text(KGML.replace('type="ortholog"', 'type="gene"'))
        for content in (KGML.replace('10,20,40,20,40,55', '10,20,40'),
                        KGML.replace('x="65"', 'x="NaN"'),
                        KGML.replace('id="2"', 'id="1"')):
            with self.assertRaises(ValueError):
                parse_kgml_text(content)

    def test_geometry_written_as_json(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "map.xml"
            source.write_text(KGML, encoding="utf8")
            output = Path(temporary) / "nested/geometry.json"
            geometry = write_geometry(source, output, "ko00010")
            self.assertEqual(json.loads(output.read_text()), geometry)
            self.assertEqual(parse_kgml(source), geometry)


@unittest.skipUnless(RSCRIPT, "Set WORKSTATION_RSCRIPT to enable real R rendering")
class CachedKGMLRendering(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="16s-ko-map-")
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        self.source = self.base / "source"
        self.source.mkdir()
        self.kgml = self.source / "map.xml"
        self.kgml.write_text(KGML, encoding="utf8")
        self.geometry = self.source / "geometry.json"
        write_geometry(self.kgml, self.geometry)
        self.stats = self.source / "saved_differential.tsv"
        self.stats.write_text(STATS, encoding="utf8")
        self.original_hashes = self.hashes()

    def hashes(self):
        return {file.name: hashlib.sha256(file.read_bytes()).hexdigest() for file in self.source.iterdir()}

    def render(self, comparison="C001", q=.05, effect=1, formats="png_svg_pdf", fail=False):
        output = self.base / f"render-{comparison}-{q}-{effect}"
        options = self.base / "options.json"
        options.write_text(json.dumps({"geometry_path": str(self.geometry), "hits_path": str(self.stats),
            "comparison": comparison, "comparison_label": "Synthetic numerator minus reference",
            "effect_label": "synthetic log2 CLR difference", "font": "sans",
            "plot_settings": {"general": {"dpi": 72, "formats": formats, "font_family": "DejaVu Sans", "show_caption": True},
                "metabolic_map": {"width": 7, "height": 5, "q_threshold": q, "effect_threshold": effect}}}), encoding="utf8")
        result = subprocess.run([RSCRIPT, "--vanilla", str(ROOT / "ko_map_render.R"), "--render", str(options), str(output)],
                                cwd=ROOT, capture_output=True, text=True, timeout=90)
        if fail:
            self.assertNotEqual(result.returncode, 0)
            return result
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertEqual(self.hashes(), self.original_hashes)
        return output

    def test_original_geometry_cached_thresholds_mixed_overlap_and_exports(self):
        output = self.render()
        for suffix in ("png", "svg", "pdf"):
            self.assertGreater((output / f"03_figures/metabolic_map.{suffix}").stat().st_size, 1000)
        self.assertEqual(struct.unpack(">II", (output / "03_figures/metabolic_map.png").read_bytes()[16:24]), (504,360))
        summary = json.loads((output / "render_summary.json").read_text())
        self.assertEqual(summary["source"], "ko_network")
        self.assertIs(summary["statistics_refit"], False)
        coverage = summary["coverage"]
        self.assertEqual(coverage["input_rows_all_comparisons"], 7)
        self.assertEqual(coverage["selected_rows"], 6)
        self.assertEqual(coverage["selected_tested_kos"], 5)
        self.assertEqual(coverage["duplicate_selected_rows"], 1)
        self.assertEqual(coverage["mapped_selected_kos"], 4)
        self.assertEqual(coverage["unmapped_selected_kos"], 1)
        self.assertEqual(coverage["selected_significant_kos"], 3)
        self.assertEqual(coverage["mapped_significant_kos"], 2)
        self.assertEqual(coverage["mixed_feature_groups"], 1)
        groups = (output / "02_mapping/feature_group_overlap.tsv").read_text()
        self.assertIn("K00001;K00002", groups)
        self.assertIn("mixed", groups)
        self.assertIn("background", groups)
        self.assertIn("K99999", (output / "02_mapping/unmapped_ko_rows.tsv").read_text())
        manifest = json.loads((output / "figure_manifest.json").read_text())
        self.assertEqual(manifest["figures"][0]["profile"], "metabolic_map")

    def test_selected_comparison_and_zero_hits_are_real_nonhighlight_plots(self):
        output = self.render(comparison="C002", formats="png_svg")
        coverage = json.loads((output / "render_summary.json").read_text())["coverage"]
        self.assertEqual(coverage["selected_tested_kos"], 1)
        self.assertEqual(coverage["highlighted_feature_groups"], 1)
        self.assertEqual(coverage["partial_hit_feature_groups"], 1)
        output = self.render(q=0, effect=100, formats="png_svg")
        coverage = json.loads((output / "render_summary.json").read_text())["coverage"]
        self.assertEqual(coverage["highlighted_feature_groups"], 0)
        self.assertEqual(coverage["mapped_significant_kos"], 0)
        svg = (output / "03_figures/metabolic_map.svg").read_text()
        self.assertIn("No KOs pass", svg)
        self.assertIn("No KO passes thresholds", svg)
        self.assertFalse((output / "03_figures/metabolic_map.pdf").exists())

    def test_conflicting_duplicate_cannot_silently_choose_direction(self):
        self.stats.write_text(STATS + "K00001\t-2\t0.01\tC001\n", encoding="utf8")
        self.original_hashes = self.hashes()
        result = self.render(fail=True)
        self.assertIn("Conflicting cached differential", result.stderr)


if __name__ == "__main__":
    unittest.main()
