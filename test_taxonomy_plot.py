"""Opt-in real R checks for cleaning and taxonomy previews from saved counts."""
import csv
import hashlib
import json
import os
import re
from pathlib import Path
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET

from abundance_tables import clean_taxonomy_table

ROOT = Path(__file__).resolve().parent
RSCRIPT = os.environ.get('WORKSTATION_RSCRIPT')


@unittest.skipUnless(RSCRIPT, 'Set WORKSTATION_RSCRIPT for real R taxonomy checks')
class TaxonomyPlots(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='16s-taxonomy-')
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        self.source = self.base / 'source'
        self.source.mkdir()
        self.ids = ['001', 'X9', 'PIGz', 'k21', 'foo']
        self.groups = ['LDC', 'CON', 'HDC', 'CON', 'LDC']
        self.features = [f'k__Bacteria; p__P{i:02}; c__C{i:02}; o__O{i:02}; f__F{i:02}; g__G{i:02}; s__S{i:02}' for i in range(1, 25)]
        self.features += ['k__Bacteria; p__P01; c__C01; o__O01; f__F01; g__uncultured_x',
                          'k__Bacteria; p__P02; c__C02; o__O02; f__F02; g__G01',
                          'D_0__Bacteria; D_1__P01; D_2__C01; D_3__O01; D_4__F01; D_5__Incertae_Sedis',
                          'k__Bacteria; p__P01; g__']
        self.values = [[100-i*3, (100-i*3)*10, 80-i*2, 20+i, 30+i] for i in range(1, 25)]
        self.values += [[5000, 4000, 3000, 2000, 1000], [4, 3, 2, 1, 0], [2, 1, 1, 1, 0], [0, 2, 3, 4, 5]]
        self.write(self.source / 'metadata.tsv', 'sample-id\tGroup\n' + ''.join(f'{s}\t{g}\n' for s, g in zip(self.ids, self.groups)))
        self.write(self.source / 'exported-table/feature-table.tsv', '#OTU ID\t' + '\t'.join(self.ids) + '\n' + ''.join(
            f'ASV{i}\t' + '\t'.join(map(str, row)) + '\n' for i, row in enumerate(self.values)))
        self.write(self.source / 'exported-taxonomy/taxonomy.tsv', 'Feature ID\tTaxon\n' + ''.join(
            f'ASV{i}\t{lineage}\n' for i, lineage in enumerate(self.features)))
        self.source_hashes = self.hashes()

    def hashes(self):
        return {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in self.source.rglob('*') if p.is_file()}

    @staticmethod
    def write(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(value, encoding='utf8')

    @staticmethod
    def read_csv(path):
        with path.open(encoding='utf-8-sig', newline='') as stream:
            return list(csv.DictReader(stream))

    def execute(self, args):
        result = subprocess.run([RSCRIPT, '--vanilla', *args], cwd=ROOT, capture_output=True, text=True, timeout=90)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def render(self, name, *, palette='user_R', target=True):
        output = self.base / name
        settings = self.base / f'{name}.json'
        self.write(settings, json.dumps({'general': {'dpi': 72, 'formats': 'png_svg'},
                                        'taxonomy': {'width': 9, 'height': 7, 'palette': palette}}))
        args = [str(ROOT / 'render.R'), str(self.source), 'sans', str(self.source), str(ROOT / 'plots'), '', str(output), str(settings)]
        if target:
            target_path = self.base / 'target.json'
            self.write(target_path, json.dumps({'profile': 'taxonomy', 'stem': 'taxonomy/genus_group_mean_unrarefied'}))
            args.append(str(target_path))
        self.execute(args)
        self.assertEqual(self.hashes(), self.source_hashes)
        return output

    def assert_clean(self, output, rank):
        names = {2:'phylum',3:'class',4:'order',5:'family',6:'genus',7:'species'}
        counts, rel, _, _ = clean_taxonomy_table((self.ids, self.features, self.values), rank)
        directory = output / f'taxonomy/cleaned_tables/L{rank}_{names[rank]}'
        for key, table in [('counts', counts), ('relative', rel)]:
            actual = {row['FeatureID']: row for row in self.read_csv(directory / f'{key}_cleaned.csv')}
            self.assertEqual(set(actual), set(table[1]))
            for label, values in zip(table[1], table[2]):
                for sample, expected in zip(self.ids, values):
                    self.assertAlmostEqual(float(actual[label][sample]), float(expected), places=12)

    def test_full_render_cleans_all_levels_and_only_collapses_requested_ranks(self):
        # A full render expects all seven saved alpha exports; these are fixtures,
        # not values computed by this test or evidence of statistical fitting.
        for metric in ('observed_features','shannon','evenness','faith_pd','chao1','simpson','goods_coverage'):
            self.write(self.source / f'exports/alpha/{metric}/alpha-diversity.tsv',
                       f'sample-id\t{metric}\n' + ''.join(f'{sid}\t0.5\n' for sid in self.ids))
        self.source_hashes = self.hashes()
        output = self.render('all', target=False)
        for rank in range(2, 8):
            self.assert_clean(output, rank)
        for rank in ('order', 'family', 'genus'):
            report = json.loads((output / f'taxonomy/{rank}_plot_report.json').read_text())
            self.assertEqual(report['top_n'], 20)
        for rank in ('phylum', 'class'):
            report = json.loads((output / f'taxonomy/{rank}_plot_report.json').read_text())
            self.assertIsNone(report['top_n'])
        self.assertEqual(len(list((output / 'taxonomy').rglob('*.png'))), 10)

    def test_top20_others_palette_metadata_and_single_preview(self):
        output = self.render('user')
        self.assert_clean(output, 6)
        report = json.loads((output / 'taxonomy/genus_plot_report.json').read_text())
        taxa = report['bottom_to_top']
        self.assertEqual(len(taxa), 21)
        self.assertEqual(taxa[-1], 'Others')
        self.assertEqual(report['colors'][taxa[0]], '#575A8F')
        self.assertEqual(report['colors']['Others'], '#857D83')
        # Check actual SVG geometry: y increases down the page. Each group's
        # largest globally ranked taxon must be at the bottom, Others at the top.
        stacks = {}
        for rect in ET.parse(output / 'taxonomy/genus_group_mean_unrarefied.svg').iter('{http://www.w3.org/2000/svg}rect'):
            color = re.search(r'fill: (#[0-9A-Fa-f]{6});', rect.get('style', ''))
            if color and color[1] in report['colors'].values() and float(rect.get('width', '0')) > 30:
                stacks.setdefault(rect.get('x'), []).append((float(rect.get('y')), color[1]))
        self.assertEqual(len(stacks), 3)
        for bars in stacks.values():
            self.assertEqual([color for _,color in sorted(bars, reverse=True)], [report['colors'][t] for t in taxa])
        clean = {row['FeatureID']: row for row in self.read_csv(output / 'taxonomy/cleaned_tables/L6_genus/relative_cleaned.csv')}
        ordered = sorted((t for t in clean if t!='Others'), key=lambda t:(-sum(float(clean[t][s]) for s in self.ids)/len(self.ids),t))
        self.assertEqual(taxa[:-1], ordered[:20])
        rows = self.read_csv(output / 'taxonomy/genus_sample_top20_plot_input.csv')
        for sid in self.ids:
            self.assertAlmostEqual(sum(float(row[sid]) for row in rows), 1)
        means = self.read_csv(output / 'taxonomy/genus_group_mean_top20_plot_input.csv')
        for row in means:
            expected = sum(float(next(r for r in rows if r['Taxon']==row['Taxon'])[sid]) for sid,g in zip(self.ids,self.groups) if g==row['Group'])/self.groups.count(row['Group'])
            self.assertAlmostEqual(float(row['Abundance']), expected)
        self.assertEqual(len(list(output.rglob('*.png'))), 1)
        original = hashlib.sha256((output / 'taxonomy/genus_group_mean_unrarefied.svg').read_bytes()).hexdigest()
        recolored = self.render('tableau', palette='Tableau 20')
        updated = json.loads((recolored / 'taxonomy/genus_plot_report.json').read_text())
        self.assertEqual(updated['colors'][taxa[0]], '#1F77B4')
        self.assertEqual((output / 'taxonomy/genus_sample_top20_plot_input.csv').read_bytes(), (recolored / 'taxonomy/genus_sample_top20_plot_input.csv').read_bytes())
        self.assertEqual(hashlib.sha256((output / 'taxonomy/genus_group_mean_unrarefied.svg').read_bytes()).hexdigest(), original)

    def test_helper_empty_rank_zero_sample_all_others_and_palette_extensions(self):
        # Exercises edge cases without requiring a full analysis or model.
        script = f'''
source({json.dumps(str(ROOT / 'taxonomy_plot.R'))})
stopifnot(taxonomy_rank_label('P__Known; g__unknown',2)=='Known')
counts <- matrix(c(0,0,4,6),nrow=2,dimnames=list(NULL,c('001','A')))
clean <- taxonomy_clean_counts(counts,c('p__Known; g__unknown','p__Known; g__'),6)
stopifnot(nrow(clean$counts)==1,all(clean$relative[,1]==0),sum(clean$relative[,2])==1)
selected <- taxonomy_top_table(clean$relative,6,list(top_mode='by_rank',top_n=3))
stopifnot(identical(rownames(selected$relative),'Others'),sum(selected$relative[,2])==1)
rel <- matrix(rep(1/25,50),nrow=25,dimnames=list(paste0('T',1:25),c('001','A')))
stopifnot(nrow(taxonomy_top_table(rel,6,list(top_mode='custom',top_n=3))$relative)==4)
stopifnot(nrow(taxonomy_top_table(rel,6,list(top_mode='all',top_n=3))$relative)==26)
registry <- jsonlite::read_json({json.dumps(str(ROOT / 'taxonomy_palettes.json'))})
for (palette in names(registry)) {{
  colors <- taxonomy_fill_colors(c(rownames(rel),'Others'),palette,registry)
  stopifnot(length(colors)==26,all(grepl('^#[A-Fa-f0-9]{{6}}$',colors)),colors['Others']=='#857D83')
}}
'''
        self.execute(['-e', script])


if __name__ == '__main__':
    unittest.main()
