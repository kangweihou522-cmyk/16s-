"""Run the real argv/parent-directory wrapper against a small fake CLI process."""
import json
import tempfile
import unittest
from pathlib import Path

import app
from run_layout import initialize_layout,artifact_path


class CommandPathTests(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup)
        self.root=Path(temporary.name);self.run=self.root/'run';self.run.mkdir();initialize_layout(self.run,'diversity')
        self.cli=self.root/'conda'
        self.cli.write_text('''#!/usr/bin/env python3
import sys,json
from pathlib import Path
args=sys.argv[sys.argv.index('qiime')+1:]
if args[:2] not in [['diversity','alpha'],['diversity','core-metrics-phylogenetic'],['diversity','beta-group-significance']]:
 sys.exit('Unknown QIIME action: '+repr(args[:2]))
Path('argv.json').write_text(json.dumps(args))
if '--output-dir' in args:
 output=Path(args[args.index('--output-dir')+1])
 if output.exists():sys.exit('Output directory must not be precreated')
 if not output.parent.is_dir():sys.exit('Output parent missing')
 output.mkdir()
else:
 option='--o-alpha-diversity' if '--o-alpha-diversity' in args else '--o-visualization'
 output=Path(args[args.index(option)+1])
 if not output.parent.is_dir():sys.exit('Output parent missing')
 output.write_bytes(b'fake CLI artifact, not actual biology')
''')
        self.cli.chmod(0o755)
        self.config={'conda':str(self.cli),'env':'fake'};self.job={'run':str(self.run),'log':str(self.run/'log.txt')}

    def argv(self):return json.loads((self.run/'argv.json').read_text())

    def test_diversity_alpha_action_is_preserved_in_executed_process(self):
        app.q(self.job,self.config,'diversity','alpha','--i-table','core-metrics-results/rarefied_table.qza','--p-metric','chao1','--o-alpha-diversity','core-metrics-results/chao1_vector.qza')
        args=self.argv();self.assertEqual(args[:2],['diversity','alpha'])
        self.assertEqual(args[args.index('--p-metric')+1],'chao1')
        self.assertTrue(artifact_path(self.run,'core-metrics-results/chao1_vector.qza').is_file())

    def test_scalar_alpha_column_is_preserved_and_missing_parents_created(self):
        app.q(self.job,self.config,'diversity','beta-group-significance','--m-metadata-file','metadata.tsv','--m-metadata-column','alpha','--o-visualization','beta/new/deep/test.qzv')
        args=self.argv();self.assertEqual(args[args.index('--m-metadata-column')+1],'alpha')
        self.assertEqual(args[args.index('--m-metadata-file')+1],'00_inputs/metadata.tsv')
        self.assertTrue(artifact_path(self.run,'beta/new/deep/test.qzv').is_file())

    def test_core_output_parent_is_created_but_output_directory_not_precreated(self):
        app.q(self.job,self.config,'diversity','core-metrics-phylogenetic','--output-dir','core-metrics-results')
        self.assertTrue(artifact_path(self.run,'core-metrics-results').is_dir())


if __name__=='__main__':unittest.main()
