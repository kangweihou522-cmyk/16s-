"""Plot controls remain presentation-only and inherit immutable task snapshots."""
import json
import tempfile
import unittest
from pathlib import Path

import app
from plot_settings import DEFAULTS, resolve_plot_settings, save_plot_settings, ui_schema


class PlotSettingsTests(unittest.TestCase):
    def test_defaults_inheritance_and_no_mutation(self):
        saved = {'plot_settings': {'general': {'dpi': 144}, 'heatmap': {'width': 11, 'show_values': True}}}
        request = {'plot_settings': {'bubble': {'max_size': 8}, 'general': {'formats': 'png_svg_pdf'},
                                    'taxonomy': {'levels': [2, 7]}}}
        before = json.dumps([saved, request, DEFAULTS])
        result = resolve_plot_settings(request, saved)
        self.assertEqual(result['general']['dpi'], 144)
        self.assertEqual(result['general']['formats'], 'png_svg_pdf')
        self.assertEqual(result['heatmap']['width'], 11)
        self.assertEqual(result['bubble']['max_size'], 8)
        self.assertIsNone(result['alpha']['point_size'])
        self.assertEqual(json.dumps([saved, request, DEFAULTS]), before)
        result['taxonomy']['levels'].append(7)
        self.assertEqual(DEFAULTS['taxonomy']['levels'], [2, 3, 4, 5, 6])
        self.assertEqual(request['plot_settings']['taxonomy']['levels'], [2, 7])
        self.assertIsNone(resolve_plot_settings({})['general']['dpi'])
        self.assertIsNone(resolve_plot_settings({'plot_settings': {'heatmap': {'width': ''}}})['heatmap']['width'])

    def test_invalid_controls_rejected(self):
        invalid = [None, {'unknown': {}}, {'heatmap': {'q_value': .05}},
                   {'heatmap': {'width': 0}}, {'general': {'dpi': float('nan')}},
                   {'general': {'dpi': 120.5}}, {'general': {'dpi': True}},
                   {'general': {'font_family': 'missing font'}}, {'bubble': {'point_alpha': 2}},
                   {'alpha': {'show_points': 'false'}}, {'taxonomy': {'levels': []}},
                   {'taxonomy': {'levels': [2, 2]}}, {'taxonomy': {'levels': [True]}},
                   {'taxonomy': {'levels': [8]}}, {'lefse_bar': {'top_n': -1}},
                   {'lefse_cladogram': {'min_point_size': 10}},
                   {'lefse_cladogram': {'max_point_size': 1}},
                   {'lefse_cladogram': {'min_point_size': 3, 'max_point_size': 3}},
                   {'beta': {'title': 'x\x00y'}}]
        for settings in invalid:
            with self.subTest(settings=settings), self.assertRaises(ValueError):
                resolve_plot_settings({'plot_settings': settings})

    def test_recorded_history_and_legacy_defaults_preserve_source(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory); app.marker(run)
            parameters = run / 'parameters.json'
            parameters.write_text(json.dumps({'config': {'samples': [], 'font': 'DejaVu Serif'}}))
            old = parameters.read_bytes()
            self.assertFalse(app.plot_settings_input(str(run))['recorded'])
            settings = resolve_plot_settings({'plot_settings': {'heatmap': {'width': 10}, 'bubble': {'show_values': True}}})
            save_plot_settings(run, settings)
            loaded = app.plot_settings_input(str(run))
            self.assertTrue(loaded['recorded'])
            self.assertEqual(loaded['plot_settings'], settings)
            self.assertEqual(parameters.read_bytes(), old)
            schema = ui_schema()
            self.assertEqual({p['key'] for p in schema}, set(DEFAULTS))
            for profile in schema:
                self.assertEqual({f['key'] for f in profile['fields']}, set(DEFAULTS[profile['key']]))


if __name__ == '__main__': unittest.main()
