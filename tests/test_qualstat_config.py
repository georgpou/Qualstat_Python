from __future__ import annotations

import math
import tempfile
import unittest
from pathlib import Path

import yaml

from tests.support import load_script_module, write_text


class QualStatConfigTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.qs = load_script_module("qualstat.py", "qualstat_config_test")

    def base_document(self):
        return {
            "input_file": "data.csv",
            "output_file": "report.log",
            "analysis_type": "abfe",
            "energy_unit": "kJ/mol",
            "bootstrap_rounds": 10,
            "random_seed": 2026,
            "significance_multiplier": 1.645,
            "statistics": {"MAD": True, "R": False},
            "cycle_analysis": False,
        }

    def write_config(self, directory: Path, document: object, name: str = "settings.yaml") -> Path:
        path = directory / name
        write_text(path, yaml.safe_dump(document, sort_keys=False))
        return path

    def test_yaml_root_must_be_a_mapping(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "settings.yaml"
            write_text(path, "- not a mapping\n")
            with self.assertRaisesRegex(ValueError, "YAML root must be a mapping"):
                self.qs.load_config(path)

    def test_required_settings_are_required(self):
        required = ("input_file", "output_file", "analysis_type", "bootstrap_rounds", "statistics")
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            for key in required:
                with self.subTest(key=key):
                    document = self.base_document()
                    document.pop(key)
                    path = self.write_config(directory, document, f"missing_{key}.yaml")
                    with self.assertRaisesRegex(ValueError, f"missing required YAML setting: {key}"):
                        self.qs.load_config(path)

    def test_invalid_analysis_type_and_bootstrap_rounds_are_rejected(self):
        cases = (
            ("analysis_type", "not-an-analysis", "analysis_type must be 'abfe' or 'rbfe'"),
            ("bootstrap_rounds", 1, "bootstrap_rounds must be an integer of at least 2"),
            ("bootstrap_rounds", 2.5, "bootstrap_rounds must be an integer of at least 2"),
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            for key, value, message in cases:
                with self.subTest(key=key, value=value):
                    document = self.base_document()
                    document[key] = value
                    path = self.write_config(directory, document, f"invalid_{key}_{value}.yaml")
                    with self.assertRaisesRegex(ValueError, message):
                        self.qs.load_config(path)

    def test_statistics_mapping_and_metric_switches_are_validated(self):
        cases = (
            ({"statistics": ["MAD"]}, "statistics must be a mapping"),
            ({"statistics": {"unknown": True}}, "unknown statistic: unknown"),
            ({"statistics": {"MAD": True, "mad": True}}, "duplicate statistic: mad"),
            ({"statistics": {"MAD": "yes"}}, "statistics.MAD must be true or false"),
            ({"statistics": {"MAD": False, "R": False}}, "at least one statistic must be true"),
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            for change, message in cases:
                with self.subTest(change=change):
                    document = self.base_document()
                    document.update(change)
                    path = self.write_config(directory, document, f"statistics_{len(list(directory.iterdir()))}.yaml")
                    with self.assertRaisesRegex(ValueError, message):
                        self.qs.load_config(path)

    def test_seed_multiplier_unit_and_cycle_flag_are_validated(self):
        cases = (
            ("random_seed", "2026", "random_seed must be an integer or null"),
            ("significance_multiplier", 0, "significance_multiplier must be a positive finite number"),
            ("significance_multiplier", -1, "significance_multiplier must be a positive finite number"),
            ("significance_multiplier", math.inf, "significance_multiplier must be a positive finite number"),
            ("energy_unit", "", "energy_unit must be a non-empty string"),
            ("cycle_analysis", "yes", "cycle_analysis must be true or false"),
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            for index, (key, value, message) in enumerate(cases):
                with self.subTest(key=key, value=value):
                    document = self.base_document()
                    document[key] = value
                    path = self.write_config(directory, document, f"invalid_{index}.yaml")
                    with self.assertRaisesRegex(ValueError, message):
                        self.qs.load_config(path)

    def test_output_cannot_overwrite_input_csv_or_yaml(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            write_text(directory / "data.csv", "header\n")
            for output_name, message in (
                ("data.csv", "output_file must differ from the input CSV and YAML configuration"),
                ("settings.yaml", "output_file must differ from the input CSV and YAML configuration"),
            ):
                with self.subTest(output_name=output_name):
                    document = self.base_document()
                    document["output_file"] = output_name
                    path = self.write_config(directory, document)
                    with self.assertRaisesRegex(ValueError, message):
                        self.qs.load_config(path)

    def test_input_and_output_paths_are_relative_to_yaml_location(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            settings_directory = directory / "nested" / "analysis"
            settings_directory.mkdir(parents=True)
            path = self.write_config(settings_directory, self.base_document())

            config = self.qs.load_config(path)

        self.assertEqual(config.input_file, settings_directory / "data.csv")
        self.assertEqual(config.output_file, settings_directory / "report.log")
        self.assertEqual(config.statistics, ("MAD",))
        self.assertEqual(config.analysis_type, "abfe")


if __name__ == "__main__":
    unittest.main()
