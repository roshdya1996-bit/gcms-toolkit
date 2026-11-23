import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from gcms_toolkit.mrm import (
    Transition,
    export_agilent_masshunter,
    export_sciex,
    export_thermo,
    filter_transitions,
    load_transitions,
    validate,
    write_transitions,
)

DATA = os.path.join(os.path.dirname(__file__), "..", "data", "mrm_transitions_gcms_example.csv")


class TestTransitionData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.transitions = load_transitions(DATA)

    def test_bundled_list_passes_validation(self):
        self.assertEqual(validate(self.transitions), [])

    def test_every_compound_has_quant_and_qual(self):
        compounds = {t.compound for t in self.transitions}
        self.assertGreaterEqual(len(compounds), 25)
        for c in compounds:
            roles = [t.role for t in self.transitions if t.compound == c]
            self.assertIn("quant", roles, c)
            self.assertIn("qual", roles, c)

    def test_product_below_precursor(self):
        for t in self.transitions:
            self.assertLess(t.product_mz, t.precursor_mz, t.compound)

    def test_blank_ce_loads_as_none(self):
        diazinon = [t for t in self.transitions if t.compound == "Diazinon"]
        self.assertTrue(all(t.ce_ev is None for t in diazinon))

    def test_known_ce_loads(self):
        chlorpyrifos = [t for t in self.transitions if t.compound == "Chlorpyrifos" and t.role == "quant"][0]
        self.assertEqual((chlorpyrifos.precursor_mz, chlorpyrifos.product_mz, chlorpyrifos.ce_ev), (314.0, 258.0, 15.0))


class TestValidation(unittest.TestCase):
    def test_missing_qualifier_flagged(self):
        ts = [Transition("X", "0-00-0", "test", 200.0, 100.0, 10.0, "quant")]
        problems = validate(ts)
        self.assertTrue(any("no qualifier" in p for p in problems))

    def test_two_quantifiers_flagged(self):
        ts = [
            Transition("X", "0-00-0", "test", 200.0, 100.0, 10.0, "quant"),
            Transition("X", "0-00-0", "test", 200.0, 150.0, 10.0, "quant"),
            Transition("X", "0-00-0", "test", 200.0, 90.0, 10.0, "qual"),
        ]
        self.assertTrue(any("exactly 1 quant" in p for p in validate(ts)))

    def test_bad_physics_rejected(self):
        with self.assertRaises(ValueError):
            Transition("X", "0-00-0", "test", 100.0, 200.0, 10.0, "quant")
        with self.assertRaises(ValueError):
            Transition("X", "0-00-0", "test", 200.0, 100.0, 500.0, "quant")


class TestFilterAndExport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.transitions = load_transitions(DATA)

    def test_filter_by_class(self):
        ops = filter_transitions(self.transitions, classes=["organophosphate"])
        self.assertTrue(ops)
        self.assertTrue(all(t.compound_class == "organophosphate" for t in ops))

    def test_filter_by_compound_case_insensitive(self):
        sub = filter_transitions(self.transitions, compounds=["malathion"])
        self.assertEqual({t.compound for t in sub}, {"Malathion"})

    def test_exporters_produce_csv_with_all_rows(self):
        n = len(self.transitions)
        for export in (export_agilent_masshunter, export_sciex, export_thermo):
            text = export(self.transitions)
            lines = [l for l in text.strip().split("\n")]
            self.assertEqual(len(lines), n + 1, export.__name__)
        self.assertIn("Compound Name", export_agilent_masshunter(self.transitions).splitlines()[0])
        self.assertIn("Q1 Mass (Da)", export_sciex(self.transitions).splitlines()[0])
        self.assertIn("Precursor (m/z)", export_thermo(self.transitions).splitlines()[0])

    def test_roundtrip_write_load(self):
        import tempfile

        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "rt.csv")
            write_transitions(self.transitions, path)
            again = load_transitions(path)
        self.assertEqual(len(again), len(self.transitions))
        self.assertEqual(validate(again), [])


if __name__ == "__main__":
    unittest.main()
