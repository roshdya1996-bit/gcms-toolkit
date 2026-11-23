import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from gcms_toolkit.retention_index import (
    Alkane,
    bracket,
    load_alkane_csv,
    retention_index,
    retention_time,
    ri_report,
)

LADDER = [
    Alkane(10, 5.0),
    Alkane(12, 10.0),
    Alkane(14, 15.0),
    Alkane(16, 20.0),
]

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")


class TestRetentionIndex(unittest.TestCase):
    def test_exact_alkane_rt_gives_100n(self):
        self.assertEqual(retention_index(10.0, LADDER), 1200.0)

    def test_midpoint_between_c12_and_c14(self):
        # halfway between C12 (10 min) and C14 (15 min): 1200 + 200 * 0.5
        self.assertAlmostEqual(retention_index(12.5, LADDER), 1300.0)

    def test_roundtrip(self):
        rt = retention_time(1340.0, LADDER)
        self.assertAlmostEqual(retention_index(rt, LADDER), 1340.0, places=9)

    def test_retention_time_value(self):
        # RI 1300 with the uniform test ladder -> 12.5 min
        self.assertAlmostEqual(retention_time(1300.0, LADDER), 12.5)

    def test_outside_ladder_raises(self):
        with self.assertRaises(ValueError):
            retention_index(30.0, LADDER)
        with self.assertRaises(ValueError):
            retention_time(500.0, LADDER)

    def test_unsorted_ladder_is_sorted_and_validated(self):
        shuffled = [LADDER[2], LADDER[0], LADDER[3], LADDER[1]]
        self.assertAlmostEqual(retention_index(12.5, shuffled), 1300.0)

    def test_decreasing_ladder_rejected(self):
        bad = [Alkane(10, 10.0), Alkane(12, 5.0)]
        with self.assertRaises(ValueError):
            retention_index(7.0, bad)

    def test_bracket(self):
        prev, nxt = bracket(11.0, LADDER)
        self.assertEqual((prev.carbons, nxt.carbons), (12, 14))

    def test_report(self):
        rows = ri_report([("b", 12.5), ("a", 7.5)], LADDER)
        self.assertEqual([r[0] for r in rows], ["a", "b"])
        self.assertAlmostEqual(rows[0][2], 1100.0)

    def test_example_csv_loads(self):
        ladder = load_alkane_csv(os.path.join(DATA_DIR, "alkanes_example.csv"))
        self.assertEqual(ladder[0].carbons, 8)
        self.assertEqual(ladder[-1].carbons, 32)


if __name__ == "__main__":
    unittest.main()
