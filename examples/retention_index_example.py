"""Compute retention indices for peaks against an alkane ladder.

Run from the repo root:

    PYTHONPATH=src python examples/retention_index_example.py
"""

import os

from gcms_toolkit import load_alkane_csv, retention_time, ri_report

here = os.path.dirname(__file__)
ladder = load_alkane_csv(os.path.join(here, "..", "data", "alkanes_example.csv"))

# Peaks from a run on the same method as the ladder injection.
peaks = [
    ("unknown-1", 9.25),
    ("unknown-2", 17.30),
]

for name, rt, ri in ri_report(peaks, ladder):
    print(f"{name:10s}  rt = {rt:6.2f} min   LRI = {ri:7.1f}")

# Reverse: literature RI -> expected retention time on this method.
lit_ri = 1680  # example literature value on a 5% phenyl column
print(f"\nLRI {lit_ri} would elute at {retention_time(lit_ri, ladder):.2f} min on this method")
