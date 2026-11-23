# gcms-toolkit

Small, open helpers for GC-MS/MS pesticide residue work:

1. **Linear retention index (LRI) tools** — compute Kovats-style linear
   retention indices from an n-alkane ladder, and invert them to predict
   retention times from literature RI values.
2. **MRM transition list management** — keep triple-quadrupole MRM
   transitions in one reviewable CSV, validate them (one quantifier, at
   least one qualifier per compound, sane m/z and CE ranges), filter by
   compound or class, and export starting templates for Agilent MassHunter,
   SCIEX, and Thermo software.

Pure Python 3.9+ standard library. No dependencies.

> **Status: early stage.** Built by an analytical chemist who runs GC-MS/MS
> pesticide residue analysis (QuEChERS, SANTE 11312/2021, ISO/IEC 17025).
> No users yet — feedback and contributions from residue labs are very
> welcome, especially verified transitions and vendor export formats.

## Quick start

```bash
# run the tests
PYTHONPATH=src python -m unittest discover -s tests -v

# retention index example
PYTHONPATH=src python examples/retention_index_example.py

# filter + export an MRM list
PYTHONPATH=src python examples/export_mrm_example.py
```

## Retention index

```python
from gcms_toolkit import Alkane, retention_index, retention_time

ladder = [Alkane(12, 10.0), Alkane(14, 15.0), Alkane(16, 20.0)]  # C, rt (min)
retention_index(16.2, ladder)   # -> 1440.0
retention_time(1440.0, ladder)  # -> 16.2
```

Use `data/alkanes_example.csv` as the format template and **replace the
times with your own alkane ladder injection** — RI interpolation is only as
good as the ladder run under the same method. Extrapolation beyond the
ladder range raises an error on purpose.

The linear (van den Dool & Kratz) form is used, which is the correct form
for temperature-programmed runs; the classic Kovats equation applies to
isothermal runs.

## MRM transition lists

The bundled list (`data/mrm_transitions_gcms_example.csv`) covers ~35
pesticides common in residue monitoring (organophosphates, pyrethroids,
triazoles, strobilurins and others), each with a quantifier and a qualifier
transition.

```python
from gcms_toolkit import load_transitions, validate, filter_transitions, EXPORTERS

ts = load_transitions("data/mrm_transitions_gcms_example.csv")
assert validate(ts) == []                      # structural + SANTE-style checks
ops = filter_transitions(ts, classes=["organophosphate"])
open("op_screen_agilent.csv", "w").write(EXPORTERS["agilent"](ops))
```

CSV columns: `compound, cas, compound_class, precursor_mz, product_mz,
ce_ev, role, rt_min, note`. `role` is `quant` or `qual`. `ce_ev` and
`rt_min` may be blank.

### Important: verify before you inject

- **Blank CE means "optimize it yourself".** Collision energies are
  instrument- and tune-dependent; even the bundled CE values are starting
  points from the cited sources, not validated settings for your
  instrument.
- The vendor exports are **starting templates** — confirm the column layout
  against your software version and set dwell times / retention windows in
  the acquisition software.
- For confirmed residue results follow SANTE 11312/2021 identification
  criteria (≥ 2 transitions, ion-ratio and retention-time tolerances).

### Data sources

- EURL-SRM (CVUA Stuttgart), *GC-MS/MS-1* method transition list:
  https://www.eurl-pesticides.eu/userfiles/file/EurlSRM/meth-GC-MSMS-1.pdf
- Table 2 of doi:10.3390/ijerph16030472 (ion pairs and collision energies
  for chlorpyrifos, myclobutanil, tebuconazole, bifenthrin,
  lambda-cyhalothrin, beta-cypermethrin, esfenvalerate, difenoconazole).
- NRL / NRC for Grapes (India), *Database on the GC-MS and LC-MS MRM
  transitions* (diazinon, profenofos):
  https://www.nrcgrapes.in/NRL/Database%20on%20the%20GC-MS%20and%20LC-MS%20MRM%20transitions.pdf

Each row's `note` column names its source. When you re-optimize a compound
on your instrument, keep the row and update the values — the git history
then documents your lab's own list.

## Contributing

Issues and pull requests welcome. Most wanted:

- verified transitions (with instrument model and source) for more
  pesticides,
- CE values measured on specific instruments,
- additional vendor export layouts (Shimadzu, Waters),
- retention index reference data on standard phases (DB-5ms, HP-5ms UI).

## License

MIT — see [LICENSE](LICENSE).
