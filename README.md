# gcms-toolkit

Small, open helpers for GC-MS/MS **food additive** analysis:

1. **Linear retention index (LRI) tools** — compute Kovats-style linear
   retention indices from an n-alkane ladder, and invert them to predict
   retention times from literature RI values. Especially useful for
   preservatives like the parabens, which share the same fragment ions on a
   5%-phenyl column and are told apart by retention.
2. **MRM transition list management** — keep triple-quadrupole MRM
   transitions in one reviewable CSV, validate them (one quantifier, at
   least one qualifier per compound, sane m/z and CE ranges), filter by
   compound or class, and export starting templates for Agilent MassHunter,
   SCIEX, and Thermo software.

Pure Python 3.9+ standard library. No dependencies.

> **Status: early stage.** Built by an analytical chemist who runs GC-MS/MS
> food-safety analysis. No users yet — feedback and contributions from food
> and residue labs are very welcome, especially verified transitions for
> more additive classes and vendor export formats.

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
isothermal runs. The bundled preservative list carries published DB-5MS
retention indices in its `note` column, so the RI tools and the MRM list
work together: measure your ladder, predict where each paraben should
elute, then set your MRM time windows.

## MRM transition lists

The bundled list (`data/mrm_transitions_food_example.csv`) covers 25
compounds in two groups relevant to food testing labs:

- **10 preservatives** — the paraben family (methyl- through benzyl-),
  benzoic acid (E210) and benzyl alcohol, with quantifier and qualifier
  transitions, collision energies **and** published DB-5MS retention
  indices.
- **15 phthalate food-contact migrants** (DMP through DNP), with
  quantifier + qualifier transitions and collision energies. Phthalates are
  packaging migrants rather than additives, but food labs test for them in
  the same GC-MS/MS workflows.

```python
from gcms_toolkit import load_transitions, validate, filter_transitions, EXPORTERS

ts = load_transitions("data/mrm_transitions_food_example.csv")
assert validate(ts) == []                      # structural + identification checks
parabens = filter_transitions(ts, classes=["paraben preservative (e218)"])
# or by name: filter_transitions(ts, compounds=["Methylparaben", "Ethylparaben"])
open("paraben_screen_agilent.csv", "w").write(EXPORTERS["agilent"](
    filter_transitions(ts, compounds=["Methylparaben", "Ethylparaben"])))
```

CSV columns: `compound, cas, compound_class, precursor_mz, product_mz,
ce_ev, role, rt_min, note`. `role` is `quant` or `qual`. `rt_min` may be
blank. Class filtering matches the exact class string, case-insensitively —
filter by compound names for finer subsets.

### Important: verify before you inject

- Every bundled transition and CE comes from the cited public sources —
  nothing is estimated — but CEs are instrument- and tune-dependent, so
  re-optimize on your own instrument before validated routine use.
- The parabens all fragment through the common m/z 121 ion on a DB-5ms-type
  column; identification leans on retention time (hence the RI tools), and
  the qualifier/quantifier ion ratio.
- The vendor exports are **starting templates** — confirm the column layout
  against your software version and set dwell times / retention windows in
  the acquisition software.
- For confirmed results follow your accreditation's identification criteria
  (e.g. SANTE 11312/2021: ≥ 2 transitions, ion-ratio and retention-time
  tolerances).

### Data sources

- Preservatives (ion pairs, CEs, DB-5MS retention indices): *Determination
  of 15 preservative allergens in cosmetics by gas chromatography–tandem
  mass spectrometry*, PMC10782272 —
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10782272/ (DB-5MS columns used;
  the same preservatives are regulated food additives E210–E219).
- Phthalates (quantification + identification ion pairs, CEs): SRM method
  for 15 phthalate esters, PMC9743266 —
  https://pmc.ncbi.nlm.nih.gov/articles/PMC9743266/

Each row's `note` column names its source. When you re-optimize a compound
on your instrument, keep the row and update the values — the git history
then documents your lab's own list.

### Planned expansion (contributions welcome)

- Antioxidants (BHA, BHT, TBHQ, gallates) — no verified public MRM table
  bundled yet.
- Synthetic sweeteners and colorants (mostly LC-MS/MS territory; an LC-MS/MS
  list format is a natural next module).
- More vendor export layouts (Shimadzu, Waters).

## Contributing

Issues and pull requests welcome. Most wanted:

- verified transitions (with instrument model and source) for more additive
  classes,
- CE values measured on specific instruments,
- retention index reference data on standard phases (DB-5ms, HP-5ms UI).

## License

MIT — see [LICENSE](LICENSE).
