"""Organize and export MRM transition lists for GC-MS/MS triple-quadrupole work.

The toolkit keeps transitions in one plain CSV (see
``data/mrm_transitions_gcms_example.csv``) so a lab can review, diff and
version them in git, then export the same list into the column layouts the
common instrument software packages import.

IMPORTANT: every bundled transition is a literature/app-note starting point.
Collision energies and even the choice of ions are instrument- and
tune-dependent. Re-optimize on your own instrument (SANTE 11312/2021 expects
identification with at least two transitions and ion-ratio tolerances) before
using a list for residue work.
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

VALID_ROLES = ("quant", "qual")

COLUMNS = [
    "compound",
    "cas",
    "compound_class",
    "precursor_mz",
    "product_mz",
    "ce_ev",
    "role",
    "rt_min",
    "note",
]


@dataclass(frozen=True)
class Transition:
    compound: str
    cas: str
    compound_class: str
    precursor_mz: float
    product_mz: float
    ce_ev: Optional[float] = None  # None = optimize on your own instrument
    role: str = "quant"  # "quant" or "qual"
    rt_min: Optional[float] = None
    note: str = ""

    def __post_init__(self) -> None:
        if self.role not in VALID_ROLES:
            raise ValueError(f"role must be one of {VALID_ROLES}, got {self.role!r}")
        if not (0 < self.product_mz < self.precursor_mz):
            raise ValueError(
                f"{self.compound}: product m/z must be positive and below the "
                f"precursor m/z ({self.precursor_mz} -> {self.product_mz})"
            )
        if self.ce_ev is not None and not (0 <= self.ce_ev <= 80):
            raise ValueError(f"{self.compound}: collision energy {self.ce_ev} eV out of range 0-80")


def load_transitions(path: str) -> List[Transition]:
    """Load transitions from the toolkit CSV format."""
    out: List[Transition] = []
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(
            line for line in fh if not line.lstrip().startswith("#")
        )
        missing = set(COLUMNS) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing columns: {sorted(missing)}")
        for row in reader:
            out.append(
                Transition(
                    compound=row["compound"].strip(),
                    cas=row["cas"].strip(),
                    compound_class=row["compound_class"].strip(),
                    precursor_mz=float(row["precursor_mz"]),
                    product_mz=float(row["product_mz"]),
                    ce_ev=float(row["ce_ev"]) if row["ce_ev"].strip() else None,
                    role=(row["role"].strip().lower() or "quant"),
                    rt_min=float(row["rt_min"]) if row["rt_min"].strip() else None,
                    note=row["note"].strip(),
                )
            )
    return out


def write_transitions(transitions: Sequence[Transition], path: str) -> None:
    """Write transitions back to the toolkit CSV format (sorted by compound)."""
    rows = sorted(transitions, key=lambda t: (t.compound.lower(), t.role != "quant", -t.product_mz))
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(COLUMNS)
        for t in rows:
            writer.writerow(
                [
                    t.compound,
                    t.cas,
                    t.compound_class,
                    _num(t.precursor_mz),
                    _num(t.product_mz),
                    "" if t.ce_ev is None else _num(t.ce_ev),
                    t.role,
                    "" if t.rt_min is None else _num(t.rt_min),
                    t.note,
                ]
            )


def _num(x: float) -> str:
    return f"{x:g}"


def validate(transitions: Sequence[Transition]) -> List[str]:
    """Sanity-check a transition list. Returns a list of problems (empty = OK).

    Checks follow common triple-quad practice and SANTE 11312/2021
    identification needs: each compound needs exactly one quantifier and at
    least one qualifier transition.
    """
    problems: List[str] = []
    by_compound: Dict[str, List[Transition]] = {}
    for t in transitions:
        by_compound.setdefault(t.compound, []).append(t)
    for name, ts in sorted(by_compound.items()):
        quants = [t for t in ts if t.role == "quant"]
        quals = [t for t in ts if t.role == "qual"]
        if len(quants) != 1:
            problems.append(f"{name}: expected exactly 1 quant transition, found {len(quants)}")
        if not quals:
            problems.append(f"{name}: no qualifier transition (SANTE needs >= 2 transitions)")
        seen = set()
        for t in ts:
            key = (t.precursor_mz, t.product_mz)
            if key in seen:
                problems.append(f"{name}: duplicate transition {key[0]:g} -> {key[1]:g}")
            seen.add(key)
    return problems


def filter_transitions(
    transitions: Sequence[Transition],
    compounds: Optional[Sequence[str]] = None,
    classes: Optional[Sequence[str]] = None,
) -> List[Transition]:
    """Subset a list by compound names and/or classes (case-insensitive)."""
    out = list(transitions)
    if compounds:
        wanted = {c.lower() for c in compounds}
        out = [t for t in out if t.compound.lower() in wanted]
    if classes:
        wanted = {c.lower() for c in classes}
        out = [t for t in out if t.compound_class.lower() in wanted]
    return out


def _csv_text(header: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return buf.getvalue()


def export_agilent_masshunter(transitions: Sequence[Transition]) -> str:
    """Export as an Agilent MassHunter dMRM/Dynamic MRM import CSV (starting template).

    Column names follow the MassHunter acquisition method CSV import. Confirm
    against your software version and adjust dwell/segment settings there.
    """
    header = [
        "Compound Name",
        "Precursor Ion",
        "Product Ion",
        "Collision Energy",
        "Retention Time (min)",
        "Qualifier",
    ]
    rows = [
        [
            t.compound,
            _num(t.precursor_mz),
            _num(t.product_mz),
            "" if t.ce_ev is None else _num(t.ce_ev),
            "" if t.rt_min is None else _num(t.rt_min),
            "FALSE" if t.role == "quant" else "TRUE",
        ]
        for t in transitions
    ]
    return _csv_text(header, rows)


def export_sciex(transitions: Sequence[Transition]) -> str:
    """Export as a SCIEX (Analyst/MultiQuant-style) scheduled MRM CSV template."""
    header = [
        "Q1 Mass (Da)",
        "Q3 Mass (Da)",
        "Compound Name",
        "CE",
        "Type",
        "Retention Time (min)",
    ]
    rows = [
        [
            _num(t.precursor_mz),
            _num(t.product_mz),
            t.compound,
            "" if t.ce_ev is None else _num(t.ce_ev),
            "Quantifier" if t.role == "quant" else "Qualifier",
            "" if t.rt_min is None else _num(t.rt_min),
        ]
        for t in transitions
    ]
    return _csv_text(header, rows)


def export_thermo(transitions: Sequence[Transition]) -> str:
    """Export as a Thermo (TraceFinder/TSQ-style) transition CSV template."""
    header = [
        "Name",
        "Precursor (m/z)",
        "Product (m/z)",
        "Collision Energy (V)",
        "Ion Type",
        "RT (min)",
    ]
    rows = [
        [
            t.compound,
            _num(t.precursor_mz),
            _num(t.product_mz),
            "" if t.ce_ev is None else _num(t.ce_ev),
            "Quant" if t.role == "quant" else "Qual",
            "" if t.rt_min is None else _num(t.rt_min),
        ]
        for t in transitions
    ]
    return _csv_text(header, rows)


EXPORTERS = {
    "agilent": export_agilent_masshunter,
    "sciex": export_sciex,
    "thermo": export_thermo,
}
