"""Linear (Kovats-style) retention index utilities for temperature-programmed GC.

The linear retention index (LRI) is the standard form for programmed-temperature
GC-MS work:

    LRI = 100 * n + 100 * (rt - rt_n) / (rt_{n+1} - rt_n)

where n is the carbon number of the n-alkane eluting just before the analyte.

References
----------
- Kovats, E. (1958). Helv. Chim. Acta 41, 1915-1932.
- van Den Dool, H. & Kratz, P. D. (1963). J. Chromatogr. 11, 463-471.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Sequence, Tuple


@dataclass(frozen=True)
class Alkane:
    """One n-alkane calibration point."""

    carbons: int
    rt_min: float

    def __post_init__(self) -> None:
        if self.carbons < 1:
            raise ValueError("carbon number must be >= 1")
        if self.rt_min <= 0:
            raise ValueError("retention time must be positive")


def _validate_ladder(ladder: Sequence[Alkane]) -> List[Alkane]:
    ladder = sorted(ladder, key=lambda a: a.rt_min)
    if len(ladder) < 2:
        raise ValueError("need at least two alkanes to interpolate")
    for prev, nxt in zip(ladder, ladder[1:]):
        if nxt.rt_min <= prev.rt_min:
            raise ValueError("alkane retention times must be strictly increasing")
        if nxt.carbons <= prev.carbons:
            raise ValueError("carbon numbers must increase with retention time")
    return list(ladder)


def bracket(rt_min: float, ladder: Sequence[Alkane]) -> Tuple[Alkane, Alkane]:
    """Return the two alkanes that bracket a retention time.

    Raises ValueError if the retention time falls outside the ladder range
    (extrapolation is deliberately not supported: RI values extrapolated
    beyond the alkane ladder are unreliable).
    """
    ladder = _validate_ladder(ladder)
    if rt_min < ladder[0].rt_min or rt_min > ladder[-1].rt_min:
        raise ValueError(
            f"rt={rt_min} min is outside the alkane ladder "
            f"({ladder[0].rt_min}-{ladder[-1].rt_min} min); extrapolation not supported"
        )
    for prev, nxt in zip(ladder, ladder[1:]):
        if prev.rt_min <= rt_min <= nxt.rt_min:
            return prev, nxt
    raise AssertionError("unreachable")


def retention_index(rt_min: float, ladder: Sequence[Alkane]) -> float:
    """Linear retention index of an analyte eluting at ``rt_min`` minutes."""
    prev, nxt = bracket(rt_min, ladder)
    if nxt.rt_min == prev.rt_min:  # guarded by validation, kept for clarity
        raise ValueError("bracketing alkanes share a retention time")
    span_rt = nxt.rt_min - prev.rt_min
    span_c = nxt.carbons - prev.carbons
    return 100.0 * (prev.carbons + span_c * (rt_min - prev.rt_min) / span_rt)


def retention_time(index: float, ladder: Sequence[Alkane]) -> float:
    """Predict the retention time (min) for a target retention index.

    Inverse of :func:`retention_index`; useful for checking identity when a
    literature RI is available for a column of the same polarity class.
    """
    ladder = _validate_ladder(ladder)
    lo = 100.0 * ladder[0].carbons
    hi = 100.0 * ladder[-1].carbons
    if index < lo or index > hi:
        raise ValueError(
            f"index {index} is outside the ladder range ({lo:.0f}-{hi:.0f})"
        )
    for prev, nxt in zip(ladder, ladder[1:]):
        ri_prev = 100.0 * prev.carbons
        ri_next = 100.0 * nxt.carbons
        if ri_prev <= index <= ri_next:
            frac = (index - ri_prev) / (ri_next - ri_prev)
            return prev.rt_min + frac * (nxt.rt_min - prev.rt_min)
    raise AssertionError("unreachable")


def load_alkane_csv(path: str) -> List[Alkane]:
    """Load an alkane ladder from a CSV with columns ``carbons,rt_min``.

    Lines starting with ``#`` are treated as comments. The file bundled in
    ``data/alkanes_example.csv`` shows the format; replace the times with your
    own ladder injection for real work.
    """
    import csv

    alkanes: List[Alkane] = []
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.reader(line for line in fh if not line.lstrip().startswith("#"))
        header = next(reader)
        if header[:2] != ["carbons", "rt_min"]:
            raise ValueError("expected header 'carbons,rt_min'")
        for row in reader:
            if not row:
                continue
            alkanes.append(Alkane(carbons=int(row[0]), rt_min=float(row[1])))
    return _validate_ladder(alkanes)


def ri_report(peaks: Iterable[Tuple[str, float]], ladder: Sequence[Alkane]) -> List[Tuple[str, float, float]]:
    """Annotate ``(name, rt_min)`` peaks with their retention index.

    Returns rows of ``(name, rt_min, ri)`` sorted by retention time. Peaks
    outside the ladder range raise ValueError.
    """
    ladder = _validate_ladder(ladder)
    out = []
    for name, rt in sorted(peaks, key=lambda p: p[1]):
        out.append((name, rt, retention_index(rt, ladder)))
    return out
