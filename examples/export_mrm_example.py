"""Filter the bundled transition list and export vendor CSV templates.

Run from the repo root:

    PYTHONPATH=src python examples/export_mrm_example.py
"""

import os

from gcms_toolkit import export_agilent_masshunter, filter_transitions, load_transitions, validate

here = os.path.dirname(__file__)
transitions = load_transitions(os.path.join(here, "..", "data", "mrm_transitions_food_example.csv"))

problems = validate(transitions)
print(f"loaded {len(transitions)} transitions, validation problems: {problems or 'none'}")

# Keep only the common parabens for a small preservative screen.
screen = filter_transitions(
    transitions,
    compounds=["Methylparaben", "Ethylparaben", "Propylparaben", "Butylparaben", "Benzoic acid"],
)
compounds = sorted({t.compound for t in screen})
print(f"\npreservative screen: {len(compounds)} compounds: {', '.join(compounds)}\n")

print(export_agilent_masshunter(screen))
