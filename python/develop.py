"""Select this checkout in JuliaCall's managed project for wrapper development.

Run after `python -m pip install -e './python[plot,test]'`. The unreleased
derivative API is not present in the registered v0.1.2 numerical core.
"""

from pathlib import Path

import juliapkg

juliapkg.add(
    "ScientificFitting",
    uuid="4ec0e560-7423-4ff9-937f-0d4da6f5f8f8",
    path=str(Path(__file__).resolve().parents[1]),
    dev=True,
)
juliapkg.resolve()
