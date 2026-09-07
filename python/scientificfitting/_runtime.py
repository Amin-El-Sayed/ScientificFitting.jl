"""Load the numerical backend lazily, without importing any plotting package."""

from functools import cache
from pathlib import Path


@cache
def _backend():
    import juliacall

    bridge = juliacall.newmodule("ScientificFittingPython")
    bridge.include(str(Path(__file__).with_name("_bridge.jl")))
    if not bridge.supports_finite_derivatives:
        raise ImportError(
            "This development wrapper needs the matching Julia checkout. "
            "Run `python python/develop.py` from the repository, then restart Python."
        )
    return bridge
