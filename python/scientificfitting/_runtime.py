"""Load the numerical backend lazily, without importing any plotting package."""

from functools import cache
from pathlib import Path


@cache
def _backend():
    import juliacall

    bridge = juliacall.newmodule("ScientificFittingPython")
    # JuliaPkg enforces this pin normally; check explicit environment overrides too.
    core_version = bridge.seval("using ScientificFitting; pkgversion(ScientificFitting)")
    if not bridge.seval("v -> v\"0.3.0\" <= v < v\"0.4.0\"")(core_version):
        raise ImportError(
            f"scientificfitting requires Julia core 0.3.x, but loaded {core_version}. "
            "Use a matching Julia environment and restart Python. For a source "
            "checkout, run `python python/develop.py` before importing the package."
        )
    bridge.include(str(Path(__file__).with_name("_bridge.jl")))
    return bridge
