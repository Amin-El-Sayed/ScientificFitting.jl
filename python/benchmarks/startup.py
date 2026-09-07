"""Measure a fresh Python process using the installed package and selected Julia environment.

Run after installation/precompilation to measure restart latency; an empty JuliaPkg
environment also includes provisioning. Timings are observations, not test budgets.
"""

import argparse
import json
from pathlib import Path
import platform
from time import perf_counter


def measure():
    seconds = {}
    start = perf_counter()
    import numpy as np
    import scientificfitting as sf
    seconds["python_import"] = perf_counter() - start

    start = perf_counter()
    import juliacall
    seconds["juliacall_import"] = perf_counter() - start
    from scientificfitting._core import Result, _model_callback
    from scientificfitting._runtime import _backend
    start = perf_counter()
    bridge = _backend()
    seconds["bridge_load"] = perf_counter() - start

    x = np.arange(6.)
    y = np.array([0.8, 2.9, 4.7, 6.8, 8.7, 10.9])
    names = ["slope", "offset"]
    callback = _model_callback(lambda x, slope, offset: slope*x + offset, names)
    start = perf_counter()
    handle = bridge.run_fit("gaussian", callback, x, y, np.array([1., 0.]),
                            {"sigma_y": np.full(6, 0.2)})
    seconds["solve"] = perf_counter() - start
    start = perf_counter()
    result = Result(handle, names, "gaussian")
    seconds["result_snapshot"] = perf_counter() - start
    start = perf_counter()
    result.report()
    result.diagnose()
    seconds["report_diagnose"] = perf_counter() - start

    # A new Python callable must reuse compiled bridge code, not recompile a solver.
    start = perf_counter()
    repeated = sf.fit_model(lambda x, slope, offset: slope*x + offset, x, y,
                           p0={"slope": 1., "offset": 0.}, sigma_y=0.2)
    seconds["new_callback_fit"] = perf_counter() - start
    design = np.column_stack([x, np.ones_like(x)])
    expected = np.linalg.lstsq(design, y, rcond=None)[0]
    covariance = np.linalg.inv(design.T @ design) * 0.2**2
    for checked in (result, repeated):
        assert checked.converged
        np.testing.assert_allclose(checked.params, expected, atol=1e-6)
        np.testing.assert_allclose(checked.covariance, covariance, rtol=1e-5)
    assert not bridge.seval('any(m -> nameof(m) in (:Makie, :CairoMakie), values(Base.loaded_modules))')
    return dict(platform=platform.platform(), python=platform.python_version(),
                julia=bridge.seval("string(VERSION)"), numpy=np.__version__,
                core=bridge.seval("string(pkgversion(ScientificFitting))"),
                core_source=bridge.seval("pathof(ScientificFitting)"),
                seconds=seconds)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = json.dumps(measure(), indent=2)
    print(report)
    if args.output:
        args.output.write_text(report + "\n", encoding="utf-8")
