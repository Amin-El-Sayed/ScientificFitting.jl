"""Measure a fresh Python process using the installed package and selected Julia environment.

Run after installation/precompilation to measure restart latency; an empty JuliaPkg
environment also includes provisioning. Timings are observations, not test budgets.
"""

import argparse
from importlib.metadata import version
import json
from math import lgamma
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
from time import perf_counter


CASES = ("gaussian", "inplace", "custom", "likelihood", "poisson", "histogram",
         "histogram_density", "unbinned", "extended_unbinned", "indexed", "multi")


def reference_case(kind):
    """Small analytic problems isolate restart overhead, not large-data throughput.

    Each returned fit call constructs a new Python model, so its second call
    checks that ordinary new callables reuse the compiled numerical bridge.
    """
    import numpy as np
    import scientificfitting as sf

    x = np.arange(6.)
    y = np.array([0.8, 2.9, 4.7, 6.8, 8.7, 10.9])
    if kind in ("gaussian", "inplace", "indexed", "multi"):
        design = np.column_stack([x, np.ones_like(x)])
        expected = np.linalg.lstsq(design, y, rcond=None)[0]
        covariance = np.linalg.inv(design.T @ design) * 0.2**2
        cost = np.sum(((y-design@expected)/0.2)**2)

        def fit():
            options = dict(p0={"slope": 1., "offset": 0.}, sigma_y=0.2)
            def model(x, slope, offset):
                return slope*x + offset
            if kind == "inplace":
                def model_inplace(out, x, slope, offset):
                    out[:] = slope*x + offset
                return sf.fit_model(model_inplace, x, y, inplace=True, **options)
            if kind == "indexed":
                coordinates = dict(zip("abcdef", x))
                return sf.fit_indexed_model(
                    lambda labels, **p: model(np.array([coordinates[i] for i in labels]), **p),
                    list(coordinates), y, **options)
            if kind == "multi":
                options["sigma_y"] = [0.2, 0.2]
                return sf.fit_multi_model([model, model], [x[:3], x[3:]], [y[:3], y[3:]], **options)
            return sf.fit_model(model, x, y, **options)

    elif kind in ("custom", "likelihood"):
        y = np.array([0.9, 1.1, 0.8, 1.2, 1.3, 0.7])
        expected, covariance = [y.mean()], [[0.2**2/len(y)]]
        cost = np.sum(((y-y.mean())/0.2)**2 + np.log(2*np.pi*0.2**2))

        def fit():
            def logprob(y, prediction, mean):
                return -0.5*(((y-prediction)/0.2)**2 + np.log(2*np.pi*0.2**2))
            if kind == "custom":
                return sf.fit_custom(lambda mean: -2*logprob(y, mean, mean).sum(),
                                     p0={"mean": 0.5}, nobs=len(y))
            return sf.fit_likelihood_model(lambda x, mean: np.full_like(x, mean), x, y,
                                           logprob=logprob, p0={"mean": 0.5})

    elif kind == "poisson":
        counts = np.array([2., 4., 5., 3., 6.])
        mean = counts.mean()
        expected, covariance = [mean], [[mean/len(counts)]]
        cost = 2*sum(mean-n*np.log(mean)+lgamma(n+1) for n in counts)

        def fit():
            return sf.fit_poisson_model(lambda x, rate: np.full_like(x, rate),
                                        np.arange(len(counts), dtype=float), counts,
                                        p0={"rate": 2.}, bounds={"rate": (0.1, 20.)})

    elif kind in ("histogram", "histogram_density"):
        # At tau=1 the integrated expectations are exactly [6, 3].
        edges, counts = np.log([1., 2., 4.]), np.array([6., 3.])
        expected, covariance = [1.], [[1/(6*np.log(2)**2)]]
        cost = 2*sum(n-n*np.log(n)+lgamma(n+1) for n in counts)

        def fit():
            options = dict(p0={"tau": 0.7}, bounds={"tau": (0.3, 4.)})
            if kind == "histogram":
                return sf.fit_histogram_model(lambda e, tau: -12*np.diff(np.exp(-e/tau)),
                                              edges, counts, **options)
            return sf.fit_histogram_density(lambda x, tau: np.exp(-x/tau)/tau,
                                            edges, counts, total_count=12,
                                            vectorized=True, rtol=1e-8, **options)

    else:
        events = np.array([0.1, 0.4, 1.2, 0.7, 2.5, 3.0])
        n = len(events)
        if kind == "unbinned":
            tau = events.mean()
            expected, covariance, cost = [tau], [[tau**2/n]], 2*n*(np.log(tau)+1)
        else:
            expected, covariance, cost = [n/4], [[n/16]], 2*n*(1-np.log(n/4))

        def fit():
            if kind == "unbinned":
                return sf.fit_unbinned_model(lambda x, tau: np.exp(-x/tau)/tau, events,
                                             p0={"tau": 1.}, bounds={"tau": (0.05, 8.)}, vectorized=True)
            return sf.fit_extended_unbinned_model(lambda x, rate: np.full_like(x, rate),
                                                  events, (0., 4.), p0={"rate": 1.},
                                                  bounds={"rate": (0.1, 20.)}, vectorized=True)
    return fit, expected, covariance, cost


def measure(kind="gaussian"):
    seconds = {}
    start = perf_counter()
    import numpy as np
    import scientificfitting as sf
    seconds["python_import"] = perf_counter() - start

    start = perf_counter()
    import juliacall
    seconds["juliacall_import"] = perf_counter() - start
    from scientificfitting._runtime import _backend
    start = perf_counter()
    bridge = _backend()
    seconds["bridge_load"] = perf_counter() - start

    fit, expected, covariance, cost = reference_case(kind)
    start = perf_counter()
    result = fit()
    seconds["first_fit"] = perf_counter() - start
    start = perf_counter()
    result.report()
    result.diagnose()
    seconds["report_diagnose"] = perf_counter() - start

    # A new Python callable must reuse compiled bridge code, not recompile a solver.
    start = perf_counter()
    repeated = fit()
    seconds["new_callback_fit"] = perf_counter() - start
    for checked in (result, repeated):
        assert checked.converged, (kind, checked.message)
        np.testing.assert_allclose(checked.params, expected, atol=1e-6)
        np.testing.assert_allclose(checked.covariance, covariance, rtol=3e-4, atol=1e-7)
        np.testing.assert_allclose(checked.statistics["cost_min"], cost, atol=1e-7)
    assert not bridge.seval('any(m -> nameof(m) in (:Makie, :CairoMakie), values(Base.loaded_modules))')
    return dict(case=kind, platform=platform.platform(), python=platform.python_version(),
                julia=bridge.seval("string(VERSION)"), numpy=np.__version__,
                wrapper=version("scientificfitting"), juliacall=version("juliacall"),
                core=bridge.seval("string(pkgversion(ScientificFitting))"),
                core_source=bridge.seval("pathof(ScientificFitting)"),
                parameters=result.params.tolist(), covariance=result.covariance.tolist(),
                cost=result.statistics["cost_min"], seconds=seconds)


def measure_all():
    """Use separate processes so one fit family cannot warm another family's code."""
    records = []
    with tempfile.TemporaryDirectory() as directory:
        for kind in CASES:
            output = Path(directory) / f"{kind}.json"
            subprocess.run([sys.executable, str(Path(__file__).resolve()), "--case", kind,
                            "--output", str(output)], check=True, stdout=sys.stderr)
            records.append(json.loads(output.read_text(encoding="utf-8")))
    return records


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--case", choices=(*CASES, "all"), default="gaussian")
    args = parser.parse_args()
    report = json.dumps(measure_all() if args.case == "all" else measure(args.case), indent=2)
    print(report)
    if args.output:
        args.output.write_text(report + "\n", encoding="utf-8")
