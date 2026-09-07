"""Compare warmed Python and Julia model callbacks under the same fit contract.

Requires the installed Python preview and SciPy for independent references.
No plotting or startup timings. Reports local observations, not speed guarantees.
"""

import argparse
from importlib.metadata import version
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln

import scientificfitting as sf
from scientificfitting._runtime import _backend


def measure(n=2000, repetitions=3):
    rng = np.random.default_rng(724)
    x = np.linspace(-2, 2, n)
    y = 1.7*x + 0.4 + rng.normal(0, 0.2, n)
    events = -1.7*np.log1p(-(np.arange(n)+0.5)/n)
    edges = np.linspace(0, 8, 41)
    counts = rng.poisson(-n*np.diff(np.exp(-edges/1.7)))
    design = np.column_stack([x, np.ones(n)])
    parameters = np.linalg.lstsq(design, y, rcond=None)[0]
    references = {"gaussian": (parameters, np.linalg.inv(design.T @ design)*0.2**2,
                               np.sum(((y-design@parameters)/0.2)**2)),
                  "unbinned": ([events.mean()], [[events.mean()**2/n]],
                               2*n*(np.log(events.mean())+1))}

    def histogram_cost(tau):
        mu = -n*np.diff(np.exp(-edges/tau))
        return 2*np.sum(mu-counts*np.log(mu)+gammaln(counts+1))

    reference = minimize_scalar(histogram_cost, bounds=(0.3, 6), method="bounded",
                                options={"xatol": 1e-12})
    assert reference.success
    tau = reference.x
    a, b = edges[:-1], edges[1:]
    ea, eb = np.exp(-a/tau), np.exp(-b/tau)
    mu = n*(ea-eb)
    first = n*(a*ea-b*eb)/tau**2
    second = n*((a*a*ea-b*b*eb)/tau**4-2*(a*ea-b*eb)/tau**3)
    curvature = 2*np.sum(second*(1-counts/mu)+counts*(first/mu)**2)
    references["histogram"] = ([tau], [[2/curvature]], reference.fun)

    bridge = _backend()
    native = bridge.seval('''
        function callback_benchmark(kind, x, y, events, edges, counts, derivatives)
            x, y, events = vector(x), vector(y), vector(events)
            edges, counts = vector(edges), vector(counts)
            mode = Symbol(derivatives)
            density(x, p) = exp(-x/p[1])/p[1]
            start = time_ns()
            result = if kind == "gaussian"
                fit_model((x,p) -> @.(p[1]*x+p[2]), x, y;
                    p0=[1., 0.], sigma_y=fill(0.2,length(x)), tol=1e-6, derivatives=mode)
            elseif kind == "unbinned"
                fit_unbinned_model(density, events; p0=[1.], bounds=([0.3], [6.]),
                    tol=1e-6, derivatives=mode)
            else
                fit_histogram_density(density, edges, counts; p0=[1.],
                    total_count=length(events), bounds=([0.3], [6.]),
                    tol=1e-6, derivatives=mode, rtol=1e-8)
            end
            return (result=result, seconds=(time_ns()-start)/1e9)
        end
        callback_benchmark
    ''')
    records = []
    for kind in references:
        calls = 0

        def density(x, tau):
            nonlocal calls
            calls += 1
            return np.exp(-x/tau)/tau

        def python_fit(vectorized=False, integrated=False):
            if kind == "gaussian":
                return sf.fit_model(lambda x, slope, offset: slope*x+offset, x, y,
                                    p0={"slope": 1., "offset": 0.}, sigma_y=0.2, tol=1e-6)
            options = dict(p0={"tau": 1.}, bounds={"tau": (0.3, 6.)}, tol=1e-6)
            if integrated:
                # The same bin expectations have an analytic integral in this example.
                return sf.fit_histogram_model(lambda e, tau: -n*np.diff(np.exp(-e/tau)),
                                              edges, counts, **options)
            options["vectorized"] = vectorized
            if kind == "unbinned":
                return sf.fit_unbinned_model(density, events, **options)
            return sf.fit_histogram_density(density, edges, counts, total_count=n, rtol=1e-8, **options)

        implementations = ("python", "julia_finite", "julia_auto") if kind == "gaussian" else (
            "python", "python_vectorized", "julia_finite", "julia_auto")
        if kind == "histogram":
            implementations += ("python_integrated",)
        for implementation in implementations:
            samples, callback_counts = [], []
            # The first run compiles this family; never include it in steady-state timing.
            for index in range(repetitions+1):
                calls = 0
                if implementation.startswith("python"):
                    start = perf_counter()
                    result = python_fit(vectorized=implementation == "python_vectorized",
                                        integrated=implementation == "python_integrated")
                    seconds = perf_counter()-start
                    covariance, cost = result.covariance, result.statistics["cost_min"]
                else:
                    timed = native(kind, x, y, events, edges, counts, implementation.removeprefix("julia_"))
                    result, seconds = timed.result, timed.seconds
                    covariance, cost = result.param_covariance, result.stats.cost_min
                assert result.converged, (kind, implementation, result.message)
                expected, expected_cov, expected_cost = references[kind]
                np.testing.assert_allclose(np.asarray(result.params), expected, atol=3e-6)
                # Compare dimensionless entries; a true zero cross-covariance has no relative error scale.
                errors = np.sqrt(np.diag(expected_cov))
                scale = np.outer(errors, errors)
                np.testing.assert_allclose(np.asarray(covariance)/scale, np.asarray(expected_cov)/scale,
                                           rtol=3e-4, atol=1e-7, err_msg=f"{kind}/{implementation}")
                np.testing.assert_allclose(cost, expected_cost, rtol=1e-10, atol=1e-7)
                if index:
                    samples.append(seconds)
                    callback_counts.append(calls)
            records.append(dict(case=kind, implementation=implementation,
                                seconds=samples, density_calls=callback_counts,
                                parameters=np.asarray(result.params).tolist(),
                                covariance=np.asarray(covariance).tolist(), cost=float(cost)))
            print(f"{kind}/{implementation}: {1000*np.median(samples):.3f} ms; "
                  f"density calls = {callback_counts}", file=sys.stderr, flush=True)
    return dict(platform=platform.platform(), python=platform.python_version(),
                numpy=np.__version__, scipy=version("scipy"), juliacall=version("juliacall"),
                julia=bridge.seval("string(VERSION)"),
                core=bridge.seval("string(pkgversion(ScientificFitting))"),
                core_source=bridge.seval("pathof(ScientificFitting)"),
                observations=n, bins=len(counts), repetitions=repetitions, records=records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observations", type=int, default=2000)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.observations < 10 or args.repetitions < 1:
        parser.error("use at least ten observations and one repetition")
    report = json.dumps(measure(args.observations, args.repetitions), indent=2)
    print(report)
    if args.output:
        args.output.write_text(report + "\n", encoding="utf-8")
