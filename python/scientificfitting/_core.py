"""Python conventions and ownership at the boundary; all inference stays in Julia."""

from collections.abc import Mapping
from functools import cache
from pathlib import Path
from types import MappingProxyType, SimpleNamespace

import numpy as np


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


def _array(values, *, ndim=1):
    array = np.array(values, dtype=np.float64, copy=True)
    if array.ndim != ndim:
        raise ValueError(f"expected a {ndim}-dimensional numeric array")
    return array


def _snapshot(values):
    array = np.array(values, dtype=np.float64, copy=True)
    array.flags.writeable = False
    return array


def _model_callback(function, names):
    def call(x, parameters):
        # Julia owns inputs during the callback; user code must not mutate them.
        values = np.asarray(x).view()
        values.flags.writeable = False
        return np.asarray(function(values, **dict(zip(names, map(float, parameters)))), dtype=np.float64)

    return call


def _parameter_callback(function, names, *, scalar=False):
    def call(parameters):
        value = function(**dict(zip(names, map(float, parameters))))
        return float(value) if scalar else np.atleast_1d(np.asarray(value, dtype=np.float64))

    return call


_OPTIONS = {
    "sigma_x", "sigma_y", "cov_x", "cov_y", "bounds", "constraints",
    "parameter_priors", "fixed_parameters", "jacobian", "x_derivative",
    "backend", "cost", "scale_covariance", "maxiters", "tol",
    "initial_guesses", "multistart", "nobs", "cost_name",
}


def _options(options, names, nobs):
    unknown = options.keys() - _OPTIONS
    if unknown:
        raise TypeError(f"unsupported fit options: {', '.join(sorted(unknown))}")
    converted = dict(options)
    for key in ("sigma_x", "sigma_y"):
        if converted.get(key) is not None:
            value = np.asarray(converted[key], dtype=float)
            converted[key] = np.full(nobs, float(value)) if value.ndim == 0 else _array(value)
    for key in ("cov_x", "cov_y"):
        if converted.get(key) is not None:
            converted[key] = _array(converted[key], ndim=2)
    if converted.get("bounds") is not None:
        # Named bounds avoid one-based/zero-based index translation in user code.
        bounds = converted["bounds"]
        if isinstance(bounds, Mapping):
            if bounds.keys() - set(names):
                raise ValueError("bounds contain unknown parameter names")
            bounds = list(zip(*(bounds.get(name, (-np.inf, np.inf)) for name in names)))
        converted["bounds"] = tuple(_array(bound) for bound in bounds)
    for key in ("fixed_parameters", "parameter_priors"):
        if converted.get(key) is not None:
            values = converted[key]
            if not isinstance(values, Mapping):
                raise TypeError(f"{key} must map parameter names to values")
            converted[key] = [
                [names.index(name) + 1, *np.atleast_1d(value).tolist()]
                for name, value in values.items()
            ]
    for key in ("jacobian", "x_derivative"):
        if converted.get(key) is not None:
            converted[key] = _model_callback(converted[key], names)
    if converted.get("constraints") is not None:
        constraints = converted["constraints"]
        if constraints.keys() - {"eq", "ineq"}:
            raise ValueError("constraints use eq == 0 and ineq <= 0")
        converted["constraints"] = {name: _parameter_callback(f, names) for name, f in constraints.items()}
    return converted


class Result:
    """Fit snapshot with NumPy values and Julia-backed predictions and profile refits.

    `params`, `stderr`, and `covariance` follow `parameter_names` order and are
    read-only snapshots. `statistics` is read-only. `report()` and `diagnose()`
    return the actual core output; plotting never invokes another fit.
    """

    def __init__(self, handle, names, kind):
        self._handle = handle
        self._kind = kind
        self.parameter_names = tuple(names)
        fields = _backend().result_values(handle)
        for name in ("params", "stderr", "covariance", "correlation"):
            setattr(self, name, _snapshot(fields[name]))
        self.converged = bool(fields["converged"])
        self.statistics = MappingProxyType(dict(fields["statistics"]))

    @property
    def values(self):
        """Named best-fit parameters; mutations cannot affect the retained fit."""
        return dict(zip(self.parameter_names, self.params))

    def report(self):
        """Core parameter/statistics report, not a reconstructed Python summary."""
        return str(_backend().result_report(self._handle, list(self.parameter_names)))

    def diagnose(self):
        """Core diagnostic dashboard with concrete follow-up actions."""
        return str(_backend().result_diagnose(self._handle))

    def predict(self, x, *, uncertainty=False):
        """Gaussian model mean, optionally `(mean, sigma)` for local mean-fit uncertainty.

        The band excludes new-observation noise and is not a prediction interval.
        Analytic model Jacobians are used when supplied; otherwise central
        differences are evaluated in the Julia core.
        """
        if self._kind != "gaussian":
            raise TypeError("predict currently supports Gaussian x-y fit results")
        prediction = _backend().prediction(self._handle, _array(x), uncertainty)
        if uncertainty:
            return _snapshot(prediction.mean), _snapshot(prediction.sigma)
        return _snapshot(prediction)

    def profile(self, parameter, **options):
        """Re-optimize nuisance parameters at each point; parameter is a name.

        Returns NumPy `values`, `delta_cost`, and the requested `threshold`.
        Options follow Julia `profile`, e.g. `values`, `npoints`, `adaptive`.
        """
        index = self.parameter_names.index(parameter) + 1
        result = _backend().run_profile(self._handle, index, options)
        return SimpleNamespace(values=_snapshot(result.values), delta_cost=_snapshot(result.delta_cost),
                               threshold=float(result.threshold))

    def contour(self, first, second, **options):
        """Two-parameter profile grid; `delta_cost[i,j]` belongs to `x[i], y[j]`.

        Matplotlib `contourf(x, y, delta_cost.T, ...)` expects the transpose.
        Failed refits remain infinite rather than being interpolated away.
        """
        indices = [self.parameter_names.index(name) + 1 for name in (first, second)]
        result = _backend().run_contour(self._handle, *indices, options)
        return SimpleNamespace(x=_snapshot(result.x_values), y=_snapshot(result.y_values),
                               delta_cost=_snapshot(result.delta_cost), levels=_snapshot(result.levels))


def _fit(kind, model, x, y, p0, options):
    if not isinstance(p0, Mapping) or not p0 or not all(isinstance(name, str) for name in p0):
        raise TypeError("p0 must be a nonempty mapping of parameter names to initial values")
    names = list(p0)
    # Mapping order defines result-array order; callbacks bind parameters by name.
    start = _array(list(p0.values()))
    x, y = _array(x), _array(y)
    keywords = _options(options, names, len(y))
    callback = _parameter_callback(model, names, scalar=True) if kind == "custom" else _model_callback(model, names)
    handle = _backend().run_fit(kind, callback, x, y, start, keywords)
    return Result(handle, names, kind)


def fit_model(model, x, y, *, p0, **options):
    """Fit a vectorized `model(x, *parameters)` to Gaussian x-y data.

    `p0={"slope": 1.0, "offset": 0.0}` sets result-array order and names.
    Callbacks receive parameters by keyword, so dictionary order cannot swap
    the meaning of model arguments. Names must match the Python signature.
    Supports scalar/pointwise sigma_x/sigma_y, dense cov_x/cov_y, bounds,
    named fixed_parameters/parameter_priors, nonlinear constraints, analytic
    jacobian/x_derivative, and Julia solver controls. All derivative paths use
    finite differences unless an analytic callback is available. Callbacks
    must be smooth and defined in a neighborhood of evaluated points.
    The core's finite-mode stopping tolerance defaults to 1e-6; explicit `tol`
    values are preserved. This is a solver criterion, not a parameter error.
    """
    return _fit("gaussian", model, x, y, p0, options)


def fit_poisson_model(model, x, counts, *, p0, **options):
    """Fit strictly positive expected counts to nonnegative integer observations."""
    return _fit("poisson", model, x, counts, p0, options)


def fit_histogram_model(expected_counts, edges, counts, *, p0, **options):
    """Fit one expected count per bin; the model must integrate over bin edges."""
    return _fit("histogram", expected_counts, edges, counts, p0, options)


def fit_custom(objective, *, p0, nobs, **options):
    """Fit a scalar `objective(*parameters)` on the normalized -2 log L scale.

    For arbitrary losses the optimum is usable, but covariance and likelihood
    summaries have no automatic statistical interpretation. `nobs` is required.
    """
    return _fit("custom", objective, [], [], p0, {**options, "nobs": nobs})
