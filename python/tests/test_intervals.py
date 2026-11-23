"""Interval methods reach the core unchanged and return its limits by parameter name."""

from juliacall import JuliaError
import numpy as np
import pytest
from scipy.optimize import brentq
from scipy.stats import chi2, norm

from scientificfitting import (
    BootstrapSamples, ConfidenceInterval, fit_custom, fit_model, fit_poisson_model,
)
from scientificfitting._core import _backend

X = np.arange(10.)
Y = 2*X + 1 + 0.5*np.array([0.3, -1.1, 0.4, 0.9, -0.2, 1.4, -0.7, 0.1, -1.6, 0.5])
START = {"slope": 1., "offset": 0.}


def line(x, slope, offset):
    return slope*x + offset


def least_squares(sigma):
    """Closed-form estimate and covariance; sigma=None scales by the residual variance."""
    design = np.column_stack([X, np.ones_like(X)])
    normal = design.T @ design / (1. if sigma is None else sigma**2)
    params = np.linalg.solve(normal, design.T @ Y / (1. if sigma is None else sigma**2))
    covariance = np.linalg.inv(normal)
    if sigma is None:
        residuals = Y - design @ params
        covariance *= residuals @ residuals / (len(Y) - len(params))
    return params, covariance


@pytest.fixture(scope="module")
def weighted():
    return fit_model(line, X, Y, p0=START, sigma_y=0.5)


@pytest.mark.parametrize("sigma", [0.5, None])
def test_linear_hessian_and_profile_match_closed_form(sigma):
    # Linear model with fixed covariance: the cost is exactly quadratic, so both
    # methods equal params +- z*stderr, including the chi2/ndf-scaled unweighted fit.
    # The tolerance covers the interpolation between profile scan points.
    params, covariance = least_squares(sigma)
    result = fit_model(line, X, Y, p0=START, sigma_y=sigma)
    for level in (0.6827, 0.9545):
        for i, name in enumerate(result.parameter_names):
            halfwidth = norm.ppf((1 + level)/2) * np.sqrt(covariance[i, i])
            for method in ("hessian", "profile"):
                interval = result.confidence_interval(name, method=method, level=level)
                assert (interval.parameter, interval.level, interval.method) == (name, level, method)
                np.testing.assert_allclose([interval.lower, interval.upper],
                                           [params[i] - halfwidth, params[i] + halfwidth],
                                           rtol=0, atol=1e-6*halfwidth)


def test_poisson_profile_and_hessian_match_their_closed_forms():
    counts = np.array([0, 0, 1, 0])
    result = fit_poisson_model(lambda x, rate: np.full_like(x, rate), np.arange(4.), counts,
                               p0={"rate": 0.3}, bounds={"rate": (0.001, 10)})
    level = 0.9
    threshold = chi2.ppf(level, 1)
    # -2 log L ratio of four Poisson bins with total count 1 around the MLE 1/4;
    # the tolerance covers interpolating the curved profile between scan points.
    delta = lambda rate: 2*(4*(rate - 0.25) - np.log(rate/0.25)) - threshold
    profile = result.confidence_interval("rate", level=level)
    np.testing.assert_allclose([profile.lower, profile.upper],
                               [brentq(delta, 0.001, 0.25), brentq(delta, 0.25, 5)], atol=5e-5)
    # Fisher information 4/rate at the MLE gives stderr sqrt(0.25/4) = 0.25.
    hessian = result.confidence_interval("rate", method="hessian", level=level)
    np.testing.assert_allclose([hessian.lower, hessian.upper],
                               0.25 + np.array([-1, 1])*norm.ppf((1 + level)/2)*0.25, rtol=1e-6)


def test_bootstrap_reproduces_gaussian_sampling_distribution(weighted):
    # For a linear model with Gaussian noise the refitted parameters are exactly
    # N(params, covariance); tolerances are four Monte Carlo standard errors.
    params, covariance = least_squares(0.5)
    stderr = np.sqrt(np.diag(covariance))
    nsamples = 400
    samples = weighted.parametric_bootstrap(nsamples=nsamples, seed=7)
    assert isinstance(samples, BootstrapSamples)
    assert samples.params.shape == (nsamples, 2)
    assert (samples.n_requested, samples.n_failed) == (nsamples, 0)
    np.testing.assert_array_less(np.abs(samples.params.mean(axis=0) - params),
                                 4*stderr/np.sqrt(nsamples))
    np.testing.assert_array_less(np.abs(samples.params.std(axis=0, ddof=1)/stderr - 1),
                                 4/np.sqrt(2*(nsamples - 1)))
    with pytest.raises(ValueError):
        samples.params[0, 0] = 0

    again = weighted.parametric_bootstrap(nsamples=nsamples, seed=np.int64(7))
    np.testing.assert_array_equal(again.params, samples.params)
    assert not np.array_equal(weighted.parametric_bootstrap(nsamples=nsamples, seed=8).params,
                              samples.params)
    assert weighted.parametric_bootstrap(nsamples=20).params.shape == (20, 2)  # Julia default rng

    # The same seed reproduces the same samples; Julia's default quantile is
    # NumPy's default linear interpolation.
    level = 0.9
    interval = weighted.confidence_interval("offset", method="bootstrap", level=level,
                                            nsamples=nsamples, seed=7)
    assert interval == ConfidenceInterval("offset", interval.lower, interval.upper, level, "bootstrap")
    np.testing.assert_allclose([interval.lower, interval.upper],
                               np.quantile(samples.params[:, 1], [(1 - level)/2, (1 + level)/2]),
                               rtol=1e-12)


def _is_julia(error, kind):
    return _backend().seval(f"e -> e isa {kind}")(error.value.exception)


@pytest.mark.parametrize("kind, call", [
    ("DomainError", lambda r: r.confidence_interval("slope", level=1.0)),
    ("DomainError", lambda r: r.confidence_interval("slope", level=0.0)),
    ("DomainError", lambda r: r.confidence_interval("slope", method="bootstrap", nsamples=0)),
    ("DomainError", lambda r: r.parametric_bootstrap(nsamples=0)),
    ("ArgumentError", lambda r: r.confidence_interval("slope", method="median")),
    # A keyword the selected method does not use is an error, never ignored.
    ("ArgumentError", lambda r: r.confidence_interval("slope", nsamples=100)),
    ("ArgumentError", lambda r: r.confidence_interval("slope", method="hessian", seed=1)),
    ("ArgumentError", lambda r: r.confidence_interval("slope", method="bootstrap", npoints=11)),
])
def test_core_errors_keep_their_julia_types(weighted, kind, call):
    with pytest.raises(JuliaError) as error:
        call(weighted)
    assert _is_julia(error, kind)


def test_unavailable_intervals_raise():
    fixed = fit_model(line, X, Y, p0=START, sigma_y=0.5, fixed_parameters={"offset": 1.})
    assert fixed.confidence_interval("slope", method="hessian").method == "hessian"
    with pytest.raises(JuliaError) as error:
        fixed.confidence_interval("offset", method="hessian")
    assert _is_julia(error, "ArgumentError")

    custom = fit_custom(lambda mu: (mu - 2)**2, p0={"mu": 0.}, nobs=10, parameter_covariance="none")
    with pytest.raises(JuliaError) as error:
        custom.confidence_interval("mu", method="hessian")
    assert _is_julia(error, "ArgumentError")
    with pytest.raises(JuliaError) as error:
        custom.parametric_bootstrap(nsamples=10, seed=1)
    assert _is_julia(error, "ArgumentError")

    with pytest.raises(ValueError):
        fixed.confidence_interval("gain")
    with pytest.raises(TypeError):
        fixed.parametric_bootstrap(nsamples=10, seed=1.5)
