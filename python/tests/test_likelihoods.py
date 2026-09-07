"""Native Python likelihoods against analytic references and independent SciPy fits."""

import numpy as np
import pytest
from scipy import optimize, stats

from scientificfitting import (
    fit_custom, fit_extended_unbinned_model, fit_histogram_density,
    fit_histogram_model, fit_indexed_model, fit_likelihood_model, fit_multi_model,
    fit_poisson_model, fit_unbinned_model,
)


def test_scipy_student_t_errors_and_profile():
    x = np.array([-1., -0.4, 0., 0.5, 1.2, 1.8, 2.4])
    y = np.array([-1.1, -0.23, 0.31, 0.94, 2.04, 2.91, 6.9])
    scales = np.array([0.12, 0.2, 0.14, 0.18, 0.11, 0.25, 0.2])
    calls = []

    def model(x, slope, offset):
        return slope*x + offset

    def logprob(y, mu, slope, offset):
        assert not y.flags.writeable and not mu.flags.writeable
        calls.append(len(y))
        return stats.t.logpdf(y, df=4, loc=mu, scale=scales)

    result = fit_likelihood_model(model, x, y, logprob=logprob, p0={"slope": 1., "offset": 0.})
    reference = optimize.minimize(
        lambda p: -2*stats.t.logpdf(y, df=4, loc=p[0]*x+p[1], scale=scales).sum(),
        [1., 0.], method="BFGS", tol=1e-7,
    )
    assert result.converged
    np.testing.assert_allclose(result.params, reference.x, atol=2e-6)
    np.testing.assert_allclose(result.statistics["cost_min"], reference.fun, atol=1e-8)
    assert np.isnan(result.statistics["pvalue"])
    assert set(calls) == {len(y)}  # no per-observation language crossings
    scan = result.profile("slope", values=result.params[0] + result.stderr[0]*np.array([-1., 0., 1.]))
    assert abs(scan.delta_cost[1]) < 1e-6
    assert "slope" in result.report()


def test_discrete_observation_probabilities():
    trials = np.array([8, 12, 10, 15, 9, 20])
    counts = np.array([2, 5, 4, 7, 3, 9])
    probability = counts.sum()/trials.sum()
    result = fit_likelihood_model(
        lambda x, probability: np.full_like(x, probability), np.arange(6), counts,
        logprob=lambda y, mu, probability: stats.binom.logpmf(y, trials, mu),
        p0={"probability": 0.3}, bounds={"probability": (0.01, 0.99)},
    )
    assert result.converged
    np.testing.assert_allclose(result.params, [probability], atol=2e-6)
    np.testing.assert_allclose(result.covariance, [[probability*(1-probability)/trials.sum()]], rtol=2e-5)


def test_scalar_unbinned_and_extended_event_models():
    data = np.array([0.12, 0.28, 0.51, 0.62, 0.75, 1.3, 1.8])

    def density(x, tau):
        assert isinstance(x, float)
        return np.exp(-x/tau)/tau

    result = fit_unbinned_model(density, data, p0={"tau": 0.5}, bounds={"tau": (0.01, 5)})
    assert result.converged
    np.testing.assert_allclose(result.params, [data.mean()], atol=2e-6)
    np.testing.assert_allclose(result.covariance, [[data.mean()**2/len(data)]], rtol=2e-5)
    assert np.isnan(result.statistics["chi2"])
    result = fit_extended_unbinned_model(lambda x, rate: rate, data, (0, 2),
                                        p0={"rate": 3.}, bounds={"rate": (0.01, 20)})
    assert result.converged
    np.testing.assert_allclose(result.params, [len(data)/2], atol=2e-5)
    np.testing.assert_allclose(result.covariance, [[len(data)/4]], rtol=2e-5)


def test_histogram_density_integrates_unequal_width_bins():
    edges = np.array([0., 0.3, 0.9, 2., 5.])
    counts = np.array([24, 35, 28, 12])
    result = fit_histogram_density(lambda x, tau: np.exp(-x/tau)/tau, edges, counts,
                                   p0={"tau": 1.}, total_count=100,
                                   bounds={"tau": (0.1, 5.)})
    reference = optimize.minimize_scalar(
        lambda tau: -2*stats.poisson.logpmf(counts, -100*np.diff(np.exp(-edges/tau))).sum(),
        bounds=(0.1, 5), method="bounded", options={"xatol": 1e-10},
    )
    assert result.converged
    np.testing.assert_allclose(result.params, [reference.x], atol=3e-6)
    np.testing.assert_allclose(result.statistics["cost_min"], reference.fun, atol=1e-8)


def test_vectorized_event_densities_and_profile_costs():
    data = np.linspace(0.01, 4., 2000)
    calls = 0

    def density(x, tau):
        nonlocal calls
        assert isinstance(x, np.ndarray) and x.ndim == 1 and not x.flags.writeable
        calls += 1
        return np.exp(-x/tau)/tau

    result = fit_unbinned_model(density, data, p0={"tau": 1.}, bounds={"tau": (0.1, 5.)},
                               vectorized=True)
    assert result.converged and calls < 1000
    np.testing.assert_allclose(result.params, [data.mean()], atol=2e-6)
    np.testing.assert_allclose(result.covariance, [[data.mean()**2/len(data)]], rtol=3e-5)
    values = data.mean()*np.array([0.9, 1., 1.1])
    scan = result.profile("tau", values=values)
    expected = 2*len(data)*(np.log(values/data.mean())+data.mean()/values-1)
    np.testing.assert_allclose(scan.delta_cost, expected, atol=2e-8)

    extended = fit_extended_unbinned_model(lambda x, rate: np.full_like(x, rate),
        [0.2, 0.5, 0.9, 1.2, 1.6], (0, 2), p0={"rate": 2.}, bounds={"rate": (0.1, 8.)},
        vectorized=True)
    assert extended.converged
    np.testing.assert_allclose(extended.params, [2.5], atol=2e-5)
    np.testing.assert_allclose(extended.covariance, [[1.25]], rtol=3e-5)


def test_vectorized_histogram_matches_integrated_bin_probabilities():
    edges, counts = np.array([0., 0.3, 0.9, 2., 5.]), [24, 35, 28, 12]
    sizes = []

    def density(x, tau):
        assert isinstance(x, np.ndarray) and not x.flags.writeable
        sizes.append(len(x))
        return np.exp(-x/tau)/tau

    options = dict(p0={"tau": 1.}, bounds={"tau": (0.1, 5.)})
    result = fit_histogram_density(density, edges, counts, total_count=100,
                                   vectorized=True, **options)
    reference = fit_histogram_model(lambda e, tau: -100*np.diff(np.exp(-e/tau)),
                                    edges, counts, **options)
    assert result.converged and sizes and min(sizes) > 1
    np.testing.assert_allclose(result.params, reference.params, atol=3e-6)
    np.testing.assert_allclose(result.covariance, reference.covariance, rtol=3e-5)
    np.testing.assert_allclose(result.statistics["cost_min"], reference.statistics["cost_min"], atol=1e-8)


def test_vectorized_density_rejects_malformed_output_and_option():
    for bad in (lambda x, tau: 1., lambda x, tau: np.ones(len(x)+1)):
        with pytest.raises(Exception, match="vectorized density|1-dimensional numeric array"):
            fit_unbinned_model(bad, [0.1, 0.5, 1.], p0={"tau": 1.}, vectorized=True)
    with pytest.raises(TypeError, match="vectorized must be a boolean"):
        fit_unbinned_model(lambda x, tau: np.exp(-x/tau)/tau, [0.1, 0.5],
                           p0={"tau": 1.}, vectorized="yes")


def test_indexed_data_and_correlated_parameter_constraint():
    indices = ["low", "middle", "high", "high"]
    coordinates = {"low": -1., "middle": 0., "high": 1.}
    y = np.array([-0.8, 0.3, 1.4, 1.3])
    design = np.column_stack(([coordinates[k] for k in indices], np.ones(4)))
    prior_mean = np.array([1.2, 0.1])
    prior_cov = np.array([[0.16, 0.02], [0.02, 0.09]])
    covariance = np.linalg.inv(design.T@design/0.2**2 + np.linalg.inv(prior_cov))
    expected = covariance @ (design.T@y/0.2**2 + np.linalg.solve(prior_cov, prior_mean))
    result = fit_indexed_model(
        lambda keys, slope, offset: np.array([slope*coordinates[k]+offset for k in keys]),
        indices, y, p0={"slope": 1., "offset": 0.}, sigma_y=0.2,
        parameter_constraints=[{"names": ["slope", "offset"], "mean": prior_mean,
                                "covariance": prior_cov}],
    )
    assert result.converged
    np.testing.assert_allclose(result.params, expected, atol=2e-6)
    np.testing.assert_allclose(result.covariance, covariance, rtol=2e-5)
    indices[:] = ["invalid"]  # retained indexed callback owns a copy
    assert np.isfinite(result.profile("slope", npoints=3, nsigma=1).delta_cost).all()


def test_multi_dataset_named_parameter_sharing():
    x = np.array([-1., 0., 0.5, 1.5])
    ya, yb = np.array([-1.5, 0.1, 1.1, 2.8]), np.array([-0.9, 0.8, 1.7, 3.4])
    line = lambda x, gain, offset: gain*x + offset
    # Deliberately scrambled global order; map local names rather than positional guesses.
    result = fit_multi_model([line, line], [x, x], [ya, yb], sigma_y=[0.2, 0.3],
        p0={"offset_b": 0., "shared_gain": 1., "offset_a": 0.},
        parameter_map=[{"gain": "shared_gain", "offset": "offset_a"},
                       {"gain": "shared_gain", "offset": "offset_b"}])
    design = np.vstack([np.column_stack([np.zeros(4), x, np.ones(4)]),
                        np.column_stack([np.ones(4), x, np.zeros(4)])])
    sigma = np.repeat([0.2, 0.3], 4)
    weighted = design / sigma[:, None]
    expected = np.linalg.lstsq(weighted, np.concatenate([ya, yb])/sigma, rcond=None)[0]
    assert result.converged
    np.testing.assert_allclose(result.params, expected, atol=2e-6)
    np.testing.assert_allclose(result.covariance, np.linalg.inv(weighted.T@weighted), rtol=2e-5)
    assert "shared_gain" in result.report()
    scan = result.profile("shared_gain", values=result.params[1] + result.stderr[1]*np.array([-1., 0., 1.]))
    np.testing.assert_allclose(scan.delta_cost, [1., 0., 1.], atol=3e-5)


def test_custom_goodness_statistic_and_validation():
    result = fit_custom(lambda mu: (mu-2)**2, p0={"mu": 0.}, nobs=5,
                        gof=lambda mu: (mu-2)**2 + 4)
    np.testing.assert_allclose(result.statistics["pvalue"], stats.chi2.sf(4, 4))
    with pytest.raises(Exception, match="one log probability"):
        fit_likelihood_model(lambda x, mu: np.full_like(x, mu), [1, 2], [3, 4],
                             logprob=lambda y, pred, mu: [0.], p0={"mu": 1.})
    with pytest.raises(ValueError, match="unknown parameter"):
        fit_custom(lambda mu: mu*mu, p0={"mu": 0.}, nobs=5, fixed_parameters={"typo": 0.})


def test_multi_dataset_rejects_complex_uncertainties_before_model_evaluation():
    def not_called(x, mu):
        raise AssertionError("invalid uncertainties must fail before evaluating the model")

    with pytest.raises(ValueError, match="real"):
        fit_multi_model([not_called], [[0., 1.]], [[1., 2.]], p0={"mu": 1.},
                        sigma_y=[np.array([0.1+0.2j, 0.1])])


def test_laplace_errors_use_derivative_free_optimization_without_hessian_errors():
    y = np.array([-1.2, -0.1, 0.2, 0.4, 0.8, 1.3, 5.0])

    def logprob(y, mu, location):
        assert isinstance(location, float) and mu.dtype == np.float64
        return stats.laplace.logpdf(y, loc=mu, scale=1)

    result = fit_likelihood_model(lambda x, location: np.full_like(x, location),
        np.arange(len(y)), y, logprob=logprob, p0={"location": 0.1},
        optimizer="nelder_mead", tol=1e-10)
    assert result.converged and result.iterations is None
    np.testing.assert_allclose(result.params, [np.median(y)], atol=1e-7)
    assert result.options["optimizer"] == "nelder_mead"
    assert result.options["parameter_covariance"] == "none"
    assert np.isnan(result.stderr).all() and np.isnan(result.covariance).all()
    codes = {finding.code for finding in result.diagnose(structured=True).findings}
    assert "covariance_not_computed" in codes and "invalid_local_covariance" not in codes
    values = np.array([-0.2, 0., 0.4, 0.7, 1.])
    scan = result.profile("location", values=values, on_failure="throw")
    reference = 2*np.abs(y[:, None]-values).sum(axis=0) - 2*np.abs(y-np.median(y)).sum()
    np.testing.assert_allclose(scan.delta_cost, reference, atol=1e-7)


def test_moving_support_is_not_clipped_or_differentiated():
    y = np.array([0.3, 0.5, 0.9, 1.2, 2.1])
    model = lambda x, location: np.full_like(x, location)
    logprob = lambda y, mu, location: stats.expon.logpdf(y, loc=mu, scale=1)
    result = fit_likelihood_model(model, np.arange(len(y)), y, logprob=logprob,
        p0={"location": 0.}, bounds={"location": (-1., 1.)},
        optimizer="nelder_mead", tol=1e-10)
    assert result.converged and result.params[0] <= y.min()
    np.testing.assert_allclose(result.params, [y.min()], atol=1e-7)
    scan = result.profile("location", values=[0., 0.2, 0.3, 0.4])
    np.testing.assert_allclose(scan.delta_cost[:3], [3., 1., 0.], atol=1e-6)
    assert np.isinf(scan.delta_cost[3])
    with pytest.raises(Exception, match="initial likelihood cost must be finite"):
        fit_likelihood_model(model, np.arange(len(y)), y, logprob=logprob,
                             p0={"location": 0.8}, optimizer="nelder_mead")


def test_derivative_free_nuisance_refits_keep_constraints_and_missing_errors():
    def cost(a, b, fixed):
        assert 0 <= a <= 2 and 0 <= b <= 2 and fixed == 0.4
        return 2*abs(a-0.7) + (b-a)**2 + fixed**2

    result = fit_custom(cost, p0={"a": 0.2, "b": 0.3, "fixed": 0.}, nobs=10,
        bounds={"a": (0., 2.), "b": (0., 2.)}, fixed_parameters={"fixed": 0.4},
        parameter_priors={"b": (0.7, 1.)}, optimizer="nelder_mead", tol=1e-10)
    assert result.converged
    np.testing.assert_allclose(result.params, [0.7, 0.7, 0.4], atol=2e-5)
    assert np.isnan(result.stderr[:2]).all() and result.stderr[2] == 0
    scan = result.profile("a", values=[0.5, 0.7, 0.9], on_failure="throw")
    np.testing.assert_allclose(scan.delta_cost, [0.42, 0., 0.42], atol=2e-6)
    stopped = fit_custom(lambda mu: (mu-1)**2, p0={"mu": 0.}, nobs=10,
                         optimizer="nelder_mead", maxiters=1)
    assert not stopped.converged


@pytest.mark.parametrize("optimizer", ["auto", "lbfgs", "nelder_mead"])
@pytest.mark.parametrize("covariance", ["none", "hessian"])
@pytest.mark.parametrize("center", [0., 1.5])
def test_likelihood_solver_and_covariance_are_orthogonal(optimizer, covariance, center):
    result = fit_custom(lambda mu: (mu-center)**2/0.04, p0={"mu": 0.2}, nobs=10,
        optimizer=optimizer, parameter_covariance=covariance, tol=1e-8, maxiters=200)
    assert result.converged
    np.testing.assert_allclose(result.params, [center], atol=2e-8)
    if covariance == "none":
        assert np.isnan(result.stderr).all()
    else:
        np.testing.assert_allclose(result.stderr, [0.2], rtol=2e-5)


def test_all_likelihood_wrappers_forward_solver_and_covariance_options():
    options = dict(p0={"mu": 1.}, optimizer="nelder_mead", parameter_covariance="none",
                   bounds={"mu": (0.1, 5.)}, tol=1e-8)
    counts = [2, 1, 3]
    constant = lambda x, mu: np.full_like(x, mu)
    fits = [
        fit_poisson_model(constant, [0, 1, 2], counts, **options),
        fit_histogram_model(lambda edges, mu: mu*np.diff(edges), [0, 1, 2, 3], counts, **options),
        fit_histogram_density(lambda x, mu: np.exp(-x/mu)/mu, [0, 1, 2, 3], counts,
                              total_count=6, **options),
        fit_unbinned_model(lambda x, mu: np.exp(-x/mu)/mu, [0.2, 0.5, 1.2], **options),
        fit_extended_unbinned_model(lambda x, mu: mu, [0.2, 0.5, 1.2], (0, 2), **options),
        fit_indexed_model(constant, [0, 1, 2], counts, **options),
        fit_multi_model([constant], [[0, 1, 2]], [counts], **options),
    ]
    for result in fits:
        assert result.converged
        assert result.options["optimizer"] == "nelder_mead"
        assert result.options["parameter_covariance"] == "none"
        assert np.isnan(result.stderr).all()
