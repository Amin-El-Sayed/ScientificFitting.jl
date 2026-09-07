"""Native Python likelihoods against analytic references and independent SciPy fits."""

import numpy as np
import pytest
from scipy import optimize, stats

from scientificfitting import (
    fit_custom, fit_extended_unbinned_model, fit_histogram_density,
    fit_indexed_model, fit_likelihood_model, fit_multi_model, fit_unbinned_model,
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
