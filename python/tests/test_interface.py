"""Execute real NumPy callbacks through the public Python interface."""

import sys

import numpy as np
import pytest

from scientificfitting import fit_custom, fit_histogram_model, fit_model, fit_poisson_model, plot_fit
from scientificfitting._core import _backend


def line(x, slope, offset):
    assert isinstance(x, np.ndarray)
    assert x.dtype == np.float64
    assert isinstance(slope, float)
    return slope * x + offset


@pytest.fixture
def calibration():
    x = np.linspace(-1, 1, 15)
    y = line(x, 1.7, 0.4) + 0.03 * np.sin(np.arange(1, 16))
    return x, y


def test_bounds_predictions_profiles_and_ownership(calibration):
    x, y = calibration
    result = fit_model(line, x, y, p0={"slope": 1.2, "offset": 0.2}, sigma_y=0.12,
                       bounds={"slope": (0, 3), "offset": (-1, 1)})
    design = np.column_stack([x, np.ones_like(x)])
    covariance = np.linalg.inv(design.T @ design) * 0.12**2
    assert result.converged
    np.testing.assert_allclose(result.params, np.linalg.lstsq(design, y, rcond=None)[0], atol=2e-6)
    np.testing.assert_allclose(result.covariance, covariance, atol=2e-8)
    mean, sigma = result.predict(x, uncertainty=True)
    np.testing.assert_allclose(mean, design @ result.params)
    np.testing.assert_allclose(sigma**2, np.diag(design @ covariance @ design.T), rtol=2e-6)
    scan = result.profile("slope", values=result.params[0] + result.stderr[0] * np.array([-1, 0, 1]))
    np.testing.assert_allclose(scan.delta_cost, [1, 0, 1], atol=2e-5)
    contour = result.contour("slope", "offset", npoints=3, nsigma=1.0)
    assert contour.delta_cost.shape == (3, 3)
    assert np.isfinite(contour.delta_cost).all()
    assert abs(contour.delta_cost[1, 1]) < 1e-6
    assert "slope" in result.report()
    assert "Fit diagnostic dashboard" in result.diagnose()
    assert isinstance(result.statistics["chi2"], float)
    x[:] = 999
    np.testing.assert_allclose(result.predict(np.linspace(-1, 1, 15)), mean)
    with pytest.raises(ValueError):
        result.params[0] = 0
    assert "matplotlib.pyplot" not in sys.modules
    assert not _backend().seval('any(m -> nameof(m) in (:Makie, :CairoMakie), values(Base.loaded_modules))')


def test_numpy_nonlinear_model_and_dense_xy_covariance():
    def decay(x, amplitude, tau, background):
        return amplitude * np.exp(-x / tau) + background

    x = np.linspace(0, 8, 32)
    y = decay(x, 3, 1.7, 0.2) + 0.03 * np.sin(np.arange(32))
    cov_x = 0.02**2 * 0.4**np.abs(np.subtract.outer(np.arange(32), np.arange(32)))
    result = fit_model(decay, x, y, p0={"amplitude": 2.8, "tau": 1.5, "background": 0.1},
                       sigma_y=0.1, cov_x=cov_x, bounds=([0.1, 0.1, -1], [10, 5, 1]))
    native = _backend().seval('''
        (x, y, cx) -> fit_model((x,p) -> @.(p[1]*exp(-x/p[2])+p[3]),
            vector(x), vector(y); p0=[2.8,1.5,0.1], sigma_y=fill(0.1,length(x)),
            cov_x=matrix(cx), bounds=([0.1,0.1,-1.0],[10.0,5.0,1.0]))
    ''')(x, y, cov_x)
    assert result.converged
    np.testing.assert_allclose(result.params, np.asarray(native.params), atol=2e-5)
    np.testing.assert_allclose(result.covariance, np.asarray(native.param_covariance), rtol=4e-3)
    assert np.isfinite(result.predict(x, uncertainty=True)[1]).all()


def test_constraints_priors_and_fixed_parameters(calibration):
    result = fit_custom(lambda a, b: (a - 2)**2 + b*b, p0={"a": 0.8, "b": 0.3}, nobs=10,
                        constraints={"eq": lambda a, b: a*a + b*b - 1}, tol=1e-10)
    np.testing.assert_allclose(result.params, [1, 0], atol=3e-5)
    assert abs(np.sum(result.params**2) - 1) < 1e-6
    x, y = calibration
    result = fit_model(line, x, y, p0={"slope": 1.2, "offset": 0.4}, sigma_y=0.12,
                       parameter_priors={"slope": (1.7, 0.5)}, fixed_parameters={"offset": 0.4})
    assert result.params[1] == 0.4
    assert result.stderr[1] == 0


def test_poisson_and_histogram():
    counts = [9, 13, 10, 12, 8, 11]
    result = fit_poisson_model(lambda x, rate: np.full_like(x, rate), range(6), counts,
                               p0={"rate": 9}, bounds={"rate": (0.1, 30)})
    np.testing.assert_allclose(result.params, [10.5], atol=1e-5)
    np.testing.assert_allclose(result.covariance, [[1.75]], rtol=2e-5)
    np.testing.assert_allclose(result.profile("rate", values=[9, 10.5, 12]).delta_cost[1], 0, atol=1e-6)
    result = fit_histogram_model(lambda edges, rate: rate * np.diff(edges), np.arange(7), counts,
                                 p0={"rate": 9}, bounds={"rate": (0.1, 30)})
    np.testing.assert_allclose(result.params, [10.5], atol=1e-5)


def test_analytic_callbacks_and_native_matplotlib(calibration, tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    calls = []

    def jacobian(x, slope, offset):
        calls.append(len(x))
        return np.column_stack([x, np.ones_like(x)])

    x, y = calibration
    result = fit_model(line, x, y, p0={"slope": 1.2, "offset": 0.2}, sigma_y=0.12,
                       jacobian=jacobian)
    before = len(calls)
    result.predict(x, uncertainty=True)
    assert len(calls) == before + 1
    with plt.rc_context({"font.size": 13}):
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
        returned_fig, returned_ax = plot_fit(result, ax=axes[0], panel=False, xlabel="x / mm", ylabel="U / V",
                                              curve_kwargs={"color": "#0072B2"})
        assert returned_fig is fig and returned_ax is axes[0]
        parameters = result.params.copy()
        axes[0].axvline(0, color="black", linestyle="--", label="reference")
        axes[0].legend(fontsize=9)
        scan = result.profile("slope", npoints=9, nsigma=2)
        axes[1].plot(scan.values, scan.delta_cost)
        axes[1].set(xlabel="slope", ylabel=r"$\Delta\chi^2$")
        fig.savefig(tmp_path / "fit.svg")
        fig.savefig(tmp_path / "fit.png", dpi=120)
        np.testing.assert_array_equal(result.params, parameters)
        plt.close(fig)


def test_errors_are_not_hidden(calibration):
    x, y = calibration
    with pytest.raises(TypeError, match="unsupported"):
        fit_model(line, x, y, p0={"a": 1, "b": 0}, pretend_option=True)
    with pytest.raises(Exception, match="deliberate callback failure"):
        def broken(x, a):
            raise ValueError("deliberate callback failure")
        fit_model(broken, x, y, p0={"a": 1})


def test_named_parameters_do_not_depend_on_dictionary_order(calibration):
    x, y = calibration
    result = fit_model(line, x, y, p0={"offset": 0.2, "slope": 1.2}, sigma_y=0.12)
    expected = np.linalg.lstsq(np.column_stack([x, np.ones_like(x)]), y, rcond=None)[0]
    np.testing.assert_allclose(result.params, expected[::-1], atol=2e-6)
