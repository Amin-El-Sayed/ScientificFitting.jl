"""Execute real NumPy callbacks through the public Python interface."""

import subprocess
import sys

import numpy as np
import pytest

from scientificfitting import fit_custom, fit_histogram_model, fit_model, fit_poisson_model, plot_fit
from scientificfitting._core import _backend, _options


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
    assert not _backend().seval('any(m -> nameof(m) in (:Makie, :CairoMakie), values(Base.loaded_modules))')


def test_import_and_report_need_no_plotting_runtime():
    # A fresh process makes this independent of other tests using Matplotlib.
    subprocess.run([sys.executable, "-c", """
import sys
import numpy as np
import scientificfitting as sf
assert 'juliacall' not in sys.modules and 'matplotlib' not in sys.modules
result = sf.fit_model(lambda x, mean: np.full_like(x, mean), [0, 1, 2], [1, 2, 3],
                      p0={'mean': 1.}, sigma_y=1.)
assert result.converged and 'mean' in result.report()
assert 'matplotlib' not in sys.modules
"""], check=True, timeout=120)


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
    assert result.converged, result.report()
    assert native.converged, str(native.message)
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


@pytest.mark.parametrize("option,ndim", [("jacobian", 2), ("x_derivative", 1)])
def test_derivative_callback_array_contract(option, ndim):
    from scientificfitting._core import _options

    x = np.arange(3.)
    expected = np.column_stack([x, np.ones_like(x)]) if ndim == 2 else np.full_like(x, 2.)

    def derivative(x, slope, offset):
        assert not x.flags.writeable
        assert slope == 2. and offset == 4.
        return expected

    callback = _options({option: derivative}, ["offset", "slope"], len(x))[option]
    np.testing.assert_array_equal(callback(x, [4., 2.]), expected)
    assert x.flags.writeable  # Borrowing must not change the owner's array flags.
    wrong_shape = np.ones(3) if ndim == 2 else np.ones((3, 1))
    callback = _options({option: lambda x, **p: wrong_shape}, ["slope"], len(x))[option]
    with pytest.raises(ValueError, match=f"{ndim}-dimensional numeric array"):
        callback(x, [2.])


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


@pytest.mark.parametrize("guesses", [
    {"offset": 0.3, "mu": -1.}, [{"offset": 0.3, "mu": -1.}],
    [-1., 0.3], [[-1., 0.3]], np.array([[-1., 0.3]]),
])
def test_initial_guess_forms_follow_p0_order(guesses):
    options = _options({"initial_guesses": guesses}, ["mu", "offset"], 10)
    np.testing.assert_array_equal(options["initial_guesses"], [[-1., 0.3]])


@pytest.mark.parametrize("guesses, message", [
    ({"mu": 1.}, "parameter names"),
    ({"mu": 1., "typo": 0.}, "parameter names"),
    ([1.], "parameter count"),
    (np.array([1.+2.j, 0.]), "real"),
])
def test_initial_guesses_reject_ambiguous_or_lossy_inputs(guesses, message):
    with pytest.raises(ValueError, match=message):
        _options({"initial_guesses": guesses}, ["mu", "offset"], 10)


def test_named_multistart_reaches_better_basin_and_respects_budget():
    def double_well(mu, offset):
        return (mu*mu-1)**2 + 0.2*mu + (offset-0.3)**2

    options = {"p0": {"mu": 1., "offset": 0.}, "nobs": 10,
               "initial_guesses": [{"offset": 0.2, "mu": -1.}]}
    local = fit_custom(double_well, **options, multistart=1)
    result = fit_custom(double_well, **options, multistart=2)
    # Stationary points solve a cubic; no competing minimizer needed as oracle.
    roots = np.roots([4., 0., -4., 0.2])
    expected = min(roots, key=lambda mu: double_well(mu, 0.3))
    assert local.converged and result.converged
    assert local.values["mu"] > 0
    np.testing.assert_allclose(result.params, [expected, 0.3], atol=2e-6)
    assert result.statistics["cost_min"] < local.statistics["cost_min"] - 0.3
