"""Independent GLS references and buffer contracts at the Python boundary."""

from dataclasses import replace

import numpy as np
import pytest
from scipy import sparse
from scipy.optimize import minimize
from scipy.sparse.linalg import spsolve

from scientificfitting import ErrorComponent, WhiteningOperator, fit_model
from scientificfitting._core import _backend
from scientificfitting._inputs import _covariance


def line(x, slope, offset):
    return slope * x + offset


def observations(n=45):
    x = np.linspace(-1, 2, n)
    y = line(x, 1.7, 0.3) + 0.07 * np.sin(np.arange(n))
    return x, y, np.column_stack([x, np.ones(n)])


def gls(design, y, covariance):
    # Solve C z = X instead of constructing the observation precision matrix.
    weighted_design = spsolve(covariance, design)
    parameter_covariance = np.linalg.inv(design.T @ weighted_design)
    return parameter_covariance @ weighted_design.T @ y, parameter_covariance


@pytest.mark.parametrize("format", ["csc", "csr", "coo", "dia", "lil", "dok"])
def test_sparse_formats_and_gls(format, monkeypatch):
    x, y, design = observations()
    c = sparse.diags([np.full(len(x)-1, 0.005), np.full(len(x), 0.04),
                      np.full(len(x)-1, 0.005)], [-1, 0, 1], format=format)
    expected, covariance = gls(design, y, c.tocsc())

    def no_dense(*args, **kwargs):
        raise AssertionError("sparse input must not be densified at the boundary")

    monkeypatch.setattr(type(c), "toarray", no_dense)
    result = fit_model(line, x, y, p0={"slope": 1.0, "offset": 0.0}, cov_y=c)
    np.testing.assert_allclose(result.params, expected, atol=2e-7)
    np.testing.assert_allclose(result.covariance, covariance, rtol=2e-6)
    assert _backend().seval("r -> issparse(r.problem.cov_y)")(result._handle)


def test_sparse_arrays_bounds_profiles_and_large_components():
    x, y, design = observations(20_000)
    covariance = sparse.csc_array(sparse.diags(
        [np.full(len(x)-1, 0.005), np.full(len(x), 0.04), np.full(len(x)-1, 0.005)], [-1, 0, 1]))
    expected, parameter_covariance = gls(design, y, covariance)
    # A component must stay sparse during both validation and optimizer setup.
    result = fit_model(line, x, y, p0={"slope": 1.0, "offset": 0.0},
                       error_components=ErrorComponent("readout", "y", "covariance", covariance),
                       bounds={"slope": (0, 3)})
    assert result.converged
    np.testing.assert_allclose(result.params, expected, atol=3e-6)
    np.testing.assert_allclose(result.covariance, parameter_covariance, rtol=3e-5)
    assert _backend().seval("r -> nnz(r.problem.error_components[1].values)")(result._handle) == 3*len(x)-2
    scan = result.profile("slope", values=result.params[0] + result.stderr[0]*np.array([-1, 0, 1]))
    np.testing.assert_allclose(scan.delta_cost, [1, 0, 1], atol=3e-4)
    before = result.predict(x)
    covariance.data[:] = 0
    np.testing.assert_array_equal(result.predict(x), before)
    assert _backend().seval("r -> r.problem.error_components[1].values[1, 1]")(result._handle) == 0.04


def test_sparse_canonicalization_and_invalid_indices():
    # Duplicate diagonal entries sum to the same covariance as a canonical CSC.
    c = sparse.csc_array(([0.01, 0.03, 0.04], [0, 0, 1], [0, 2, 3]), shape=(2, 2))
    payload = _covariance(c)
    np.testing.assert_allclose(payload["data"], [0.04, 0.04])
    assert len(c.data) == 3  # The caller's sparse buffers are not canonicalized in place.
    with pytest.raises(ValueError, match="real"):
        _covariance(c.astype(complex).todok())
    c.indices[0] = 2
    with pytest.raises(ValueError):
        _covariance(c)


@pytest.mark.parametrize("inplace", [False, True])
def test_matrix_free_whitening_and_jacobian_views(inplace):
    x, y, design = observations()
    sigma, rho = 0.2, 0.65
    covariance = sigma**2 * rho**np.abs(np.subtract.outer(np.arange(len(x)), np.arange(len(x))))
    expected, parameter_covariance = gls(design, y, sparse.csc_array(covariance))
    views = []

    def whiten(out, residual):
        assert isinstance(out, np.ndarray) and out.flags.writeable
        assert not residual.flags.writeable
        views.append(out.strides)
        out[0] = residual[0] / sigma
        out[1:] = (residual[1:] - rho*residual[:-1]) / (sigma*np.sqrt(1-rho*rho))

    def allocating(residual):
        out = np.empty_like(residual)
        whiten(out, residual)
        return out

    operator = WhiteningOperator(whiten if inplace else allocating,
        2*len(x)*np.log(sigma) + (len(x)-1)*np.log1p(-rho*rho),
        marginal_sigma=sigma, inplace=inplace)
    result = fit_model(line, x, y, p0={"slope": 1.0, "offset": 0.0}, whitening=operator,
                       jacobian=lambda x, slope, offset: np.column_stack([x, np.ones_like(x)]),
                       bounds={"slope": (0, 3)})
    np.testing.assert_allclose(result.params, expected, atol=2e-6)
    np.testing.assert_allclose(result.covariance, parameter_covariance, rtol=3e-5)
    residual = y-design @ expected
    chi2 = residual @ np.linalg.solve(covariance, residual)
    np.testing.assert_allclose(result.statistics["chi2"], chi2, rtol=1e-7)
    # The determinant, not just W r, contributes to the normalized likelihood.
    minus2loglik = chi2 + np.linalg.slogdet(covariance)[1] + len(x)*np.log(2*np.pi)
    np.testing.assert_allclose(result.statistics["aic"], minus2loglik+4, atol=1e-6)
    mean, uncertainty = result.predict(x, uncertainty=True)
    np.testing.assert_allclose(uncertainty**2, np.einsum("ij,jk,ik->i", design, parameter_covariance, design), rtol=3e-5)
    assert views


@pytest.mark.parametrize("analytic", [False, True])
def test_inplace_model_buffers_fixed_parameters_and_prediction(analytic):
    x, y, design = observations()
    calls = []

    def model(out, x, slope, offset):
        assert out.dtype == np.float64 and out.flags.writeable
        assert not x.flags.writeable
        np.multiply(x, slope, out=out)
        out += offset
        calls.append("model")

    def jacobian(out, x, slope, offset):
        assert out.shape == (len(x), 2) and out.flags.writeable
        out[:, 0], out[:, 1] = x, 1
        calls.append("jacobian")

    result = fit_model(model, x, y, p0={"slope": 1.0, "offset": 0.3}, sigma_y=0.2,
                       inplace=True, jacobian=jacobian if analytic else None,
                       fixed_parameters={"offset": 0.3})
    expected = x @ (y-0.3) / (x @ x)
    np.testing.assert_allclose(result.params, [expected, 0.3], atol=1e-6)
    np.testing.assert_allclose(result.covariance, [[0.04/(x @ x), 0], [0, 0]], atol=1e-8)
    mean, uncertainty = result.predict(x, uncertainty=True)
    np.testing.assert_allclose(mean, design @ result.params)
    np.testing.assert_allclose(uncertainty, np.abs(x)*result.stderr[0], atol=1e-8)
    scan = result.profile("slope", values=result.params[0] + result.stderr[0]*np.array([-1, 0, 1]))
    np.testing.assert_allclose(scan.delta_cost, [1, 0, 1], atol=1e-5)
    assert ("jacobian" in calls) == analytic


def test_error_components_parameter_dependent_normalization_and_toggle():
    x, y, design = observations()
    sources = [ErrorComponent("resolution", "y", "absolute", 0.12),
               ErrorComponent("calibration", "y", "relative", 0.03),
               ErrorComponent("gain", "y", "model_relative", 0.04),
               ErrorComponent("position", "x", "absolute", 0.02)]

    def objective(p):
        mu = design @ p
        variance = 0.12**2 + (0.03*y)**2 + (0.04*mu)**2 + (0.02*p[0])**2
        return np.sum((y-mu)**2/variance + np.log(2*np.pi*variance))

    expected = minimize(objective, [1.0, 0.0], method="BFGS", tol=1e-8)
    result = fit_model(line, x, y, p0={"slope": 1.0, "offset": 0.0}, error_components=sources)
    np.testing.assert_allclose(result.params, expected.x, atol=3e-6)
    np.testing.assert_allclose(result.statistics["cost_min"], objective(result.params), atol=1e-8)
    disabled = [sources[0], *(replace(source, active=False) for source in sources[1:])]
    simpler = fit_model(line, x, y, p0={"slope": 1.0, "offset": 0.0}, error_components=disabled)
    np.testing.assert_allclose(simpler.params, np.linalg.lstsq(design, y, rcond=None)[0], atol=1e-6)


def test_invalid_callback_contracts_fail_before_optimization():
    x, y, _ = observations()
    arguments = dict(p0={"slope": 1.0, "offset": 0.0})
    with pytest.raises(Exception, match="return None"):
        fit_model(lambda out, x, slope, offset: out.fill(0) or out, x, y, inplace=True, **arguments)
    with pytest.raises(Exception, match="finite"):
        fit_model(lambda out, x, slope, offset: out.__setitem__(0, slope), x, y, inplace=True, **arguments)
    with pytest.raises(Exception, match="one real value"):
        fit_model(line, x, y, whitening=WhiteningOperator(lambda r: r[:-1], 0), **arguments)
    with pytest.raises(Exception, match="whitening"):
        fit_model(line, x, y, sigma_y=0.1, whitening=WhiteningOperator(lambda r: r, 0), **arguments)
    with pytest.raises(ValueError, match="real"):
        fit_model(line, x.astype(complex), y, **arguments)
    with pytest.raises(ValueError, match="real"):
        fit_model(line, x, y, sigma_y=np.full(len(x), 0.1+0.01j), **arguments)
    with pytest.raises(Exception, match="must be real"):
        fit_model(lambda x, slope, offset: slope*x + offset + 1j, x, y, **arguments)
