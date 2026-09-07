"""Reports and scan geometry must survive the boundary without reinterpretation."""

from dataclasses import FrozenInstanceError
from pathlib import Path
import re

import numpy as np
import pytest
from scipy.optimize import brentq

from scientificfitting import (
    DiagnosticReport, FitReport, ProfileResult, ProfileInterval, ProfileMatrixResult,
    fit_custom, fit_model, fit_poisson_model,
)
from scientificfitting._core import _backend


@pytest.fixture
def calibration():
    calls = []

    def line(x, slope, offset):
        calls.append(1)
        return slope*x + offset

    x = np.linspace(-1, 2, 17)
    y = line(x, 1.7, 0.3) + 0.05*np.sin(np.arange(len(x)))
    result = fit_model(line, x, y, p0={"slope": 1., "offset": 0.}, sigma_y=0.1)
    return result, calls


def test_result_and_report_snapshots(calibration):
    result, calls = calibration
    before = len(calls)
    report = result.report(structured=True)
    assert isinstance(report, FitReport)
    assert result.report() == report.text == str(report)
    assert report.parameters["slope"].value == result.values["slope"]
    assert report.parameters["slope"].uncertainty == result.stderr[0]
    assert report.parameters["slope"].fixed is False
    assert result.backend == report.backend == "lsqfit"
    assert isinstance(result.message, str)
    assert result.iterations is None or isinstance(result.iterations, int)
    assert _backend().seval("r -> ismissing(r.iterations)")(result._handle) == (result.iterations is None)
    assert result.options["backend"] == "auto"
    assert result.statistics["cost"] == "chi2"
    np.testing.assert_allclose(result.residuals, result.y-result.model_y)
    np.testing.assert_allclose(result.weighted_residuals, result.residuals/0.1)
    np.testing.assert_allclose(result.jacobian, np.asarray(result._handle.jacobian))
    expected_nll = result.statistics["chi2"] + len(result.y)*np.log(2*np.pi*0.1**2)
    assert result.statistics["minus2loglik_min"] == pytest.approx(expected_nll)
    assert report.statistics == result.statistics
    for key in ("warnings", "active_bounds", "findings"):
        assert report.numerical_diagnostics[key] == result.numerical_diagnostics[key]
    for key in ("covariance_condition", "hessian_condition"):
        np.testing.assert_allclose(report.numerical_diagnostics[key], result.numerical_diagnostics[key])
    assert len(calls) == before  # Reports and stored fields must not reevaluate a model.
    with pytest.raises(ValueError):
        result.residuals[:] = 0
    with pytest.raises(TypeError):
        report.parameters["slope"] = None
    with pytest.raises(FrozenInstanceError):
        report.parameters["slope"].value = 0


def test_likelihood_and_fixed_parameter_metadata():
    result = fit_custom(lambda offset, gain: (gain-2)**2, p0={"offset": 1., "gain": 1.5},
                        fixed_parameters={"offset": 1.}, nobs=20)
    assert result.residuals is None and result.jacobian is None
    report = result.report(structured=True)
    assert report.parameters["offset"].fixed
    assert report.parameters["offset"].uncertainty == 0
    assert report.parameters["gain"].value == pytest.approx(2, abs=1e-6)
    assert np.isnan(report.statistics["pvalue"])  # No fabricated chi-square reference.
    with pytest.raises(ValueError, match="errors"):
        result.report(errors="guess")


@pytest.mark.parametrize("inplace", [False, True])
@pytest.mark.parametrize("analytic", [False, True])
def test_solver_limit_and_actual_status(inplace, analytic):
    def model(x, log_rate):
        return np.full(len(x), np.exp(log_rate))

    def model_inplace(out, x, log_rate):
        out[:] = np.exp(log_rate)

    def jacobian(x, log_rate):
        return np.full((len(x), 1), np.exp(log_rate))

    def jacobian_inplace(out, x, log_rate):
        out[:] = np.exp(log_rate)

    y = np.array([1.9, 2.1, 2.0, 1.8, 2.2, 2.1, 1.9, 2.0])
    result = fit_model(model_inplace if inplace else model, np.arange(8), y,
                       p0={"log_rate": -2.0}, sigma_y=1.0, maxiters=1, inplace=inplace,
                       jacobian=(jacobian_inplace if inplace else jacobian) if analytic else None)
    assert not result.converged and not result.report(structured=True).converged
    assert result.params[0] == pytest.approx(-0.7474443456, abs=1e-8)
    assert "did not converge" in result.message
    assert "optimizer_not_converged" in {f.code for f in result.diagnose(structured=True).findings}


def test_nonconverged_nuisance_fit_is_not_a_profile_minimum():
    result = fit_custom(lambda a, b: a*a+(b*b-a-1)**2,
                        p0={"a": 0., "b": 1.}, nobs=20, maxiters=1)
    assert result.converged  # The original parameters are exactly at the minimum.
    scan = result.profile("a", values=[-0.5, 0., 0.5])
    assert np.isinf(scan.delta_cost[2])
    assert np.isnan(scan.interval().lower) and np.isnan(scan.interval().upper)
    assert "profile_refit_failed" in {f.code for f in scan.diagnostics.findings}
    with pytest.raises(Exception, match="did not converge"):
        result.profile("a", values=[-0.5, 0., 0.5], on_failure="throw")


def test_diagnostics_are_exact_core_findings_and_actions():
    x = np.linspace(-2, 2, 21)
    result = fit_model(lambda x, slope, offset: slope*x+offset, x, 0.5*x*x+1,
                       p0={"slope": 0.1, "offset": 0.0}, sigma_y=0.01)
    report = result.diagnose(structured=True)
    assert isinstance(report, DiagnosticReport)
    assert report.status == "stop"
    assert report.severity_counts["critical"] > 0
    assert str(report) == result.diagnose()
    assert report.text == _backend().seval("r -> diagnose_text(r)")(result._handle)
    assert report.dashboard_text == _backend().seval("r -> diagnostic_dashboard_text(r)")(result._handle)
    assert all(f.code and f.evidence and f.recommendation for f in report.findings)
    limited = result.diagnose(structured=True, max_actions=0)
    assert limited.findings == report.findings
    assert limited.next_actions == ()
    with pytest.raises(Exception, match="non-negative"):
        result.diagnose(max_actions=-1)
    with pytest.raises(TypeError):
        report.severity_counts["critical"] = 0


def test_profile_interval_reuses_scan_and_does_not_hide_missing_crossings(calibration):
    result, calls = calibration
    center, sigma = result.params[0], result.stderr[0]
    scan = result.profile("slope", values=center+sigma*np.array([-2, -1, 0, 1, 2]))
    assert isinstance(scan, ProfileResult)
    np.testing.assert_allclose(scan.delta_cost, [4, 1, 0, 1, 4], atol=2e-5)
    np.testing.assert_allclose(scan.cost_values, scan.delta_cost+result.statistics["cost_min"])
    before = len(calls)
    interval = scan.interval()
    assert isinstance(interval, ProfileInterval) and interval.profile_result is scan
    np.testing.assert_allclose([interval.lower, interval.upper], [center-sigma, center+sigma], atol=1e-6)
    assert len(calls) == before
    assert scan.diagnostics.status == "ok"
    narrow = result.profile("slope", values=center+sigma*np.array([-0.2, 0, 0.2]))
    missing = narrow.interval()
    assert np.isnan(missing.lower) and np.isnan(missing.upper)
    assert "profile_threshold_not_bracketed" in {f.code for f in narrow.diagnostics.findings}
    report = result.report(structured=True, errors="profile", profile_npoints=5, profile_nsigma=0.2)
    assert all(np.isnan(p.uncertainty_minus) and np.isnan(p.uncertainty_plus)
               for p in report.parameters.values())
    assert "NaN" in report.text  # Requested profile errors cannot silently become local ones.


def test_poisson_asymmetric_interval_and_failed_grid_point():
    counts = np.array([0, 0, 1, 0])
    result = fit_poisson_model(lambda x, rate: np.full_like(x, rate), np.arange(4), counts,
                              p0={"rate": 0.3}, bounds={"rate": (0.001, 10)})
    interval = result.profile_interval("rate", npoints=61, nsigma=4.0)
    delta = lambda rate: 2*(4*(rate-0.25)-np.log(rate/0.25))-1
    expected = [brentq(delta, 0.001, 0.25), brentq(delta, 0.25, 2)]
    np.testing.assert_allclose([interval.lower, interval.upper], expected, atol=2e-4)
    assert interval.uncertainty_minus < interval.uncertainty_plus
    assert "profile_not_parabolic" in {f.code for f in interval.profile_result.diagnostics.findings}
    failed = result.profile("rate", values=[-0.5, 0.25, 0.5])
    assert np.isinf(failed.cost_values[0]) and np.isinf(failed.delta_cost[0])
    assert "profile_refit_failed" in {f.code for f in failed.diagnostics.findings}
    with pytest.raises(Exception):
        result.profile("rate", values=[-0.5, 0.25, 0.5], on_failure="throw")


def test_named_profile_matrix_geometry_and_order():
    center = np.array([1., -1.2, 0.7])
    covariance = np.array([[0.09, 0.04, -0.015], [0.04, 0.16, 0.03], [-0.015, 0.03, 0.25]])
    precision = np.linalg.inv(covariance)
    calls = []

    def objective(a, b, c):
        calls.append(1)
        d = np.array([a, b, c])-center
        return d @ precision @ d

    result = fit_custom(objective, p0={"a": 0., "b": 0., "c": 0.}, nobs=20)
    matrix = result.profile_matrix(["c", "a", "b"], npoints_profile=7, npoints_contour=5, nsigma=2.5)
    assert isinstance(matrix, ProfileMatrixResult)
    assert matrix.parameters == ("c", "a", "b")
    np.testing.assert_allclose(matrix.local_covariance, covariance[np.ix_([2, 0, 1], [2, 0, 1])], atol=2e-6)
    assert set(matrix.profiles) == {"a", "b", "c"}
    assert set(matrix.contours) == {("c", "a"), ("c", "b"), ("a", "b")}
    for names, scan in matrix.contours.items():
        indices = [["a", "b", "c"].index(name) for name in names]
        local_precision = np.linalg.inv(covariance[np.ix_(indices, indices)])
        u, v = np.meshgrid(scan.x-center[indices[0]], scan.y-center[indices[1]], indexing="ij")
        expected = local_precision[0, 0]*u*u + 2*local_precision[0, 1]*u*v + local_precision[1, 1]*v*v
        np.testing.assert_allclose(scan.delta_cost, expected, atol=1e-5)
        np.testing.assert_allclose(scan.best_values, center[indices], atol=2e-6)
        assert scan.parameters == names and scan.diagnostics.status == "ok"
    before = len(calls)
    rows = matrix.triage(include_ok=True)
    assert len(rows) == 6 and matrix.triage() == ()
    assert all(row.status == matrix.panel_status[row.parameters] for row in rows)
    assert len(calls) == before
    assert not _backend().seval('any(m -> nameof(m) in (:Makie, :CairoMakie), values(Base.loaded_modules))')
    with pytest.raises(TypeError):
        matrix.contours["c", "a"] = None
    with pytest.raises(Exception, match="unique"):
        result.profile_matrix(["a", "a"])


def test_python_guide_executes_in_document_order(tmp_path, monkeypatch):
    """Run the published cells themselves, not a separately maintained facsimile."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("MPLBACKEND", "Agg")
    page = Path(__file__).resolve().parents[2] / "docs" / "src" / "python.md"
    cells = re.findall(r"^```python\n(.*?)^```", page.read_text(), flags=re.M | re.S)
    assert len(cells) >= 7
    namespace = {"__name__": "__main__"}
    for i, code in enumerate(cells):
        exec(compile(code, f"{page}:cell-{i+1}", "exec"), namespace)
    assert (tmp_path / "calibration.pdf").is_file()
    assert (tmp_path / "student_t_errors.pdf").is_file()
    assert namespace["robust_result"].converged
    assert np.isnan(namespace["robust_result"].statistics["pvalue"])
    assert isinstance(namespace["matrix"], ProfileMatrixResult)
    assert namespace["interval"].lower < 0.25 < namespace["interval"].upper
    namespace["plt"].close("all")
