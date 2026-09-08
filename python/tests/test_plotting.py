"""Numerical and layout contracts for native Matplotlib, including failure cases."""

from dataclasses import replace

import numpy as np
import pytest

from scientificfitting import (add_report, fit_custom, fit_model, plot_contour, plot_diagnostics,
                              plot_fit, plot_profile, plot_profile_matrix, plot_residuals, plot_style)
from scientificfitting._runtime import _backend


@pytest.fixture
def plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    yield plt
    plt.close("all")


@pytest.fixture(scope="module")
def line_fit():
    x = np.array([0., 0.3, 0.7, 1.1, 1.6, 2., 2.5, 3.1, 3.8, 4.6, 5.5, 6.6])
    y = np.array([0.18, 0.36, 0.89, 1.20, 1.62, 2.19, 2.54, 3.31, 3.85, 4.82, 5.50, 6.72])
    return fit_model(lambda x, slope, offset: slope*x+offset, x, y,
                     p0={"slope": 1., "offset": 0.}, sigma_y=0.12)


@pytest.fixture(scope="module")
def matrix():
    center = np.array([1., -1.2, 0.7])
    covariance = np.array([[0.09, 0.04, -0.015], [0.04, 0.16, 0.03], [-0.015, 0.03, 0.25]])
    precision = np.linalg.inv(covariance)
    calls = []

    def objective(a, b, c):
        calls.append(1)
        delta = np.array([a, b, c])-center
        return delta @ precision @ delta

    fit = fit_custom(objective, p0={"a": 0., "b": 0., "c": 0.}, nobs=20)
    return fit.profile_matrix(["c", "a", "b"], npoints_profile=9, npoints_contour=9, nsigma=3), calls


def assert_inside_canvas(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for artist in [*fig.axes, *fig.legends]:
        box = artist.get_tightbbox(renderer)
        if box is not None:
            assert box.x0 >= -1 and box.y0 >= -1, box
            assert box.x1 <= fig.bbox.width+1 and box.y1 <= fig.bbox.height+1, box


@pytest.mark.parametrize("style", ["sans", "tex"])
@pytest.mark.parametrize("panel", [True, False])
def test_styles_and_panels_are_orthogonal(plt, line_fit, style, panel, tmp_path):
    original = dict(plt.rcParams)
    with plt.rc_context(plot_style(style)):
        fig, ax = plot_fit(line_fit, panel=panel, title="Calibration", xlabel="x / mm", ylabel="U / V")
        assert len(fig.legends) == int(panel)
        assert_inside_canvas(fig)
        assert ax.get_window_extent().width >= 400
        ax.axvline(2.5, color="black", linestyle=":", label="reference")
        fig.savefig(tmp_path / f"{style}-{panel}.svg")
        fig.savefig(tmp_path / f"{style}-{panel}.pdf")
        fig.savefig(tmp_path / f"{style}-{panel}.png", dpi=120)
    assert dict(plt.rcParams) == original


def test_wider_figure_gives_space_to_the_data_not_the_report(plt, line_fit):
    widths, reports = [], []
    for width in (10, 13):
        fig, ax = plot_fit(line_fit, figsize=(width, 4.8))
        assert_inside_canvas(fig)
        np.testing.assert_allclose(fig.get_size_inches(), (width, 4.8))
        widths.append(ax.get_window_extent().width)
        reports.append(fig.legends[0].get_window_extent().width)
    assert widths[1]-widths[0] == pytest.approx(300., abs=2.)
    assert reports[0] == pytest.approx(reports[1], abs=1.)


def test_nonlinear_multiplot_and_long_report_share_native_layout(plt, tmp_path):
    t = np.array([0., 0.3, 0.7, 1.1, 1.6, 2., 2.5, 3.1, 3.8, 4.6, 5.5, 6.6])
    y = np.array([3.21, 2.61, 2.27, 1.69, 1.44, 1.05, 0.97, 0.66, 0.61, 0.35, 0.38, 0.20])
    sigma = np.array([.13, .12, .12, .10, .10, .09, .09, .08, .08, .07, .07, .07])
    fit = fit_model(lambda t, amplitude, tau, background: amplitude*np.exp(-t/tau)+background,
                    t, y, p0={"amplitude": 2.8, "tau": 1.5, "background": .1}, sigma_y=sigma,
                    bounds={"amplitude": (0, 10), "tau": (.1, 10), "background": (0, 1)})
    assert fit.converged
    # The default local range crosses B=0; automatic scans must respect that
    # bound instead of manufacturing a failed profile point outside the model.
    profile = fit.profile("background", npoints=9, nsigma=3, on_failure="throw")
    assert profile.values[0] == 0
    assert np.isfinite(profile.cost_values).all()
    contour = fit.contour("tau", "background", npoints=5, on_failure="throw")
    assert contour.y[0] == 0
    assert np.isfinite(contour.cost_values).all()
    for style in ("sans", "tex"):
        with plt.rc_context(plot_style(style)):
            fig, axes = plt.subplots(2, 1, figsize=(6.4, 6), sharex=True, layout="constrained")
            plot_fit(fit, ax=axes[0], panel=False, xlabel="", ylabel="U / V")
            axes[0].set_title("Decay with detector background")
            axes[0].axhline(fit.values["background"], color="black", label="fitted background")
            plot_residuals(fit, ax=axes[1], xlabel="t / s")
            panel = add_report(fig, fit, ax=axes[0], expand=True,
                               parameter_labels={"tau": r"detector response lifetime $\tau$ / s"},
                               statistics=("chi2", "chi2_ndf", "pvalue", "aic", "bic"))
            assert_inside_canvas(fig)
            for ax in axes:
                assert ax.get_window_extent().width > 450
                assert not ax.get_tightbbox(fig.canvas.get_renderer()).overlaps(panel.get_window_extent())
            fig.savefig(tmp_path / f"nonlinear-{style}.png", dpi=120)
            fig.savefig(tmp_path / f"nonlinear-{style}.pdf")


def test_completed_asymmetric_report_does_not_scan_or_invent_gof(plt):
    calls = []

    def cost(rate):
        calls.append(1)
        return 2*(4*rate-np.log(rate))

    result = fit_custom(cost, p0={"rate": .3}, bounds={"rate": (.001, 10)}, nobs=4)
    report = result.report(structured=True, errors="profile", profile_npoints=61, profile_nsigma=4)
    before = len(calls)
    fig, _ = plt.subplots(layout="constrained")
    panel = add_report(fig, report, position="bottom", expand=True)
    labels = [text.get_text() for text in panel.get_texts()]
    assert any("^{+" in text for text in labels)
    assert any("unavailable" in text for text in labels)
    assert len(calls) == before
    assert_inside_canvas(fig)


def test_native_artist_customization_and_explicit_existing_axes(plt, line_fit):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    before = axes[0].get_subplotspec()
    plot_fit(line_fit, ax=axes[0], panel=False, band=False,
             curve_kwargs={"color": "crimson", "linewidth": 2.5},
             point_kwargs={"markersize": 2, "ecolor": "black"})
    assert axes[0].get_subplotspec() == before and len(fig.axes) == 2
    assert axes[0].lines[0].get_color() == "crimson"
    assert axes[0].lines[0].get_linewidth() == 2.5
    panel = add_report(fig, line_fit, ax=axes[0], expand=True,
                       parameters=["slope"], parameter_labels={"slope": r"$m$"}, statistics=["aic"])
    labels = [text.get_text() for text in panel.get_texts()]
    assert any(r"$m$ =" in text for text in labels)
    assert not any("offset =" in text for text in labels)
    assert any("AIC =" in text for text in labels)
    assert_inside_canvas(fig)
    panel.remove()
    assert not fig.legends
    for layout in (None, "tight"):
        plain, _ = plt.subplots(layout=layout)
        with pytest.raises(ValueError, match="constrained"):
            add_report(plain, line_fit)


def test_profile_and_matrix_consume_completed_scans(plt, matrix, tmp_path):
    result, calls = matrix
    before = len(calls)
    for style in ("sans", "tex"):
        with plt.rc_context(plot_style(style)):
            fig, axes = plot_profile_matrix(result, parameter_labels={"c": r"$\gamma$"})
            assert axes.shape == (3, 3)
            assert axes[2, 0].get_xlabel() == r"$\gamma$"
            assert axes[2, 0].get_ylabel() == "b"
            assert axes[1, 1].get_title() == "a"
            assert axes[2, 2].get_title() == "b"
            assert_inside_canvas(fig)
            assert min(axes[i, i].get_window_extent().width for i in range(3)) > 140
            labels = [t.get_text() for t in fig.legends[0].get_texts()]
            assert any("2 parameters" in text for text in labels)
            assert any("1 parameter" in text for text in labels)
            assert any("parabolic" in text for text in labels)
            fig.savefig(tmp_path / f"matrix-{style}.png", dpi=120)
            fig.savefig(tmp_path / f"matrix-{style}.pdf")
    scan = result.contours["c", "a"]
    fig, ax = plot_contour(scan, region_kwargs={"cmap": "Blues_r"})
    assert_inside_canvas(fig)
    precision = np.linalg.inv(scan.local_covariance)
    for line, level in zip(ax.lines[:2], scan.levels):
        delta = line.get_xydata()-scan.best_values
        np.testing.assert_allclose(np.einsum("ni,ij,nj->n", delta, precision, delta), level, atol=1e-12)
    _, ax = plot_profile(result.profiles["c"])
    assert len(ax.lines[0].get_xdata()) == len(result.profiles["c"].values)
    assert len(ax.lines[1].get_xdata()) == 256
    assert len(calls) == before
    assert not _backend().seval('any(m -> nameof(m) in (:Makie, :CairoMakie), values(Base.loaded_modules))')


def test_failed_costs_and_invalid_covariance_are_not_invented(plt, matrix):
    result, _ = matrix
    profile = result.profiles["c"]
    costs = profile.delta_cost.copy()
    costs[2] = np.inf
    costs[4] = -0.5
    fig, ax = plot_profile(replace(profile, delta_cost=costs), delta_max=5)
    assert np.isnan(ax.lines[0].get_ydata()[2])
    assert ax.get_ylim()[0] <= -0.5
    assert any(line.get_label() == "failed refit" for line in ax.lines)
    scan = result.contours["c", "a"]
    costs = np.full_like(scan.delta_cost, np.inf)
    fig, ax = plot_contour(replace(scan, delta_cost=costs, local_covariance=np.full((2, 2), np.nan)))
    assert any("No finite" in text.get_text() for text in ax.texts)
    assert len(ax.lines) == 1
    assert_inside_canvas(fig)
    with pytest.raises(ValueError, match="levels"):
        plot_contour(scan, region_kwargs={"levels": [1, 4]})


def test_derivative_free_profiles_render_without_invented_curvature(plt, tmp_path):
    result = fit_custom(lambda a, b: 2*abs(a-0.7) + 3*abs(b-1.3),
        p0={"a": 0.2, "b": 1.}, nobs=10, optimizer="nelder_mead", tol=1e-10)
    profile = result.profile("a", values=np.linspace(-1, 2, 9))
    fig, ax = plot_profile(profile)
    assert not any("parabola" in line.get_label() for line in ax.lines)
    assert_inside_canvas(fig)
    scan = result.contour("a", "b", xvalues=np.linspace(-1, 2, 9),
                          yvalues=np.linspace(0, 3, 9), levels=[1., 3.])
    fig, ax = plot_contour(scan)
    assert len(ax.lines) == 1  # Minimum only; there is no local covariance ellipse.
    assert "local covariance unavailable" in [t.get_text() for t in fig.legends[0].get_texts()]
    assert_inside_canvas(fig)
    matrix = result.profile_matrix(npoints_profile=5, npoints_contour=5, nsigma=10,
                                   contour_levels=[1., 3.])
    fig, axes = plot_profile_matrix(matrix)
    assert "local correlation unavailable" in [t.get_text() for t in axes[0, 1].texts]
    assert_inside_canvas(fig)
    fig.savefig(tmp_path / "derivative-free-profiles.pdf")


def test_residuals_use_core_whitening_and_exclude_auxiliary_terms(plt):
    calls = []

    def line(x, slope, offset):
        calls.append(1)
        return slope*x+offset

    x = np.arange(7.)
    covariance = 0.04*0.5**np.abs(np.subtract.outer(x, x))
    y = np.array([0.6, 1.4, 2.8, 3.5, 4.7, 5.6, 6.8])
    result = fit_model(line, x, y, p0={"slope": 1., "offset": 0.}, cov_y=covariance,
                       parameter_priors={"slope": (1., 0.3)})
    before = len(calls)
    fig, ax = plot_residuals(result)
    plotted = ax.containers[0].lines[0].get_ydata()
    expected = np.linalg.solve(np.linalg.cholesky(covariance), result.residuals)
    np.testing.assert_allclose(plotted, expected, atol=1e-12)
    assert len(plotted) == len(x) and len(result.weighted_residuals) == len(x)+1
    assert len(calls) == before
    fig, axes = plot_diagnostics(result, xlabel="time / s")
    assert len(axes) == 3
    assert_inside_canvas(fig)
    np.testing.assert_allclose(axes[0].containers[0].lines[0].get_ydata(), result.residuals)
    np.testing.assert_allclose(axes[2].containers[0].lines[0].get_ydata(), y/result.model_y)


def test_zero_model_keeps_other_diagnostics_available(plt):
    result = fit_model(lambda x, level: np.full_like(x, level), np.arange(4), [0.1, -0.1, 0.1, -0.1],
                       p0={"level": 0.}, fixed_parameters={"level": 0.}, sigma_y=0.1)
    with pytest.raises(Exception, match="zero"):
        plot_residuals(result, kind="ratio")
    fig, axes = plot_diagnostics(result)
    assert len(axes) == 3
    assert any("unavailable" in text.get_text() for text in axes[2].texts)
    assert len(axes[0].containers) == 1 and len(axes[1].containers) == 1
    assert_inside_canvas(fig)
