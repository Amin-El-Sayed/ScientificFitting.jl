"""Execute the worked Python sources; compare their inference and rendered data."""

from pathlib import Path
import re
import runpy

import numpy as np
import pytest
from scipy import optimize, stats

from scientificfitting import plot_style


EXAMPLES = Path(__file__).resolve().parents[2] / "examples" / "python"


def test_documented_python_cells_execute_in_order(tmp_path, monkeypatch):
    """Run the actual page, not separately maintained copies of its snippets."""
    source = EXAMPLES.parents[1] / "docs" / "src" / "python.md"
    monkeypatch.chdir(tmp_path)
    namespace = {"__name__": "__documentation__"}
    for cell in re.findall(r"^```python\n(.*?)^```", source.read_text(), re.M | re.S):
        exec(compile(cell, str(source), "exec"), namespace)
    np.testing.assert_allclose(namespace["laplace_result"].params, [0.4], atol=1e-7)
    assert np.isnan(namespace["laplace_result"].stderr).all()
    np.testing.assert_allclose(namespace["waiting_fit"].params, [namespace["waiting_times"].mean()], atol=2e-6)


@pytest.fixture(scope="module")
def counts_example():
    source = runpy.run_path(str(EXAMPLES / "likelihood_workflows.py"))
    return source, source["fit_examples"]()


@pytest.fixture(scope="module")
def multi_example():
    source = runpy.run_path(str(EXAMPLES / "multi_dataset_calibration.py"))
    return source, source["fit_examples"]()


def test_count_examples_match_independent_likelihoods(counts_example):
    source, results = counts_example
    for result, model, x, y, scale in zip(
        results, (source["decay"], source["spectrum"]),
        (source["time_min"], source["edges"]), (source["decay_counts"], source["bin_counts"]),
        (np.array([40., .2, 3.]), np.array([210., 3., 1., 2.])),
    ):
        def objective(q):
            parameters = dict(zip(result.parameter_names, q*scale))
            return -2*stats.poisson.logpmf(y, model(x, **parameters)).sum()

        # Independent SciPy minimization in well-scaled coordinates.
        reference = optimize.minimize(objective, result.params/scale, method="BFGS",
                                      options={"gtol": 1e-7})
        assert result.converged
        assert np.linalg.norm(reference.jac, ord=np.inf) < 2e-4
        np.testing.assert_allclose(result.params, reference.x*scale, rtol=2e-4, atol=2e-5)
        np.testing.assert_allclose(result.statistics["cost_min"], reference.fun, atol=2e-7)
        mean = model(x, **result.values)
        deviance = 2*np.sum(mean-y + np.where(y > 0, y*np.log(np.maximum(y, 1)/mean), 0))
        assert result.statistics["chi2"] == pytest.approx(deviance, abs=1e-10)

    # Full observed Hessian: do not drop the model's second derivatives.
    count_fit = results[0]
    t, n = source["time_min"], source["decay_counts"]
    signal, rate, background = count_fit.params
    e = np.exp(-rate*t)
    mean = signal*e + background
    jac = np.column_stack([e, -signal*t*e, np.ones_like(t)])
    precision = jac.T @ ((n/mean**2)[:, None]*jac)
    precision[0, 1] += np.sum((1-n/mean)*(-t*e))
    precision[1, 0] = precision[0, 1]
    precision[1, 1] += np.sum((1-n/mean)*signal*t*t*e)
    np.testing.assert_allclose(count_fit.covariance, np.linalg.inv(precision), rtol=5e-4, atol=1e-7)


def test_multi_example_matches_block_linear_algebra(multi_example):
    source, results = multi_example
    for result, maps in zip(results, (source["shared_map"], source["partial_map"])):
        blocks = []
        for x, mapping in zip(source["xs"], maps):
            block = np.zeros((len(x), len(result.params)))
            block[:, result.parameter_names.index(mapping["gain"])] = x
            block[:, result.parameter_names.index(mapping["offset"])] = 1
            blocks.append(block)
        sigma, y = np.concatenate(source["sigmas"]), np.concatenate(source["ys"])
        design = np.vstack(blocks)/sigma[:, None]
        expected = np.linalg.lstsq(design, y/sigma, rcond=None)[0]
        covariance = np.linalg.inv(design.T @ design)
        residuals = y/sigma-design @ expected
        assert result.converged
        np.testing.assert_allclose(result.params, expected, atol=2e-6)
        np.testing.assert_allclose(result.covariance, covariance, rtol=2e-5, atol=2e-9)
        assert result.statistics["chi2"] == pytest.approx(residuals @ residuals, abs=2e-8)
    gap, error = source["gain_difference"](results[1])
    assert gap/error == pytest.approx(5.58, abs=.01)
    assert results[0].statistics["chi2"] - results[1].statistics["chi2"] == pytest.approx(31.1659, abs=1e-4)


@pytest.mark.parametrize("style", ["sans", "tex"])
@pytest.mark.parametrize("panel", [True, False])
def test_worked_figures_native_layout_and_actual_data(counts_example, multi_example, style, panel, tmp_path):
    import matplotlib.pyplot as plt

    counts, (decay, spectrum) = counts_example
    multi, (shared, partial) = multi_example
    with plt.rc_context(plot_style(style)):
        figures = [counts["plot_decay"](decay, panel=panel),
                   counts["plot_spectrum"](spectrum, panel=panel),
                   multi["plot_comparison"](shared, partial, panel=panel)]
        # Quantile edges must be horizontal/vertical, not diagonally interpolated.
        vertices = figures[0][1].collections[0].get_paths()[0].vertices
        changes = np.diff(vertices, axis=0)
        assert np.all((changes[:, 0] == 0) | (changes[:, 1] == 0))
        np.testing.assert_array_equal(figures[0][1].lines[-1].get_ydata(), counts["decay_counts"])
        histogram = figures[1][1].patches[-1].get_data()
        np.testing.assert_array_equal(histogram.values, counts["bin_counts"])
        np.testing.assert_array_equal(histogram.edges, counts["edges"])
        for ax, result, maps in zip(figures[2][1][1:], (shared, partial),
                                    (multi["shared_map"], multi["partial_map"])):
            for i, artist in enumerate(ax.lines[:3]):
                model = multi["line"](multi["xs"][i], **multi["channel_parameters"](result, maps[i]))
                np.testing.assert_allclose(artist.get_ydata(), (multi["ys"][i]-model)/multi["sigmas"][i])
        for i, (fig, _) in enumerate(figures):
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            assert len(fig.legends) == int(panel)
            for ax in fig.axes:
                assert ax.get_window_extent().width >= 450
                if panel:
                    assert not ax.get_tightbbox(renderer).overlaps(fig.legends[0].get_window_extent())
            for artist in [*fig.axes, *fig.legends]:
                box = artist.get_tightbbox(renderer)
                assert box.x0 >= -1 and box.y0 >= -1
                assert box.x1 <= fig.bbox.width+1 and box.y1 <= fig.bbox.height+1
            if panel and i < 2:
                labels = [t.get_text() for t in fig.legends[0].get_texts()]
                assert any(r"$D/\mathrm{ndf}$" in t for t in labels)
                assert not any(r"\chi^2" in t for t in labels)
            fig.savefig(tmp_path / f"example-{i}-{style}-{panel}.png", dpi=120)
            fig.savefig(tmp_path / f"example-{i}-{style}-{panel}.pdf")
            plt.close(fig)


def test_multi_plot_resizing_changes_graph_width_not_report(multi_example):
    import matplotlib.pyplot as plt

    source, (shared, partial) = multi_example
    widths, legend_widths = [], []
    for width in (11, 14):
        fig, axes = source["plot_comparison"](shared, partial, figsize=(width, 7.2))
        fig.canvas.draw()
        np.testing.assert_allclose(fig.get_size_inches(), (width, 7.2))
        widths.append(axes[0].get_window_extent().width)
        legend_widths.append(fig.legends[0].get_window_extent().width)
        plt.close(fig)
    assert widths[1]-widths[0] == pytest.approx(300, abs=2)
    assert legend_widths[0] == pytest.approx(legend_widths[1], abs=1)
