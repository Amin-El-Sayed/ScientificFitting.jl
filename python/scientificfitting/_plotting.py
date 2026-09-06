"""Matplotlib rendering only: prediction/error propagation belongs to the Julia core."""

import numpy as np
import matplotlib.pyplot as plt

from ._core import _backend


def plot_fit(result, *, ax=None, x=None, band=True, data=True,
             xlabel="x", ylabel="y", title=None,
             curve_kwargs=None, point_kwargs=None, band_kwargs=None):
    """Draw an existing Gaussian fit and return native `(Figure, Axes)` objects.

    `ax` embeds the plot in an existing subplot without changing its layout or
    global rcParams. `x` selects the finite prediction grid; the default spans
    the measured range (no implicit extrapolation). `band` is the local 1-sigma
    uncertainty of the model mean, not future measurement noise. Disabling it
    also avoids evaluating a model Jacobian. `data=False` hides observations.

    `curve_kwargs`, `point_kwargs`, `band_kwargs` are normal Matplotlib keyword
    dictionaries for plot, errorbar, and fill_between. Further artists, legends,
    mathtext, styles, and savefig remain standard Matplotlib operations. This
    initial adapter deliberately does not reproduce the Julia report-panel UI;
    use `result.report()` for the complete numerical output.
    """
    if result._kind != "gaussian":
        raise TypeError("plot_fit currently supports Gaussian x-y fit results")
    if ax is None:
        _, ax = plt.subplots(layout="constrained")
    problem = result._handle.problem
    observed_x, observed_y = np.asarray(problem.x), np.asarray(problem.y)
    grid = np.linspace(observed_x.min(), observed_x.max(), 400) if x is None else np.asarray(x)
    if band:
        mean, sigma = result.predict(grid, uncertainty=True)
    else:
        mean = result.predict(grid)
    curve, = ax.plot(grid, mean, **{"label": "fit", **(curve_kwargs or {})})
    if band:
        ax.fill_between(grid, mean - sigma, mean + sigma, **{
            "color": curve.get_color(), "alpha": 0.18, "linewidth": 0,
            "label": r"local $1\sigma$ mean-fit band", **(band_kwargs or {}),
        })
    if data:
        errors = _backend().plot_errors(result._handle)
        xerr, yerr = (None if value is None else np.asarray(value) for value in errors)
        ax.errorbar(observed_x, observed_y, xerr=xerr, yerr=yerr, **{
            "fmt": "o", "color": "black", "markersize": 3,
            "elinewidth": 0.8, "capsize": 2, "capthick": 0.8,
            "label": "data", **(point_kwargs or {}),
        })
    ax.set(xlabel=xlabel, ylabel=ylabel)
    if title is not None:
        ax.set_title(title)
    return ax.figure, ax
