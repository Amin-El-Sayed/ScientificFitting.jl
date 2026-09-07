"""Matplotlib rendering only: prediction/error propagation belongs to the Julia core."""

import numpy as np

from ._runtime import _backend
from ._results import FitReport


BLUE = "#0072B2"
INK = "#252a30"


def plot_style(name="sans"):
    """Return ordinary rcParams for use with ``plt.rc_context(plot_style(...))``.

    ``sans`` uses bundled DejaVu Sans; ``tex`` uses bundled STIX fonts and
    mathtext, not an external TeX installation. Neither changes panel visibility
    or global settings. Override any returned setting using normal Matplotlib.
    """
    if name not in ("sans", "tex"):
        raise ValueError('style must be "sans" or "tex"')
    return {"font.family": "DejaVu Sans" if name == "sans" else "STIXGeneral",
            "mathtext.fontset": "dejavusans" if name == "sans" else "stix",
            "text.usetex": False, "font.size": 11, "axes.labelsize": 12,
            "axes.titlesize": 14, "xtick.labelsize": 11, "ytick.labelsize": 11,
            "legend.fontsize": 11, "axes.linewidth": 0.9,
            "axes.edgecolor": INK, "text.color": INK,
            "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
            "axes.spines.top": name == "tex", "axes.spines.right": name == "tex"}


def _legend(fig, handles, labels, *, position="bottom", expand=False, **kwargs):
    """Let constrained layout reserve a measured, native outside legend."""
    from matplotlib.layout_engine import ConstrainedLayoutEngine

    if not isinstance(fig.get_layout_engine(), ConstrainedLayoutEngine):
        raise ValueError('outside legends require a figure with layout="constrained"')
    if position not in ("right", "bottom"):
        raise ValueError('legend position must be "right" or "bottom"')
    legend = fig.legend(handles, labels, **{
        "loc": "outside right upper" if position == "right" else "outside lower left",
        "frameon": False, "alignment": "left", **kwargs})
    if expand:
        # Only auto-sized, newly created figures grow. No draw callbacks,
        # hand-set axes positions, or changes to a caller's existing layout.
        box = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width, height = fig.get_size_inches()
        pad = 0.15
        if position == "right":
            fig.set_size_inches(width + box.width + pad, max(height, box.height + 2*pad))
        else:
            fig.set_size_inches(max(width, box.width + 2*pad), height + box.height + pad)
    return legend


def _number(value, sigdigits):
    if not np.isfinite(value):
        return r"\mathrm{unavailable}"
    parts = f"{value:.{sigdigits}g}".split("e")
    return parts[0] if len(parts) == 1 else rf"{parts[0]}\times10^{{{int(parts[1])}}}"


def add_report(fig, result, *, ax=None, parameters=None, parameter_labels=None,
               statistics=("chi2_ndf", "pvalue"), statistic_labels=None, position="right", sigdigits=4,
               expand=False, **legend_kwargs):
    """Add a native Legend containing fit estimates; return the editable Legend.

    ``result`` is a Result or a completed FitReport, including a report with
    profile errors. This never fits or scans. ``ax`` optionally contributes its
    data/curve legend entries above the report. Names and statistic keys select
    fields; labels may use Matplotlib mathtext. Missing values remain explicitly
    unavailable. The full numerical/text report remains on ``result.report()``.
    ``statistic_labels`` overrides displayed field labels without changing their
    values, e.g. label chi2_ndf as deviance/ndf for a Poisson fit.

    Use a figure created with ``layout="constrained"``. ``expand=True`` adds
    room based on the legend's measured extent; otherwise the figure size is
    unchanged. Subsequent resize/save/customization are ordinary Matplotlib.
    Panel visibility and plot_style are independent.
    """
    from matplotlib.lines import Line2D

    if not isinstance(sigdigits, int) or isinstance(sigdigits, bool) or sigdigits < 1:
        raise ValueError("sigdigits must be a positive integer")
    report = result if isinstance(result, FitReport) else result.report(structured=True)
    names = list(report.parameters if parameters is None else parameters)
    labels = parameter_labels or {}
    handles, texts = ([], []) if ax is None else ax.get_legend_handles_labels()
    if texts:
        handles.append(Line2D([], [], visible=False))
        texts.append("")
    for name in names:
        p = report.parameters[name]
        value = _number(p.value, sigdigits)
        if p.fixed:
            error = r"\quad\mathrm{(fixed)}"
        elif p.uncertainty_minus == p.uncertainty_plus:
            error = rf"\ \pm\ {_number(p.uncertainty, sigdigits)}"
        else:
            error = rf"^{{+{_number(p.uncertainty_plus, sigdigits)}}}_{{-{_number(p.uncertainty_minus, sigdigits)}}}"
        handles.append(Line2D([], [], visible=False))
        texts.append(f"{labels.get(name, name)} = ${value}{error}$")
    stat_labels = {"chi2": r"$\chi^2$", "chi2_ndf": r"$\chi^2/\mathrm{ndf}$",
                   "pvalue": r"$p$ (goodness of fit)", "aic": "AIC", "bic": "BIC"}
    stat_labels.update(statistic_labels or {})
    for key in statistics:
        value = report.statistics[key]
        rendered = f"${_number(value, sigdigits)}$" if isinstance(value, (float, int)) else str(value)
        handles.append(Line2D([], [], visible=False))
        texts.append(f"{stat_labels.get(key, key)} = {rendered}")
    status = "converged" if report.converged else "not converged"
    return _legend(fig, handles, texts, position=position, expand=expand,
                   **{"title": f"Fit result ({status})", **legend_kwargs})


def plot_fit(result, *, ax=None, x=None, band=True, data=True,
             xlabel="x", ylabel="y", title=None,
             panel=True, panel_kwargs=None, figsize=None,
             curve_kwargs=None, point_kwargs=None, band_kwargs=None):
    """Draw an existing Gaussian fit and return native `(Figure, Axes)` objects.

    `ax` embeds the plot in an existing subplot. Use `panel=False` to leave the
    caller's layout alone; otherwise the figure needs constrained layout for
    an outside report. `x` selects the finite prediction grid; the default spans
    the measured range (no implicit extrapolation). `band` is the local 1-sigma
    uncertainty of the model mean, not future measurement noise. Disabling it
    also avoids evaluating a model Jacobian. `data=False` hides observations.

    `curve_kwargs`, `point_kwargs`, `band_kwargs` are normal Matplotlib keyword
    dictionaries for plot, errorbar, and fill_between. Further artists, legends,
    mathtext, styles, and savefig remain standard Matplotlib operations. Panel
    visibility is independent of style; panel_kwargs go to add_report. Without
    an explicit figsize, a new figure grows by the measured panel size rather
    than squeezing its data axes. Explicit figsize uses inches and stays exact.
    """
    if result._kind != "gaussian":
        raise TypeError("plot_fit currently supports Gaussian x-y fit results")
    own_figure = ax is None
    if own_figure:
        import matplotlib.pyplot as plt

        _, ax = plt.subplots(figsize=figsize, layout="constrained")
    elif figsize is not None:
        raise ValueError("figsize belongs to a new figure; resize ax.figure directly")
    observed_x, observed_y = result.x, result.y
    grid = np.linspace(observed_x.min(), observed_x.max(), 400) if x is None else np.asarray(x)
    if band:
        mean, sigma = result.predict(grid, uncertainty=True)
    else:
        mean = result.predict(grid)
    curve, = ax.plot(grid, mean, **{"label": "fit", "color": BLUE, **(curve_kwargs or {})})
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
    if panel:
        add_report(ax.figure, result, ax=ax,
                   **{"expand": own_figure and figsize is None, **(panel_kwargs or {})})
    return ax.figure, ax
