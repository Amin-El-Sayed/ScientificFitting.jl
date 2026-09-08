"""Native Matplotlib diagnostics consuming completed Julia result snapshots."""

import numpy as np

from ._runtime import _backend
from ._plotting import BLUE, INK, _legend


def _axis(ax, figsize):
    if ax is None:
        import matplotlib.pyplot as plt

        _, ax = plt.subplots(figsize=figsize, layout="constrained")
        return ax, True
    if figsize is not None:
        raise ValueError("resize an existing ax.figure directly instead of passing figsize")
    return ax, False


def _status(ax, diagnosis, mode):
    if mode not in ("issues", "all", "none"):
        raise ValueError('panel_status must be "issues", "all", or "none"')
    if mode == "none" or (mode == "issues" and diagnosis.status == "ok"):
        return
    ax.text(0.03, 0.97, diagnosis.status, transform=ax.transAxes, va="top",
            color="#a33a2b" if diagnosis.status == "stop" else INK,
            bbox={"facecolor": ax.get_facecolor(), "edgecolor": "none", "alpha": 0.85, "pad": 1})


def _level_label(level):
    sigma = "1" if np.isclose(level, 2.30, atol=0.015) else "2" if np.isclose(level, 6.18, atol=0.015) else None
    name = rf"${sigma}\sigma$ profile region" if sigma else "profile region"
    return rf"{name} (2 parameters, $\Delta$cost = {level:g})"


def _profile(ax, scan, *, local, delta_max, line_kwargs, local_kwargs, threshold_kwargs):
    values = np.where(np.isfinite(scan.delta_cost), scan.delta_cost, np.nan)
    ax.plot(scan.values, values, **{"color": BLUE, "label": "profile cost", **(line_kwargs or {})})
    if local and np.isfinite(scan.local_stderr) and scan.local_stderr > 0:
        # Only the analytic reference gets a dense grid; do not smooth refit data.
        grid = np.linspace(scan.values.min(), scan.values.max(), 256)
        ax.plot(grid, ((grid-scan.best_value)/scan.local_stderr)**2,
                **{"color": INK, "linestyle": "--", "label": "local covariance parabola",
                   **(local_kwargs or {})})
    ax.axhline(scan.threshold, **{
        "color": INK, "linestyle": ":", "linewidth": 1,
        "label": rf"interval threshold ($\Delta$cost = {scan.threshold:g}, 1 parameter)",
        **(threshold_kwargs or {})})
    failed = ~np.isfinite(scan.delta_cost)
    if failed.any():
        ax.plot(scan.values[failed], np.full(failed.sum(), 0.025), "x", color=INK,
                transform=ax.get_xaxis_transform(), label="failed refit")
    if delta_max is not None:
        if not np.isfinite(delta_max) or delta_max <= 0:
            raise ValueError("delta_max must be finite and positive")
        finite = values[np.isfinite(values)]
        ax.set_ylim(min(0., finite.min(initial=0.)), delta_max)
    ax.set(xlabel=scan.parameter, ylabel=r"$\Delta$cost")


def plot_profile(scan, *, ax=None, local=True, delta_max=None, panel_status="issues",
                 show_legend=True, figsize=None, title="Parameter profile",
                 line_kwargs=None, local_kwargs=None, threshold_kwargs=None, legend_kwargs=None):
    """Plot a completed ProfileResult; return (Figure, Axes), without more fits.

    The solid profile is compared with the dashed local covariance parabola.
    The dotted line is this scan's interval threshold, not a posterior density.
    Failed fits are gaps with crosses at the bottom; negative costs are not
    reset to zero. delta_max only limits the displayed y range. All *_kwargs
    dictionaries are Matplotlib styling options. Existing axes are not resized.
    """
    ax, own = _axis(ax, figsize)
    _profile(ax, scan, local=local, delta_max=delta_max, line_kwargs=line_kwargs,
             local_kwargs=local_kwargs, threshold_kwargs=threshold_kwargs)
    ax.set_title(title)
    _status(ax, scan.diagnostics, panel_status)
    if show_legend:
        handles, labels = ax.get_legend_handles_labels()
        if own:
            _legend(ax.figure, handles, labels, expand=figsize is None, **(legend_kwargs or {}))
        else:
            ax.legend(**{"frameon": False, **(legend_kwargs or {})})
    return ax.figure, ax


def _contour(ax, scan, *, local, regions, region_kwargs, contour_kwargs, local_kwargs):
    from matplotlib.lines import Line2D

    levels = np.unique(scan.levels)
    if not len(levels) or not np.all(np.isfinite(levels) & (levels > 0)):
        raise ValueError("contour levels must be finite and positive")
    for options in (region_kwargs, contour_kwargs):
        if options and "levels" in options:
            raise ValueError("contour levels belong to the computed scan, not styling kwargs")
    costs = np.ma.masked_invalid(scan.delta_cost.T)
    handles, labels = [], []
    if costs.count():
        if regions:
            # Non-overlapping fill intervals avoid repeatedly stacking alpha.
            colors = [(0., 0.447, 0.698, 0.26/(i+1)) for i in range(len(levels))]
            options = {"colors": colors, "corner_mask": False, **(region_kwargs or {})}
            if region_kwargs and "cmap" in region_kwargs:
                options.pop("colors", None)
            artist = ax.contourf(scan.x, scan.y, costs,
                                 levels=np.r_[min(0., float(costs.min())), levels], **options)
        else:
            options = {"colors": BLUE, **(contour_kwargs or {})}
            if contour_kwargs and "cmap" in contour_kwargs:
                options.pop("colors", None)
            artist = ax.contour(scan.x, scan.y, costs, levels=levels, **options)
        handles, _ = artist.legend_elements()
        labels = [_level_label(level) for level in levels]
    else:
        ax.text(0.5, 0.5, "No finite profile costs", ha="center", transform=ax.transAxes)

    if local:
        covariance = scan.local_covariance
        valid = covariance.shape == (2, 2) and np.isfinite(covariance).all()
        eigenvalues, vectors = np.linalg.eigh(covariance) if valid else (np.array([np.nan]), None)
        if valid and (eigenvalues > 0).all():
            angle = np.linspace(0, 2*np.pi, 256)
            circle = np.vstack([np.cos(angle), np.sin(angle)])
            for i, level in enumerate(levels):
                # Drawing a covariance ellipse is a coordinate transform, not a refit.
                ellipse = scan.best_values[:, None] + (vectors*np.sqrt(level*eigenvalues)) @ circle
                line, = ax.plot(*ellipse, **{"color": INK, "linestyle": "--",
                                            "linewidth": 1.2, **(local_kwargs or {})})
                if i == 0:
                    handles.append(line)
                    labels.append("local covariance (parabolic approximation)")
        else:
            handles.append(Line2D([], [], visible=False))
            labels.append("local covariance unavailable")
    minimum, = ax.plot(*scan.best_values, marker="+", color=INK, linestyle="none", markersize=8)
    handles.append(minimum)
    labels.append("fit minimum")
    ax.set(xlim=(scan.x.min(), scan.x.max()), ylim=(scan.y.min(), scan.y.max()),
           xlabel=scan.parameters[0], ylabel=scan.parameters[1])
    return handles, labels


def plot_contour(scan, *, ax=None, local=True, regions=True, panel_status="issues",
                 show_legend=True, figsize=None, title="Profile regions and local covariance",
                 region_kwargs=None, contour_kwargs=None, local_kwargs=None, legend_kwargs=None):
    """Draw completed two-parameter contours; return native (Figure, Axes).

    Filled regions use the scan's delta-cost thresholds (usually 2.30/6.18 for
    asymptotic joint 1/2-sigma regions). Dashed ellipses show the local parabolic
    approximation at those same levels. They are not posterior samples.
    regions=False draws profile contour lines instead. Failed grid points are
    masked, never replaced with large finite costs. Styling kwargs go directly
    to contourf, contour, plot, and legend; levels remain those of the scan.
    """
    ax, own = _axis(ax, figsize)
    handles, labels = _contour(ax, scan, local=local, regions=regions, region_kwargs=region_kwargs,
                               contour_kwargs=contour_kwargs, local_kwargs=local_kwargs)
    ax.set_title(title)
    _status(ax, scan.diagnostics, panel_status)
    if show_legend:
        if own:
            _legend(ax.figure, handles, labels, expand=figsize is None, **(legend_kwargs or {}))
        else:
            ax.legend(handles, labels, **{"frameon": False, **(legend_kwargs or {})})
    return ax.figure, ax


def plot_profile_matrix(matrix, *, axes=None, figsize=None, parameter_labels=None,
                        local=True, panel_status="issues", delta_max=6.5, show_legend=True,
                        line_kwargs=None, local_kwargs=None, region_kwargs=None, legend_kwargs=None):
    """Render a completed ProfileMatrixResult; return (Figure, n-by-n Axes).

    Diagonal: profile versus local parabola. Lower triangle: joint profile
    regions versus covariance ellipses. Upper triangle: local correlations.
    Input selection order is preserved. No new fits or scans are run; compute
    them explicitly with result.profile_matrix(...). Existing axes must have
    the matching square shape. Returned axes/artists stay fully editable.
    """
    names, labels = matrix.parameters, parameter_labels or {}
    n = len(names)
    own = axes is None
    if own:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(n, n, figsize=figsize or (2.7*n, 2.5*n),
                                 squeeze=False, layout="constrained")
    else:
        axes = np.asarray(axes, dtype=object)
        if axes.shape != (n, n):
            raise ValueError(f"axes must have shape {(n, n)}")
        fig = axes[0, 0].figure
        if figsize is not None or any(ax.figure is not fig for ax in axes.flat):
            raise ValueError("existing axes must share a figure; resize that figure directly")
    handles, texts = [], []
    for row, yname in enumerate(names):
        for col, xname in enumerate(names):
            ax = axes[row, col]
            if row == col:
                scan = matrix.profiles[xname]
                _profile(ax, scan, local=local, delta_max=delta_max, line_kwargs=line_kwargs,
                         local_kwargs=local_kwargs, threshold_kwargs=None)
                current_handles, current_texts = ax.get_legend_handles_labels()
            elif row > col:
                scan = matrix.contours[xname, yname]
                current_handles, current_texts = _contour(ax, scan, local=local, regions=True,
                    region_kwargs=region_kwargs, contour_kwargs=None, local_kwargs=local_kwargs)
            else:
                ax.set_axis_off()
                value = matrix.local_correlation[row, col]
                text = rf"local $\rho$ = {value:.2f}" if np.isfinite(value) else "local correlation unavailable"
                ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center", wrap=True)
                continue
            _status(ax, scan.diagnostics, panel_status)
            ax.set_xlabel(labels.get(xname, xname) if row == n-1 else "")
            ax.set_ylabel(r"$\Delta$cost" if row == col else labels.get(yname, yname))
            ax.tick_params(labelbottom=row == n-1)
            if row == col:
                ax.set_title(labels.get(xname, xname))
            for handle, text in zip(current_handles, current_texts):
                if text not in texts:
                    handles.append(handle)
                    texts.append(text)
    if show_legend:
        _legend(fig, handles, texts, expand=own and figsize is None,
                **{"ncols": 2, **(legend_kwargs or {})})
    return fig, axes


def _residual_data(result, kind):
    if result._kind != "gaussian":
        raise TypeError("x-y residual plots require a Gaussian fit; use profiles for a general likelihood")
    if kind not in ("residual", "pull", "ratio"):
        raise ValueError('kind must be "residual", "pull", or "ratio"')
    return _backend().diagnostic_data(result._handle, kind)


def _residual(ax, values, kind, point_kwargs, reference_kwargs):
    x, y, errors, title, ylabel, reference = values
    if kind == "pull":
        ax.axhspan(-2, 2, color=BLUE, alpha=0.07, linewidth=0)
        ax.axhspan(-1, 1, color=BLUE, alpha=0.10, linewidth=0)
    ax.axhline(reference, **{"color": INK, "linestyle": "--", "linewidth": 1,
                            **(reference_kwargs or {})})
    ax.errorbar(np.asarray(x), np.asarray(y), yerr=None if errors is None else np.asarray(errors),
                **{"fmt": "o", "color": INK, "markersize": 3, "elinewidth": 0.8,
                   "capsize": 2, "capthick": 0.8, **(point_kwargs or {})})
    ax.set(title=str(title), ylabel=str(ylabel))


def plot_residuals(result, *, kind="pull", ax=None, figsize=None, xlabel="x",
                   point_kwargs=None, reference_kwargs=None):
    """Plot Gaussian residuals, whitened residuals, or data/fit ratios.

    Returns (Figure, Axes). Values and error bars come from the same core helper
    as Makie. Pulls exclude auxiliary parameter terms; with correlated data
    they are whitened coordinates, not independent pointwise error-bar ratios.
    Shading at +/-1 and +/-2 is a reference guide, not fitted interval coverage.
    A ratio is undefined at zero model predictions and raises a clear error.
    """
    values = _residual_data(result, kind)
    ax, _ = _axis(ax, figsize)
    _residual(ax, values, kind, point_kwargs, reference_kwargs)
    ax.set_xlabel(xlabel)
    return ax.figure, ax


def plot_diagnostics(result, *, kinds=("residual", "pull", "ratio"), axes=None,
                     figsize=None, xlabel="x", point_kwargs=None, reference_kwargs=None):
    """Make a Gaussian diagnostic stack; return (Figure, 1-D array of Axes).

    kinds selects independent panels. Existing axes must match their number.
    At zero predictions, the ratio panel explicitly says it is unavailable;
    useful residual/pull panels still render. No refit is performed.
    """
    kinds = tuple(kinds)
    if not kinds:
        raise ValueError("kinds must contain at least one diagnostic")
    if result._kind != "gaussian":
        raise TypeError("x-y diagnostics require a Gaussian fit; use profiles for a general likelihood")
    values = [None if kind == "ratio" and np.any(result.model_y == 0) else _residual_data(result, kind)
              for kind in kinds]
    if axes is None:
        import matplotlib.pyplot as plt

        fig, grid = plt.subplots(len(kinds), 1, figsize=figsize or (6.4, 2.3*len(kinds)),
                                 sharex=True, squeeze=False, layout="constrained")
        axes = grid[:, 0]
    else:
        axes = np.asarray(axes, dtype=object).reshape(-1)
        if len(axes) != len(kinds) or figsize is not None:
            raise ValueError("provide one existing axis per kind and resize their figure directly")
        fig = axes[0].figure
        if any(ax.figure is not fig for ax in axes):
            raise ValueError("diagnostic axes must share a figure")
    for ax, kind, value in zip(axes, kinds, values):
        if value is None:
            ax.text(0.5, 0.5, "Ratio unavailable: model is zero at a measurement",
                    ha="center", va="center", transform=ax.transAxes, wrap=True)
            ax.set(ylabel="data / fit", title="Ratio")
        else:
            _residual(ax, value, kind, point_kwargs, reference_kwargs)
    axes[-1].set_xlabel(xlabel)
    return fig, axes
