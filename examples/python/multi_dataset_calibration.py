"""Can three channels share one gain? Joint fits and ordinary Matplotlib.

Fixed teaching data from the Multi-Dataset gallery, with known independent
Gaussian measurement errors. No generated data and no fits inside plot code.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2
from scientificfitting import add_report, fit_multi_model, plot_style


xs = [np.array([0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.]),
      np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]),
      np.array([0., 1.25, 2.5, 3.75, 5., 6.25, 7.5, 8.75, 10.])]
ys = [np.array([0.744750, 2.444135, 4.420060, 6.098325, 8.141240, 9.763075,
                11.660295, 13.287080, 15.417610, 17.051490, 19.047875]),
      np.array([0.391920, 2.378570, 4.007500, 5.927490, 7.877840, 9.433180,
                11.431380, 13.147100, 15.144640, 16.759710]),
      np.array([0.159600, 2.417162, 5.013387, 7.207425, 9.690625,
                12.237350, 14.323312, 16.937275, 19.023650])]
sigmas = [np.array([0.075, 0.083, 0.091, 0.099, 0.107, 0.115, 0.123, 0.131, 0.139, 0.147, 0.155]),
          np.array([0.088, 0.094, 0.100, 0.106, 0.112, 0.118, 0.124, 0.130, 0.136, 0.142]),
          np.array([0.080, 0.09125, 0.1025, 0.11375, 0.125, 0.13625, 0.1475, 0.15875, 0.170])]
shared_map = [{"gain": "gain", "offset": f"offset_{channel}"} for channel in "abc"]
partial_map = [{"gain": "gain_ab" if channel in "ab" else "gain_c", "offset": f"offset_{channel}"}
               for channel in "abc"]
COLORS, MARKERS = ("#0072B2", "#c33c40", "#237a57"), ("o", "s", "D")


def line(x, gain, offset):
    """Channel response in volts at a dimensionless reference setting."""
    return gain*x + offset


def fit_examples():
    """Compare all-shared and partially shared gains on exactly the same data."""
    shared = fit_multi_model([line]*3, xs, ys, sigma_y=sigmas,
                            p0={"gain": 1.85, "offset_a": 0.7, "offset_b": -0.4, "offset_c": 0.1},
                            parameter_map=shared_map)
    partial = fit_multi_model([line]*3, xs, ys, sigma_y=sigmas,
                             p0={"gain_ab": 1.82, "offset_a": 0.7, "offset_b": -0.4,
                                 "gain_c": 1.90, "offset_c": 0.1}, parameter_map=partial_map)
    return shared, partial


def channel_parameters(result, mapping):
    """Resolve local model names from the fit's explicit sharing map."""
    return {local: result.values[global_name] for local, global_name in mapping.items()}


def gain_difference(result):
    """Propagate the full parameter covariance for g_C - g_AB."""
    gradient = np.array([{"gain_ab": -1., "gain_c": 1.}.get(name, 0.) for name in result.parameter_names])
    return gradient @ result.params, np.sqrt(gradient @ result.covariance @ gradient)


def plot_comparison(shared, partial, *, panel=True, figsize=None):
    """Return editable main/pull axes; native constrained layout owns spacing."""
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=figsize or (6.8, 7.2),
                             height_ratios=(3, 1, 1), layout="constrained")
    grid = np.linspace(0., 10., 300)
    for i, (x, y, sigma, color, marker) in enumerate(zip(xs, ys, sigmas, COLORS, MARKERS)):
        selected = [partial.parameter_names.index(name) for name in partial_map[i].values()]
        jacobian = np.column_stack([grid, np.ones_like(grid)])
        covariance = partial.covariance[np.ix_(selected, selected)]
        # The linear model's mean band uses the joint fit covariance, not the scatter of its data.
        error = np.sqrt(np.einsum("ni,ij,nj->n", jacobian, covariance, jacobian))
        mean = line(grid, **channel_parameters(partial, partial_map[i]))
        axes[0].fill_between(grid, mean-error, mean+error, color=color, alpha=0.16, linewidth=0,
                              label=r"local $1\sigma$ mean-fit band" if i == 0 else None)
        axes[0].plot(grid, mean, color=color, linewidth=1.6)
        axes[0].plot(grid, line(grid, **channel_parameters(shared, shared_map[i])),
                      color=color, linewidth=1., linestyle="--")
        axes[0].errorbar(x, y, yerr=sigma, fmt=marker, color=color, markersize=3,
                         ecolor="#252a30", elinewidth=0.8, capsize=2, capthick=0.8,
                         label=f"channel {'ABC'[i]}")
        for ax, result, mapping in zip(axes[1:], (shared, partial), (shared_map, partial_map)):
            pull = (y-line(x, **channel_parameters(result, mapping[i])))/sigma
            ax.plot(x, pull, marker=marker, color=color, markersize=3, linewidth=0.8)
    axes[0].plot([], [], color="#252a30", linewidth=1.6, label="partially shared gains")
    axes[0].plot([], [], color="#252a30", linewidth=1., linestyle="--", label="all gains shared")
    axes[0].set(title="Three-channel calibration transfer", ylabel="response / V")
    for ax, title in zip(axes[1:], ("Pulls: all gains shared", "Pulls: partially shared gains")):
        ax.axhspan(-1, 1, color="#0072B2", alpha=0.08, linewidth=0)
        ax.axhline(0, color="#252a30", linewidth=0.8)
        ax.set(ylabel=r"$r_i$", ylim=(-3.5, 3.5))
        ax.set_title(title, fontsize=11, loc="left")
    axes[-1].set_xlabel("reference setting x")
    if panel:
        # Select scientific essentials for the figure; the complete report is still available.
        add_report(fig, partial, ax=axes[0], expand=figsize is None,
                   parameters=["gain_ab", "gain_c"],
                   parameter_labels={"gain_ab": r"$g_{AB}$ / V", "gain_c": r"$g_C$ / V"})
    else:
        axes[0].legend(frameon=False, fontsize=9, loc="upper left")
    return fig, axes


def main(output=Path(__file__).resolve().parents[1] / "output"):
    """Print actual comparison results, then export both typography/panel choices."""
    output.mkdir(exist_ok=True, parents=True)
    shared, partial = fit_examples()
    for name, result in (("All gains shared", shared), ("Partially shared gains", partial)):
        print(f"\n{name}\n{result.report()}\n{result.diagnose()}")
    gap, error = gain_difference(partial)
    delta = shared.statistics["chi2"] - partial.statistics["chi2"]
    # One additional linear parameter with known Gaussian errors: a nested chi-square test.
    print(f"gain C - gain AB = {gap:.6g} +/- {error:.3g} V")
    print(f"nested comparison: delta chi2 = {delta:.6g}, p = {chi2.sf(delta, 1):.6g}")
    for style in ("sans", "tex"):
        with plt.rc_context(plot_style(style)):
            for panel in (True, False):
                fig, _ = plot_comparison(shared, partial, panel=panel)
                for suffix in ("png", "pdf"):
                    fig.savefig(output / f"python_multi_{style}_panel_{panel}.{suffix}", dpi=160)
                plt.close(fig)


if __name__ == "__main__":
    main()
