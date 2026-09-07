"""Poisson decay and unequal-bin spectroscopy, with native Matplotlib exports.

The fixed teaching records match the Poisson and Histograms gallery page.
Install the development wrapper and SciPy as described in docs/src/python.md.
Neither the observations nor their uncertainty are invented by the renderer.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf
from scipy.stats import poisson
from scientificfitting import add_report, fit_histogram_model, fit_poisson_model, plot_style


# Independent 10 s counting windows, starting one minute apart.
time_min = np.array([0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18.])
decay_counts = np.array([48, 37, 35, 27, 27, 17, 22, 13, 16, 8, 13, 5, 11, 4, 7, 2, 6, 1, 5])
edges = np.array([0., 0.4, 0.9, 1.5, 2.2, 3., 4., 5.2, 6.6, 8.2, 10.])
bin_counts = np.array([0, 3, 9, 24, 47, 69, 51, 24, 8, 4])
BLUE, INK = "#0072B2", "#252a30"


def decay(t, signal, decay_constant, background):
    """Expected counts per 10 s window; t is in minutes, not seconds."""
    return signal * np.exp(-decay_constant*t) + background


def spectrum(edges, peak_yield, centroid, width, background_density):
    """Integrate a Gaussian peak plus uniform background over each bin."""
    peak_probability = 0.5*np.diff(erf((edges-centroid)/(np.sqrt(2)*width)))
    return peak_yield*peak_probability + background_density*np.diff(edges)


def fit_examples():
    """Fit both fixed records once, using three explicit candidate starts."""
    count_fit = fit_poisson_model(
        decay, time_min, decay_counts,
        p0={"signal": 40., "decay_constant": 0.15, "background": 3.},
        bounds={"signal": (1e-6, 200.), "decay_constant": (1e-6, 2.), "background": (1e-6, 50.)},
        initial_guesses=[
            {"signal": 70., "decay_constant": 0.30, "background": 2.},
            {"signal": 25., "decay_constant": 0.08, "background": 5.},
        ], multistart=3,
    )
    spectrum_fit = fit_histogram_model(
        spectrum, edges, bin_counts,
        p0={"peak_yield": 210., "centroid": 3.8, "width": 1., "background_density": 1.},
        bounds={"peak_yield": (1e-6, 1000.), "centroid": (0., 10.),
                "width": (0.05, 5.), "background_density": (1e-6, 100.)},
        initial_guesses=[
            {"peak_yield": 300., "centroid": 4.2, "width": 1.5, "background_density": 0.5},
            {"peak_yield": 150., "centroid": 3.2, "width": 0.7, "background_density": 2.},
        ], multistart=3,
    )
    return count_fit, spectrum_fit


def plot_decay(result, *, panel=True, figsize=None):
    """Show observations and conditional predictive quantiles, not a mean-fit band."""
    fig, ax = plt.subplots(figsize=figsize, layout="constrained")
    grid = np.linspace(time_min[0], time_min[-1], 2001)
    mean = decay(grid, **result.values)
    lower, upper = poisson.ppf([0.16, 0.84], mean[:, None]).T
    # Quantiles are integer-valued. Step geometry avoids slanted transitions.
    ax.fill_between(grid, lower, upper, step="post", color=BLUE, alpha=0.18,
                    linewidth=0, label="Poisson 16-84% quantiles")
    ax.plot(grid, mean, color=BLUE, linewidth=1.6, label="expected counts")
    ax.plot(time_min, decay_counts, "o", color=INK, markersize=3, label="observed counts")
    ax.set(xlabel=r"$t$ / min", ylabel="counts per 10 s window", ylim=(0, None),
           title="Radioactive decay with detector background")
    if panel:
        add_report(fig, result, ax=ax, expand=figsize is None,
                   parameter_labels={"signal": r"$S_0$", "decay_constant": r"$\lambda$ / min$^{-1}$",
                                     "background": r"$B$"},
                   statistic_labels={"chi2_ndf": r"$D/\mathrm{ndf}$", "pvalue": r"$p$ (asymptotic)"})
    else:
        ax.legend(frameon=False)
    return fig, ax


def plot_spectrum(result, *, panel=True, figsize=None):
    """Keep the real bin edges; heights are counts per bin, not a density."""
    fig, ax = plt.subplots(figsize=figsize, layout="constrained")
    mean = spectrum(edges, **result.values)
    lower, upper = poisson.ppf([0.16, 0.84], mean[:, None]).T
    ax.stairs(upper, edges, baseline=lower, fill=True, color=BLUE, alpha=0.18,
              linewidth=0, label="Poisson 16-84% quantiles")
    ax.stairs(mean, edges, baseline=None, color=BLUE, linewidth=1.6, label="integrated model")
    ax.stairs(bin_counts, edges, baseline=None, color=INK, linewidth=1., label="observed counts")
    ax.set(xlabel="pulse amplitude / V", ylabel="counts per bin", ylim=(0, None),
           title="Pulse spectrum with unequal bin widths")
    if panel:
        add_report(fig, result, ax=ax, expand=figsize is None,
                   parameter_labels={"peak_yield": r"$N_{\mathrm{peak}}$", "centroid": r"$\mu$ / V",
                                     "width": r"$s$ / V", "background_density": r"$b$ / V$^{-1}$"},
                   statistic_labels={"chi2_ndf": r"$D/\mathrm{ndf}$", "pvalue": r"$p$ (asymptotic)"})
    else:
        ax.legend(frameon=False)
    return fig, ax


def main(output=Path(__file__).resolve().parents[1] / "output"):
    """Print actual core reports and export each completed fit without refitting."""
    output.mkdir(exist_ok=True, parents=True)
    for name, result, render in zip(("counts", "spectrum"), fit_examples(), (plot_decay, plot_spectrum)):
        print(f"\n{name}\n{result.report()}\n{result.diagnose()}")
        for style in ("sans", "tex"):
            with plt.rc_context(plot_style(style)):
                for panel in (True, False):
                    fig, _ = render(result, panel=panel)
                    for suffix in ("png", "pdf"):
                        fig.savefig(output / f"python_{name}_{style}_panel_{panel}.{suffix}", dpi=160)
                    plt.close(fig)


if __name__ == "__main__":
    main()
