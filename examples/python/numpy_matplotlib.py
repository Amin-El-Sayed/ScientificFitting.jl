"""Run the development wrapper: NumPy model, Julia inference, native Matplotlib.

From the repository: pip install -e './python[plot]', then python python/develop.py.
The fixed arrays are an illustrative decay dataset, not a measurement record.
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scientificfitting import (add_report, fit_model, plot_diagnostics, plot_fit,
                              plot_profile_matrix, plot_style)


def decay(t, amplitude, tau, background):
    """Exponential signal above a constant detector background."""
    return amplitude * np.exp(-t / tau) + background


# Uneven sampling and visible pointwise uncertainty, without generated data.
t = np.array([0.0, 0.3, 0.7, 1.1, 1.6, 2.0, 2.5, 3.1, 3.8, 4.6, 5.5, 6.6])
signal = np.array([3.21, 2.61, 2.27, 1.69, 1.44, 1.05, 0.97, 0.66, 0.61, 0.35, 0.38, 0.20])
sigma = np.array([0.13, 0.12, 0.12, 0.10, 0.10, 0.09, 0.09, 0.08, 0.08, 0.07, 0.07, 0.07])

result = fit_model(
    decay, t, signal,
    p0={"amplitude": 2.8, "tau": 1.5, "background": 0.1},
    sigma_y=sigma,
    bounds={"amplitude": (0.0, 10.0), "tau": (0.1, 10.0), "background": (0.0, 1.0)},
)
print(result.report())
print(result.diagnose())

# Scans refit nuisance parameters in Julia. Compute once, not once per style.
matrix = result.profile_matrix(
    ["amplitude", "tau", "background"], npoints_profile=31, npoints_contour=21, nsigma=3,
)
print(matrix.diagnostics)
output = Path(__file__).resolve().parents[1] / "output"
output.mkdir(exist_ok=True)
labels = {"amplitude": r"$A$ / V", "tau": r"$\tau$ / s", "background": r"$B$ / V"}


def save_figure(fig, name):
    """Export the same native figure as a bitmap and a vector PDF."""
    fig.savefig(output / f"{name}.png", dpi=160)
    fig.savefig(output / f"{name}.pdf")
    plt.close(fig)


for style in ("sans", "tex"):
    with plt.rc_context(plot_style(style)):
        # Panel visibility is independent of typography. Add artists first so
        # the report's legend includes them; native Matplotlib owns the layout.
        for panel in (True, False):
            fig, ax = plot_fit(result, panel=False, xlabel="t / s", ylabel="U / V",
                               title="Decay with detector background")
            ax.axhline(result.values["background"], color="black", linestyle="--", linewidth=1,
                       label="fitted background")
            if panel:
                add_report(fig, result, ax=ax, parameter_labels=labels, expand=True)
            else:
                ax.legend(frameon=False)
            save_figure(fig, f"python_decay_{style}_panel_{panel}")

        fig, axes = plot_diagnostics(result, kinds=("residual", "pull"), xlabel="t / s")
        save_figure(fig, f"python_diagnostics_{style}")
        fig, axes = plot_profile_matrix(matrix, parameter_labels=labels)
        save_figure(fig, f"python_profiles_{style}")
