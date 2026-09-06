"""Run the development wrapper: NumPy model, Julia inference, native Matplotlib.

From the repository: pip install -e './python[plot]', then python python/develop.py.
The fixed arrays are an illustrative decay dataset, not a measurement record.
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scientificfitting import fit_model, plot_fit


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

# Normal Matplotlib customization. No Makie, no new fit to edit the figure.
with plt.rc_context({"font.size": 12, "axes.labelsize": 13, "axes.titlesize": 14}):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    plot_fit(result, ax=axes[0], xlabel="t / s", ylabel="U / V", title="Decay with background",
             curve_kwargs={"color": "#0072B2"})
    axes[0].axhline(result.values["background"], color="black", linestyle="--", linewidth=1,
                   label="fitted background")
    axes[0].legend(fontsize=9, frameon=False)

    # Profiles still refit nuisance parameters in Julia, using the NumPy model.
    scan = result.profile("tau", npoints=31, nsigma=3)
    axes[1].plot(scan.values, scan.delta_cost, color="#0072B2")
    axes[1].axhline(1, color="black", linestyle="--", linewidth=1, label=r"$\Delta C=1$")
    axes[1].set(xlabel=r"$\tau$ / s", ylabel=r"$\Delta C$", title="Lifetime profile")
    axes[1].legend(fontsize=10, frameon=False)

    output = Path(__file__).resolve().parents[1] / "output"
    output.mkdir(exist_ok=True)
    fig.savefig(output / "python_numpy_matplotlib.png", dpi=160)
    fig.savefig(output / "python_numpy_matplotlib.pdf")
