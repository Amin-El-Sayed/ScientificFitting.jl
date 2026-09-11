# Python Interface

NumPy model functions, the Julia numerical core, editable Matplotlib figures.
Python 3.10+ is required; Makie is neither installed nor loaded.

```bash
python -m pip install 'scientificfitting[plot]'
# Optional: SciPy distributions and sparse covariance matrices.
python -m pip install scipy
```

Omit `[plot]` for fitting and reports only. JuliaCall provisions Julia and the
compatible 0.2.x core on first use, with network access and compilation.
This does **not** remove the Julia runtime's size or startup cost.
Importing `scientificfitting` alone does not start Julia.

[PyPI](https://pypi.org/project/scientificfitting/) ·
[Conda-forge recipe status](https://github.com/conda-forge/staged-recipes/pull/34795).
The pip command also works inside a Conda environment.

## Fit, Inspect, Plot

Code cells run in order in one Python session.

```python
import numpy as np
import matplotlib.pyplot as plt
from scientificfitting import fit_model, plot_fit

def line(x, slope, offset):
    return slope * x + offset

x = np.array([0., 1., 2., 3.])
y = np.array([0.1, 1.2, 1.9, 3.2])
result = fit_model(line, x, y, p0={"slope": 1., "offset": 0.}, sigma_y=0.2)
print(result.report())
print(result.diagnose())

# A normal Matplotlib axis: adding a marker does not repeat the fit.
fig, ax = plot_fit(result, xlabel="x / mm", ylabel="U / V")
ax.axvline(1.5, color="black", linestyle="--")
fig.savefig("calibration.pdf")
plt.close(fig)
```

**Callback contract:** `model(x, **parameters)` returns one prediction per
coordinate. `x` is read-only. `p0` keys bind parameters by name and set result-array
order. Use deterministic models; derivatives use finite differences, including
profile refits. Optional `jacobian` and `x_derivative` callbacks differentiate the
unweighted model. `tol` controls optimization, not statistical uncertainty.

```python
# Read values directly; never parse the terminal report.
slope_value = result.values["slope"]
report = result.report(structured=True)
slope = report.parameters["slope"]
print(slope.value, slope.uncertainty, slope.fixed)

diagnosis = result.diagnose(structured=True, max_actions=3)
for finding in diagnosis.findings:
    print(finding.code, finding.recommendation)

# Local uncertainty of the fitted mean, excluding observation noise.
prediction, uncertainty = result.predict(x, uncertainty=True)
```

Arrays are independent, read-only snapshots. `result.converged` records solver
convergence, not model validity; an unknown iteration count is `None`.
`diagnosis.status == "ok"` means the implemented checks found no issue.

## Choose The Error Model

| Data / uncertainty | Entry point or option |
|---|---|
| Gaussian x/y errors | `fit_model(..., sigma_x=..., sigma_y=...)` |
| Dense or SciPy sparse covariance | `fit_model(..., cov_x=..., cov_y=...)` |
| Named error sources | `ErrorComponent`; see `help(ErrorComponent)` |
| Non-Gaussian measurement errors | `fit_likelihood_model(..., logprob=...)` |
| Poisson observations | `fit_poisson_model(...)` |
| Bin counts | `fit_histogram_model(...)` for expected counts; `fit_histogram_density(...)` to integrate a density |
| Independent events | `fit_unbinned_model(...)`; `fit_extended_unbinned_model(...)` also fits event yields |
| Dependent observations / custom likelihood | `fit_custom(...)` with the joint `-2 log L` and `nobs` |

Student-t errors allow heavier tails than Gaussian errors. Here the scales are
known separately for each observation:

```python
from scipy import stats
from scientificfitting import fit_likelihood_model

locations = np.array([-1., -0.4, 0., 0.5, 1.2, 1.8, 2.4])
readings = np.array([-1.1, -0.23, 0.31, 0.94, 2.04, 2.91, 6.9])
scales = np.array([0.12, 0.20, 0.14, 0.18, 0.11, 0.25, 0.20])

def measurement_logprob(y, prediction, slope, offset):
    # One normalized log density per independent observation, not a sum.
    return stats.t.logpdf(y, df=4, loc=prediction, scale=scales)

robust_result = fit_likelihood_model(
    line, locations, readings, logprob=measurement_logprob,
    p0={"slope": 1., "offset": 0.},
)
print(robust_result.report())
```

For `df=4`, standard deviation = ``\sqrt{2}\times`` scale. The last point has
less influence, but is not discarded. Choose the distribution from the error
process, not to conceal a wrong mean model. Discrete observations need a log
**mass** (`logpmf`); dependent observations need a joint likelihood. No universal
goodness-of-fit p-value is assigned to this custom distribution.

```python
from scientificfitting import add_report, plot_style

with plt.rc_context(plot_style("sans")):
    fig, ax = plt.subplots(layout="constrained")
    grid = np.linspace(locations.min(), locations.max(), 300)
    ax.plot(grid, line(grid, **robust_result.values), color="#0072B2", label="fitted mean")
    # Error-distribution quantiles, NOT uncertainty of the fitted line.
    ax.errorbar(locations, readings, yerr=stats.t.ppf(0.84, df=4)*scales,
                fmt="o", color="black", markersize=3, elinewidth=0.8, capsize=2,
                label="data; 16-84% error range")
    ax.set(xlabel="reference setting", ylabel="response / V", title="Student-t errors")
    add_report(fig, robust_result, ax=ax, statistics=("cost_min",),
               statistic_labels={"cost_min": r"$-2\log L$"}, expand=True)
    fig.savefig("student_t_errors.pdf")
    plt.close(fig)
```

### Non-Smooth Likelihoods

A Laplace location fit estimates the sample median. Its cost has corners;
use a derivative-free solver rather than interpreting a Hessian there:

```python
readings_laplace = np.array([-1.2, -0.1, 0.2, 0.4, 0.8, 1.3, 5.0])
def constant(x, location):
    return np.full_like(x, location)

laplace_result = fit_likelihood_model(
    constant, np.arange(len(readings_laplace)), readings_laplace,
    logprob=lambda y, mu, location: stats.laplace.logpdf(y, loc=mu, scale=1),
    p0={"location": 0.1}, optimizer="nelder_mead", tol=1e-10,
)
print(laplace_result.report())
# No local covariance scale: supply the actual trial values.
scan = laplace_result.profile("location", values=[-0.2, 0., 0.4, 0.7, 1.])
```

Nelder-Mead supports bounds, fixed values and Gaussian parameter terms, but not
nonlinear constraints. It defaults to `parameter_covariance="none"` (`NaN`
free errors, not zero); `maxiters` limits function evaluations. Select
`"hessian"` only for a locally smooth cost. Zero probability returns `-np.inf`;
start at a finite likelihood. Moving support can invalidate standard profile
thresholds: see the [support-boundary calculation](likelihood_models.md#A-Moving-Support-Boundary).

## Profiles And Contours

A profile fixes one parameter and refits the others. Four exposures with one
observed event give an asymmetric rate interval:

```python
from scientificfitting import fit_poisson_model

counts = np.array([0, 0, 1, 0])
count_fit = fit_poisson_model(
    lambda x, rate: np.full_like(x, rate), np.arange(4), counts,
    p0={"rate": 0.3}, bounds={"rate": (0.001, 10)},
)
interval = count_fit.profile_interval("rate", npoints=61, nsigma=4)
print(interval.lower, interval.upper)  # expected counts per exposure
```

`delta_cost=1` has an **asymptotic**, not exact low-count, 68.27% interpretation.
Two-parameter 68.27%/95.45% regions instead use `[2.30, 6.18]`. Missing crossings
remain `NaN`; failed refits remain gaps. See [Profiles and Contours](profiles_contours.md).

```python
from scientificfitting import (plot_contour, plot_diagnostics, plot_profile,
                              plot_profile_matrix)

# Explicit selection and grid sizes bound the number of nuisance refits.
matrix = result.profile_matrix(
    ["slope", "offset"], npoints_profile=9, npoints_contour=7, nsigma=3,
)
pair = matrix.contours["slope", "offset"]
actions = matrix.triage()  # panels needing attention; no new fits

# Rendering completed scans does not minimize again.
fig, ax = plot_profile(interval.profile_result, delta_max=5)
plt.close(fig)
fig, ax = plot_contour(pair)
plt.close(fig)
fig, axes = plot_profile_matrix(matrix)
fig.savefig("calibration_profiles.pdf")
plt.close(fig)

# Reassess stored scans without changing their confidence thresholds.
profile_review = interval.profile_result.diagnose(structured=True)
pair_review = pair.diagnose(tolerance=0.5, structured=True)
```

`tolerance` above measures departure from the local parabola/ellipse, **not**
confidence level. Arrays store `delta_cost[i, j]` at `(x[i], y[j])`; direct
Matplotlib calls need `pair.delta_cost.T`. `result.report(errors="profile")`
computes asymmetric errors with additional fits; its covariance remains local.

## Customize With Matplotlib

Style (`"sans"` / `"tex"`) and report visibility (`panel=True` / `False`) are
independent. `"tex"` uses bundled STIX fonts and MathText, not external TeX.
`rc_context` leaves global settings unchanged.

```python
with plt.rc_context(plot_style("tex")):
    fig, ax = plot_fit(result, panel=False, xlabel="x / mm", ylabel="U / V")
    # Add labeled artists before constructing the report legend.
    ax.axvline(1.5, color="black", linestyle=":", label="reference position")
    panel = add_report(
        fig, result, ax=ax, expand=True,
        parameter_labels={"slope": r"$m$", "offset": r"$b$"},
        statistics=("chi2_ndf", "pvalue"),
    )
    fig.savefig("calibration_with_reference.pdf")
    plt.close(fig)

fig, axes = plot_diagnostics(result, kinds=("residual", "pull"), xlabel="x / mm")
fig.savefig("calibration_diagnostics.pdf")
plt.close(fig)
```

| Change | Control |
|---|---|
| Use an existing subplot | `plot_fit(result, ax=ax, panel=False)` |
| More room for the graph | `figsize=(width, height)` in inches |
| Report below the axes | `add_report(..., position="bottom")` |
| Edit/remove the report | The returned object is a Matplotlib `Legend`; `panel.remove()` removes it |
| Style individual artists | `curve_kwargs`, `point_kwargs`, `band_kwargs` |
| Show asymmetric estimates | Pass a completed `FitReport` to `add_report` |

Use `layout="constrained"` for outside reports on your figures. `expand=True`
allows canvas growth; an explicitly fixed tiny canvas can still be too small.
`plot_fit` and x-y residual helpers target Gaussian results; other fit families
use ordinary Matplotlib plus `add_report`. Correlated pulls are whitened
residuals; their unit bands are reference guides, not coverage intervals.

## Large Datasets And In-Place Models

SciPy sparse `cov_x`/`cov_y` inputs are not densified. Static y-covariance reuses
its factorization; parameter-dependent covariance does not. Sparse factors can
still fill in. A known whitening operator avoids storing the matrix:

```python
from scientificfitting import WhiteningOperator

# Stationary AR(1): C[i,j] = sigma**2 * rho**abs(i-j), equal sample spacing.
sigma, rho, n = 0.2, 0.6, len(y)
def whiten(out, residual):
    out[0] = residual[0] / sigma
    # Remove the predictable neighbor contribution, then standardize.
    out[1:] = (residual[1:] - rho * residual[:-1]) / (sigma * np.sqrt(1 - rho**2))

noise = WhiteningOperator(
    whiten, 2*n*np.log(sigma) + (n-1)*np.log1p(-rho**2),  # log(det(C))
    marginal_sigma=sigma, inplace=True,  # marginal_sigma is only plot metadata
)
correlated = fit_model(line, x, y, p0={"slope": 1., "offset": 0.}, whitening=noise)

def line_inplace(out, x, slope, offset):
    np.multiply(x, slope, out=out)
    out += offset

inplace_result = fit_model(line_inplace, x, y, p0={"slope": 1., "offset": 0.},
                           sigma_y=0.2, inplace=True)
```

Whitening replaces other observation errors. In-place callbacks fill every
entry and return `None`; never retain the borrowed arrays. An in-place Jacobian
fills an `(n, k)` matrix, columns in `p0` order. This avoids output copies, not all
allocations. Without `x_derivative`, the model must also be defined near the
measured coordinates for finite differencing.

## Event Densities And Complete Workflows

An unbinned fit models the observations themselves. For uncensored positive
waiting times, evaluate the exponential density in one NumPy call:

```python
from scientificfitting import fit_unbinned_model

waiting_times = np.array([0.12, 0.28, 0.51, 0.62, 0.75, 1.3, 1.8])
def waiting_pdf(t, tau):
    return np.exp(-t/tau)/tau

waiting_fit = fit_unbinned_model(
    waiting_pdf, waiting_times, p0={"tau": 0.5}, bounds={"tau": (0.01, 5.)},
    vectorized=True,  # one array call, not one call per event
)
```

`tau` estimates the mean waiting time. A detection threshold requires a
correspondingly normalized truncated density. `vectorized=True` also supports
extended likelihoods and adaptive histogram integration; `rtol` controls the
integral accuracy. For known bin integrals, use `fit_histogram_model` directly.

| Complete script in `examples/python/` | Demonstrates |
|---|---|
| `numpy_matplotlib.py` | Nonlinear decay, editable reports, residuals, profiles, both styles |
| `likelihood_workflows.py` | Poisson decay and unequal-width bins; expectations integrated per bin |
| `multi_dataset_calibration.py` | Named `parameter_map`, shared gains, full covariance of their difference, nested-model comparison |

The scripts fit and scan once, then render. Poisson quantile bands describe count
fluctuations conditional on the fitted mean, not parameter uncertainty.
Backgrounds near zero need profiles; local symmetric errors can cross the bound.

## Development Setup

Only contributors need a local Julia checkout. From the repository root, in a
separate Python environment:

```bash
python -m pip install -e './python[plot,test]'
python python/develop.py
python examples/python/numpy_matplotlib.py
python -m pytest python/tests
```

`develop.py` persistently selects this checkout in its JuliaPkg environment.
Use a **fresh environment** for registry-installation checks. Installed-wheel
checks live in `python/tests/check_install.py`; register the Julia core before
publishing a Python release that depends on it. See
[Python startup measurements](performance.md#Python-Startup) for installation
size, latency, and reproducible benchmarks rather than timing promises.
