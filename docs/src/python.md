# Python Interface (Development Preview)

The `python/` directory contains an unreleased wrapper, not yet a PyPI or conda
release. It uses NumPy models, the same Julia numerical core, and optional
native Matplotlib plots. It neither installs nor loads Makie. From this
checkout, in a Python 3.10+ virtual environment:

```bash
python -m pip install -e './python[plot,test]'
python python/develop.py
python examples/python/numpy_matplotlib.py
python examples/python/likelihood_workflows.py
python examples/python/multi_dataset_calibration.py
python -m pytest python/tests
```

The development step selects this checkout because the wrapper requires the
new **0.2.x core**, not the registered v0.1.2 API. JuliaCall manages Julia and
its packages automatically; first use needs network access and compilation.
This avoids a manual Julia installation, but does **not** remove the Julia
runtime's disk footprint or startup cost.

```python
import numpy as np
from scientificfitting import fit_model, plot_fit

def line(x, slope, offset):
    return slope * x + offset

x = np.array([0., 1., 2., 3.])
y = np.array([0.1, 1.2, 1.9, 3.2])
result = fit_model(line, x, y, p0={"slope": 1., "offset": 0.}, sigma_y=0.2)
print(result.report())

fig, ax = plot_fit(result, xlabel="x / mm", ylabel="U / V")
ax.axvline(1.5, color="black", linestyle="--")
fig.savefig("calibration.pdf")
```

Models receive a read-only NumPy coordinate array and named float parameters.
The keys of `p0` must match the model's parameter names; mapping order sets
the order of result arrays, not how parameters bind to the model. Models must
be deterministic and, for the default solvers, smooth around evaluated points. Finite differences are
used consistently, including when bounds or profiles trigger a different
solver. Optional analytic `jacobian` and `x_derivative` callbacks use the same
Python argument convention. `jacobian(x, **parameters)` returns an `(n, k)`
matrix with columns in `p0` key order; `x_derivative(x, **parameters)` returns
an `(n,)` vector. Both differentiate the unweighted model; the core handles
uncertainty propagation and whitening. Do not flatten the Jacobian matrix.
The default solver `tol=1e-6` accounts for differenced-gradient noise;
`tol` remains configurable and is not a bound on parameter error.
`maxiters` limits each solver run. `result.converged` and the report reflect
the actual solver status, including when the iteration limit is too small.
All fit families accept `initial_guesses` as additional named dictionaries
or numerical vectors in `p0` order. Set `multistart` explicitly: it is the
total candidate budget **including `p0`**, not the number of extra runs.
For two distinct additional guesses, use `multistart=3`. With its default of
one, only `p0` is used. The core selects the lowest converged result; trying
several starts is useful for local minima, but does not prove global optimality.

The preview covers Gaussian fits with x/y errors, dense or SciPy sparse
covariance, named `ErrorComponent` sources, and matrix-free `WhiteningOperator`,
Poisson counts, expected-count and integrated-density histograms, unbinned
and extended-unbinned likelihoods, indexed and multi-dataset fits, custom
costs, and user-defined observation distributions. Parameter bounds, named
fixed values and Gaussian priors, correlated parameter constraints, nonlinear
constraints, profiles, contours, and core reports reuse the Julia engine.
`result.predict(x, uncertainty=True)` returns a
Gaussian model mean and its local standard uncertainty, without observation
noise. `plot_fit` accepts an existing Matplotlib `ax` and ordinary artist
keyword dictionaries; it does not change global plotting settings or refit.

The plotting facade now includes Gaussian fit/residual panels, editable
reports, and profiles/contours for every supported fit family. Like the Julia
renderer, `plot_fit` and x-y residual helpers target Gaussian results;
likelihood observations and multiple datasets can be drawn using ordinary
Matplotlib with `add_report` for the fit estimates. Worked examples of those
compositions are included below. This preview is
not a claim of v0.2 release readiness. Platform results and a clean installation
against the registered 0.2 core must be verified before publishing a wheel.

For non-Gaussian regression, `fit_likelihood_model` accepts a vectorized
`logprob(y, prediction, **parameters)` returning one normalized log density
or log probability mass per independent observation. For example, SciPy's
`stats.t.logpdf(y, df=4, loc=prediction, scale=scales)` describes Student-t
errors; the scale is not their standard deviation. Capture known per-point
scales in the callback. Unlike `fit_custom`, no manual likelihood summation or
observation count is needed. Both callbacks are batched. For event or histogram
density fits, use `vectorized=True` with NumPy-compatible densities to avoid
one Python call per point. Dependent observations require a joint likelihood.

## Choose The Observation Distribution

A few large residuals can be plausible under a heavy-tailed error model without
requiring a different mean response. This controlled example uses a Student-t
distribution with four degrees of freedom and known, point-specific scales:

```python
from scipy import stats
from scientificfitting import fit_likelihood_model

locations = np.array([-1., -0.4, 0., 0.5, 1.2, 1.8, 2.4])
readings = np.array([-1.1, -0.23, 0.31, 0.94, 2.04, 2.91, 6.9])
scales = np.array([0.12, 0.20, 0.14, 0.18, 0.11, 0.25, 0.20])

def measurement_logprob(y, prediction, slope, offset):
    # One normalized log density per observation, not a summed objective.
    return stats.t.logpdf(y, df=4, loc=prediction, scale=scales)

robust_result = fit_likelihood_model(
    line, locations, readings, logprob=measurement_logprob,
    p0={"slope": 1., "offset": 0.},
)
print(robust_result.report())
```

For this distribution, the standard deviation is ``\sqrt{2}`` times the scale.
Its negative log density grows only logarithmically for large residuals, rather
than quadratically as in a Gaussian model. The last observation therefore has
less influence, but it is not removed. This assumption needs a physical or
empirical justification; it is not a way to hide a wrong response model.
No generic goodness-of-fit p-value is invented for this custom distribution.
For count data, use a log **mass**, such as `stats.binom.logpmf`, instead of a
continuous density. See the [Likelihood Fitting API](api_fitting.md).

These likelihood results can be drawn with the same editable Matplotlib tools:

```python
import matplotlib.pyplot as plt
from scientificfitting import add_report, plot_style

with plt.rc_context(plot_style("sans")):
    fig, ax = plt.subplots(layout="constrained")
    grid = np.linspace(locations.min(), locations.max(), 300)
    ax.plot(grid, line(grid, **robust_result.values), color="#0072B2", label="fitted mean")
    ax.errorbar(locations, readings, yerr=stats.t.ppf(0.84, df=4)*scales,
                fmt="o", color="black", markersize=3, elinewidth=0.8, capsize=2,
                label="data; 16-84% error range")
    ax.set(xlabel="reference setting", ylabel="response / V", title="Student-t measurement errors")
    add_report(fig, robust_result, ax=ax, statistics=("cost_min",),
               statistic_labels={"cost_min": r"$-2\log L$"}, expand=True)
    fig.savefig("student_t_errors.pdf")
    plt.close(fig)
```

The bars show the specified error distribution's central range, not the
uncertainty of the fitted line. For the latter, use the fitted parameter
likelihood and its local approximation or profiles, as appropriate.

### Non-Smooth Errors And Hard Support

Solver choice and local error calculation are separate controls. All likelihood
helpers accept `optimizer="auto"` (LBFGS, or IPNewton for nonlinear constraints),
`"lbfgs"`, `"ipnewton"`, or `"nelder_mead"`. NLopt's Nelder-Mead needs no
derivatives and respects bounds, fixed values and Gaussian parameter terms.
It rejects nonlinear constraints rather than dropping them. These are local
searches over continuous parameters, even when observations are discrete.

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
scan = laplace_result.profile("location", values=[-0.2, 0., 0.4, 0.7, 1.])
```

The fitted location is the sample median, 0.4. The cost has corners at the
observations, so Nelder-Mead defaults to `parameter_covariance="none"`: free
errors are `NaN`, not zero. Explicit `"hessian"` is available for locally smooth
costs, and `"none"` also works with the gradient solvers. Profiles preserve
both settings; specify actual scan values when no local error scale is available.
For Nelder-Mead, `maxiters` limits objective evaluations, `tol` sets absolute and
relative parameter tolerances, and `result.iterations` is `None` because the
backend exposes no iteration count. Exhausting the budget is not convergence.
Scale parameters appropriately; a tiny `tol` is not a statistical precision
claim. Equal objective values alone do not stop the simplex search.

For a distribution with a moving support boundary, return `-np.inf` for zero
probability and start at a finite likelihood. For example,
`stats.expon.logpdf(y, loc=mu, scale=1)` describes a location plus a positive
exponential delay; its fitted location is the smallest observation. No density
floor is inserted. Finding that minimum does **not** validate the usual
chi-square profile thresholds: the [exact support-boundary example](likelihood_models.md#A-Moving-Support-Boundary)
shows why a nominal 68% threshold can have only 39% coverage.

## Counts, Histograms, And Shared Parameters

Two complete scripts use the same fixed arrays as the Julia gallery. Both print
actual core reports and export PNG/PDF in `sans` and `tex`, with panels independently
on or off. Each fit runs once, not once per output style.

| Script in `examples/python/` | Scientific calculation | Native Matplotlib composition |
|---|---|---|
| `likelihood_workflows.py` | Poisson decay with background; a Gaussian peak integrated over unequal bins | Count observations, expected counts, and conditional Poisson quantile regions |
| `multi_dataset_calibration.py` | One shared gain versus a separate gain for channel C | Three datasets, mean-fit bands, two pull panels, and a selected-parameter report |

For the count example, each observation is a separate **10 s** exposure, one
minute apart. The signal and background are counts per exposure, while the decay
constant is in inverse minutes. The Poisson 16th/84th percentiles are computed
at the fitted mean. Their steps reflect discrete counts; they neither include
parameter uncertainty nor guarantee exactly 68% probability at low means.
The fitted background is weakly determined: its local symmetric error extends
below zero, despite a nonnegative physical bound. Use a background profile
before interpreting that number as an interval; convergence does not make
the local approximation reliable at a boundary.
Matplotlib's [step-filled regions](https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.fill_between.html)
and [stairs](https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.stairs.html)
preserve vertical transitions. Histogram heights are **counts per bin**, and the
expected counts integrate the model over the actual, unequal bin widths.

In the channel example, `parameter_map` maps local model names to global fit
names, for example `{"gain": "gain_ab", "offset": "offset_a"}`. Reusing
`gain_ab` in two mappings shares that parameter, while `gain_c` remains separate.
The script propagates the **full** parameter covariance for the gain difference.
Its nested-model test has one additional parameter, identical observations, and
the same known Gaussian errors. These conditions justify that comparison;
different likelihoods or uncertainty assumptions cannot be interchanged silently.

## Batched Event And Histogram Densities

An unbinned fit models the distribution of the observations themselves, not
the scatter around a response curve. For uncensored positive waiting times,
an exponential model has density ``f(t;\tau)=\exp(-t/\tau)/\tau``:

```python
from scientificfitting import fit_unbinned_model

waiting_times = np.array([0.12, 0.28, 0.51, 0.62, 0.75, 1.3, 1.8])

def waiting_pdf(t, tau):
    return np.exp(-t/tau)/tau

# One NumPy call evaluates all events at each trial parameter value.
waiting_fit = fit_unbinned_model(
    waiting_pdf, waiting_times, p0={"tau": 0.5}, bounds={"tau": (0.01, 5.)},
    vectorized=True,
)
```

Here the fitted `tau` equals the sample mean. If a detection threshold removes
short waiting times, normalize a truncated density instead; the formula above
does not describe that observation process unchanged.

`vectorized=True` also works with `fit_histogram_density` and
`fit_extended_unbinned_model`. Callbacks receive read-only, one-dimensional
NumPy arrays and return one density/intensity per entry. For integrals, Julia
uses [QuadGK's batched adaptive quadrature](https://juliamath.github.io/QuadGK.jl/stable/quadgk-examples/#Batched-integrand-evaluation);
`rtol` still controls integration, separately for each histogram bin. This is
not midpoint sampling or a fixed coarse grid. The scalar default remains
available for callbacks using `math` functions or scalar conditionals.
If an analytic bin integral is available, pass its expected counts directly to
`fit_histogram_model` and avoid quadrature altogether.

## Large Gaussian Fits

SciPy sparse matrices and arrays are copied as canonical CSC buffers, not
converted to dense arrays. Install `./python[sparse]` to add SciPy if needed.
The same `cov_x`/`cov_y` arguments accept either representation. Static sparse
y covariance reuses its factorization during optimization; parameter-dependent
effective covariance must be refactorized at each parameter point. Sparse factorization
can still produce fill-in: input sparsity alone does not guarantee linear
memory or runtime. A known whitening operator can avoid the matrix entirely:

```python
from scientificfitting import WhiteningOperator

sigma = 0.2  # Known marginal standard deviation of stationary AR(1) noise.
rho = 0.6    # Known correlation between neighboring observations.
n = len(y)

def whiten(out, residual):
    out[0] = residual[0] / sigma
    out[1:] = (residual[1:] - rho * residual[:-1]) / (sigma * np.sqrt(1 - rho**2))

noise = WhiteningOperator(
    whiten, 2*n*np.log(sigma) + (n-1)*np.log1p(-rho**2),
    marginal_sigma=sigma, inplace=True,
)
correlated = fit_model(line, x, y, p0={"slope": 1., "offset": 0.}, whitening=noise)
```

This operator describes ``C_{ij}=\sigma^2\rho^{|i-j|}``, with observations
ordered at equal intervals. Subtracting ``\rho r_{i-1}`` removes the predictable
part of each residual; dividing by its innovation standard deviation gives
unit-variance residuals. The second constructor argument is ``\log\det C``,
not the determinant of the whitening operator. It is required for the
normalized likelihood. The operator replaces all other observation errors;
`marginal_sigma` is plotting metadata, not an additional error source.

For a model that fills an existing buffer, use `inplace=True` on the fit:

```python
def line_inplace(out, x, slope, offset):
    np.multiply(x, slope, out=out)
    out += offset

result = fit_model(line_inplace, x, y, p0={"slope": 1., "offset": 0.},
                   sigma_y=0.2, inplace=True)
```

An optional in-place `jacobian(out, x, **parameters)` fills an
``n\times k`` array in `p0` order. Both callbacks must fill every entry and
return `None`; the supplied arrays borrow Julia memory and must not be kept.
`x_derivative(x, **parameters)` still returns an array. In-place callbacks
avoid model-output copies; they do not make the complete fit allocation-free.

## Results Without Text Parsing

Result arrays are independent read-only NumPy snapshots. Their order follows
`result.parameter_names`; use `result.values["slope"]` for named access.
`statistics` contains the complete core summary, including the cost name and
`minus2loglik_min`. `backend`, `message`, and `iterations` preserve solver
status; an unknown iteration count is `None`, not the configured limit.

```python
report = result.report(structured=True)
slope = report.parameters["slope"]
slope.value, slope.uncertainty, slope.fixed

diagnosis = result.diagnose(structured=True, max_actions=3)
[(finding.code, finding.evidence, finding.recommendation)
 for finding in diagnosis.findings]
```

Without `structured=True`, both methods return the core's actual text.
The structured diagnosis includes `status`, `severity_counts`, all findings,
and deduplicated `next_actions`. Limiting actions does not remove findings.
`"ok"` means the implemented checks found no warning or critical issue, not
that the physical model is established.

Gaussian results also expose `x`, `y`, `model_y`, `residuals`,
`weighted_residuals`, and `jacobian`. The Jacobian differentiates the weighted
residual vector, including auxiliary parameter terms where present. Whitened
residuals are not pointwise error-bar pulls when the observations are correlated.
General likelihood results use `None` for these fields because an arbitrary
objective need not define an x-y residual. Numerical condition estimates,
warnings, and active-bound parameter names are in `numerical_diagnostics`.

## Profiles And Asymmetric Intervals

A profile fixes one parameter and refits the others. Its cost increase answers
which values remain compatible with the data, without requiring a parabolic
shape. Consider four example exposures with counts `[0, 0, 1, 0]`:

```python
from scientificfitting import fit_poisson_model

counts = np.array([0, 0, 1, 0])
count_fit = fit_poisson_model(
    lambda x, rate: np.full_like(x, rate), np.arange(4), counts,
    p0={"rate": 0.3}, bounds={"rate": (0.001, 10)},
)
interval = count_fit.profile_interval("rate", npoints=61, nsigma=4)
interval.lower, interval.upper
```

Here `rate` is the expected count **per exposure**, not per second unless
each exposure lasts one second. The interval is asymmetric because the
likelihood with so few counts is not a parabola. The default threshold is
`delta_cost=1`; its usual 68.27% interpretation is an asymptotic approximation,
not an exact low-count coverage guarantee. See [Profiles and Contours](profiles_contours.md).

`interval.profile_result` retains every scan value, absolute cost, cost
increase, and diagnosis against the local covariance approximation. If a scan
already exists, `scan.interval()` extracts its crossings without fitting again.
Unbracketed endpoints remain `NaN`; failed scan points remain infinite and
produce findings. The wrapper never fills a failed region by interpolation.
Automatic profile and contour grids respect declared parameter bounds. Explicit
grids remain unchanged, so deliberately testing a forbidden value still records
a failure (or raises with `on_failure="throw"`). A physical bound is not silently
reported as a threshold crossing.

`result.report(errors="profile")` can use asymmetric profile errors in the
parameter report. It performs additional fits. Missing crossings remain
`NaN` there too; the reported covariance matrix is still the **local**
covariance, not a new global covariance inferred from those intervals.

## Several Parameters

```python
# Reuse the completed line fit; choose parameters explicitly to bound scan work.
matrix = result.profile_matrix(
    ["slope", "offset"], npoints_profile=9, npoints_contour=7, nsigma=3,
)
pair = matrix.contours["slope", "offset"]
actions = matrix.triage()  # Core-ordered panels needing attention; no new fits.
```

`parameters` preserves selection order. `profiles` is keyed by parameter
name; each contour key is `(x_parameter, y_parameter)`. The pair above stores
`x`, `y`, `cost_values`, `delta_cost`, `levels`, `best_values`,
`local_covariance`, and `diagnostics`. Its array element `delta_cost[i,j]`
belongs to `x[i], y[j]`; Matplotlib expects the transpose:

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(layout="constrained")
ax.contour(pair.x, pair.y, pair.delta_cost.T, levels=pair.levels)
ax.set(xlabel="slope", ylabel="offset")
plt.close(fig)
```

The matrix also retains `local_stderr`, `local_correlation`,
`panel_status`, and the combined `diagnostics`. `triage(include_ok=True)`
includes successful panels. All findings and their ordering come from the
Julia core. Consuming these results never repeats the scans.

To reassess a stored scan without fitting again:

```python
profile_review = interval.profile_result.diagnose(structured=True)
pair_review = pair.diagnose(tolerance=0.5, structured=True)
```

`tolerance` is the allowed absolute difference in delta-cost from the local
parabola/ellipse (defaults: 0.25 for profiles, 0.5 for contours). It is **not**
a confidence level: changing it does not change any scan value, crossing, or
contour level, and does not repair the uncertainty approximation. These methods
return a new core report; the original `diagnostics` remains available. Failed
refits and missing crossings are still reported, independently of this
shape-comparison setting. Without `structured=True`, the returned value is the
core's dashboard text. `max_actions` only limits its short action list.

For two profiled parameters the default levels `[2.30, 6.18]` differ from
the one-parameter thresholds `[1, 4]`. They represent the usual asymptotic
68.27% and 95.45% joint regions, not marginal errors or posterior probabilities.

## Native Matplotlib, Independent Controls

`plot_style("sans")` and `plot_style("tex")` return ordinary Matplotlib rcParams.
The latter uses bundled STIX fonts and MathText, not an external TeX installation.
Both support a report panel; `panel=False` disables it without changing the font,
colors, uncertainty calculation, or fit. No global style is changed.

```python
from scientificfitting import add_report, plot_style

with plt.rc_context(plot_style("tex")):
    fig, ax = plot_fit(result, panel=False, xlabel="x / mm", ylabel="U / V")
    ax.axvline(1.5, color="black", linestyle=":", label="reference position")
    panel = add_report(
        fig, result, ax=ax, expand=True,
        parameter_labels={"slope": r"$m$", "offset": r"$b$"},
        statistics=("chi2_ndf", "pvalue"),
    )
    fig.savefig("calibration_with_reference.pdf")
    plt.close(fig)
```

Add labeled artists **before** constructing the report if they should appear in
its legend. The returned `panel` is a normal Matplotlib `Legend`; it can be edited
or removed with `panel.remove()`. `parameters` selects the displayed estimates,
`parameter_labels` changes their labels, and `statistics` selects report fields.
`statistic_labels` changes field labels without changing their values: for a
Poisson fit, the core's `chi2_ndf` field contains **deviance/ndf**, not a Gaussian
residual sum. The count script labels it accordingly.
Passing a completed `FitReport` also supports previously computed asymmetric
profile errors without launching another scan. Unavailable values stay visible
as unavailable, not zero. Reports show the actual optimizer convergence flag;
convergence alone does not validate the model.

Outside reports use Matplotlib's
[constrained layout](https://matplotlib.org/stable/users/explain/axes/constrainedlayout_guide.html).
A default new `plot_fit` figure grows by the measured report extent, rather
than squeezing the graph. An explicit `figsize=(width, height)` is respected
in inches; increasing width gives the data axes more room while the report text
keeps its size. Very small explicit canvases can still be too small for the
chosen text. `add_report(..., expand=True)` opts into content-sized expansion;
`position="bottom"` places it below instead of beside the graph.

To use an existing subplot, pass `ax=...` and `panel=False`; its layout is left
alone. For an outside report, create the parent figure with
`plt.subplots(..., layout="constrained")`. Artist dictionaries
`curve_kwargs`, `point_kwargs`, and `band_kwargs` are passed to Matplotlib's
`plot`, `errorbar`, and `fill_between`, respectively. Adding markers, changing
limits, or exporting PDF/SVG/PNG remains ordinary Matplotlib code.

## Plot Diagnostics Without Repeating Fits

```python
from scientificfitting import (plot_contour, plot_diagnostics, plot_profile,
                              plot_profile_matrix, plot_residuals)

fig, axes = plot_diagnostics(result, kinds=("residual", "pull"), xlabel="x / mm")
fig.savefig("calibration_diagnostics.pdf")
plt.close(fig)

# These consume the completed scans above; no additional minimizations.
fig, ax = plot_profile(interval.profile_result, delta_max=5)
plt.close(fig)
fig, ax = plot_contour(pair)
plt.close(fig)
fig, axes = plot_profile_matrix(matrix)
fig.savefig("calibration_profiles.pdf")
plt.close(fig)
```

`plot_residuals(..., kind="residual" | "pull" | "ratio")` returns one axis;
`plot_diagnostics` returns a selectable stack. Both use the same numerical
helper as the Julia renderer. Pulls use only stored **observation** residuals,
not appended parameter-prior terms. With covariance they are whitened
coordinates, not individual measurement discrepancies divided by marginal
errors. The shaded unit bands are reference guides, not calibrated coverage
intervals for fitted residuals. Ratios are undefined at zero predictions:
the single-plot call raises an error; the stack marks that panel unavailable
while retaining the useful residual panels.

Profile plots compare the refitted cost with the local covariance parabola.
Contour plots distinguish filled profile regions from dashed covariance
ellipses, with explicit one- versus two-parameter thresholds in the legend.
Matrix plots put profiles on the diagonal, joint regions below it, and local
correlations above it. Failed costs remain gaps; invalid local covariance is
marked unavailable. Negative cost differences, which can reveal a better
minimum, are not reset to zero. `panel_status="none"` hides the compact status
label but does not alter results; detailed findings remain on `scan.diagnostics`.

The complete nonlinear decay example in `examples/python/numpy_matplotlib.py`
exports both styles, with and without panels, plus residual and three-parameter
profile plots. It fits once and computes its scans once before rendering.

## Installation And Packaging

The wheel contains Python code, a small Julia bridge, dependency metadata, and
the MIT license. It does not contain a second numerical implementation, a Julia
runtime, or Makie. [JuliaPkg](https://github.com/JuliaPy/pyjuliapkg) selects a
compatible Julia executable or downloads one and installs the numerical core.
The core is constrained to `~0.2.0`; an explicitly overridden Julia environment
is checked before the bridge loads. Package import alone does not start Julia.
Compatible Julia 0.2.x bugfix releases do not require a new Python wheel;
the wrapper's version need not advance with every core patch.

On macOS ARM64 with Python 3.12.4, JuliaCall 0.9.35, and Julia 1.12.7, a fresh
Python process using an already installed environment needed **14-23 s** to
import the backend and obtain its first fit result. This covers one small
analytic case per fit family plus the in-place interface, each in its own
process. A new fit in the same process took **1.1-3.3 ms**. Shared callback types
and a small Julia precompile workload reduce compilation, without starting
Julia at Python import. See [Performance](performance.md#Python-Startup) for
the cases, reproducible probe, and timing boundaries. These are local
observations, not latency guarantees for larger models or other machines.

First-ever use also downloads Julia/dependencies and builds caches. A clean
provisioning test before these precompile changes took about five minutes; that
is not a new cold-install timing. The current wheel is about 30 KB, but `du`
reported 798 MiB for Julia and about 460 MiB for packages/caches, excluding Python
and optional plotting dependencies. Roughly 96 MiB is ScientificFitting's own
compiled cache; building it took about 53 s. Precompilation moves work to setup
and uses disk space; it does not remove the runtime's installation cost.

The same sources build with standard tools:

```bash
python -m pip install build
python -m build python --outdir python/dist
conda build python/recipe --override-channels -c conda-forge --no-anaconda-upload
```

Use a conda-forge build environment with
[Conda's libmamba dependency solver](https://docs.conda.io/projects/conda/en/26.7.x/user-guide/concepts/conda-performance.html),
as in the Python CI workflow.
The Conda recipe reads the Python version and requirements from `pyproject.toml`;
only the dependency name changes from `juliacall` to conda-forge's `pyjuliacall`.
Matplotlib and SciPy remain optional. No post-install scripts modify the user's
Julia installation. This recipe is a local build, **not** an existing conda-forge
release. Its installed-package reference fits also pass on macOS ARM64 with
Conda Python 3.14.7 and NumPy 2.5.3, reusing the managed Julia runtime above.

`python/tests/check_install.py` checks a base wheel in a fresh environment without
Matplotlib or SciPy: fitted coefficients/covariance against linear algebra,
profile costs, Poisson estimates, and text diagnostics. Before registration,
`--source /path/to/checkout` selects the pending Julia core without providing a
Julia executable. The Python CI workflow runs installed-wheel tests on Linux,
macOS, and Windows; a workflow definition is not evidence that every platform
has passed. Its reports retain first-fit and repeat-fit times.
They also record Pkg's resolved source, tracking mode, and tree hash. The check
requires the requested checkout with `--source`; without it, retained local-path
or repository overrides are rejected rather than counted as registry installs.
Wheel checks compare every packaged Python/Julia runtime file with its source,
including the wheel rebuilt from the source archive.

Release order matters: register the Julia 0.2 core first, then repeat the wheel
installation check **in a new environment without** `--source`, and only then
publish Python artifacts. A source selection persists in its JuliaPkg environment;
omitting the argument on a later run does not turn it into a registry test.
There is deliberately no fallback to the incompatible 0.1 core.
