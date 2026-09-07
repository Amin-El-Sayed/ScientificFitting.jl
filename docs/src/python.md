# Python Interface (Development Preview)

The `python/` directory contains an unreleased wrapper, not yet a PyPI or conda
release. It uses NumPy models, the same Julia numerical core, and optional
native Matplotlib plots. It neither installs nor loads Makie. From this
checkout, in a Python 3.10+ virtual environment:

```bash
python -m pip install -e './python[plot,test]'
python python/develop.py
python examples/python/numpy_matplotlib.py
python -m pytest python/tests
```

The development step selects this checkout because the registered v0.1.2 core
does not contain the new finite-derivative API. JuliaCall manages Julia and
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
ax.legend()
fig.savefig("calibration.pdf")
```

Models receive a read-only NumPy coordinate array and named float parameters.
The keys of `p0` must match the model's parameter names; mapping order sets
the order of result arrays, not how parameters bind to the model. Models must
be deterministic and smooth around evaluated points. Finite differences are
used consistently, including when bounds or profiles trigger a different
solver. Optional analytic `jacobian` and `x_derivative` callbacks use the same
Python argument convention.
The default solver `tol=1e-6` accounts for differenced-gradient noise;
`tol` remains configurable and is not a bound on parameter error.
`maxiters` limits each solver run. `result.converged` and the report reflect
the actual solver status, including when the iteration limit is too small.

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

Not yet covered by the Python facade: the full diagnostic/report-panel
plotting suite. The Julia APIs remain available; this preview is
not a claim of v0.2 feature parity. Cross-platform clean-install checks and
release dependency pins are also still required before publishing a wheel.

For non-Gaussian regression, `fit_likelihood_model` accepts a vectorized
`logprob(y, prediction, **parameters)` returning one normalized log density
or log probability mass per independent observation. For example, SciPy's
`stats.t.logpdf(y, df=4, loc=prediction, scale=scales)` describes Student-t
errors; the scale is not their standard deviation. Capture known per-point
scales in the callback. Unlike `fit_custom`, no manual likelihood summation or
observation count is needed. Both callbacks are batched; scalar-density
quadrature helpers cross the language boundary once per evaluation point and
can be slower. Dependent observations require a joint likelihood.

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

For two profiled parameters the default levels `[2.30, 6.18]` differ from
the one-parameter thresholds `[1, 4]`. They represent the usual asymptotic
68.27% and 95.45% joint regions, not marginal errors or posterior probabilities.
