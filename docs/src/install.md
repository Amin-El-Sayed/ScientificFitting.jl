# Installation

ScientificFitting declares support for Julia 1.10 and later. CI runs the core
and full-package gates on Julia 1.10 and Julia 1.12. Julia 1.10 is the
compatibility floor.

!!! note "First public release"
    Version 0.1 is an early work-in-progress release. Please report unexpected
    behavior with a minimal reproducer and the complete diagnostic output.

## Install The First Release

Install the registered package from Julia's General registry:

```julia
using Pkg
Pkg.add("ScientificFitting")
```

Then load the numerical core:

```julia
using ScientificFitting
```

This loads fitting, likelihoods, diagnostics, profiles, contours, and text
reports. It does **not** load Makie.

## Work From A Checkout

From a terminal, instantiate the numerical core in the repository root:

```bash
cd /path/to/ScientificFitting
julia --project=. --startup-file=no -e 'using Pkg; Pkg.instantiate()'
```

Start a Julia session in that environment:

```bash
julia --project=.
```

The checked-out repository keeps plotting and documentation dependencies in a
separate environment. Instantiate it when running the gallery or building the
site:

```bash
julia --project=docs --startup-file=no -e 'using Pkg; Pkg.instantiate()'
julia --project=docs examples/gallery/01_quickstart_linear.jl
```

The example writes its figure to the ignored `examples/output/` directory.

For static PNG, PDF, and SVG plots, install CairoMakie in the same environment:

```julia
using Pkg
Pkg.add(["ScientificFitting", "CairoMakie"])

using ScientificFitting
using CairoMakie
```

Keeping CairoMakie optional is deliberate. A batch analysis can fit data,
create reports, and run diagnostics without compiling a graphics stack.

## First Use And Compilation

The first `using ScientificFitting` in a new environment compiles the numerical core.
The first `using CairoMakie` and first rendered figure take longer because Julia
also compiles Makie's layout, text, and rendering methods. Later sessions reuse
the precompile cache unless Julia, package versions, preferences, or the target
environment change.

Do not use the full package test suite to check an installation; it is a slow
release gate. A core-only check is enough:

```bash
julia --project=. --startup-file=no -e 'using ScientificFitting; println("ScientificFitting core ready")'
```

For plotting, run the tracked quickstart example shown above. It exercises the
same API used by the first tutorial and confirms CairoMakie export.

## Python Interface (Development Preview)

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

Not yet covered by the Python facade: complete structured diagnostics,
profile matrices, and the full diagnostic/
report-panel plotting suite. The Julia APIs remain available; this preview is
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

### Large Gaussian Fits

SciPy sparse matrices and arrays are copied as canonical CSC buffers, not
converted to dense arrays. Install `./python[sparse]` to add SciPy if needed.
The same `cov_x`/`cov_y` arguments accept either representation. Static sparse
y covariance is factorized once per fit; parameter-dependent effective
covariance must be refactorized at each parameter point. Sparse factorization
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

## Troubleshooting

| symptom | first check |
| --- | --- |
| `using ScientificFitting` is slow once | Let precompilation finish; this is not fit runtime. |
| Every fresh session recompiles | Reuse the same project and depot; check whether Julia or package versions keep changing. |
| `plot_fit` says the extension is unavailable | Add and load `CairoMakie` before calling plotting functions. |
| PDF or SVG export fails | Verify a minimal CairoMakie figure in the same environment; inspect backend and font errors first. |
| A fit is unexpectedly slow | Check for dense covariance, bounds, constraints, priors, parameter-dependent covariance, or pointwise x-derivatives. These select more general numerical paths. |
| Package versions will not resolve | Confirm Julia is at least 1.10 and instantiate a clean environment rather than mixing incompatible manifests. |

Continue with the [Quickstart](quickstart.md). For package internals and scaling
limits, see [Backend Design](backend_design.md) and [Performance](performance.md).
