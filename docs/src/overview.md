# Reference Map

Use this page when the scientific model is clear but the JuFitter entry point
is not. It maps a measurement workflow to the relevant fit, uncertainty,
diagnostic, and output APIs. Exact signatures and defaults are in the
[API Reference](api.md); statistical assumptions are derived in
[Statistical Foundations](statistical_foundations.md).

## Choose The Fit From The Sampling Process

The data-generating process determines the objective. Plot style and optimizer
choice do not.

| Observation model | Entry point | Result |
|---|---|---|
| Numeric ``x`` and ``y`` with Gaussian measurement uncertainty | `fit_model` | `FitResult` |
| The same Gaussian x-y workflow, fitted and plotted in one call | `fitplot` | named tuple `(result, figure)` |
| Independent integer counts with one expected count per observation | `fit_poisson_model` | `LikelihoodFitResult` |
| Histogram counts with a model that returns expected bin counts | `fit_histogram_model` | `LikelihoodFitResult` |
| Histogram counts with a normalized density integrated over each bin | `fit_histogram_density` | `LikelihoodFitResult` |
| Unbinned independent events from a normalized density | `fit_unbinned_model` | `LikelihoodFitResult` |
| Unbinned events whose total expected rate also depends on the parameters | `fit_extended_unbinned_model` | `LikelihoodFitResult` |
| Gaussian observations addressed by labels or other non-numeric indices | `fit_indexed_model` | `LikelihoodFitResult` |
| Several datasets sharing some or all parameters | `fit_multi_model` | `LikelihoodFitResult` |
| A scalar objective with an application-specific interpretation | `fit_custom` | `LikelihoodFitResult` |

For a reusable low-level object, construct `FitProblem` or
`LikelihoodFitProblem` and call `fit(problem)`. Start with a high-level entry
point unless the problem itself must be stored, modified, or repeatedly
refitted.

## Describe Measurement Uncertainty Once

Each physical uncertainty source should enter the model exactly once. Use the
simplest representation that preserves the correlations present in the
measurement.

| Measurement situation | Representation |
|---|---|
| Independent y uncertainties | `sigma_y` |
| Complete correlated y uncertainty | `cov_y` |
| Independent x uncertainties propagated through ``\partial f/\partial x`` | `sigma_x` |
| Complete correlated x uncertainty | `cov_x` |
| Several named contributions, such as readout noise and gain uncertainty | `ErrorComponent` |
| A complete static covariance with exploitable structure | `WhiteningOperator` |

`sigma_y` and `cov_y` are alternatives, as are `sigma_x` and `cov_x`.
`ErrorComponent`s are actual additive covariance contributions: their names
make sources auditable, and `active=false` excludes a source from a particular
fit without deleting its definition. A `WhiteningOperator` is the complete
observation covariance model, not one more contribution, so it is exclusive
with `sigma_*`, `cov_*`, and `error_components`.

Dense covariance requires approximately ``O(n^2)`` storage and ``O(n^3)``
factorization. For a long series with known structure, represent the complete
static whitening operation directly instead of materializing a dense matrix.
JuFitter can validate the operator interface and finite output; the scientist
remains responsible for the identity ``W^\mathsf{T}W=C^{-1}`` and the supplied
``\log\det C``.

## Separate Parameter Knowledge From Data Uncertainty

Observation uncertainty describes the measured data. Parameter controls encode
physical domains or external information about model parameters.

| Intent | Control | Statistical effect |
|---|---|---|
| Parameter is known exactly for this fit | `FixedParameter` | Removed from the optimizer. |
| Independent external measurement of one parameter | `ParameterPrior` | Adds a Gaussian or split-normal likelihood term. |
| Correlated external measurements of several parameters | `ParameterConstraint` | Adds one multivariate Gaussian term. |
| Parameter is physically restricted to an interval | `bounds` | Restricts the search domain; it is not extra data. |
| Parameters obey a general relation | `ConstraintSpec` | Enforces nonlinear equality or inequality constraints. |

`ConstraintSpec` callbacks receive the complete parameter vector in `p0`
order, including fixed entries. They therefore express the physical relation
directly; JuFitter handles reduced optimizer coordinates internally.

An uncertainty stored in `FixedParameter(index, value, sigma)` appears in the
report and in that parameter's diagonal covariance entry, with zero fitted
cross-covariances. It does not change the objective, make the parameter free,
or propagate uncertainty into free fitted parameters. Use a prior or correlated
parameter constraint when external uncertainty must affect the fitted result.

## Know Which Result You Have

Both result types expose fitted parameters, local covariance, solver status,
fit statistics, and diagnostics:

```julia
result.params
result.param_stderr
result.param_covariance
result.stats
result.diagnostics
result.converged
```

`FitResult` additionally stores the validated x-y problem, fitted model values,
residuals, weighted residuals, and Jacobian. That information makes
`plot_fit(result)` and residual plots reproducible without fitting again.

`LikelihoodFitResult` deliberately has no universal curve representation. A
Poisson count model, an unbinned density, and a custom objective need different
visualizations, so `plot_fit` does not accept `LikelihoodFitResult`. Use the
fitted parameters in a workflow-specific Makie figure; the Poisson and
multi-dataset gallery pages show this pattern.

## Inspect Before Reporting

Given a result returned by one of the entry points above, the shortest
Makie-free inspection path is:

```julia
println(report_text(result))
println(diagnostic_dashboard_text(result))
```

The report records numerical values. The dashboard converts structured
findings into a status and concrete next actions. Inspect at least these
questions in order:

1. Did the optimizer converge to a finite result?
2. Is the goodness-of-fit statistic defined for this sampling model?
3. Do residuals or pulls contain location-dependent structure?
4. Are bounds active or the local covariance poorly conditioned?
5. Are symmetric local errors adequate for the conclusion being reported?

`diagnose(result)` returns the full machine-readable findings.
`diagnostic_dashboard(result)` provides the condensed status, severity counts,
and deduplicated actions used by the text dashboard.

## Escalate From Local Errors To Profiles

`param_covariance` is a local quadratic approximation around the fitted
minimum. It is efficient and often sufficient, but active bounds, strong
correlations, weak data, nonlinear response, or an asymmetric likelihood can
make symmetric errors misleading.

Use a one-parameter profile for an asymmetric interval:

```julia
interval = profile_interval(result, 1; adaptive=true)
diagnose(interval.profile_result)
```

Use a profile matrix when several parameters interact:

```julia
matrix = profile_matrix(result; parameters=[1, 2, 3], adaptive=true)
profile_matrix_triage(matrix)
```

The diagonal scans test the local parabolic approximation. Pairwise contours
test whether the local covariance ellipses describe the actual refitted cost
surface. Default profile and contour thresholds are Wilks-based asymptotic
regions, not exact finite-sample coverage guarantees.

Load CairoMakie only when a figure is needed:

```julia
using CairoMakie

plot_profile_matrix(matrix)
```

Passing the precomputed matrix renders the stored scans and does not rerun the
profile fits.

## Choose The Output Surface

Fitting, text reports, diagnostics, profiles, and profile-matrix computation do
not require Makie. Plotting is an optional CairoMakie extension.

For a Gaussian x-y result:

```julia
using CairoMakie

fig = plot_fit(
    result;
    theme=:screen,
    show_stats=true,
    show_legend=true,
)
```

`theme=:screen` and `:article` select the output context;
`appearance=:light` or `:dark` selects the color scheme independently. Explicit
Makie keyword arguments override the style only for the supplied element.

Add scientific annotations without reconstructing or refitting the model:

```julia
ax = fit_axis(fig)
add_vline!(ax, threshold; label="threshold")
add_curve!(ax, reference_model; xspan=(0, 10), label="reference")
```

For a one-call Gaussian workflow, `fitplot(...; report=:plot | :console |
:both | :none)` controls whether statistics appear in the figure, terminal,
both, or neither.

## Continue From Here

| Need | Page |
|---|---|
| First successful fit and plot | [Quickstart](quickstart.md) |
| Complete scientific workflows | [Gallery](gallery.md) |
| A suspicious fit and a practical next action | [Fitting for Practitioners](fitting_for_practitioners.md) |
| Plot styles, panels, and custom Makie composition | [Plotting and Customization](plotting_design.md) |
| Statistical assumptions and derivations | [Statistical Foundations](statistical_foundations.md) and its focused method chapters |
| Every keyword, field, default, and failure mode | [API Reference](api.md) |
| Internal validation, solver, and result flow | [Backend Design](backend_design.md) |
