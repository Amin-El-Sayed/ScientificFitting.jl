# Fitting for Practitioners

Use this page to decide what to fit, which uncertainty model belongs to the
experiment, and whether the result is trustworthy; the
[mathematical derivations](statistics.md) and the complete
[API reference](api.md) stand behind it. Choose the probability model before
choosing an optimizer.

## 1. Identify What Was Observed

Start from the measurement process, not from the curve shape. The
[entry-point table](api.md#Choose-An-Entry-Point) in the API reference maps
each data and sampling model to its fitting function and the contract the
model must satisfy.

If the table row is unclear, write the measurement equation first. For a
Gaussian experiment it is often

```math
y_i = f(x_i,p) + \epsilon_i,
\qquad
\epsilon \sim \mathcal{N}(0,V).
```

## 2. Build The Uncertainty Model

### Independent y uncertainty

When ``V`` is diagonal, each residual is divided by its own standard
uncertainty ([Gaussian Least Squares](statistics.md#Gaussian-Least-Squares)).
The uncertainty scale is part of the model being fitted: two residuals of
``+0.2`` and ``-0.2`` contribute ``\chi^2=8`` when ``\sigma=0.1``, but only
``\chi^2=2`` when ``\sigma=0.2``.

Use `sigma_y` only when the entries are standard uncertainties and the
off-diagonal covariances are negligible. Heteroskedastic values are expected.

### Correlated uncertainty

For correlated Gaussian measurements, ScientificFitting evaluates the
generalized chi-square through whitened residuals, without forming ``V^{-1}``
explicitly
([Correlated Measurements And Whitening](statistics.md#Correlated-Measurements-And-Whitening)).
Positive correlation makes a common residual pattern more plausible than an
opposing one; a diagonal approximation erases that experimental information.

Use `cov_y` for a moderate dense covariance (roughly ``O(n^2)`` memory,
``O(n^3)`` factorization). For a large time series or detector vector with
known structure, use a
[structured `WhiteningOperator`](statistics.md#Structured-Whitening), verified
against a small dense reference before use at scale.

### External systematic uncertainty

Not every systematic effect belongs in `cov_y`. A shared gain
``g=1.000\pm0.015`` that the model can contain explicitly is fitted as a
nuisance parameter with a Gaussian `ParameterPrior`, which propagates the
shared uncertainty into every parameter that depends on it.

Reserve `FixedParameter` for a quantity treated as exact: fixing an uncertain
calibration constant hides its contribution, and an uncertainty attached to
`FixedParameter` is report metadata outside the fitted covariance
([External Parameter Information](statistics.md#External-Parameter-Information)).
Do not reuse information from the fitted dataset as a prior; that counts the
same evidence twice.

### X uncertainty

For a smooth one-dimensional model and small x uncertainty, ScientificFitting
propagates x errors into a parameter-dependent effective covariance
([Uncertainty In X](statistics.md#Uncertainty-In-X)). `cov_x` and `cov_y`
supply the two covariance matrices; `sigma_x` and `sigma_y` are their diagonal
special cases. Because the effective covariance depends on the fitted
parameters, `cost=:auto` uses the full Gaussian likelihood cost, including its
log-determinant term.

The linearization requires a model nearly linear over each x-error interval;
large x errors, a sharp threshold, or strong curvature need an explicit
errors-in-variables or latent-x model, and a small `sigma_x` keyword cannot
make the linearization exact.

## 3. Fit Once, Then Inspect One Result

The explicit problem API makes the statistical inputs reviewable:

```julia
problem = FitProblem(
    model,
    x,
    y;
    p0=[1.0, 0.0],
    sigma_y=sigma_y,
)

result = fit(problem)
```

`fit_model(model, x, y; ...)` is the shorter wrapper around the same
construction. Both return one result object that feeds reports, diagnostics,
profiles, contours, and optional Makie plots:

```julia
println(report_text(result))
println(diagnostic_dashboard_text(result))

# Plotting remains optional for headless or server-side work.
using CairoMakie
fig = plot_fit(result; xlabel="x", ylabel="y")
```

The two low-level problem types separate observation contracts, not levels of
sophistication:

| problem type | what it stores | typical constructors |
| --- | --- | --- |
| `FitProblem` | x-y observations, model predictions, and a Gaussian residual covariance | `fit_model`, direct `FitProblem` |
| `LikelihoodFitProblem` | a complete ``-2\log L`` objective, optional goodness-of-fit statistic, and observation count | Poisson, histogram, unbinned, indexed, multi-dataset, and `fit_custom` helpers |

Direct construction is useful when a problem must be stored or inspected
explicitly.

Do not refit merely to change a label or a plot style; plot extensions operate
on the existing `FitResult` (see
[Plotting and Customization](plotting_design.md)).

## 4. Read The Fit In A Defensible Order

Read a result in this order:

1. **Validity:** did the optimizer converge, are the parameters finite, and is
   the fitted point a minimum with usable local curvature?
2. **Residual structure:** are discrepancies random, or do they form runs,
   trends, oscillations, or isolated extreme points?
3. **Goodness of fit:** is the observed cost plausible under the stated
   probability model?
4. **Parameter geometry:** are parameters strongly correlated, at bounds, or
   described poorly by a local quadratic approximation?
5. **Scientific interpretation:** are the fitted values physically meaningful,
   and does the uncertainty include every relevant source?

### Chi-square depends on ndf

Degrees-of-freedom counting — auxiliary Gaussian terms, correlated
constraints, fixed parameters — is defined in
[Degrees Of Freedom](statistics.md#Degrees-Of-Freedom). The natural spread of
``\chi^2/\mathrm{ndf}`` shrinks with sample size, so there is no universal
acceptable interval for the ratio; the sample-size-aware statement is the
upper-tail p-value ([Goodness Of Fit](statistics.md#Goodness-Of-Fit)).
Non-Gaussian likelihoods may lack an asymptotic chi-square interpretation; use
likelihood-specific diagnostics, simulation, or a parametric bootstrap.

### Residuals and pulls locate the failure

[Residuals And Pulls](statistics.md#Residuals-And-Pulls) defines pulls and the
whitened residuals used for correlated data. Whitening tests statistical
scale; the data-space residuals show where a physical pattern occurs.

Runs of one sign, large neighboring correlations, or a sinusoidal residual
pattern are often more informative than the largest single pull: they point to
missing physics, a time-dependent baseline, an omitted resonance, or an
incorrect covariance model.

## 5. Use Diagnostics As Triage

Starting from an existing result:

```julia
details = diagnose(result)
dashboard = diagnostic_dashboard(result)
```

`diagnose` returns structured findings — severity, numeric evidence,
recommended action — and `diagnostic_dashboard` summarizes them.

- **Optimizer did not converge:** the reported point is not established as a
  minimum. First rescale parameters, improve starting values, simplify the
  model, or use multistart.
- **Non-positive ndf:** reduced chi-square and its p-value do not test fit
  quality. Add independent observations or reduce the number of free
  parameters.
- **Large pulls or a long same-sign run:** a data region disagrees coherently
  with the model. Inspect the reported point/x interval and the raw measurement
  there.
- **Large chi-square or small p-value:** model, uncertainty, correlation, or
  optimizer assumptions conflict with the data. Inspect residual structure
  before inflating uncertainties.
- **Unusually small chi-square or p-value near one:** the data are less variable
  than the uncertainty model predicts. Check smoothing, averaging, duplicated
  information, and ignored correlations.
- **Strong parameter correlation:** the data constrain a combination better
  than individual parameters. Reparameterize or inspect a two-parameter
  contour.
- **Ill-conditioned covariance or Hessian:** local symmetric errors are
  numerically or statistically fragile. Rescale and inspect profiles/contours
  before reporting intervals.
- **Active bound:** the local Gaussian approximation is truncated. Decide
  whether the bound is physical, then use a profile interval.

The dashboard status values are:

- `ok - no immediate issue`: none of the current checks produced a warning or
  critical finding;
- `review - inspect diagnostics`: at least one warning needs interpretation;
- `critical - fix before use`: a critical defect must be corrected before the
  result is used for conclusions.

`ok` does not certify the model: no finite checklist can detect missing
physics that leaves the tested residuals apparently benign.

## 6. Replace Local Errors When The Cost Is Not Parabolic

The reported parameter covariance is a
[local curvature approximation](statistics.md#Local-Parameter-Covariance) on
the ``-2\log L`` [cost scale](statistics.md#The-Cost-Convention). It gives
compact symmetric errors, but it can fail near bounds, with weak data, in a
nonlinear model, or for asymmetric likelihoods.

### One parameter: profile interval

A profile fixes one parameter and refits every remaining nuisance parameter;
the refit, not a frozen slice, defines the uncertainty
(see [Profiles And Contours](statistics.md#Profiles-And-Contours)):

```julia
interval = profile_interval(result, 1; threshold=1.0)
prof = interval.profile_result
println(diagnose(prof; local_sigma=result.param_stderr[1]))

fig = plot_profile(
    prof;
    local_sigma=result.param_stderr[1],
    threshold_label="68.3% profile threshold",
)
```

Under regular one-parameter likelihood assumptions, ``\Delta C=1`` corresponds
approximately to a 68.3% interval. Read the profile shape, not only the
crossings:

- a symmetric parabola supports the local standard error;
- different left and right crossings require an asymmetric interval;
- one missing crossing means the scan range is too narrow or the parameter is
  not bounded on that side;
- clipping at a physical bound calls for a one-sided interpretation;
- a second minimum means the local covariance describes only one basin.

Adaptive refinement is enabled by default in `profile_interval`. A manual scan
can request the same threshold-focused refinement:

```julia
prof = profile(result, 1; adaptive=true, threshold=1.0)
```

### Two parameters: contour geometry

A contour fixes two parameters on a grid and refits all remaining nuisance
parameters:

```julia
cont = ScientificFitting.contour(
    result,
    1,
    2;
    adaptive=true,
    levels=[2.30, 6.18],
)

println(diagnose(
    cont;
    local_covariance=result.param_covariance,
    local_center=result.params[[1, 2]],
))

fig = plot_contour(
    cont;
    local_covariance=result.param_covariance,
    local_center=result.params[[1, 2]],
    xlabel="parameter 1",
    ylabel="parameter 2",
)
```

For two parameters under regular likelihood assumptions, ``\Delta C=2.30`` and
``6.18`` are the approximate 68.3% and 95.4% joint-confidence thresholds. In
the resulting plot, the filled regions show the profiled cost; the dashed line
overlay shows the local parabolic covariance approximation.

Use the geometry as a decision tool:

- matching ellipses support the local covariance approximation;
- a tilted ellipse shows correlation but can still be locally Gaussian;
- a curved or banana-shaped region means symmetric marginal errors hide the
  joint geometry;
- a contour cut by a bound requires bounded or one-sided intervals;
- an open contour means the scan or the data do not close the confidence
  region.

Adaptive contour refinement concentrates expensive refits near requested
thresholds; it does not repair an inadequate scan range, failed refits, or an
unidentifiable model.

## 7. Report What The Analysis Actually Supports

Before reporting a fit, work through the
[reporting checklist](statistics.md#Reporting-Checklist) in the statistics
reference.

For complete worked examples, follow the gallery progression from
[Linear Calibration](quickstart.md) through
[XY Uncertainties](gallery/xy_uncertainties.md),
[Full Covariance](gallery/full_covariance.md), and
[Constraints and Profiles](gallery/constraints_profiles.md).
