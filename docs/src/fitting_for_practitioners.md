# Fitting for Practitioners

Use this page when you have measurements in front of you and need to decide
what to fit, which uncertainty model belongs to the experiment, and whether the
result is trustworthy. It is a decision guide, not a substitute for the
[mathematical derivations](statistical_foundations.md) or the complete
[API reference](api.md).

The central rule is simple: choose the probability model before choosing an
optimizer. A solver can locate a minimum. It cannot make an inappropriate cost
function scientifically meaningful.

## 1. Identify What Was Observed

Start from the measurement process, not from the curve shape.

- **Continuous y values with known, independent standard uncertainties:** use
  Gaussian residuals with `fit_model(...; sigma_y=...)`. Verify that the quoted
  values are standard deviations and the points are genuinely independent.
- **Continuous x and y values with small measurement uncertainties:** use
  locally propagated Gaussian uncertainty with
  `fit_model(...; sigma_x=..., sigma_y=...)`. Verify that the model is smooth
  and approximately linear across each x uncertainty.
- **A Gaussian data vector with shared noise:** use `cov_y=...` or a structured
  `WhiteningOperator(...)`. The covariance must represent the actual shared
  readout, baseline, or calibration mechanism.
- **Counts at known x positions:** use `fit_poisson_model(...)`. The model must
  predict strictly positive expected counts at every fitted point.
- **Counts collected into bins:** use `fit_histogram_model(...)` or
  `fit_histogram_density(...)`. Integrate the model over each bin instead of
  sampling only at its center.
- **Individual samples from a distribution:** use `fit_unbinned_model(...)` and
  verify that the density is normalized on the observation domain used by the
  experiment.
- **Individual samples when the total event yield also carries information:**
  use `fit_extended_unbinned_model(...)`. The model must describe both shape
  and expected event count.

If the row is unclear, write the measurement equation first. For a Gaussian
experiment it is often

```math
y_i = f(x_i,p) + \epsilon_i,
\qquad
\epsilon \sim \mathcal{N}(0,V).
```

That one line states the model prediction, the random quantity, and the
covariance that defines which residual patterns are plausible.

## 2. Build The Uncertainty Model

### Independent y uncertainty

When ``V`` is diagonal, each residual is divided by its own standard
uncertainty:

```math
\chi^2(p)
=
\sum_i
\left(
\frac{y_i-f(x_i,p)}{\sigma_{y,i}}
\right)^2.
```

The uncertainty scale changes the scientific interpretation even when the
points and curve are unchanged. Two residuals of ``+0.2`` and ``-0.2``
contribute ``\chi^2=8`` when ``\sigma=0.1``, but only ``\chi^2=2`` when
``\sigma=0.2``. An uncertainty is therefore not a plotting decoration. It is
part of the model being fitted.

Use `sigma_y` only when the entries are standard uncertainties and the
off-diagonal covariances are negligible. Heteroskedastic values are expected:
every point may have a different uncertainty.

### Correlated uncertainty

For correlated Gaussian measurements, JuFitter evaluates

```math
\chi^2(p)=r(p)^T V^{-1}r(p),
\qquad
r(p)=y-f(x,p).
```

Consider two measurements with standard uncertainty ``0.1`` and correlation
``\rho=0.5``:

```math
V = 0.1^2
\begin{pmatrix}
1 & 0.5\\
0.5 & 1
\end{pmatrix}.
```

A common residual ``r=(0.1,0.1)`` then gives ``\chi^2=4/3``. The opposing
pattern ``r=(0.1,-0.1)`` gives ``\chi^2=4``. Positive correlation says that two
points moving together is more plausible than the same points moving in
opposite directions. Replacing this matrix by its diagonal would assign
``\chi^2=2`` to both patterns and erase that experimental information.

Use `cov_y` for a moderate dense covariance. JuFitter factorizes the matrix and
solves linear systems; it does not form ``V^{-1}`` explicitly. Dense covariance
still requires roughly ``O(n^2)`` memory and ``O(n^3)`` factorization. For a
large time series or detector vector with known structure, use a verified
`WhiteningOperator` rather than materializing a dense matrix.
Verify a custom whitening operation and its covariance log-determinant against
a small dense reference before applying it to a large dataset.

### External systematic uncertainty

Not every systematic effect belongs in `cov_y`. Suppose all measurements use a
gain ``g=1.000\pm0.015``. If the model can contain ``g`` explicitly, fit it as a
nuisance parameter with a Gaussian `ParameterPrior`. The fit then propagates
the shared gain uncertainty into every parameter that depends on it.

Use a `FixedParameter` without uncertainty only for a quantity that is treated
as exact in the stated analysis. Fixing an uncertain calibration constant hides
its contribution instead of propagating it. Likewise, do not reuse information
from the fitted dataset as a prior; that would count the same evidence twice.
An uncertainty attached to `FixedParameter` is retained for reporting; it does
not propagate into the fitted covariance. Use a free nuisance parameter with a
prior when that uncertainty must affect the fit result.

### X uncertainty

For a smooth one-dimensional model and small x uncertainty, JuFitter uses the
first-order effective variance

```math
\sigma_{\mathrm{eff},i}^2(p)
=
\sigma_{y,i}^2
+
\left(
\frac{\partial f(x_i,p)}{\partial x}
\right)^2
\sigma_{x,i}^2.
```

For correlated x and y measurements, the same linearization is written in
matrix form:

```math
V_{\mathrm{eff}}(p)
=
V_y + D(p)V_xD(p)^T,
\qquad
D_{ij}(p)=\delta_{ij}
\frac{\partial f(x_i,p)}{\partial x_i}.
```

This has a direct experimental meaning. A common calibration shift in two x
values is represented by an off-diagonal entry in ``V_x``; multiplication by
``D`` converts that shared horizontal displacement into the corresponding
shared displacement of the model predictions. `cov_x` and `cov_y` supply the
two matrices. `sigma_x` and `sigma_y` are their diagonal special cases.

Because this covariance depends on the fitted parameters, `cost=:auto` uses the
full Gaussian likelihood cost on the ``-2\log L`` scale, including its
log-determinant term. This is more than replacing an error bar after the fit.

The approximation is appropriate when the model is locally smooth and nearly
linear over each x-error interval. Large x errors, a sharp threshold, or a
strongly curved model may require an explicit errors-in-variables or latent-x
model instead. A small `sigma_x` keyword cannot make the linearization exact.

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

Do not refit merely to change a label, add a threshold, or switch the side
panel. Plot extensions operate on the existing `FitResult`; see
[Plotting and Customization](plotting_design.md).

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

A smooth curve only answers part of item 2. It does not establish the other
four.

### Chi-square depends on ndf

For a regular Gaussian fit with ``N`` independent observations and ``k`` free
parameters,

```math
\mathrm{ndf} \approx N-k.
```

Independent Gaussian priors or correlated parameter constraints add genuine
external observations to this count. Fixed parameters are not free parameters.
When the assumptions are valid,

```math
E[\chi^2]=\mathrm{ndf},
\qquad
\operatorname{sd}(\chi^2)=\sqrt{2\,\mathrm{ndf}}.
```

This is why there is no universal acceptable interval such as
``0.5<\chi^2/\mathrm{ndf}<1.5``. At ``\mathrm{ndf}=10``, the natural standard
deviation of the ratio is about ``0.45``. At ``\mathrm{ndf}=1000``, it is only
about ``0.045``. The same ratio can therefore be ordinary in a small fit and
decisive evidence of mismatch in a large one.

JuFitter reports the upper-tail probability

```math
p = P\!\left(\chi^2_{\nu}\geq\chi^2_{\mathrm{observed}}\right).
```

A very small value means that residuals this large would be rare if the model
and uncertainty assumptions were correct. A value very close to one can signal
overestimated uncertainties, ignored correlations, smoothing, or selection.
The p-value is not the probability that the model is true, and it should not be
used as a pass/fail switch without inspecting residuals and the measurement
process.

For non-Gaussian likelihoods, an asymptotic chi-square interpretation may be
unavailable or inappropriate. Use likelihood-specific diagnostics, simulation,
or a parametric bootstrap rather than forcing a Gaussian p-value onto the
problem.

### Residuals and pulls locate the failure

For independent Gaussian data, a pull is

```math
u_i=\frac{y_i-f(x_i,\hat p)}{\sigma_i}.
```

Plausible pulls fluctuate around zero with a scale near one. With correlated
data, JuFitter uses whitened residuals so the covariance is accounted for.
Inspect the original residuals as well: whitening tests statistical scale,
while the data-space residuals show where a physical pattern occurs.

Runs of one sign, large neighboring correlations, or a sinusoidal residual
pattern are often more informative than the largest single pull. They point to
missing physics, a time-dependent baseline, an omitted resonance, or an
incorrect covariance model.

## 5. Use Diagnostics As Triage

Starting from an existing result:

```julia
details = diagnose(result)
dashboard = diagnostic_dashboard(result)
```

`diagnose` returns structured findings. Each finding carries a severity,
numeric evidence, and a recommended action. `diagnostic_dashboard` summarizes
the same findings; it does not invent a second set of statistical rules.

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

The dashboard status is deliberately operational:

- `ok - no immediate issue`: none of the current checks produced a warning or
  critical finding;
- `review - inspect diagnostics`: at least one warning needs interpretation;
- `critical - fix before use`: a critical defect must be corrected before the
  result is used for conclusions.

`ok` does not certify the model. No finite checklist can detect missing physics
that leaves the tested residuals apparently benign.

## 6. Replace Local Errors When The Cost Is Not Parabolic

The reported parameter covariance is a local curvature approximation. In one
parameter it assumes

```math
\Delta C(p_i)
\approx
\left(
\frac{p_i-\hat p_i}{\sigma_i}
\right)^2,
```

where JuFitter uses a cost convention compatible with ``-2\log L``. This gives
compact symmetric errors, but it can fail near bounds, with weak data, in a
nonlinear model, or for asymmetric likelihoods.

### One parameter: profile interval

A profile fixes one parameter and refits every remaining nuisance parameter:

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

Adaptive refinement is already enabled by default in `profile_interval`. A
manual scan can request the same threshold-focused refinement:

```julia
prof = profile(result, 1; adaptive=true, threshold=1.0)
```

### Two parameters: contour geometry

A contour fixes two parameters on a grid and refits all remaining nuisance
parameters:

```julia
cont = JuFitter.contour(
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
thresholds. It improves resolution but does not repair an inadequate scan
range, failed refits, or an unidentifiable model.

## 7. Report What The Analysis Actually Supports

Before reporting a fit, record:

- the measured quantities and units;
- the model formula and meaning of every fitted parameter;
- the uncertainty or likelihood model, including correlations and external
  constraints;
- the cost convention and number of free parameters;
- convergence and diagnostic status;
- residual or pull evidence;
- local symmetric errors or profile-based asymmetric intervals, as justified;
- important bounds, assumptions, and known limitations.

For complete worked examples, follow the gallery progression from
[Linear Calibration](gallery/linear_calibration.md) through
[XY Uncertainties](gallery/xy_uncertainties.md),
[Full Covariance](gallery/full_covariance.md), and
[Constraints and Profiles](gallery/constraints_profiles.md). The
[Statistical Foundations](statistical_foundations.md) page derives the cost
functions and uncertainty conventions used here.
