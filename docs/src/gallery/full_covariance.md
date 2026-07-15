# Full Covariance

This workflow fits an exponential decay when neighboring measurements share
readout noise. The essential point is that correlated uncertainty changes the
statistical problem. A dense covariance matrix is not a prettier way to draw
error bars; it tells the fit which residual patterns are plausible together.

```@raw html
<img class="jufitter-plot jufitter-plot-light" data-jufitter-plot-group="full-covariance" data-jufitter-plot-style="lab" src="../assets/gallery/full_covariance_decay_lab_light.png" alt="Full covariance exponential fit in lab style">
<img class="jufitter-plot jufitter-plot-dark" data-jufitter-plot-group="full-covariance" data-jufitter-plot-style="lab" src="../assets/gallery/full_covariance_decay_lab_dark.png" alt="Full covariance exponential fit in lab dark style">
<img class="jufitter-plot jufitter-plot-light" data-jufitter-plot-group="full-covariance" data-jufitter-plot-style="modern" src="../assets/gallery/full_covariance_decay_modern_light.png" alt="Full covariance exponential fit in modern style">
<img class="jufitter-plot jufitter-plot-dark" data-jufitter-plot-group="full-covariance" data-jufitter-plot-style="modern" src="../assets/gallery/full_covariance_decay_modern_dark.png" alt="Full covariance exponential fit in modern dark style">
<img class="jufitter-plot jufitter-plot-light" data-jufitter-plot-group="full-covariance" data-jufitter-plot-style="article" src="../assets/gallery/full_covariance_decay_article_light.png" alt="Full covariance exponential fit in article style">
<img class="jufitter-plot jufitter-plot-dark" data-jufitter-plot-group="full-covariance" data-jufitter-plot-style="article" src="../assets/gallery/full_covariance_decay_article_dark.png" alt="Full covariance exponential fit in article dark style">
```

## Question

A decaying signal is sampled repeatedly with the same readout chain. Neighboring
points are not independent because baseline drift, electronics, or smoothing in
the acquisition system affects several samples at once. The scientific question
is:

```math
y(t) = A e^{-\lambda t} + C,
```

with parameter uncertainties that account for correlated residuals.

If the correlations are ignored, the fit can look artificially precise. Many
small residuals in the same direction do not provide as much independent
information as the same number of uncorrelated residuals.

## Data and Covariance

The documentation example uses a controlled decay dataset with a known
correlation length. Think of each point as containing two kinds of uncertainty:

- a local statistical part, such as readout noise that changes independently
  from sample to sample,
- a shared part, such as baseline drift, temperature drift, filtering, or
  electronics noise that affects nearby samples in similar directions.

If two points are close in acquisition order, they see more of the same shared
disturbance. If they are far apart, the shared disturbance has mostly changed.
The compact model used here is

```math
V_{ij} = \sigma^2 \exp\!\left(-\frac{|i-j|}{\ell}\right).
```

The parameter ``\ell`` is a correlation length in sample-index units. At
``|i-j|=0`` the covariance is the variance ``\sigma^2``. At separations much
larger than ``\ell``, the covariance becomes small and the points behave nearly
independently.

This is a useful first model when the acquisition chain has a finite memory:
baseline estimates, smoothing filters, thermal drift, or slowly varying
electronics offsets all tend to create nearby residuals with the same sign. It
is not the right model for every experiment. A constant calibration scale error,
for example, is a systematic effect that should usually be propagated to the
final result or represented by a separate nuisance parameter, not hidden inside
this short-range covariance matrix.

A quick sketch of the idea is:

```math
\begin{array}{c|cccc}
|i-j| & 0 & 1 & 2 & \ldots \\
\hline
V_{ij}/\sigma^2 & 1 & e^{-1/\ell} & e^{-2/\ell} & \ldots
\end{array}
```

## Model and Cost

With independent y errors, the weighted chi-square is

```math
\chi^2 = \sum_i \left(\frac{y_i - f(t_i,p)}{\sigma_i}\right)^2.
```

With a full covariance matrix ``V``, the Gaussian cost uses the quadratic form

```math
\chi^2 = r^\mathsf{T} V^{-1} r,
\qquad
r_i = y_i - f(t_i,p).
```

JuFitter does not form an explicit inverse for this calculation. It factors the
covariance and evaluates the whitened residuals, which is numerically safer:

```math
V = L L^\mathsf{T},
\qquad
\chi^2 = \lVert L^{-1} r \rVert^2.
```

The dense path is appropriate for small and medium datasets where the full
matrix is meaningful and affordable. It has ``O(n^2)`` memory cost and an
``O(n^3)`` factorization cost, so huge correlated datasets should eventually use
structured covariance or custom whitening operators rather than dense matrices.

## Fit

This is the complete code for the documentation example:

```julia
using JuFitter
using LinearAlgebra

# Measured decay samples. The values are listed explicitly because the fit
# should read like an analysis notebook, not like a data simulator.
t = [0.0, 0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 0.8333,
     0.9524, 1.0714, 1.1905, 1.3095, 1.4286, 1.5476, 1.6667, 1.7857,
     1.9048, 2.0238, 2.1429, 2.2619, 2.3810, 2.5]
y = [2.25155, 2.00932, 1.79547, 1.60660, 1.44018, 1.29422, 1.16695,
     1.05651, 0.96079, 0.87740, 0.80381, 0.73758, 0.67658, 0.61929,
     0.56494, 0.51354, 0.46572, 0.42253, 0.38510, 0.35431, 0.33052,
     0.31346]

# Neighboring samples share a readout chain, so the uncertainty model is a
# covariance matrix instead of independent error bars.
n = length(t)
base_sigma = 0.055
correlation_length = 2.3

cov_y = [
    base_sigma^2 * exp(-abs(i - j) / correlation_length)
    for i in 1:n, j in 1:n
]

decay_model(t, p) = @. p[1] * exp(p[2] * t) + p[3]

result = fit_model(
    decay_model,
    t,
    y;
    p0=[1.5, -0.7, 0.0],
    cov_y=cov_y,
)

amplitude, decay_rate, offset = result.params
sigma_amplitude, sigma_decay_rate, sigma_offset = result.param_stderr

println("A = ", amplitude, " +/- ", sigma_amplitude)
println("lambda = ", -decay_rate, " +/- ", sigma_decay_rate)
println("C = ", offset, " +/- ", sigma_offset)
println(diagnostic_dashboard_text(result))
```

```@raw html
<div class="jufitter-cell-output">
<div class="jufitter-cell-output-label">Real output (abridged)</div>
<pre>A = 2.0987750476512956 +/- 0.09109589472907365
lambda = 0.9956123089992535 +/- 0.10600822134107704
C = 0.14759627451357205 +/- 0.08708023506916575

Fit diagnostic dashboard
status = review - inspect diagnostics
critical = 0, warning = 5, info = 0
5 warning(s). Inspect before trusting uncertainties or conclusions.

Next actions:
  1. Use a covariance model, inspect acquisition order/time dependence, or fit a model with the missing systematic component.
  2. The uncertainties may be too large, correlations may be ignored, or the data may not be independent.
  3. Inspect this acquisition interval for drift, missing model structure, a calibration offset, or correlated uncertainty.
  4. Look for missing curvature, drift, hysteresis, time dependence, or an incorrect independent variable transformation.
  5. Uncertainties may be overestimated, correlations may be ignored, or the data may have been smoothed/averaged.</pre>
</div>
```

The gallery figure is generated from the same result. It shows a 1-sigma
prediction band, not only the model uncertainty. The side report lists the
parameters and goodness-of-fit values computed from the full covariance model.

## Diagnostics

For a full-covariance fit, inspect more than the parameter table:

- Verify that `cov_y` is symmetric and positive definite. Invalid covariance
  matrices are scientific input errors, not style warnings.
- Compare the fitted parameters with a diagonal-error fit only as a diagnostic.
  Agreement in central values does not mean the uncertainties are equivalent.
- Inspect residuals in acquisition order. Long runs of same-sign residuals are
  less surprising under a correlated covariance model than under an independent
  model.
- Check whether the covariance model is justified by the instrument or
  acquisition process. A mathematically valid matrix is not automatically a
  physically valid uncertainty model.

The diagnostic dashboard should be treated as a first pass. If it reports a
review condition, inspect the residual structure and the covariance assumptions
before publishing the fit.

## Interpretation

The fitted amplitude ``A`` and offset ``C`` describe the scale and baseline of
the signal. The decay rate is reported as ``\lambda = -p_2`` because the model
uses ``A\exp(p_2 t)+C`` with a negative fitted exponent.

For the dataset shown here, the fit gives approximately

```math
A = 2.10 \pm 0.09,\qquad
\lambda = 1.00 \pm 0.11,\qquad
C = 0.15 \pm 0.09.
```

The parameter covariance now includes the assumed point-to-point correlation.
That does not make the uncertainty automatically larger in every parameter, but
it changes which residual patterns count as independent evidence. This is the
main reason to use the full matrix instead of independent error bars.

For the controlled dataset shown here, the dashboard usually asks for review
because the residual pattern is smooth. That is scientifically useful: a
covariance matrix can account for known correlation, but it does not prove that
the exponential model captures every systematic feature of the data.

## What Can Go Wrong

Do not create a dense covariance matrix for huge data just because the API
accepts one. Dense covariance is exact but expensive. Large time series,
detector channels, and spatial grids need structured models if performance and
memory matter.

Do not add a small diagonal "jitter" silently to force a bad covariance matrix
to factor. If regularization is scientifically justified, it should be explicit
and reported.

Do not interpret a good p-value as proof that the covariance model is correct.
A wrong covariance can hide model mismatch or make uncertainties look more
credible than they are.

Next useful pages: [Damped Oscillator](resonance_decay.md),
[Fitting for Practitioners](@ref), and [Statistical Foundations](@ref).
