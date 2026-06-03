# Full Covariance

This workflow fits an exponential decay when neighboring measurements share
readout noise. The essential point is that correlated uncertainty changes the
statistical problem. A dense covariance matrix is not a prettier way to draw
error bars; it tells the fit which residual patterns are plausible together.

```@raw html
<img class="jufitter-plot jufitter-plot-light" src="../assets/gallery/full_covariance_decay_light.png" alt="Full covariance exponential fit">
<img class="jufitter-plot jufitter-plot-dark" src="../assets/gallery/full_covariance_decay_dark.png" alt="Full covariance exponential fit in dark mode">
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

## Data And Covariance

The documentation example uses a controlled decay dataset with a known
correlation length. The covariance between two samples ``i`` and ``j`` is

```math
V_{ij} = \sigma^2 \exp\!\left(-\frac{|i-j|}{\ell}\right).
```

This is not meant to be the only valid covariance model. It is a compact way to
represent the common laboratory situation where nearby samples share readout or
baseline noise.

## Model And Cost

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

t = collect(range(0.0, 2.5; length=22))
decay_model(t, p) = @. p[1] * exp(p[2] * t) + p[3]

n = length(t)
base_sigma = 0.055
correlation_length = 2.3

cov_y = [
    base_sigma^2 * exp(-abs(i - j) / correlation_length)
    for i in 1:n, j in 1:n
]

y = decay_model(t, [2.0, -1.12, 0.24]) .+
    0.75 .* base_sigma .* (sin.(1.8 .* t) .+ 0.28 .* cos.(4.1 .* t))

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

Next useful pages: [XY Uncertainties](@ref),
[Fitting for Practitioners](@ref), and [Statistical Foundations](@ref).
