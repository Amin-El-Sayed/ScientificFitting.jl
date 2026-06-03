# Fitting for Practitioners

This page is the practical entry point before the full statistical reference.
It is written for engineers, physicists, and applied scientists who want to know
which model and uncertainty options to use without starting from likelihood
theory.

## What a Fit Actually Does

A fit chooses parameters ``p`` so that a model ``f(x, p)`` describes measured
data ``(x_i, y_i)`` as well as possible. The important question is not only
where the best curve lies, but also how strongly the data constrain the
parameters.

For independent Gaussian y uncertainties, the standard cost is

```math
\chi^2(p) = \sum_i \left(\frac{y_i - f(x_i, p)}{\sigma_{y,i}}\right)^2.
```

Use this when each point has its own standard uncertainty and the points are not
correlated.

## Which Uncertainty Input Should I Use?

Use `sigma_y` when every point has an independent vertical uncertainty. This is
the right default for most lab measurements, calibrations, and sensor curves.

Use `sigma_x` when the x coordinate is measured with relevant uncertainty. For
smooth models, JuFitter propagates this through the local slope of the model.
This is useful for calibration curves, threshold fits, and steep resonances.

Use `cov_y` when y uncertainties are correlated. This is required when several
points share a calibration constant, a baseline correction, a normalization, or
any other systematic effect.

Use likelihood fits for counts, histograms, and unbinned samples. Least squares
is not the right statistical object for small Poisson counts or distribution
fits.

## What the Plot Band Means

`band=:confidence` shows uncertainty of the fitted mean curve from parameter
covariance. It answers: where could the true model curve plausibly be?

`band=:prediction` combines parameter uncertainty with observation uncertainty.
It answers: where would a new measurement plausibly land, given the fitted model
and the measurement errors?

For dense and precise data, confidence bands can be visually very narrow. That
is not a plotting bug; it means the parameters are strongly constrained. If the
data still scatter outside the band, inspect residuals and goodness-of-fit
statistics instead of increasing the band for visual comfort.

## First Diagnostics

The first number to inspect is not the fitted parameter, but the normalized
residual

```math
r_i=\frac{y_i-f(x_i,p)}{\sigma_i}.
```

For an appropriate model and correctly estimated independent Gaussian
uncertainties, these residuals should look like random draws with width about
one. The chi-square is the sum of their squares:

```math
\chi^2=\sum_i r_i^2.
```

The number of degrees of freedom is approximately

```math
\mathrm{ndf}=N_\mathrm{data}-N_\mathrm{free\ parameters}.
```

If the model and uncertainties are right, ``\chi^2`` should fluctuate around
``\mathrm{ndf}`` with a natural width of about ``\sqrt{2\,\mathrm{ndf}}``.
That is why ``\chi^2/\mathrm{ndf}`` should usually be near one.

Practical rules of thumb:

- ``\chi^2/\mathrm{ndf}\approx 1``: statistically plausible if the residuals
  also look structureless.
- ``0.5 \lesssim \chi^2/\mathrm{ndf} \lesssim 1.5``: often acceptable for
  moderate or large datasets, but still inspect residuals.
- ``\chi^2/\mathrm{ndf} \gg 1``: missing physics, underestimated errors,
  outliers, wrong correlations, or a failed optimizer are likely.
- ``\chi^2/\mathrm{ndf} \ll 1``: uncertainties may be overestimated, points may
  be strongly correlated but treated as independent, or the data may be
  over-smoothed.
- ``P(\chi^2)<0.01``: the result is statistically suspicious under the stated
  assumptions.
- ``P(\chi^2)>0.99``: also suspicious; the data are too good for the assigned
  errors.

The p-value is not the probability that the model is true. It is the
probability, assuming the model and uncertainty model are true, to observe a
chi-square at least as extreme as the measured one. With correct assumptions,
p-values are uniformly distributed across repeated experiments.

Always inspect:

- parameter values and standard errors,
- ``\chi^2/\mathrm{ndf}`` and ``P(\chi^2)`` when the Gaussian assumptions are
  valid,
- residual or pull plots,
- active bounds or failed optimizer convergence,
- whether a covariance matrix is ill-conditioned.

A beautiful curve is not enough. A good fit has interpretable parameters,
credible uncertainties, and residuals that do not show obvious missing physics.

## Quick Diagnosis In The Lab

When a fit looks wrong, start with:

```julia
result = fit_model(model, x, y; p0, sigma_y)
diagnose(result)
```

`diagnose` returns a structured `DiagnosticReport`. In the REPL or a notebook it
prints findings with:

- severity: `INFO`, `WARNING`, or `CRITICAL`,
- evidence: the number that triggered the finding,
- action: the next thing to inspect or change.

The current checks cover:

- optimizer non-convergence,
- non-positive degrees of freedom,
- unavailable goodness-of-fit statistics,
- ill-conditioned covariance or Hessian matrices,
- strong parameter correlations,
- active parameter bounds,
- very large or very small ``\chi^2/\mathrm{ndf}``,
- suspiciously small or large ``P(\chi^2)``,
- large pulls and extreme pulls,
- structured residual signs,
- lag-1 autocorrelation in pulls.

The intent is not to replace scientific judgment. The intent is to make the
first troubleshooting step fast: if the fit is statistically implausible,
overconstrained, bounded, locally degenerate, or visibly structured, JuFitter
should say that directly and suggest what to try next.

## Why Profile And Contour Diagnostics Matter

The parameter covariance matrix is a local approximation. Near the minimum it
assumes the cost function is quadratic:

```math
\Delta \chi^2(\theta_i) \approx
\left(\frac{\theta_i-\hat{\theta}_i}{\sigma_i}\right)^2 .
```

That approximation is convenient because it gives symmetric one-number standard
errors. It can fail when the model is nonlinear, parameters are constrained,
parameters are strongly correlated, the data are weak, or the likelihood is
asymmetric.

A profile scan checks one parameter at a time. JuFitter fixes that parameter,
refits all remaining free parameters, and plots the actual increase in cost. If
the profile follows the local parabola, the local standard error is plausible.
If it is skewed or non-parabolic, use profile intervals instead of symmetric
errors.

```julia
prof = profile(result, 1)
diagnose(prof; local_sigma=result.param_stderr[1])
```

This reports when the scan range is too narrow to bracket the chosen threshold
or when the actual profile disagrees with the local parabolic approximation.
For interval work, the default `profile_interval` path uses adaptive refinement
near threshold crossings. You can request the same behavior explicitly:

```julia
prof = profile(result, 1; adaptive=true, threshold=1.0)
```

Adaptive refinement keeps the broad scan range coarse, then adds points where
the profile crosses the selected ``\Delta C`` threshold. That is usually a
better use of refits than making the entire grid dense.

A contour scan checks pairs of parameters. JuFitter fixes two parameters on a
grid, refits the rest, and plots contours of constant ``\Delta\chi^2``. If the
actual contour follows the local covariance ellipse, the correlation estimate is
locally trustworthy. If the contour bends, clips against a bound, or becomes
banana-shaped, the parameters are not well described by a local Gaussian
approximation.

```julia
cont = contour(result, 1, 2)
diagnose(cont; local_covariance=result.param_covariance,
              local_center=result.params[[1, 2]])
```

This reports when requested contour levels are outside the scan range or when
the profiled contour disagrees with the local covariance ellipse.

For expensive contour scans, use:

```julia
cont = contour(result, 1, 2; adaptive=true, levels=[2.30, 6.18])
```

The adaptive contour pass refines grid cells whose corner values bracket a
requested contour level. It does not make the result magically exact; it makes
the expensive refits concentrate near the contour geometry that will actually
be interpreted.
