# Damped Oscillator: When A Good-Looking Fit Is Wrong

A dense nonlinear fit can look almost perfect while its residuals reject the
model. This workflow uses a laboratory recording of a freely decaying
mechanical oscillator to answer two questions:

1. What damping rate and oscillation frequency describe the record?
2. Is a constant-frequency damped oscillator an adequate model?

The second question changes the conclusion.

## Question

The measurement asks whether a standard constant-frequency damped oscillator is
an adequate description of the recorded motion, or whether the data require a
small additional frequency drift.

```@raw html
<img class="jufitter-plot jufitter-plot-light" data-jufitter-plot-group="damped-oscillator" data-jufitter-plot-style="workbench" src="../assets/gallery/damped_oscillator_decay_workbench_light.png" alt="Damped oscillator model comparison in workbench style">
<img class="jufitter-plot jufitter-plot-dark" data-jufitter-plot-group="damped-oscillator" data-jufitter-plot-style="workbench" src="../assets/gallery/damped_oscillator_decay_workbench_dark.png" alt="Damped oscillator model comparison in workbench dark style">
<img class="jufitter-plot jufitter-plot-light" data-jufitter-plot-group="damped-oscillator" data-jufitter-plot-style="showcase" src="../assets/gallery/damped_oscillator_decay_showcase_light.png" alt="Damped oscillator model comparison in showcase style">
<img class="jufitter-plot jufitter-plot-dark" data-jufitter-plot-group="damped-oscillator" data-jufitter-plot-style="showcase" src="../assets/gallery/damped_oscillator_decay_showcase_dark.png" alt="Damped oscillator model comparison in showcase dark style">
<img class="jufitter-plot jufitter-plot-light" data-jufitter-plot-group="damped-oscillator" data-jufitter-plot-style="publication" src="../assets/gallery/damped_oscillator_decay_publication_light.png" alt="Damped oscillator model comparison in publication style">
<img class="jufitter-plot jufitter-plot-dark" data-jufitter-plot-group="damped-oscillator" data-jufitter-plot-style="publication" src="../assets/gallery/damped_oscillator_decay_publication_dark.png" alt="Damped oscillator model comparison in publication dark style">
```

The main panel alone barely distinguishes the two models. The pull panels do:
the constant-frequency model leaves coherent deviations, while a weak frequency
drift removes most of that structure.

## The Measurement

The distributed CSV contains 300 angle measurements between approximately
20 s and 80 s:

| Column | Meaning | Unit |
| --- | --- | --- |
| `time_s` | acquisition timestamp | s |
| `phi_rad` | measured angular displacement | rad |
| `sigma_phi_rad` | conservative angle repeatability scale | rad |

The analysis below uses half of the conservative angle repeatability scale as
the pointwise statistical uncertainty and assigns a 0.5 ms standard timestamp
uncertainty from the acquisition timing resolution. Both x and y uncertainty
therefore enter the fit. The scale choice is not hidden in the plotting code:
it is part of the statistical model and should be changed if repeated
measurements or instrument specifications justify a different value.

The record has dense sampling, periodic parameters, a slowly changing envelope,
and residual structure that a plot of the fitted curve can hide.

## Model: Start With The Physical Baseline

For an underdamped linear oscillator with constant coefficients,

```math
\phi(t)
= A_\mathrm{ref}
  e^{-\lambda \tau}
  \cos\!\left(\omega_\mathrm{ref}\tau+\phi_\mathrm{ref}\right),
\qquad
\tau=t-t_\mathrm{ref}.
```

The parameters are:

- ``A_\mathrm{ref}``: amplitude at the reference time,
- ``\lambda``: exponential damping rate,
- ``\omega_\mathrm{ref}``: angular frequency,
- ``\phi_\mathrm{ref}``: phase at the reference time.

The time coordinate is centered at the middle of the record. This is not a
cosmetic rewrite. Without centering, phase and frequency must compensate for a
large arbitrary time origin and become more strongly correlated.

The damping time is the derived quantity

```math
\tau_d = \frac{1}{\lambda}.
```

## Diagnostics: Fit and Diagnose the Baseline

The fit uses a Gaussian likelihood with supplied angle uncertainty and
effective-variance propagation of timestamp uncertainty. Multiple initial
guesses are important because phase-periodic models have repeated local minima.

```julia
using JuFitter

function load_damped_oscillator(path)
    rows = readlines(path)[2:end]
    parsed = [parse.(Float64, split(row, ",")) for row in rows if !isempty(strip(row))]
    return (
        time=[row[1] for row in parsed],
        angle=[row[2] for row in parsed],
        sigma_angle=[row[3] for row in parsed],
    )
end

data = load_damped_oscillator(
    "examples/data/damped_oscillator/pohl_wheel_free_decay.csv",
)
time = data.time
angle = data.angle
sigma_angle = 0.5 .* data.sigma_angle
sigma_time = fill(0.0005, length(time))
time_reference = (minimum(time) + maximum(time)) / 2

constant_frequency_model(t, p) = @. p[1] * exp(-p[4] * (t - time_reference)) *
                                     cos(p[2] * (t - time_reference) + p[3])

constant_result = fit_model(
    constant_frequency_model,
    time,
    angle;
    p0=[1.6, 3.26, 0.0, 0.0035],
    sigma_y=sigma_angle,
    sigma_x=sigma_time,
    bounds=([0.0, 2.0, -20.0, 0.0], [5.0, 5.0, 20.0, 0.05]),
    initial_guesses=[
        [1.6, 3.26, 0.0, 0.0035],
        [1.8, 3.20, 2.0, 0.0020],
        [1.5, 3.35, -2.0, 0.0060],
    ],
    maxiters=3000,
)

println(report_text(
    constant_result;
    parameter_names=["A_ref", "omega_ref", "phi_ref", "lambda"],
))
println(diagnostic_dashboard_text(constant_result))
```

```@raw html
<div class="jufitter-cell-output">
<div class="jufitter-cell-output-label">Real output (abridged)</div>
<pre>Fit report
backend = optimization
converged = true
iterations = 111
message = Success

Parameters:
  A_ref = 1.60613 +/- 0.000457675
  omega_ref = 3.26011 +/- 1.7327e-5
  phi_ref = -0.763199 +/- 0.0003003
  lambda = 0.00348387 +/- 1.639e-5

Statistics:
  cost = gaussian_nll
  cost_min = -885.658
  nll_min = -885.658
  chi2 = 1659.13
  ndf = 296
  chi2/ndf = 5.60516
  pvalue = 4.3964999999999994e-188
  AIC = -877.658
  BIC = -862.843

Fit diagnostic dashboard
status = critical - fix before use
critical = 3, warning = 1, info = 0
3 critical issue(s), 1 warning(s). Fix the issue before using this result for conclusions.

Next actions:
  1. Inspect the corresponding data point, uncertainty, units, and possible outlier handling before trusting the fit.
  2. Under the stated assumptions this fit is statistically implausible. Inspect residuals and the uncertainty model.
  3. Look for missing physics, underestimated uncertainties, wrong correlations, outliers, or a failed optimizer.
  4. Use a covariance model, inspect acquisition order/time dependence, or fit a model with the missing systematic component.
</pre>
</div>
```

The fit converges and gives plausible parameter values, but convergence answers
only whether the optimizer found a minimum. It does not validate the model.

For this fit,

```math
\frac{\chi^2}{\mathrm{ndf}} = 5.61,
\qquad
P(\chi^2) = 4.4\times 10^{-188}.
```

Under the stated independent Gaussian uncertainty model, residuals this
incompatible with the fit would be extraordinarily unlikely. JuFitter therefore
returns a critical dashboard status, not `ok`.

The first pull panel explains why. The deviations change coherently over time
instead of scattering without structure around zero. Adding more digits to the
reported damping rate would not repair this.

## Test A Specific Missing Effect

A slowly changing oscillation frequency produces an accumulating phase error.
The smallest useful extension adds a linear frequency drift:

```math
\phi(t)
= A_\mathrm{ref}e^{-\lambda\tau}
  \cos\!\left(
    \omega_\mathrm{ref}\tau
    + \frac{1}{2}\beta\tau^2
    + \phi_\mathrm{ref}
  \right),
```

where

```math
\omega(t)=\frac{\mathrm{d}}{\mathrm{d}t}
\left(
  \omega_\mathrm{ref}\tau+\frac{1}{2}\beta\tau^2+\phi_\mathrm{ref}
\right)
=\omega_\mathrm{ref}+\beta\tau.
```

The new parameter ``\beta`` has units ``\mathrm{rad\,s^{-2}}`` and measures the
rate of change of angular frequency.

```julia
frequency_drift_model(t, p) = @. p[1] * exp(-p[4] * (t - time_reference)) *
                                  cos(
    p[2] * (t - time_reference) +
    0.5 * p[5] * (t - time_reference)^2 +
    p[3],
)

drift_result = fit_model(
    frequency_drift_model,
    time,
    angle;
    p0=[1.6, 3.26, 0.0, 0.0035, 0.0],
    sigma_y=sigma_angle,
    sigma_x=sigma_time,
    bounds=([0.0, 2.0, -20.0, 0.0, -0.01], [5.0, 5.0, 20.0, 0.05, 0.01]),
    initial_guesses=[
        [1.6, 3.26, 0.0, 0.0035, 0.0],
        [1.8, 3.20, 2.0, 0.0020, 0.0001],
        [1.5, 3.35, -2.0, 0.0060, -0.0001],
    ],
    maxiters=4000,
)

println(report_text(
    drift_result;
    parameter_names=["A_ref", "omega_ref", "phi_ref", "lambda", "beta"],
))
println(diagnostic_dashboard_text(drift_result))
```

```@raw html
<div class="jufitter-cell-output">
<div class="jufitter-cell-output-label">Real output (abridged)</div>
<pre>Fit report
backend = optimization
converged = true
iterations = 80
message = Success

Parameters:
  A_ref = 1.60605 +/- 0.000457643
  omega_ref = 3.26016 +/- 1.73882e-5
  phi_ref = -0.77562 +/- 0.000448643
  lambda = 0.00348309 +/- 1.63859e-5
  beta = 8.32827e-5 +/- 2.23515e-6

Statistics:
  cost = gaussian_nll
  cost_min = -2274.67
  nll_min = -2274.67
  chi2 = 270.124
  ndf = 295
  chi2/ndf = 0.915675
  pvalue = 0.847734
  AIC = -2264.67
  BIC = -2246.15

Fit diagnostic dashboard
status = ok - no immediate issue
critical = 0, warning = 0, info = 0
No major diagnostic issues detected by the current checks.
No next action required by the current diagnostic checks.</pre>
</div>
```

## Read The Comparison, Not Just The Better Curve

The fitted drift and damping parameters are

```math
\beta
= (8.33 \pm 0.22)\times 10^{-5}\ \mathrm{rad\,s^{-2}},
```

```math
\lambda
= (3.483 \pm 0.016)\times 10^{-3}\ \mathrm{s^{-1}},
\qquad
\tau_d
= (287.1 \pm 1.4)\ \mathrm{s}.
```

Across the recorded interval, the fitted angular frequency changes by

```math
\Delta\omega
= \beta(t_\max-t_\min)
= (4.98 \pm 0.13)\times 10^{-3}\ \mathrm{rad\,s^{-1}}.
```

That change is only about 0.15% of ``\omega_\mathrm{ref}``, yet its phase effect
accumulates over many cycles and becomes obvious in the pulls.

The two fits use the same observations and likelihood, so their AIC values may
be compared. The drift model improves AIC by approximately 1387 despite adding
only one parameter. The constant-frequency model is therefore inadequate for
this record.

The improved model also passes the automatic first-pass diagnostics:

```math
\frac{\chi^2_\mathrm{drift}}{\mathrm{ndf}}=0.916,
\qquad
P(\chi^2_\mathrm{drift})=0.848.
```

That does not prove the model is true, but it means the first residual and
goodness-of-fit checks no longer reject the stated model and uncertainty scale.

## What The Band Means

The main panel shows the drift model's **local 1σ prediction band**. It combines:

```math
\sigma_\mathrm{pred}^2(t)
= J_p(t)\,\mathrm{Cov}(p)\,J_p^\mathsf{T}(t)
+ \sigma_\phi^2
+ \left(\frac{\partial\phi}{\partial t}\sigma_t\right)^2.
```

The band is narrow compared with the full oscillation amplitude and is
therefore difficult to judge in the main panel. The pull panels display the
same uncertainty scale directly: the darker region is ``\pm1\sigma`` and the
lighter region is ``\pm2\sigma``.

This band uses the local covariance matrix. It is conditional on the fitted
model and does not include uncertainty about whether frequency drift is the
correct physical explanation.

## Complete Reproducible Figure

The tracked script contains the complete fit, diagnostics, derived quantities,
prediction-band propagation, and light/dark Makie figure:

```bash
julia --project=docs examples/gallery/08_damped_oscillator_decay.jl
```

It writes the documentation assets and prints both diagnostic dashboards.

## What To Do Before Reporting A Physical Result

For a critical analysis, the next work is experimental rather than cosmetic:

1. Determine whether neighboring angle samples share acquisition or filtering
   correlations.
2. Validate the supplied angle uncertainty against repeated measurements or
   stationary segments.
3. Test whether the frequency drift repeats in independent decay records.
4. Compare the drift model with a physically motivated nonlinear-oscillator or
   friction model.
5. Refit after defining the correct covariance model, then reassess pulls and
   parameter intervals.

The correct conclusion from this page is not “the more flexible model wins.”
It is: **the baseline model fails; frequency drift explains the dominant
structure; the uncertainty model still requires investigation.**

Next, use [Full Covariance](full_covariance.md) to see how correlated
measurements change fit interpretation, or
[Constraints and Profiles](constraints_profiles.md) when local covariance
errors may be unreliable.
