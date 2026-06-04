# Damped Oscillator: When A Good-Looking Fit Is Wrong

A dense nonlinear fit can look almost perfect while its residuals reject the
model. This workflow uses a laboratory recording of a freely decaying
mechanical oscillator to answer two questions:

1. What damping rate and oscillation frequency describe the record?
2. Is a constant-frequency damped oscillator an adequate model?

The second question changes the conclusion.

```@raw html
<img class="jufitter-plot jufitter-plot-light" src="../assets/gallery/damped_oscillator_decay_light.png" alt="Damped oscillator model comparison with a prediction band and pull panels">
<img class="jufitter-plot jufitter-plot-dark" src="../assets/gallery/damped_oscillator_decay_dark.png" alt="Damped oscillator model comparison with a prediction band and pull panels in dark mode">
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
| `sigma_phi_rad` | supplied standard uncertainty of the angle | rad |

The example additionally assigns a 0.5 ms standard timestamp uncertainty from
the acquisition timing resolution. Both x and y uncertainty therefore enter the
fit.

This is not a synthetic perfect-data exercise. The record has dense sampling,
periodic parameters, a slowly changing envelope, and residual structure that a
plot of the fitted curve can hide.

## Start With The Physical Baseline

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

## Fit And Diagnose The Baseline

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
sigma_angle = data.sigma_angle
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

The fit converges and gives plausible parameter values, but convergence answers
only whether the optimizer found a minimum. It does not validate the model.

For this fit,

```math
\frac{\chi^2}{\mathrm{ndf}} = 1.564,
\qquad
P(\chi^2) = 1.77\times 10^{-9}.
```

Under the stated independent Gaussian uncertainty model, residuals this
incompatible with the fit would be extraordinarily unlikely. JuFitter therefore
returns diagnostic status `stop`, not `ok`.

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

The new parameter ``\beta`` has units rad s``^{-2}`` and measures the rate of
change of angular frequency.

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

## Read The Comparison, Not Just The Better Curve

The fitted drift and damping parameters are

```math
\beta
= (8.33 \pm 0.42)\times 10^{-5}\ \mathrm{rad\,s^{-2}},
```

```math
\lambda
= (3.484 \pm 0.032)\times 10^{-3}\ \mathrm{s^{-1}},
\qquad
\tau_d
= (287.1 \pm 2.6)\ \mathrm{s}.
```

Across the recorded interval, the fitted angular frequency changes by

```math
\Delta\omega
= \beta(t_\max-t_\min)
= (4.98 \pm 0.25)\times 10^{-3}\ \mathrm{rad\,s^{-1}}.
```

That change is only about 0.15% of ``\omega_\mathrm{ref}``, yet its phase effect
accumulates over many cycles and becomes obvious in the pulls.

The two fits use the same observations and likelihood, so their AIC values may
be compared. The drift model improves AIC by approximately 388 despite adding
only one parameter. The constant-frequency model is therefore inadequate for
this record.

That still does **not** make the drift model publication-ready:

```math
\frac{\chi^2_\mathrm{drift}}{\mathrm{ndf}}=0.248,
\qquad
P(\chi^2_\mathrm{drift})\approx 1.
```

The remaining pulls are substantially narrower than the expected unit scale.
Possible explanations include:

- conservative angle uncertainties,
- correlations between neighboring samples,
- quantization or preprocessing that changes the residual distribution,
- a flexible model absorbing a systematic effect that should instead be
  described by the uncertainty model.

JuFitter reports status `review` for this model. The data support frequency
drift as a useful description, but the local parameter uncertainties should not
be treated as final until the uncertainty model has been audited.

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
julia --project=. examples/gallery/08_damped_oscillator_decay.jl
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
[Constraints And Profiles](constraints_profiles.md) when local covariance
errors may be unreliable.
