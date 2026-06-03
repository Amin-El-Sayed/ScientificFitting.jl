# XY Uncertainties

This controlled calibration workflow shows what changes when the independent
variable has uncertainty too. Horizontal error bars are not only a plotting
feature: if the model is steep enough, uncertainty in ``x`` contributes to the
statistical cost and to the fitted parameter errors.

```@raw html
<img class="jufitter-plot jufitter-plot-light" src="../assets/gallery/xy_uncertainties_light.png" alt="XY uncertainty fit">
<img class="jufitter-plot jufitter-plot-dark" src="../assets/gallery/xy_uncertainties_dark.png" alt="XY uncertainty fit in dark mode">
```

## Question

A sensor is calibrated by measuring pairs ``(x_\mathrm{meas}, y_\mathrm{meas})``.
Both instruments have finite resolution. The scientific question is:

```math
y = m x + b,
```

with realistic uncertainty on both ``m`` and ``b``. If ``x`` uncertainty is
ignored, the fit treats the measured abscissa as exact and usually reports
parameter errors that are too small.

## Data

The example uses a controlled dataset rather than a perfect synthetic line. The
measured ``x`` values are shifted by a smooth readout error, and ``y`` contains
a different structured measurement error. This keeps the example reproducible
while still producing residuals a scientist should inspect.

The uncertainties are:

- ``\sigma_x = 0.16`` for every measured ``x`` point,
- ``\sigma_y = 0.10`` for every measured ``y`` point.

The visible horizontal and vertical error bars correspond to these standard
uncertainties. The 1-sigma prediction band in the plot includes the fitted
parameter uncertainty and the observation noise used by the plotting routine.

## Model And Cost

For a model ``f(x,p)``, an uncertainty in ``x`` changes the vertical residual
through the local model slope:

```math
f(x + \delta x, p) \approx f(x,p)
+ \frac{\partial f}{\partial x}(x,p)\,\delta x.
```

The effective vertical variance is therefore approximated as

```math
\sigma_{\mathrm{eff},i}^2
= \sigma_{y,i}^2
+ \left(\frac{\partial f}{\partial x}(x_i,p)\sigma_{x,i}\right)^2.
```

For a straight line this becomes

```math
\sigma_{\mathrm{eff},i}^2 = \sigma_y^2 + (m\sigma_x)^2.
```

JuFitter uses this effective variance when `sigma_x` is supplied. This is a
local first-order approximation. It is appropriate for smooth models and
moderate x errors; it is not a full errors-in-variables model.

## Fit

This is the complete code for the documentation example:

```julia
using JuFitter

x_true = collect(range(0.0, 4.0; length=18))
line_model(x, p) = @. p[1] * x + p[2]

sigma_x = fill(0.16, length(x_true))
sigma_y = fill(0.10, length(x_true))

x_measured = x_true .+ sigma_x .* cos.(2.2 .* x_true)
y_measured = line_model(x_measured, [0.9, 1.2]) .+
             sigma_y .* sin.(3.1 .* x_measured)

result = fit_model(
    line_model,
    x_measured,
    y_measured;
    p0=[0.5, 0.5],
    sigma_y=sigma_y,
    sigma_x=sigma_x,
)

slope, intercept = result.params
sigma_slope, sigma_intercept = result.param_stderr

println("m = ", slope, " +/- ", sigma_slope)
println("b = ", intercept, " +/- ", sigma_intercept)
println(diagnostic_dashboard_text(result))
```

The documentation asset generator uses the same data and result, then renders a
light and dark Makie plot with x/y error bars, the fitted line, a 1-sigma
prediction band, legend, and side report panel.

## Diagnostics

This workflow is a good place to check whether a result is numerically
reasonable before moving to more complicated models:

- Compare the result with and without `sigma_x`. The central value should not
  jump wildly for this dataset, but the parameter uncertainty should increase
  when x errors are included.
- Inspect residuals for structure. Smooth residuals mean that the line may be
  an incomplete model, even when the error bars were propagated correctly.
- Check ``\chi^2/\mathrm{ndf}``. Values far above one suggest underestimated
  uncertainties or model mismatch; values far below one suggest overestimated
  uncertainties or correlated residuals.
- For steep nonlinear models, run profiles or a more explicit measurement-error
  model before trusting local symmetric errors.

The diagnostic dashboard is deliberately part of the complete code. For this
controlled example it reports `status = review`, because the residuals are
smooth rather than randomly scattered. That is the correct lesson: x-uncertainty
propagation fixes one part of the uncertainty model, but it does not explain
every structured effect in the data.

## Interpretation

The fitted slope is the calibration sensitivity, and the intercept is the
offset. Because ``\sigma_x`` contributes through the slope, the uncertainty of
``m`` and ``b`` depends on the fitted model itself. This is why x errors cannot
be fixed later by drawing larger horizontal error bars on a finished plot.

For the dataset shown here, the fitted line is approximately

```math
y = (0.854 \pm 0.032)x + (1.290 \pm 0.076).
```

The central values are close to the calibration used to construct the
controlled example, but the dashboard still asks for review because the
residuals contain a smooth component. A serious report would mention that
limitation or improve the measurement model.

The reported covariance is still local. It describes the curvature of the cost
near the best fit after the effective-variance approximation has been applied.
If the approximation is poor, the covariance can be precise but statistically
misleading.

## What Can Go Wrong

Do not use `sigma_x` as a visual-only option. It changes the cost function. If
you only want horizontal error bars in a custom plot, keep that separate from
the fit.

Do not use effective variance blindly for discontinuous or kinked models. A
first-order derivative approximation is not meaningful at a threshold, clipping
point, or sharp regime switch.

Do not assume this solves all x-error problems. Large x errors, latent true
abscissae, and calibration transfer problems may require a full measurement
model with nuisance parameters or a structured covariance description.

Next useful pages: [Photoelectric Work Function](@ref), [Full Covariance](@ref),
and [Statistical Foundations](@ref).
