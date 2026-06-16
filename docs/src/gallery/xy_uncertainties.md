# XY Uncertainties

This controlled calibration workflow shows what changes when the independent
variable has uncertainty too. Horizontal error bars are not only a plotting
feature: if the model is steep enough, uncertainty in ``x`` contributes to the
statistical cost and to the fitted parameter errors.

```@raw html
<img class="jufitter-plot jufitter-plot-light" data-jufitter-plot-group="xy-uncertainties" data-jufitter-plot-style="lab" src="../assets/gallery/xy_uncertainties_lab_light.png" alt="XY uncertainty fit in lab style">
<img class="jufitter-plot jufitter-plot-dark" data-jufitter-plot-group="xy-uncertainties" data-jufitter-plot-style="lab" src="../assets/gallery/xy_uncertainties_lab_dark.png" alt="XY uncertainty fit in lab dark style">
<img class="jufitter-plot jufitter-plot-light" data-jufitter-plot-group="xy-uncertainties" data-jufitter-plot-style="modern" src="../assets/gallery/xy_uncertainties_modern_light.png" alt="XY uncertainty fit in modern style">
<img class="jufitter-plot jufitter-plot-dark" data-jufitter-plot-group="xy-uncertainties" data-jufitter-plot-style="modern" src="../assets/gallery/xy_uncertainties_modern_dark.png" alt="XY uncertainty fit in modern dark style">
<img class="jufitter-plot jufitter-plot-light" data-jufitter-plot-group="xy-uncertainties" data-jufitter-plot-style="article" src="../assets/gallery/xy_uncertainties_article_light.png" alt="XY uncertainty fit in article style">
<img class="jufitter-plot jufitter-plot-dark" data-jufitter-plot-group="xy-uncertainties" data-jufitter-plot-style="article" src="../assets/gallery/xy_uncertainties_article_dark.png" alt="XY uncertainty fit in article dark style">
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

## Model and Cost

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

The size of the effect is easy to estimate before fitting. If ``m\approx0.85``
and ``\sigma_x=0.16``, the x-resolution contributes about
``m\sigma_x\approx0.14`` in vertical units. That is already larger than
``\sigma_y=0.10``. In this situation drawing horizontal error bars but fitting
as if x were exact would understate the parameter uncertainty.

JuFitter uses this effective variance when `sigma_x` is supplied. This is a
local first-order approximation. It is appropriate for smooth models and
moderate x errors; it is not a full errors-in-variables model.

## Fit

This is the complete code for the documentation example:

```julia
using JuFitter

# Both coordinates are measured. sigma_x and sigma_y are standard
# uncertainties, not visual-only error-bar lengths.
x_measured = [0.1600, 0.3743, 0.5522, 0.7087, 0.8645, 1.0403, 1.2519,
              1.5053, 1.7958, 2.1091, 2.4246, 2.7213, 2.9831, 3.2032,
              3.3854, 3.5437, 3.6982, 3.8702]
y_measured = [1.3916, 1.6286, 1.7960, 1.9189, 2.0226, 2.1280, 2.2593,
              2.4549, 2.7506, 3.1234, 3.4764, 3.7327, 3.9024, 4.0345,
              4.1591, 4.2893, 4.4392, 4.6294]
sigma_x = fill(0.16, length(x_measured))
sigma_y = fill(0.10, length(x_measured))

line_model(x, p) = @. p[1] * x + p[2]

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

```@raw html
<div class="jufitter-cell-output">
<div class="jufitter-cell-output-label">Real output (abridged)</div>
<pre>m = 0.8542596063464186 +/- 0.03220459924722076
b = 1.289983732168427 +/- 0.07606775503455655

Fit diagnostic dashboard
status = review - inspect diagnostics
critical = 0, warning = 3, info = 0
3 warning(s). Inspect before trusting uncertainties or conclusions.

Next actions:
  1. Use a covariance model, inspect acquisition order/time dependence, or fit a model with the missing systematic component.
  2. The uncertainties may be too large, correlations may be ignored, or the data may not be independent.
  3. Uncertainties may be overestimated, correlations may be ignored, or the data may have been smoothed/averaged.</pre>
</div>
```

The documentation asset generator uses these same arrays and result, then
renders a light/dark Makie plot for each supported documentation style. The
visible band is a 1-sigma prediction band; it is not a profile interval and not
a confidence band for the mean line alone.

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

Next useful pages: [Full Covariance](@ref),
[Damped Oscillator](resonance_decay.md), and
[Statistical Foundations](@ref).
