# How JuFitter Works

JuFitter is built around one idea: keep the scientific objects explicit, then
let fitting, diagnostics, reports, and plots read from those objects. A fit is
not only an optimizer call. It is a chain from data and uncertainty assumptions
to a cost function, a numerical solver, a `FitResult`, and user-facing output.

```@raw html
<div class="jufitter-fit-flow" aria-label="JuFitter fit pipeline">
  <div class="jufitter-fit-track">
    <section class="jufitter-fit-stage api">
      <div class="jufitter-fit-stage-title">User API</div>
      <div class="jufitter-fit-node public">Julia notebook or script</div>
      <div class="jufitter-fit-node optional">arrays, CSV, DataFrames</div>
      <div class="jufitter-fit-node public">fit_model, fitplot, FitProblem</div>
    </section>
    <div class="jufitter-fit-arrow">→</div>
    <section class="jufitter-fit-stage object">
      <div class="jufitter-fit-stage-title">FitProblem</div>
      <div class="jufitter-fit-node object">data <span>x, y, counts, bins, samples</span></div>
      <div class="jufitter-fit-node object">model <span>Julia function f(x, p)</span></div>
      <div class="jufitter-fit-node object">uncertainties <span>σy, σx, covariance, components</span></div>
      <div class="jufitter-fit-node object">parameter control <span>p0, fixed values, bounds, priors</span></div>
    </section>
    <div class="jufitter-fit-arrow">→</div>
    <section class="jufitter-fit-stage check">
      <div class="jufitter-fit-stage-title">Validation</div>
      <div class="jufitter-fit-node check">shape checks <span>matching lengths and parameter maps</span></div>
      <div class="jufitter-fit-node check">finite values <span>data, errors, model output</span></div>
      <div class="jufitter-fit-node check">uncertainty checks <span>positive σ, positive-definite covariance</span></div>
    </section>
    <div class="jufitter-fit-arrow">→</div>
    <section class="jufitter-fit-stage stats">
      <div class="jufitter-fit-stage-title">Cost Construction</div>
      <div class="jufitter-fit-branch">
        <div class="jufitter-fit-node stats">Gaussian χ² <span>diagonal weights or whitened covariance</span></div>
        <div class="jufitter-fit-node stats">likelihood cost <span>Poisson, histogram, unbinned, custom</span></div>
        <div class="jufitter-fit-node stats">constraint terms <span>bounds, fixed parameters, Gaussian priors</span></div>
      </div>
      <div class="jufitter-fit-merge">merge → C(p)</div>
    </section>
    <div class="jufitter-fit-arrow">→</div>
    <section class="jufitter-fit-stage solver">
      <div class="jufitter-fit-stage-title">Solver Dispatch</div>
      <div class="jufitter-fit-node solver">LsqFit fast path <span>static unconstrained least squares</span></div>
      <div class="jufitter-fit-node solver">Optimization.jl path <span>bounds, constraints, likelihoods, effective variance</span></div>
      <div class="jufitter-fit-node solver">profile / contour scans <span>hold parameters fixed and refit nuisance parameters</span></div>
    </section>
    <div class="jufitter-fit-arrow">→</div>
    <section class="jufitter-fit-stage result">
      <div class="jufitter-fit-stage-title">FitResult</div>
      <div class="jufitter-fit-node result">minimum and parameters</div>
      <div class="jufitter-fit-node result">covariance and residuals</div>
      <div class="jufitter-fit-node result">ndf, p-values, AIC/BIC, diagnostics</div>
    </section>
    <div class="jufitter-fit-arrow">→</div>
    <section class="jufitter-fit-stage output">
      <div class="jufitter-fit-stage-title">Output</div>
      <div class="jufitter-fit-node public">plot_fit / fitplot <span>Makie figure and style contract</span></div>
      <div class="jufitter-fit-node public">report_text <span>console, plot panel, or both</span></div>
      <div class="jufitter-fit-node public">diagnostics <span>what to inspect next</span></div>
      <div class="jufitter-fit-node optional">Makie extensions <span>markers, bands, thresholds</span></div>
    </section>
  </div>
</div>
<div class="jufitter-fit-legend">
  <span><b class="api">blue</b> public API and user-facing output</span>
  <span><b class="object">teal</b> scientific problem objects</span>
  <span><b class="check">gray</b> validation and safety checks</span>
  <span><b class="stats">amber</b> statistical cost construction</span>
  <span><b class="solver">red</b> numerical backend selection</span>
  <span><b class="optional">dashed</b> optional input or extension</span>
</div>
```

## The Four Inputs

Every ordinary fit starts with four concepts:

- **Data:** measured `x` and `y`, counts, histogram bins, or indexed
  observations.
- **Model:** a Julia function that maps data coordinates and parameters to
  predictions.
- **Uncertainty model:** `sigma_y`, `sigma_x`, dense covariance, named error
  components, Poisson counts, histogram likelihoods, or custom objectives.
- **Parameter control:** starting values, bounds, fixed parameters, priors, and
  Gaussian parameter constraints.

These inputs become a `FitProblem`. Convenience functions such as `fitplot` and
`fit_model` build that object for common cases, but the statistical meaning is
the same.

## What Happens Internally

For Gaussian fits, JuFitter builds residuals and turns uncertainty assumptions
into a weighted cost. Independent errors divide residuals by their standard
deviation. Dense covariance is handled by factorization and whitening:

```math
V = L L^\mathsf{T}, \qquad
\chi^2 = \lVert L^{-1}(y-f(x,p)) \rVert^2.
```

For likelihood fits, JuFitter minimizes the appropriate negative
log-likelihood or deviance. Poisson and histogram workflows do not invent
Gaussian error bars for low counts.

The backend is selected by the problem:

- simple static least-squares fits use the fast `LsqFit` path,
- bounds, constraints, effective variance, priors, and likelihoods use the
  `Optimization.jl` path,
- profiles and contours refit while scanning parameters to test whether local
  covariance errors are trustworthy.

## What A `FitResult` Contains

The result stores the numerical minimum and the statistical interpretation:

- best-fit parameters and local covariance,
- fitted model values and residuals,
- chi-square, likelihood, p-value, AIC, BIC, and degrees of freedom where
  meaningful,
- optimizer status and diagnostics,
- enough problem metadata for plots, reports, profiles, contours, and
  downstream analysis.

The important design rule is that output functions consume the result. They do
not re-encode the fit.

## Output Is Switchable

The same fit can produce nothing, terminal output, a plot-side report, or both.
Use these controls deliberately:

- `report=:none` suppresses automatic reporting.
- `report=:console` prints a text report.
- `report=:plot` shows a right-side plot panel.
- `report=:both` does both.
- `show_stats=false` removes the statistics panel from a plot.
- `show_legend=false` removes the legend.
- `stats_position=:right` keeps results outside the data axis.
- `stats_position=:inside` uses a compact in-axis box when space is limited.

Diagnostics are separate. `diagnostic_dashboard(result)` and
`diagnostic_dashboard_text(result)` are for deciding what to inspect next; they
are not mandatory output.

## Plots Stay Extensible

`plot_fit(result)` returns a Makie `Figure`. You can add experiment-specific
content without fitting again:

```julia
fig = plot_fit(result; report=:plot, show_legend=true)
ax = fit_axis(fig)

add_vline!(ax, threshold; color=:gray40, linestyle=:dash, label="threshold")
add_curve!(ax, reference_model; color=:black, linestyle=:dot)
add_points!(ax, x_special, y_special; marker=:star5, color=:gray25)
```

Use this for thresholds, extrapolations, accepted regions, literature values,
or derived-quantity markers. The fit remains a `FitResult`; the extra visual
elements remain Makie objects.

## Why Profiles and Contours Exist

Local covariance assumes the cost is approximately parabolic near the minimum.
That is often good for well-constrained linear problems. It can fail for weak
data, bounds, nonlinear parameters, or asymmetric likelihoods.

Profiles and contours answer a practical question: if one or two parameters are
moved away from the minimum, how much worse can the best refitted model become?
If the profile curve is not parabolic, or a contour does not resemble the local
covariance ellipse, report profile intervals or contour regions instead of
pretending that symmetric local errors are enough.

Next useful pages: [Quickstart](@ref), [Gallery](@ref), and
[Statistical Foundations](@ref).
