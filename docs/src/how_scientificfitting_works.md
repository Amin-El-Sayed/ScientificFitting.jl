# How ScientificFitting Works

ScientificFitting is built around one rule: statistical assumptions become explicit
problem objects before a solver is selected. Reports, diagnostics, profiles,
and plots then read the fitted result instead of reconstructing the analysis.

```@raw html
<div class="scientificfitting-fit-flow" aria-label="ScientificFitting fit pipeline" data-flow-direction="top-to-bottom">
  <div class="scientificfitting-fit-track">
    <section class="scientificfitting-fit-stage api">
      <div class="scientificfitting-fit-stage-head">
        <span class="scientificfitting-fit-step">01</span>
        <div><div class="scientificfitting-fit-stage-title">Define the scientific problem</div><span>Public API</span></div>
      </div>
      <div class="scientificfitting-fit-stage-body">
        <div class="scientificfitting-fit-branch four">
          <div class="scientificfitting-fit-node">observations <span>x/y values, counts, bins, or samples</span></div>
          <div class="scientificfitting-fit-node">model <span>Julia function and starting parameters</span></div>
          <div class="scientificfitting-fit-node">uncertainty or sampling law <span>standard deviations, covariance, whitening, or likelihood</span></div>
          <div class="scientificfitting-fit-node">parameter control <span>fixed values, bounds, priors, and constraints</span></div>
        </div>
        <div class="scientificfitting-fit-merge">normalized problem <span><code>FitProblem</code> or <code>LikelihoodFitProblem</code></span></div>
      </div>
    </section>
    <div class="scientificfitting-fit-arrow" aria-hidden="true">↓</div>
    <section class="scientificfitting-fit-stage check">
      <div class="scientificfitting-fit-stage-head">
        <span class="scientificfitting-fit-step">02</span>
        <div><div class="scientificfitting-fit-stage-title">Validate before optimization</div><span>Scientific input checks</span></div>
      </div>
      <div class="scientificfitting-fit-stage-body">
        <div class="scientificfitting-fit-branch">
          <div class="scientificfitting-fit-node">dimensions and mappings <span>matching observations, model output, and parameter indices</span></div>
          <div class="scientificfitting-fit-node">finite physical inputs <span>data, starts, bounds, errors, and model predictions</span></div>
          <div class="scientificfitting-fit-node">valid uncertainty structures <span>positive σ and factorizable covariance or operator output</span></div>
        </div>
        <div class="scientificfitting-fit-stop">Invalid scientific input stops here with an actionable error.</div>
      </div>
    </section>
    <div class="scientificfitting-fit-arrow" aria-hidden="true">↓</div>
    <section class="scientificfitting-fit-stage stats">
      <div class="scientificfitting-fit-stage-head">
        <span class="scientificfitting-fit-step">03</span>
        <div><div class="scientificfitting-fit-stage-title">Construct one objective</div><span>Statistical model</span></div>
      </div>
      <div class="scientificfitting-fit-stage-body">
        <div class="scientificfitting-fit-branch">
          <div class="scientificfitting-fit-node">Gaussian residual cost <span>diagonal weights, covariance factorization, or whitening operator</span></div>
          <div class="scientificfitting-fit-node">likelihood cost <span>Poisson, histogram, unbinned, extended, indexed, or custom</span></div>
          <div class="scientificfitting-fit-node">additive parameter information <span>Gaussian priors and correlated parameter constraints</span></div>
        </div>
        <div class="scientificfitting-fit-merge">objective <span><code>C(p)</code>; bounds and fixed parameters restrict the parameter space rather than adding a penalty</span></div>
      </div>
    </section>
    <div class="scientificfitting-fit-arrow" aria-hidden="true">↓</div>
    <section class="scientificfitting-fit-stage solver">
      <div class="scientificfitting-fit-stage-head">
        <span class="scientificfitting-fit-step">04</span>
        <div><div class="scientificfitting-fit-stage-title">Dispatch a compatible solver</div><span>Numerical backend</span></div>
      </div>
      <div class="scientificfitting-fit-stage-body">
        <div class="scientificfitting-fit-branch two">
          <div class="scientificfitting-fit-node"><code>LsqFit</code> fast path <span>static, unconstrained Gaussian least squares</span></div>
          <div class="scientificfitting-fit-node"><code>Optimization.jl</code> path <span>bounds, parameter-dependent uncertainty, priors, constraints, and likelihoods</span></div>
        </div>
        <div class="scientificfitting-fit-note">Explicitly incompatible backend requests fail instead of silently dropping statistical terms.</div>
      </div>
    </section>
    <div class="scientificfitting-fit-arrow" aria-hidden="true">↓</div>
    <section class="scientificfitting-fit-stage result">
      <div class="scientificfitting-fit-stage-head">
        <span class="scientificfitting-fit-step">05</span>
        <div><div class="scientificfitting-fit-stage-title">Build the fitted result</div><span>Single source of truth</span></div>
      </div>
      <div class="scientificfitting-fit-stage-body">
        <div class="scientificfitting-fit-branch">
          <div class="scientificfitting-fit-node">minimum and parameters <span>solver status and best-fit values</span></div>
          <div class="scientificfitting-fit-node">local uncertainty <span>Jacobian/Hessian covariance and correlations</span></div>
          <div class="scientificfitting-fit-node">fit statistics <span>residuals, cost, ndf, p-value, AIC/BIC where meaningful</span></div>
        </div>
        <div class="scientificfitting-fit-merge"><code>FitResult</code> or <code>LikelihoodFitResult</code></div>
      </div>
    </section>
    <div class="scientificfitting-fit-arrow" aria-hidden="true">↓</div>
    <section class="scientificfitting-fit-stage output">
      <div class="scientificfitting-fit-stage-head">
        <span class="scientificfitting-fit-step">06</span>
        <div><div class="scientificfitting-fit-stage-title">Inspect and communicate</div><span>Post-fit tools</span></div>
      </div>
      <div class="scientificfitting-fit-stage-body">
        <div class="scientificfitting-fit-branch four">
          <div class="scientificfitting-fit-node"><code>report_text</code> <span>reproducible numerical summary</span></div>
          <div class="scientificfitting-fit-node"><code>diagnose</code> <span>structured findings and next actions</span></div>
          <div class="scientificfitting-fit-node"><code>profile</code> / <code>contour</code> <span>controlled refits away from the minimum</span></div>
          <div class="scientificfitting-fit-node optional"><code>plot_fit</code> and Makie tools <span>optional CairoMakie extension</span></div>
        </div>
      </div>
    </section>
  </div>
</div>
```

## The Four Inputs

Every ordinary fit starts with four concepts:

- **Data:** measured `x` and `y`, counts, histogram bins, or indexed
  observations.
- **Model:** a Julia function that maps data coordinates and parameters to
  predictions.
- **Uncertainty model:** `sigma_y`, `sigma_x`, dense/sparse covariance,
  matrix-free static whitening, named error components, Poisson counts,
  histogram likelihoods, or custom objectives.
- **Parameter control:** starting values, bounds, fixed parameters, priors, and
  Gaussian parameter constraints.

Gaussian x-y workflows become a `FitProblem`, which retains x-y observations,
model predictions, residuals, and a Gaussian uncertainty model. Count,
histogram, sample, indexed, and multi-dataset likelihood workflows become a
`LikelihoodFitProblem`, which retains a scalar ``-2\log L`` objective, an
optional goodness-of-fit statistic, and the number of observations; a generic
likelihood need not have a y residual or a natural fit curve.

With the measured inputs `model`, `x`, `y`, and `sigma_y` from the
[Quickstart](@ref), the explicit core path is:

```julia
problem = FitProblem(model, x, y; p0=[1.0, 0.0], sigma_y=sigma_y)
result = fit(problem)
summary = report_text(result)
```

`fit_model(...)` performs the first two lines; after `using CairoMakie`,
`fitplot(...)` adds the plotting workflow.

## What Happens Internally

For Gaussian fits, the uncertainty model — pointwise ``\sigma_i``, dense or
sparse covariance, or a matrix-free `WhiteningOperator` — becomes one whitened
residual cost; the derivation and worked examples are in
[Correlated Measurements And Whitening](statistics.md#Correlated-Measurements-And-Whitening).

Likelihood fits minimize the appropriate ``-2\log L`` objective or deviance on
the scale fixed by [The Cost Convention](statistics.md#The-Cost-Convention);
Poisson and histogram workflows do not invent Gaussian error bars for low
counts.

Backend dispatch (stage 04) follows from the problem; an explicitly requested
incompatible backend fails before optimization.

## What A `FitResult` Contains

`FitResult` and `LikelihoodFitResult` store the numerical minimum and its
statistical interpretation:

- best-fit parameters and local covariance,
- fitted model values and residuals,
- chi-square, the ``-2\log L`` minimum, p-value, AIC, BIC, and degrees of
  freedom where meaningful,
- optimizer status and diagnostics,
- enough problem metadata for plots, reports, profiles, contours, and
  downstream analysis.

## Output Is Switchable

`fitplot` keeps output surfaces independent:

- `print_report=true` prints a text report.
- `show_panel=false` removes the statistics panel from a plot.
- `show_legend=false` removes the legend.
- `stats_position=:right` keeps results outside the data axis.
- `stats_position=:inside` uses a compact in-axis box when space is limited.

`diagnostic_dashboard(result)` and `diagnostic_dashboard_text(result)` are
separate and suggest what to inspect next.

## Plots Stay Extensible

After `using CairoMakie`, `plot_fit(result)` returns a Makie `Figure`; the
fragment below adds experiment-specific content without refitting and assumes
the named values and reference function exist:

```julia
fig = plot_fit(result; show_panel=true, show_legend=true)
ax = fit_axis(fig)

add_vline!(ax, threshold; color=:gray40, linestyle=:dash, label="threshold")
add_curve!(ax, reference_model; color=:black, linestyle=:dot)
add_points!(ax, x_special, y_special; marker=:star5, color=:gray25)
```

Use this for thresholds, extrapolations, accepted regions, literature values,
or derived-quantity markers.

## Why Profiles and Contours Exist

Local covariance is a parabolic approximation at the minimum and can fail for
weak data, bounds, nonlinear parameters, or asymmetric likelihoods;
[Profiles And Contours](statistics.md#Profiles-And-Contours) derives the
nuisance-parameter refit and its thresholds. When a profile is not parabolic,
or a contour does not resemble the local covariance ellipse, report profile
intervals or contour regions instead of symmetric local errors.

Next useful pages: [Quickstart](@ref), the [worked examples](gallery.md), and
the [Statistics Reference](statistics.md).
