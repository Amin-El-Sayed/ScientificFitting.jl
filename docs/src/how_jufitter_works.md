# How JuFitter Works

JuFitter is built around one rule: statistical assumptions become explicit
problem objects before a solver is selected. Reports, diagnostics, profiles,
and plots then read the fitted result instead of reconstructing the analysis.
This makes the one-line interface and the low-level API two views of the same
pipeline.

```@raw html
<div class="jufitter-fit-flow" aria-label="JuFitter fit pipeline" data-flow-direction="top-to-bottom">
  <div class="jufitter-fit-track">
    <section class="jufitter-fit-stage api">
      <div class="jufitter-fit-stage-head">
        <span class="jufitter-fit-step">01</span>
        <div><div class="jufitter-fit-stage-title">Define the scientific problem</div><span>Public API</span></div>
      </div>
      <div class="jufitter-fit-stage-body">
        <div class="jufitter-fit-branch four">
          <div class="jufitter-fit-node">observations <span>x/y values, counts, bins, or samples</span></div>
          <div class="jufitter-fit-node">model <span>Julia function and starting parameters</span></div>
          <div class="jufitter-fit-node">uncertainty or sampling law <span>standard deviations, covariance, whitening, or likelihood</span></div>
          <div class="jufitter-fit-node">parameter control <span>fixed values, bounds, priors, and constraints</span></div>
        </div>
        <div class="jufitter-fit-merge">normalized problem <span><code>FitProblem</code> or <code>LikelihoodFitProblem</code></span></div>
      </div>
    </section>
    <div class="jufitter-fit-arrow" aria-hidden="true">↓</div>
    <section class="jufitter-fit-stage check">
      <div class="jufitter-fit-stage-head">
        <span class="jufitter-fit-step">02</span>
        <div><div class="jufitter-fit-stage-title">Validate before optimization</div><span>Scientific input checks</span></div>
      </div>
      <div class="jufitter-fit-stage-body">
        <div class="jufitter-fit-branch">
          <div class="jufitter-fit-node">dimensions and mappings <span>matching observations, model output, and parameter indices</span></div>
          <div class="jufitter-fit-node">finite physical inputs <span>data, starts, bounds, errors, and model predictions</span></div>
          <div class="jufitter-fit-node">valid uncertainty structures <span>positive σ and factorizable covariance or operator output</span></div>
        </div>
        <div class="jufitter-fit-stop">Invalid scientific input stops here with an actionable error.</div>
      </div>
    </section>
    <div class="jufitter-fit-arrow" aria-hidden="true">↓</div>
    <section class="jufitter-fit-stage stats">
      <div class="jufitter-fit-stage-head">
        <span class="jufitter-fit-step">03</span>
        <div><div class="jufitter-fit-stage-title">Construct one objective</div><span>Statistical model</span></div>
      </div>
      <div class="jufitter-fit-stage-body">
        <div class="jufitter-fit-branch">
          <div class="jufitter-fit-node">Gaussian residual cost <span>diagonal weights, covariance factorization, or whitening operator</span></div>
          <div class="jufitter-fit-node">likelihood cost <span>Poisson, histogram, unbinned, extended, indexed, or custom</span></div>
          <div class="jufitter-fit-node">additive parameter information <span>Gaussian priors and correlated parameter constraints</span></div>
        </div>
        <div class="jufitter-fit-merge">objective <span><code>C(p)</code>; bounds and fixed parameters restrict the parameter space rather than adding a penalty</span></div>
      </div>
    </section>
    <div class="jufitter-fit-arrow" aria-hidden="true">↓</div>
    <section class="jufitter-fit-stage solver">
      <div class="jufitter-fit-stage-head">
        <span class="jufitter-fit-step">04</span>
        <div><div class="jufitter-fit-stage-title">Dispatch a compatible solver</div><span>Numerical backend</span></div>
      </div>
      <div class="jufitter-fit-stage-body">
        <div class="jufitter-fit-branch two">
          <div class="jufitter-fit-node"><code>LsqFit</code> fast path <span>static, unconstrained Gaussian least squares</span></div>
          <div class="jufitter-fit-node"><code>Optimization.jl</code> path <span>bounds, parameter-dependent uncertainty, priors, constraints, and likelihoods</span></div>
        </div>
        <div class="jufitter-fit-note">Explicitly incompatible backend requests fail instead of silently dropping statistical terms.</div>
      </div>
    </section>
    <div class="jufitter-fit-arrow" aria-hidden="true">↓</div>
    <section class="jufitter-fit-stage result">
      <div class="jufitter-fit-stage-head">
        <span class="jufitter-fit-step">05</span>
        <div><div class="jufitter-fit-stage-title">Build the fitted result</div><span>Single source of truth</span></div>
      </div>
      <div class="jufitter-fit-stage-body">
        <div class="jufitter-fit-branch">
          <div class="jufitter-fit-node">minimum and parameters <span>solver status and best-fit values</span></div>
          <div class="jufitter-fit-node">local uncertainty <span>Jacobian/Hessian covariance and correlations</span></div>
          <div class="jufitter-fit-node">fit statistics <span>residuals, cost, ndf, p-value, AIC/BIC where meaningful</span></div>
        </div>
        <div class="jufitter-fit-merge"><code>FitResult</code> or <code>LikelihoodFitResult</code></div>
      </div>
    </section>
    <div class="jufitter-fit-arrow" aria-hidden="true">↓</div>
    <section class="jufitter-fit-stage output">
      <div class="jufitter-fit-stage-head">
        <span class="jufitter-fit-step">06</span>
        <div><div class="jufitter-fit-stage-title">Inspect and communicate</div><span>Post-fit tools</span></div>
      </div>
      <div class="jufitter-fit-stage-body">
        <div class="jufitter-fit-branch four">
          <div class="jufitter-fit-node"><code>report_text</code> <span>reproducible numerical summary</span></div>
          <div class="jufitter-fit-node"><code>diagnose</code> <span>structured findings and next actions</span></div>
          <div class="jufitter-fit-node"><code>profile</code> / <code>contour</code> <span>controlled refits away from the minimum</span></div>
          <div class="jufitter-fit-node optional"><code>plot_fit</code> and Makie tools <span>optional CairoMakie extension</span></div>
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

Gaussian x-y workflows become a `FitProblem`; count, histogram, sample, indexed,
and multi-dataset likelihood workflows become a `LikelihoodFitProblem`.
Convenience functions build these objects for common cases, but the normalized
problem is what the solver receives.

The distinction is the information available to the core. A `FitProblem`
retains x-y observations, model predictions, residuals, and a Gaussian
uncertainty model. A `LikelihoodFitProblem` retains an already defined scalar
``-2\log L`` objective, an optional goodness-of-fit statistic, and the number of
observations; a generic likelihood need not have a y residual or even a natural
fit curve. Both support the same parameter controls, diagnostics, local
covariance, profiles, and contours.

The explicit core path is short. Assuming `model`, `x`, `y`, and `sigma_y` are
the measured inputs from the [Quickstart](@ref):

```julia
problem = FitProblem(model, x, y; p0=[1.0, 0.0], sigma_y=sigma_y)
result = fit(problem)
summary = report_text(result)
```

`fit_model(...)` performs the first two lines. After `using CairoMakie`,
`fitplot(...)` combines the same fit with the plotting workflow.

## What Happens Internally

For Gaussian fits, JuFitter builds residuals and turns uncertainty assumptions
into a weighted cost. Independent errors divide residuals by their standard
deviation. Dense covariance is handled by factorization and whitening:

```math
V = L L^\mathsf{T}, \qquad
\chi^2 = \lVert L^{-1}(y-f(x,p)) \rVert^2.
```

Here ``V=LL^T`` is a Cholesky factorization and
``z=L^{-1}(y-f(x,p))`` is the whitened residual vector. Multiplication by
``L^{-1}`` removes the scale and correlation encoded by ``V``; under a correct
Gaussian model, ``\operatorname{Cov}(z)=I`` and ``\chi^2=z^Tz``. Whitening is
therefore a coordinate transformation of the residuals, not filtering or
smoothing of the data.

For scale, a residual of ``0.2`` with an independent standard deviation of
``0.1`` contributes ``(0.2/0.1)^2=4`` to ``\chi^2``. If several measurements
share the same offset or gain fluctuation, treating those errors independently
would count the same information repeatedly; the covariance or whitening
operator represents that shared motion.

For large structured covariance, a `WhiteningOperator` supplies the equivalent
``L^{-1}`` action without storing the dense matrix. The statistical cost is the
same; only the representation and asymptotic scaling change.

For likelihood fits, JuFitter minimizes the appropriate ``-2\log L`` objective
or deviance. This common scale makes likelihood-ratio thresholds, Hessian
covariance, and information criteria use one convention throughout the
package. Poisson and histogram workflows do not invent Gaussian error bars for
low counts.

After validation and objective construction, the backend is selected by the
problem:

- simple static least-squares fits use the fast `LsqFit` path,
- bounds, constraints, effective variance, priors, and likelihoods use the
  `Optimization.jl` path,
- incompatible explicit backend choices fail before optimization rather than
  silently ignoring bounds, priors, constraints, or uncertainty terms.

Profiles and contours are post-fit analyses, not a third primary solver. They
reuse the normalized problem, hold one or two parameters at requested values,
and repeatedly refit the remaining nuisance parameters.

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

The important design rule is that output functions consume this result. They do
not carry a second copy of the model, parameters, or fit statistics.

## Output Is Switchable

`fitplot` keeps output surfaces independent. Use these controls deliberately:

- `print_report=true` prints a text report.
- `show_panel=false` removes the statistics panel from a plot.
- `show_legend=false` removes the legend.
- `stats_position=:right` keeps results outside the data axis.
- `stats_position=:inside` uses a compact in-axis box when space is limited.

Diagnostics are separate. `diagnostic_dashboard(result)` and
`diagnostic_dashboard_text(result)` are for deciding what to inspect next; they
are not mandatory output.

## Plots Stay Extensible

After `using CairoMakie`, `plot_fit(result)` returns a Makie `Figure`. You can
add experiment-specific content without fitting again. The following is a
fragment that assumes the named physical values and reference function already
exist:

```julia
fig = plot_fit(result; show_panel=true, show_legend=true)
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
moved away from the minimum, how much worse can the best refitted model become
after all nuisance parameters are optimized again? If the profile is not
parabolic, or a contour does not resemble the local covariance ellipse, report
profile intervals or contour regions instead of treating symmetric local errors
as the final uncertainty statement.

Next useful pages: [Quickstart](@ref), [Gallery](@ref), and
[Statistical Foundations](@ref).
