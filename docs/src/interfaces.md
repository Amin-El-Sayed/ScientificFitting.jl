# Packages And Interfaces

Keep your probability model, choose a minimizer, and reuse the fitted object.
These are independent interfaces, not separate workflows for different sciences.

!!! note "v0.3 development API"
    The distribution-object, BuildConstructors and NativeMinuit adapters below
    are on the development branch, not in v0.2. Optional packages must be loaded
    to activate their extensions. The [Python API](python.md) retains NumPy
    callbacks; it does not yet wrap these Julia-specific objects.

## Probability Models

| Package | Interface in ScientificFitting |
|:---|:---|
| [Distributions.jl](https://juliastats.org/Distributions.jl/stable/) | Pass a fixed distribution as `error`, or a parameterized factory to `fit_distribution`. Continuous, discrete and multivariate observations are supported. |
| [NumericalDistributions.jl](https://github.com/mmikhasenko/NumericalDistributions.jl) | Reuse numerical normalization through the same factory interface, including PDF-only components in mixtures and products. Bin probabilities use a CDF or adaptive quadrature. |
| [DistributionsHEP.jl](https://github.com/JuliaHEP/DistributionsHEP.jl) | Reuse compatible shapes and native `ExtendedMixtureModel` objects. Extended fits retain component yields for both event samples and histograms. |

An ordinary factory is sufficient when parameters fit naturally in a vector:

```@example interfaces
using ScientificFitting, Distributions, OptimizationOptimJL

observations = [4.8, 5.1, 5.6, 5.4, 5.0]
factory = p -> Normal(p[1], exp(p[2]))  # positive scale via log parameter
result = fit_distribution(factory, observations; p0=[5.0, log(0.4)])
println(report_text(result))
@assert result.converged # hide
@assert isapprox(result.params[1], sum(observations)/length(observations); atol=1e-5) # hide
```

`fit_distribution` models the observations themselves. A fixed `error` distribution
instead describes additive noise around predictions:

```@example interfaces
x, y = [0., 1., 2., 3.], [0.1, 1.2, 1.8, 3.4]
line(x, p) = @. p[1]*x + p[2]
regression = fit_likelihood_model(line, x, y;
    error=Normal(0, 0.2), p0=[1., 0.])
println("Slope, intercept: ", round.(regression.params; digits=4))
@assert regression.converged # hide
@assert isapprox(regression.params, [1.05, 0.05]; atol=1e-5) # hide
```

For multivariate samples, set `obsdim` explicitly; each vector
observation uses its joint density. Between-observation dependence requires a
joint likelihood, not a product of marginal densities.

The factory is evaluated once per objective evaluation, not once per observation.
It must preserve dual-number types for automatic differentiation. Select
`derivatives=:finite` for models that cannot do that. A PDF that already
underflows to zero cannot be repaired by taking its logarithm; prefer a stable
upstream `logpdf`. See [Distribution Objects](api_fitting.md#Distribution-Objects)
for support, truncation and bin-boundary conventions.

## Named Model Construction

[BuildConstructors.jl](https://github.com/RUB-EP1/BuildConstructors.jl) separates
model construction from parameter metadata. Its optional adapter is independent
of NativeMinuit: names, starts, bounds, fixed values and shared parameters enter
the same fitting API with either solver.

The following small, controlled histogram illustrates composition, not a
detector measurement. A normalized Gaussian component and uniform background
share the observation window; their coefficients are expected counts.

```@example interfaces
using BuildConstructors, DistributionsHEP

@with_parameters(Peak; center::P, width::P, window::Tuple{Float64,Float64}, begin
    truncated(Normal(center, width), window...)
end)

@with_parameters(Spectrum; peak, background, signal_yield::P, background_yield::P, begin
    ExtendedMixtureModel(
        [build_model(peak, pars), background],  # build the nested component
        [signal_yield, background_yield],
    )
end)

constructor = ConstructorOfSpectrum(
    ConstructorOfPeak(
        AdvancedParameter("center", 5.0; boundaries=(4.0, 6.0)),
        AdvancedParameter("width", 0.4; fixed=true), (0.0, 10.0)),
    Uniform(0.0, 10.0),
    AdvancedParameter("signal_yield", 70.0; boundaries=(0.0, 200.0)),
    AdvancedParameter("background_yield", 40.0; boundaries=(0.0, 200.0)),
)
edges = collect(0.0:0.5:10.0)
counts = [2, 1, 3, 1, 2, 2, 1, 1, 3, 14, 33, 22, 7, 2, 1, 2, 3, 1, 1, 2]
nothing # hide
```

`::P` marks a parameter descriptor. Repeated names share one fitted parameter;
conflicting metadata is rejected. The adapter snapshots the constructor rather
than modifying it during fitting. Here `fixed=true` treats the width as exactly
known. For an uncertain calibration use an explicit
[parameter constraint](api_fitting.md#Constraints-And-Uncertainty-Objects);
BuildConstructors' `uncertainty` metadata does not add a prior or measurement term.

For bin ``i``, the expected count is an integral, not the density at its center:

```math
\mu_i=N_s\int_{e_i}^{e_{i+1}}f_s(x)\,dx
     +N_b\int_{e_i}^{e_{i+1}}f_b(x)\,dx.
```

The uniform component contributes ``N_b/20`` to each bin. Component yields
already set the normalization: do not add `total_count` or a second Poisson
term. A non-extended distribution instead requires an explicit `total_count`
for a histogram fit. For individual events use
`fit_distribution(constructor, events; solver=...)`.

## Minimizers

| Package | Role and current integration |
|:---|:---|
| [LsqFit.jl](https://julianlsolvers.github.io/LsqFit.jl/latest/) | Levenberg-Marquardt for compatible Gaussian residual problems; used by the least-squares path. |
| [Optimization.jl](https://docs.sciml.ai/Optimization/stable/) / [Optim.jl](https://docs.sciml.ai/Optimization/stable/optimization_packages/optim/) | Solver interface / Julia algorithms for scalar objectives. Select `OptimizationSolver(algorithm)`; bounds and nonlinear constraints require a compatible algorithm. |
| [NativeMinuit.jl](https://github.com/fkguo/NativeMinuit.jl) | Optional Julia-native MIGRAD adapter: `NativeMinuitSolver()`. Requires Julia 1.11+ for NativeMinuit 0.7; the core still supports Julia 1.10. |
| [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl) | Provides the bounded Nelder-Mead path, `optimizer=:nelder_mead`. This derivative-free choice does not imply the objective has a meaningful Hessian. |
| [NonlinearSolve.jl](https://docs.sciml.ai/NonlinearSolve/stable/solvers/nonlinear_least_squares_solvers/) | Julia residual-based solvers; not currently integrated. |
| [Minuit2.jl](https://github.com/JuliaHEP/Minuit2.jl) | Julia bindings to C++ Minuit2; distinct from NativeMinuit and not currently integrated. |

Change the solver without rebuilding the statistical model:

```@example interfaces
import NativeMinuit  # activates the optional extension

optim = fit_distribution(constructor, edges, counts;
    solver=OptimizationSolver(LBFGS()), tol=1e-7)
minuit = fit_distribution(constructor, edges, counts;
    solver=NativeMinuitSolver(), tol=1e-3)

println(report_text(minuit))
@assert optim.converged && minuit.converged # hide
@assert isapprox(optim.params, minuit.params; rtol=1e-4) # hide
@assert isapprox(optim.param_covariance, minuit.param_covariance; rtol=1e-3) # hide
@assert isapprox(optim.stats.cost_min, minuit.stats.cost_min; atol=1e-6) # hide
```

Both minimize the same ``-2\log L`` with the same fixed values and bounds.
Their tolerances have different meanings: Optim uses its convergence criteria,
whereas MIGRAD uses an estimated-distance-to-minimum (EDM) criterion. Equal
numbers would not imply equal accuracy; this example checks the fitted values,
cost and covariance against each other during the documentation build.
Omit `tol` for solver-specific defaults: `1e-10` for Optim with automatic
derivatives, `1e-6` with finite differences, and MIGRAD's native `0.1`.
The example tightens MIGRAD's tolerance to compare costs within `1e-6`.

The adapter sets `errordef=1`. Reported covariance comes from ScientificFitting's
local curvature calculation, not native MINOS intervals. Native failure details
remain available in `result.solver_result.raw`. See
[Solver Adapters](api_fitting.md#Solver-Adapters) for the extension contract.

## Reuse Results And Profile Scans

```@example interfaces
model = fitted_model(minuit)  # a native ExtendedMixtureModel, not a replacement type
values = parameter_values(minuit)  # named values, including fixed parameters
println("Center: ", round(values.center; digits=4))
println("Expected count: ", round(DistributionsHEP.total_yield(model); digits=3))

# Re-optimize the yields at each trial center; retain the chosen solver.
scan = profile(minuit, 1; nsigma=2.5, npoints=25, on_failure=:throw)
interval = profile_interval(scan)
println("Profile interval: ", round(interval.lower; digits=3),
    " to ", round(interval.upper; digits=3))
comparison = profile(optim, 1; values=scan.values, on_failure=:throw) # hide
@assert isapprox(scan.delta_cost, comparison.delta_cost; atol=1e-5) # hide
@assert isfinite(interval.lower) && isfinite(interval.upper) # hide
```

`model(x)` evaluates the fitted intensity; `MixtureModel(model)` gives the
normalized distribution. Reconstruction and plotting do not refit.
`update!(constructor, values)` explicitly replaces constructor starts if wanted.

Profile crossings at ``\Delta(-2\log L)=1`` have an approximate one-parameter
68.3% interpretation under regular likelihood assumptions. A boundary or an
unidentified component can invalidate it. Compare the scan with the local
parabola; increase the grid resolution before quoting more digits.

```@example interfaces
using CairoMakie  # only needed for the figure
plot_options = (local_sigma=minuit.param_stderr[1],
    title="Component-location likelihood", xlabel="center",
    ylabel="Delta (-2 log L)")
fig = plot_profile(scan; plot_options...)
nothing # hide
```

```@setup interfaces
for style in (:sans, :tex), appearance in (:light, :dark)
    figure = plot_profile(scan; theme=style, appearance, plot_options...)
    save("interface_profile_$(style)_$(appearance).svg", figure)
end
```

```@raw html
<img class="scientificfitting-plot scientificfitting-plot-light" data-scientificfitting-plot-group="interface-profile" data-scientificfitting-plot-style="sans" src="interface_profile_sans_light.svg" alt="Refitted component-location profile compared with the local covariance parabola, sans style">
<img class="scientificfitting-plot scientificfitting-plot-dark" data-scientificfitting-plot-group="interface-profile" data-scientificfitting-plot-style="sans" src="interface_profile_sans_dark.svg" alt="Refitted component-location profile compared with the local covariance parabola, dark sans style">
<img class="scientificfitting-plot scientificfitting-plot-light" data-scientificfitting-plot-group="interface-profile" data-scientificfitting-plot-style="tex" src="interface_profile_tex_light.svg" alt="Refitted component-location profile compared with the local covariance parabola, tex style">
<img class="scientificfitting-plot scientificfitting-plot-dark" data-scientificfitting-plot-group="interface-profile" data-scientificfitting-plot-style="tex" src="interface_profile_tex_dark.svg" alt="Refitted component-location profile compared with the local covariance parabola, dark tex style">
```

## Related Workflows

These packages offer different modeling workflows, not interchangeable minimizers:

| Package | When to consider it |
|:---|:---|
| [RooFit](https://root.cern.ch/manual/roofit/) / [RooFitLite.jl](https://github.com/JuliaHEP/RooFitLite.jl) | Compositional probability modeling in ROOT / a RooFit-style Julia interface. ScientificFitting does not convert their model graphs. A custom likelihood callback is possible if you supply a compatible scalar objective; that is not a native adapter. |
| [GLM.jl](https://juliastats.org/GLM.jl/stable/) | Linear/generalized linear models with formulas, tables and link functions. |
| [Turing.jl](https://turinglang.org/docs/) | Probabilistic models and posterior inference. Posterior credible intervals and profile confidence intervals answer different questions; a posterior sampler is not a drop-in fitting backend. |

Custom density and objective callbacks remain available alongside the adapters.
When transferring a likelihood to a probabilistic model, do not include the same
prior twice. For kafe2's influence and software attribution, see
[Citation and License](citation.md).
