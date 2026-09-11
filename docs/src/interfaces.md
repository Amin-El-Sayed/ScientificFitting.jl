# Packages And Interfaces

Fit distribution objects, attach named parameters, and select a solver.
The [LHCb mass spectrum](gallery/lhcb_mass_spectrum.md) combines these interfaces
on collision data.

!!! note "v0.3 development API"
    The distribution-object, BuildConstructors and NativeMinuit adapters below
    are on the development branch, not in v0.2. Optional packages must be loaded
    to activate their extensions. The [Python API](python.md) retains NumPy
    callbacks; it does not yet wrap these Julia-specific objects.

## Probability Models

| Package | Interface in ScientificFitting |
|:---|:---|
| [Distributions.jl](https://juliastats.org/Distributions.jl/stable/) | Pass a fixed distribution as `error`, or a function `p -> distribution` to `fit_distribution`. Supports continuous, discrete and multivariate observations. |
| [NumericalDistributions.jl](https://github.com/mmikhasenko/NumericalDistributions.jl) | Numerically normalized densities, also inside mixtures and products. Bin probabilities use a CDF or adaptive quadrature when no CDF is available. |
| [DistributionsHEP.jl](https://github.com/JuliaHEP/DistributionsHEP.jl) | Reuse compatible shapes and native `ExtendedMixtureModel` objects. Extended fits retain component yields for both event samples and histograms. |

Supply a function that constructs the distribution from parameter vector `p`.
Here five illustrative observations determine a Gaussian center and scale:

```@example interfaces
using ScientificFitting, Distributions, OptimizationOptimJL

observations = [4.8, 5.1, 5.6, 5.4, 5.0]
normal_model(p) = Normal(p[1], exp(p[2]))  # p[2] = log(scale), so scale > 0
result = fit_distribution(normal_model, observations; p0=[5.0, log(0.4)])
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

The model-construction function (sometimes called a *factory*) runs once per
objective evaluation, not once per observation. Preserve numeric types for
automatic differentiation; avoid converting parameters to `Float64`. Select
`derivatives=:finite` for models that cannot do that. A PDF that already
underflows to zero cannot be repaired by taking its logarithm; prefer a stable
upstream `logpdf`. See [Distribution Objects](api_fitting.md#Distribution-Objects)
for support, truncation and bin-boundary conventions.

## Named Model Construction

[BuildConstructors.jl](https://github.com/RUB-EP1/BuildConstructors.jl) attaches
names, starts, bounds and fixed values to model parameters. Reusing a name
shares that parameter between components.

**Data:** an illustrative count array, in 20 equal bins on ``[0,10]``.
**Model:** a Gaussian peak with fixed width ``0.4``, plus a uniform background.
Fit the center and expected signal/background counts ``N_s,N_b``.
Both densities are normalized within the observation window:

```math
\mu_i=N_s\int_{e_i}^{e_{i+1}}f_s(x)\,dx
     +N_b\int_{e_i}^{e_{i+1}}f_b(x)\,dx.
```

The uniform component contributes ``N_b/20`` to each bin. The fit uses these
integrals as Poisson means, not densities evaluated at bin centers.

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

optim = fit_distribution(constructor, edges, counts;
    solver=OptimizationSolver(LBFGS()), tol=1e-7)
println(report_text(optim))
@assert optim.converged # hide
```

`::P` marks a parameter descriptor; conflicting metadata for shared names is
rejected. Fitting leaves the constructor unchanged. `fixed=true` treats the
width as exactly known. For an uncertain calibration use an explicit
[parameter constraint](api_fitting.md#Constraints-And-Uncertainty-Objects);
BuildConstructors' `uncertainty` metadata does not add a prior or measurement term.

Component yields set the normalization: omit `total_count` and any additional
Poisson count term. A non-extended distribution requires an explicit `total_count`
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
import NativeMinuit  # load the extension without importing NativeMinuit.profile

minuit = fit_distribution(constructor, edges, counts;
    solver=NativeMinuitSolver(), tol=1e-3)

for (label, fit) in (("Optim", optim), ("MIGRAD", minuit))
    println(label, ": center = ", round(parameter_values(fit).center; digits=5),
        "; -2 log L = ", round(fit.stats.cost_min; digits=5))
end
@assert optim.converged && minuit.converged # hide
@assert isapprox(optim.params, minuit.params; rtol=1e-4) # hide
@assert isapprox(optim.param_covariance, minuit.param_covariance; rtol=1e-3) # hide
@assert isapprox(optim.stats.cost_min, minuit.stats.cost_min; atol=1e-6) # hide
```

Julia's [`import`](https://docs.julialang.org/en/v1/manual/modules/#Standalone-using-and-import)
loads the package without bringing its exports into scope. Both packages export
`profile`; `import NativeMinuit` avoids that name conflict.

Both solvers minimize the same ``-2\log L`` with identical bounds and fixed values.
Optim and MIGRAD use different stopping criteria; equal tolerances do not mean
equal accuracy. Omit `tol` for solver defaults: `1e-10` for Optim with automatic
derivatives, `1e-6` with finite differences, and MIGRAD's native `0.1`.
This example tightens MIGRAD's tolerance and checks costs within `1e-6`, plus
parameters and covariance, during the documentation build.

The adapter sets `errordef=1`. Reported covariance comes from ScientificFitting's
local curvature calculation, not native MINOS intervals. Native failure details
remain available in `result.solver_result.raw`. See
[Solver Adapters](api_fitting.md#Solver-Adapters) for the extension contract.

## Reuse Results And Profile Scans

```@example interfaces
model = fitted_model(minuit)  # returns DistributionsHEP.ExtendedMixtureModel
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
