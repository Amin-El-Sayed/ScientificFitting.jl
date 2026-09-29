using ScientificFitting
using BuildConstructors
using BuildConstructors: AbstractConstructor
using Distributions
using DistributionsHEP
import NativeMinuit
using Test

struct ShiftedPeakConstructor{M, S} <: AbstractConstructor
    mean::M
    scale::S
    offset::Float64
end
BuildConstructors.build_model(c::ShiftedPeakConstructor, pars) =
    Normal(BuildConstructors.value(c.mean; pars) + c.offset, BuildConstructors.value(c.scale; pars))

struct TwoPeakConstructor{L, R, A, B} <: AbstractConstructor
    signal::L
    background::R
    signal_yield::A
    background_yield::B
end
BuildConstructors.build_model(c::TwoPeakConstructor, pars) = ExtendedMixtureModel(
    [build_model(c.signal, pars), build_model(c.background, pars)],
    [BuildConstructors.value(c.signal_yield; pars), BuildConstructors.value(c.background_yield; pars)],
)

@testset "Named shared model, extended likelihood, interchangeable solvers" begin
    mu = AdvancedParameter("mu", 0.; boundaries=(-1., 1.))
    sigma = AdvancedParameter("sigma", 0.6; fixed=true)
    c = TwoPeakConstructor(ShiftedPeakConstructor(mu, sigma, 0.), ShiftedPeakConstructor(mu, sigma, 3.),
        AdvancedParameter("signal_yield", 4.; boundaries=(0.01, 50.)),
        AdvancedParameter("background_yield", 4.; boundaries=(0.01, 50.)))
    data = [-0.4, -0.2, 0., 0.2, 0.5, 2.9, 3.2, 3.4]
    @test_logs (:warn, r"Shared parameters") validate_parameters(c)
    optim = @test_logs (:warn, r"Shared parameters") fit_distribution(c, data;
        solver=OptimizationSolver(ScientificFitting.OptimizationOptimJL.LBFGS()), tol=1e-7)
    minuit = @test_logs (:warn, r"Shared parameters") fit_distribution(c, data;
        solver=NativeMinuitSolver(steps=[0.1, 0.1, 0.5, 0.5]), tol=1e-7)
    @test optim.converged && minuit.converged
    @test optim.params ≈ minuit.params atol=2e-4
    @test optim.param_covariance ≈ minuit.param_covariance rtol=1e-3
    @test optim.stats.cost_min ≈ minuit.stats.cost_min atol=1e-7
    @test sum(minuit.params[3:4]) ≈ length(data) atol=1e-3
    @test minuit.problem.parameter_names == ["mu", "sigma", "signal_yield", "background_yield"]
    @test minuit.solver_result.parameter_indices == [1, 3, 4]
    @test collect(minuit.solver_result.raw.parameters) == ["mu", "signal_yield", "background_yield"]
    @test parameter_values(c) == (mu=0., sigma=0.6, signal_yield=4., background_yield=4.)
    model = fitted_model(minuit)
    @test minuit.stats.cost_min ≈ 2extended_negative_log_likelihood(model, data) rtol=1e-10
    values = [minuit.params[1] - 0.1, minuit.params[1], minuit.params[1] + 0.1]
    a = profile(optim, 1; values, on_failure=:throw)
    b = profile(minuit, 1; values, on_failure=:throw)
    @test a.delta_cost ≈ b.delta_cost atol=1e-5
    refit = ScientificFitting._refit_with_fixed(minuit, [FixedParameter(1, values[1])])
    @test refit.options.solver === minuit.options.solver
    @test collect(refit.solver_result.raw.parameters) == ["signal_yield", "background_yield"]

    bad = TwoPeakConstructor(c.signal,
        ShiftedPeakConstructor(AdvancedParameter("mu", 0.2), sigma, 3.), c.signal_yield, c.background_yield)
    @test_throws ArgumentError fit_distribution(bad, data)
    @test_throws ArgumentError fit_custom(p -> sum(abs2, p); p0=[1., 1.], nobs=4, parameter_names=["mu", "mu"])

    edges, counts = [-2.5, -0.5, 0.5, 2.5, 3.5, 5.5], [2, 14, 3, 9, 2]
    binned_optim = @test_logs (:warn, r"Shared parameters") fit_distribution(c, edges, counts;
        solver=OptimizationSolver(ScientificFitting.OptimizationOptimJL.LBFGS()), tol=1e-7)
    binned_minuit = @test_logs (:warn, r"Shared parameters") fit_distribution(c, edges, counts;
        solver=NativeMinuitSolver(steps=[0.1, 0.1, 0.5, 0.5]), tol=1e-7)
    @test binned_optim.converged && binned_minuit.converged
    @test binned_optim.params ≈ binned_minuit.params atol=2e-4
    @test binned_optim.param_covariance ≈ binned_minuit.param_covariance rtol=1e-3
    @test binned_optim.stats.cost_min ≈ binned_minuit.stats.cost_min atol=1e-7
    values = [binned_optim.params[1] - 0.1, binned_optim.params[1], binned_optim.params[1] + 0.1]
    @test profile(binned_optim, 1; values, on_failure=:throw).delta_cost ≈
        profile(binned_minuit, 1; values, on_failure=:throw).delta_cost atol=1e-5
end

@testset "NativeMinuit defaults support binned nuisance refits" begin
    # At the generic 1e-10 tolerance the peak fit converged, but MIGRAD rejected
    # a nuisance refit on EDM. Use the native criterion, not an implicit retry.
    edges = collect(0.0:0.5:10.0)
    counts = [2, 1, 3, 1, 2, 2, 1, 1, 3, 14, 33, 22, 7, 2, 1, 2, 3, 1, 1, 2]
    model(p) = ExtendedMixtureModel(
        [truncated(Normal(p[1], 0.4), 0.0, 10.0), Uniform(0.0, 10.0)], p[2:3])
    settings = (; p0=[5., 70., 40.], bounds=([4., 0., 0.], [6., 200., 200.]))
    r = fit_distribution(model, edges, counts; settings..., solver=NativeMinuitSolver())
    reference = fit_distribution(model, edges, counts; settings..., tol=1e-7)
    @test r.converged && reference.converged
    @test r.options.tol == NativeMinuit.Minuit(p -> sum(abs2, p), zeros(3)).tol
    @test maximum(abs.((r.params .- reference.params) ./ reference.param_stderr)) < 0.01
    @test r.param_covariance ≈ reference.param_covariance rtol=1e-3
    @test r.stats.cost_min ≈ reference.stats.cost_min atol=2e-4
    scan = profile(r, 1; nsigma=2.5, npoints=25, on_failure=:throw)
    comparison = profile(reference, 1; values=scan.values, on_failure=:throw)
    # Native EDM thresholds are in cost units, not parameter-error units.
    @test scan.delta_cost ≈ comparison.delta_cost atol=2e-3
    interval, ref_interval = profile_interval(scan), profile_interval(comparison)
    @test interval.lower ≈ ref_interval.lower atol=1e-4
    @test interval.upper ≈ ref_interval.upper atol=1e-4
    refit = ScientificFitting._refit_with_fixed(r, [FixedParameter(1, scan.values[1])])
    @test refit.converged
    @test refit.options.tol == refit.solver_result.raw.tol == r.options.tol
end
