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
end
