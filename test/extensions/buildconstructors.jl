using ScientificFitting
using Distributions
using BuildConstructors
using BuildConstructors: AbstractConstructor
using Statistics
using Test

struct TestNormalConstructor{M, S} <: AbstractConstructor
    mean::M
    scale::S
end
BuildConstructors.build_model(c::TestNormalConstructor, pars) =
    Normal(BuildConstructors.value(c.mean; pars), BuildConstructors.value(c.scale; pars))

@testset "BuildConstructors works without Minuit" begin
    @test Base.get_extension(ScientificFitting, :ScientificFittingNativeMinuitExt) === nothing
    mu = AdvancedParameter("mu", 0.; boundaries=(-2., 2.), uncertainty=100.)
    sigma = AdvancedParameter("sigma", 0.5; boundaries=(0.1, 2.), uncertainty=100., fixed=true)
    constructor = TestNormalConstructor(mu, sigma)
    original = parameter_metadata(constructor)
    data = [-0.4, -0.1, 0.3, 0.5, 0.8]
    result = fit_distribution(constructor, data)
    @test result.converged
    @test result.params ≈ [mean(data), 0.5] atol=1e-6
    @test result.param_stderr ≈ [0.5/sqrt(5), 0.] atol=1e-6
    @test isempty(result.problem.parameter_priors)
    @test result.problem.parameter_names == ["mu", "sigma"]
    @test result.solver_result.parameter_indices == [1]
    @test parameter_metadata(constructor) == original
    rebuilt = fitted_model(result)
    @test params(rebuilt) == Tuple(result.params)
    @test parameter_values(constructor) == (mu=0., sigma=0.5)
    refit = ScientificFitting._refit_with_fixed(result, [FixedParameter(1, 0.2)])
    @test refit.params == [0.2, 0.5]
    @test params(fitted_model(refit)) == (0.2, 0.5)

    # Updating the user's constructor must not change an existing fit's objective.
    cost = result.problem.objective(result.params)
    BuildConstructors.update!(constructor, (mu=1., sigma=1.2))
    @test result.problem.objective(result.params) == cost
    @test params(fitted_model(result)) == Tuple(result.params)
    @test parameter_values(constructor) == (mu=1., sigma=1.2)

    running = TestNormalConstructor(Running("mu"), BuildConstructors.Fixed(0.5))
    @test_throws ArgumentError fit_distribution(running, data)
    explicit = fit_distribution(running, data; p0=(mu=0.1,))
    @test explicit.params[1] ≈ mean(data) atol=1e-6
    @test_throws ArgumentError fit_distribution(running, data; p0=(unknown=1.,))
    @test_throws ArgumentError fit_distribution(constructor, data; bounds=([-1., 0.1], [1., 2.]))
    @test_throws ArgumentError fit_distribution(constructor, data; p0=(sigma=0.7,))
    @test parameter_values(result) == (mu=result.params[1], sigma=0.5)
    BuildConstructors.update!(constructor, parameter_values(result))
    @test parameter_values(constructor) == parameter_values(result)
    custom = fit_custom(p -> sum(abs2, p); p0=[1.], nobs=3)
    @test_throws ArgumentError parameter_values(custom)
    @test_throws ArgumentError fitted_model(custom)

    edges, counts = [-1., 0., 1., 2.], [5, 10, 3]
    binned = fit_distribution(running, edges, counts; p0=(mu=0.1,), total_count=20.)
    direct = fit_distribution(p -> Normal(p[1], 0.5), edges, counts; p0=[0.1], total_count=20.)
    @test binned.converged && direct.converged
    @test binned.params ≈ direct.params atol=1e-8
    @test binned.problem.parameter_names == ["mu"]
    @test params(fitted_model(binned)) == (binned.params[1], 0.5)
    @test_throws ArgumentError fit_distribution(running, edges, counts;
        p0=(mu=0.1,), total_count=20., parameter_names=["replacement"])
end
