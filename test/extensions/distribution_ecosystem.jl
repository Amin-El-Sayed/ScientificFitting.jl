using ScientificFitting
using Distributions
using DistributionsHEP
using NumericalDistributions
using ForwardDiff
using Test

@testset "HEP distributions and numerical normalization" begin
    # The normalization depends on p; its derivatives must not be discarded.
    builds = Ref(0)
    function numerical(p)
        builds[] += 1
        return NumericallyIntegrable(x -> exp(p[1]*x), (0., 1.))
    end
    data = [0.12, 0.37, 0.69, 0.84]
    fitted = fit_distribution(numerical, data; p0=[1.], tol=1e-7)
    before = builds[]
    cost = fitted.problem.objective([1.])
    @test builds[] == before + 1
    reference(p) = -2 * (p[1]*sum(data) - length(data)*log(expm1(p[1])/p[1]))
    @test cost ≈ reference([1.]) rtol=1e-8
    @test ForwardDiff.gradient(fitted.problem.objective, [1.]) ≈ ForwardDiff.gradient(reference, [1.]) rtol=1e-6
    @test ForwardDiff.hessian(fitted.problem.objective, [1.]) ≈ ForwardDiff.hessian(reference, [1.]) rtol=1e-5
    @test fitted.converged

    # Shape and event yield are separate pieces of an extended likelihood.
    parts = [Normal(-1., 0.5), Normal(1., 0.5)]
    events = [-1.4, -0.9, -0.7, 0.8, 1.2]
    factory(p) = ExtendedMixtureModel(parts, exp.(p))
    result = fit_distribution(factory, events; p0=log.([2., 2.]), tol=1e-7)
    @test result.converged
    @test sum(exp.(result.params)) ≈ length(events) atol=1e-5
    @test result.stats.cost_min ≈ 2extended_negative_log_likelihood(factory(result.params), events) rtol=1e-10
    actual(p) = ScientificFitting._distribution_cost(factory(p), events)
    expected(p) = 2extended_negative_log_likelihood(factory(p), events)
    @test ForwardDiff.gradient(actual, [0.1, 0.2]) ≈ ForwardDiff.gradient(expected, [0.1, 0.2])
    @test ForwardDiff.hessian(actual, [0.1, 0.2]) ≈ ForwardDiff.hessian(expected, [0.1, 0.2])

    tail = ExtendedMixtureModel([Normal()], [2.])
    @test isfinite(ScientificFitting._distribution_cost(tail, [40.]))
    zero_component = ExtendedMixtureModel(parts, [0., 2.])
    @test ScientificFitting._distribution_cost(zero_component, events) ≈ 2extended_negative_log_likelihood(zero_component, events)
    @test_throws ArgumentError ScientificFitting._distribution_cost(ExtendedMixtureModel(parts, [-1., 2.]), events)
    @test_throws ArgumentError ScientificFitting._distribution_cost(ExtendedMixtureModel(parts, [0., 0.]), events)
end
