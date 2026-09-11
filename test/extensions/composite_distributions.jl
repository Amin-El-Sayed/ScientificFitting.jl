using ScientificFitting
using Distributions
using DistributionsHEP
using NumericalDistributions
using ForwardDiff
using Test

@testset "Numerical densities inside native mixtures and products" begin
    data = [0.12, 0.37, 0.69, 0.84]
    builds = Ref(0)
    function numerical(p)
        builds[] += 1
        return NumericallyIntegrable(x -> exp(p[1]*x), (0., 1.))
    end
    mixture(p) = MixtureModel([numerical(p), Uniform()], [p[2], 1-p[2]])
    cost(p) = ScientificFitting._distribution_cost(mixture(p), data)
    # The exponential normalization is parameter-dependent, including inside a mixture.
    reference(p) = -2sum(log(p[2]*exp(p[1]*x)*p[1]/expm1(p[1]) + 1-p[2]) for x in data)
    p = [1.2, 0.35]
    @test cost(p) ≈ reference(p) rtol=1e-8
    before = builds[]
    cost(p)
    @test builds[] == before + 1
    @test ForwardDiff.gradient(cost, p) ≈ ForwardDiff.gradient(reference, p) rtol=1e-6
    @test ForwardDiff.hessian(cost, p) ≈ ForwardDiff.hessian(reference, p) rtol=1e-5

    # Both public product constructors must retain independent factors, not refit them.
    samples = [0.1 -0.2 0.5 0.4; data']
    tuple_product(p) = product_distribution(Normal(p[2], 1.), numerical(p))
    vector_product(p) = product_distribution([Normal(p[2], 1.), numerical(p)])
    expected(p) = -2sum(logpdf(Normal(p[2], 1.), x[1]) +
        p[1]*x[2] - log(expm1(p[1])/p[1]) for x in eachcol(samples))
    for factory in (tuple_product, vector_product)
        actual(p) = ScientificFitting._distribution_cost(factory(p), samples)
        @test actual(p) ≈ expected(p) rtol=1e-8
        @test ForwardDiff.gradient(actual, p) ≈ ForwardDiff.gradient(expected, p) rtol=1e-6
        @test ForwardDiff.hessian(actual, p) ≈ ForwardDiff.hessian(expected, p) rtol=1e-5
    end

    background = product_distribution(Normal(-1., 1.5), Uniform())
    nested(p) = MixtureModel([tuple_product(p), background], [0.4, 0.6])
    nested_reference(p) = -2sum(log(
        0.4*pdf(Normal(p[2], 1.), x[1])*exp(p[1]*x[2])*p[1]/expm1(p[1]) +
        0.6*pdf(background, x)) for x in eachcol(samples))
    nested_cost(p) = ScientificFitting._distribution_cost(nested(p), samples)
    @test nested_cost(p) ≈ nested_reference(p) rtol=1e-8
    @test ForwardDiff.gradient(nested_cost, p) ≈ ForwardDiff.gradient(nested_reference, p) rtol=1e-6
    @test ForwardDiff.hessian(nested_cost, p) ≈ ForwardDiff.hessian(nested_reference, p) rtol=1e-5

    extended(p) = ExtendedMixtureModel([tuple_product(p), background], [5., 7.])
    extended_cost(p) = ScientificFitting._distribution_cost(extended(p), samples)
    extended_reference(p) = 24 - 2sum(log(
        5pdf(Normal(p[2], 1.), x[1])*exp(p[1]*x[2])*p[1]/expm1(p[1]) +
        7pdf(background, x)) for x in eachcol(samples))
    @test extended_cost(p) ≈ extended_reference(p) rtol=1e-8
    @test ForwardDiff.gradient(extended_cost, p) ≈ ForwardDiff.gradient(extended_reference, p) rtol=1e-6
    @test ForwardDiff.hessian(extended_cost, p) ≈ ForwardDiff.hessian(extended_reference, p) rtol=1e-5

    # Fixed measurement errors use the same preparation, including joint residuals.
    @test -2ScientificFitting._error_loglikelihood(mixture(p), data) ≈ reference(p) rtol=1e-8
    @test ScientificFitting._error_loglikelihood(tuple_product(p), samples[:, 1]) ≈
        logpdf(Normal(p[2], 1.), samples[1, 1]) + p[1]*data[1] - log(expm1(p[1])/p[1])
    native = MixtureModel([Normal(), Normal(1., 2.)], [0.4, 0.6])
    @test ScientificFitting._with_logpdf(native, 0.1) === native
    @test ScientificFitting._with_logpdf(background, samples[:, 1]) === background
    original = mixture(p)
    prepared = ScientificFitting._with_logpdf(original, data[1])
    @test components(prepared)[1].distribution === components(original)[1]
    @test !applicable(logpdf, components(original)[1], data[1])

    # Zero weights must stay differentiable, also below another mixture/product.
    for point in ([1.2, 0.], [1.2, 1.])
        @test ForwardDiff.gradient(cost, point) ≈ ForwardDiff.gradient(reference, point) rtol=1e-6
        @test ForwardDiff.hessian(cost, point) ≈ ForwardDiff.hessian(reference, point) rtol=1e-5
    end
    zero_nested(p) = MixtureModel([tuple_product(p), background], [p[2], 1-p[2]])
    zero_reference(p) = -2sum(log(
        p[2]*pdf(Normal(p[2], 1.), x[1])*exp(p[1]*x[2])*p[1]/expm1(p[1]) +
        (1-p[2])*pdf(background, x)) for x in eachcol(samples))
    zero_cost(p) = ScientificFitting._distribution_cost(zero_nested(p), samples)
    @test ForwardDiff.gradient(zero_cost, [1.2, 0.]) ≈ ForwardDiff.gradient(zero_reference, [1.2, 0.]) rtol=1e-6
    @test ForwardDiff.hessian(zero_cost, [1.2, 0.]) ≈ ForwardDiff.hessian(zero_reference, [1.2, 0.]) rtol=1e-5
    # A squared coefficient has zero gradient but nonzero curvature at zero.
    squared_cost(p) = cost([p[1], p[2]^2])
    squared_reference(p) = reference([p[1], p[2]^2])
    @test ForwardDiff.hessian(squared_cost, [1.2, 0.]) ≈
        ForwardDiff.hessian(squared_reference, [1.2, 0.]) rtol=1e-5
end
