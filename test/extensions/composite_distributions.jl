using ScientificFitting
using Distributions
using DistributionsHEP
using NumericalDistributions
using ForwardDiff
using Test

struct BatchObservedNormal{T} <: ContinuousUnivariateDistribution
    location::T
    batch_sizes::Vector{Int}
end
Distributions.logpdf(d::BatchObservedNormal, x::Real) = logpdf(Normal(d.location, 1.), x)
function Distributions.logpdf(d::BatchObservedNormal, x::AbstractVector{<:Real})
    push!(d.batch_sizes, length(x))
    return logpdf.(Ref(Normal(d.location, 1.)), x)
end

@testset "Event reductions bound scratch batches without losing observations" begin
    sizes = Int[]
    data = collect(range(-2., 2.; length=10_003))
    model(p) = MixtureModel([BatchObservedNormal(p[1], sizes), Normal(-1., 2.)], [p[2], 1-p[2]])
    actual(p) = ScientificFitting._distribution_cost(model(p), data)
    reference(p) = -2sum(log(p[2]*pdf(Normal(p[1], 1.), x) +
        (1-p[2])*pdf(Normal(-1., 2.), x)) for x in data)
    for point in ([.4, .3], [.4, 0.])
        empty!(sizes)
        @test actual(point) ≈ reference(point) rtol=1e-12
        point[2] == 0 || @test sum(sizes) == length(data)
        @test ForwardDiff.gradient(actual, point) ≈ ForwardDiff.gradient(reference, point) rtol=1e-10
        @test ForwardDiff.hessian(actual, point) ≈ ForwardDiff.hessian(reference, point) rtol=1e-9
        @test maximum(sizes) <= 4096
    end

    # Joint events are columns: blocking must never split their coordinates.
    joint = vcat(data', (data ./ 2)')
    joint_model(p) = MixtureModel([
        product_distribution(Normal(p[1], 1.), Normal()),
        product_distribution(Uniform(-3., 3.), Normal(1., 2.))], [.3, .7])
    joint_cost(p) = ScientificFitting._distribution_cost(joint_model(p), joint)
    joint_reference(p) = -2sum(logpdf(joint_model(p), x) for x in eachcol(joint))
    @test joint_cost([.4]) ≈ joint_reference([.4]) rtol=1e-12
    @test ForwardDiff.gradient(joint_cost, [.4]) ≈ ForwardDiff.gradient(joint_reference, [.4]) rtol=1e-10
    @test ForwardDiff.hessian(joint_cost, [.4]) ≈ ForwardDiff.hessian(joint_reference, [.4]) rtol=1e-9
end

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

@testset "Mixture batches preserve densities, derivatives and linear storage" begin
    points = [-5., -2., -0.1, 0.4, 2., 5.]
    peak(p) = truncated(Normal(p[1], p[2]), -1., 1.)
    background = Normal(2., 3.)
    model(p) = MixtureModel([peak(p), background], [p[3], 1-p[3]])
    actual(p) = ScientificFitting._distribution_cost(model(p), points)
    reference(p) = -2sum(log(p[3]*pdf(peak(p), x) + (1-p[3])*pdf(background, x)) for x in points)
    # Some observations lie outside the signal window. The background still
    # gives them positive probability, also at a zero signal yield.
    for p in ([0.2, 0.7, 0.4], [0.2, 0.7, 0.])
        @test actual(p) ≈ reference(p) rtol=1e-12
        @test ForwardDiff.gradient(actual, p) ≈ ForwardDiff.gradient(reference, p) atol=1e-10
        @test ForwardDiff.hessian(actual, p) ≈ ForwardDiff.hessian(reference, p) atol=1e-9
    end
    # At the first observation the variable peak vanishes, and both weights
    # are fixed. Later observations must still retain its parameter derivatives.
    fixed_weights(p) = MixtureModel([peak(p), background], [0.4, 0.6])
    fixed_cost(p) = ScientificFitting._distribution_cost(fixed_weights(p), points)
    fixed_reference(p) = reference([p[1], p[2], 0.4])
    @test ForwardDiff.hessian(fixed_cost, [0.2, 0.7]) ≈
        ForwardDiff.hessian(fixed_reference, [0.2, 0.7]) atol=1e-9

    # Both components vanish at the first point; a later component restores
    # support. Nested and flat representations describe the same probability.
    inner = MixtureModel([Uniform(0., 1.), Uniform(2., 3.)], [0.25, 0.75])
    nested = MixtureModel([inner, Normal()], [0.4, 0.6])
    flat = MixtureModel([Uniform(0., 1.), Uniform(2., 3.), Normal()], [0.1, 0.3, 0.6])
    @test ScientificFitting._distribution_cost(nested, points) ≈ -2loglikelihood(flat, points)
    @test ScientificFitting._distribution_cost(inner, [-1.]) == Inf
    tail = MixtureModel([Normal(), Normal(1000.)], [0., 1.])
    @test ScientificFitting._distribution_cost(tail, [0.]) ≈ -2logpdf(Normal(1000.), 0.)
    tiny_weight = MixtureModel([Normal(), Normal(1000.)], [1e-300, 1.])
    @test ScientificFitting._distribution_cost(tiny_weight, [0., 500., 1000.]) ≈
        -2loglikelihood(tiny_weight, [0., 500., 1000.])
    discrete = MixtureModel([Poisson(0.5), Poisson(3.)], [0.3, 0.7])
    @test ScientificFitting._distribution_cost(discrete, [0, 1, 5]) ≈ -2loglikelihood(discrete, [0, 1, 5])
    # Even a mutated upstream probability vector must not become a signed density.
    invalid = MixtureModel([Normal(), Uniform()], [0.3, 0.7])
    probs(invalid)[1] = -0.1
    @test_throws ArgumentError ScientificFitting._distribution_cost(invalid, [0.1, 0.3])

    error = MixtureModel([Normal(), Uniform(-10., 10.)], [0.9, 0.1])
    residual_cost(p) = -2ScientificFitting._error_loglikelihood(error, points .- p[1])
    residual_reference(p) = -2sum(log(0.9pdf(Normal(), x-p[1]) + 0.005) for x in points)
    @test ForwardDiff.gradient(residual_cost, [0.2]) ≈ ForwardDiff.gradient(residual_reference, [0.2])
    @test ForwardDiff.hessian(residual_cost, [0.2]) ≈ ForwardDiff.hessian(residual_reference, [0.2])

    data = collect(range(-4., 4.; length=10_000))
    mixture = MixtureModel([Normal(), Uniform(-4., 4.)], [0.6, 0.4])
    evaluate(data) = ScientificFitting._distribution_cost(mixture, data)
    evaluate(data) # Warm specialization before checking allocation volume.
    @test (@allocated evaluate(data)) < 100length(data) + 100_000
    # A concrete component type already has a cheap native loop. Keep it rather
    # than introducing event-sized scratch arrays for every mixture.
    homogeneous = MixtureModel([Normal(), Normal(2.)], [0.6, 0.4])
    homogeneous_cost(data) = ScientificFitting._distribution_cost(homogeneous, data)
    native_cost(data) = -2loglikelihood(homogeneous, data)
    homogeneous_cost(data); native_cost(data)
    @test (@allocated homogeneous_cost(data)) < (@allocated native_cost(data)) + 4096
end
