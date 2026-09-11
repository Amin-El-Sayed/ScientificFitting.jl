using ScientificFitting
using Distributions
using ForwardDiff
using Test

struct PDFOnlySlope{T} <: ContinuousUnivariateDistribution
    rate::T
end
Distributions.pdf(d::PDFOnlySlope, x::Real) =
    0 <= x <= 1 ? exp(d.rate*x)*d.rate/expm1(d.rate) : zero(x)

@testset "Distribution histogram likelihood" begin
    edges, counts = [-2., -0.8, 0.2, 1.5, 3.], [5, 15, 18, 7]
    builds = Ref(0)
    factory(p) = (builds[] += 1; Normal(p[1], 1.))
    r = fit_distribution(factory, edges, counts; p0=[0.], total_count=50.)
    reference = fit_histogram_model((e, p) -> 50diff(cdf.(Normal(p[1], 1.), e)), edges, counts; p0=[0.])
    @test r.converged
    @test r.params ≈ reference.params atol=1e-6
    @test r.param_covariance ≈ reference.param_covariance rtol=1e-6
    @test r.stats.cost_min ≈ reference.stats.cost_min atol=1e-8
    @test r.stats.chi2 ≈ reference.stats.chi2 atol=1e-8
    @test params(fitted_model(r)) == (r.params[1], 1.)
    before = builds[]
    expected = r.problem.objective([0.2])
    @test builds[] == before + 1
    edges[1], counts[1] = -20., 100
    @test r.problem.objective([0.2]) == expected

    # The Bernoulli MLE and curvature check discrete endpoint semantics as well.
    binomial = fit_distribution(p -> Bernoulli(p[1]), [-0.5, 0.5, 1.5], [5, 7];
        p0=[0.5], total_count=12., bounds=([0.01], [0.99]))
    @test binomial.converged
    @test binomial.params[1] ≈ 7/12 atol=1e-6
    @test binomial.param_covariance[1, 1] ≈ (7/12)*(5/12)/12 rtol=1e-5

    bins = ScientificFitting.DistributionHistogram([-0.4, 0.1, 1.4, 2.], [2., 7., 3.], 15., :cdf, 1e-9)
    quadrature = ScientificFitting.DistributionHistogram(bins.edges, bins.counts, 15., :quadgk, 1e-9)
    normal(p) = Normal(p[1], exp(p[2]))
    cdf_cost(p) = ScientificFitting._distribution_cost(normal(p), bins)
    quad_cost(p) = ScientificFitting._distribution_cost(normal(p), quadrature)
    p = [0.2, 0.1]
    @test cdf_cost(p) ≈ quad_cost(p) atol=1e-8
    @test ForwardDiff.gradient(cdf_cost, p) ≈ ForwardDiff.gradient(quad_cost, p) rtol=1e-7
    @test ForwardDiff.hessian(cdf_cost, p) ≈ ForwardDiff.hessian(quad_cost, p) rtol=1e-6

    # Log bin probabilities must survive ordinary probability underflow.
    tail = ScientificFitting.DistributionHistogram([40., 40.1], [1.], 1., :cdf, 1e-8)
    tail_cost(p) = ScientificFitting._distribution_cost(Normal(p[1], 1.), tail)
    @test tail_cost([0.]) ≈ -2logdiffcdf(Normal(), 40.1, 40.)
    @test all(isfinite, ForwardDiff.gradient(tail_cost, [0.]))
    @test all(isfinite, ForwardDiff.hessian(tail_cost, [0.]))
    # An inactive peak must not set the scale and underflow the active tail.
    inactive = MixtureModel([Normal(40., 1.), Normal()], [0., 1.])
    @test ScientificFitting._distribution_cost(inactive, tail) ≈ tail_cost([0.])
    @test ScientificFitting._distribution_cost(inactive, [40.]) ≈ -2logpdf(Normal(), 40.)

    empty_bins = fit_distribution(p -> Uniform(0., p[1]), [-2., -1., 0., 0.5, 1., 2.], [0, 0, 5, 7, 0];
        p0=[1.], total_count=12., fixed_parameters=1 => 1.)
    @test isfinite(empty_bins.stats.cost_min)
    @test isnan(empty_bins.stats.chi2) && isnan(empty_bins.stats.pvalue)
    @test ScientificFitting._poisson_minus2loglik_terms([0.], [0.]) == 0
    @test ScientificFitting._poisson_minus2loglik_terms([1.], [0.]) == Inf
    zero_rate = fit_poisson_model((x, p) -> [0., exp(p[1])], [1., 2.], [0, 4]; p0=[1.])
    @test zero_rate.params[1] ≈ log(4) atol=1e-6
    @test isnan(zero_rate.stats.pvalue)

    @test_throws ArgumentError fit_distribution(normal, [0., 1.], [2]; p0=p)
    @test_throws ArgumentError fit_distribution(normal, [0., 1.], [2]; p0=p, total_count=-1.)
    @test_throws ArgumentError fit_distribution(normal, [0., 1.], [2]; p0=p, total_count=3., rtol=0.)
    @test_throws ArgumentError fit_distribution(normal, [0., 1.], [2]; p0=p, total_count=3., integration=:guess)
    @test_throws ArgumentError fit_distribution(normal, [1., 0.], [2]; p0=p, total_count=3.)
    @test_throws ArgumentError fit_distribution(normal, [0., 1.], [2.5]; p0=p, total_count=3.)
    @test_throws ArgumentError fit_distribution(p -> Bernoulli(p[1]), [-0.5, 0.5, 1.5], [2, 3];
        p0=[0.5], total_count=5., integration=:quadgk)
end

@testset "Truncated distribution bins preserve boundary derivatives" begin
    bins = ScientificFitting.DistributionHistogram([0., 0.5], [1.], 1., :cdf, 1e-9)
    factory(p) = truncated(Normal(p[1], exp(p[2])), 0., 10.)
    logmass(p) = ScientificFitting._bin_logmass(factory(p), 0., 0.5, bins)
    # This parameter point previously made the lower-edge CDF derivative NaN.
    reference(p) = log(ScientificFitting.quadgk(x -> pdf(factory(p), x), 0., 0.5; rtol=1e-9)[1])
    point = [5.358013439054251, log(0.4)]
    @test logmass(point) ≈ reference(point) atol=1e-10
    @test ForwardDiff.gradient(logmass, point) ≈ ForwardDiff.gradient(reference, point) rtol=1e-8
    @test ForwardDiff.hessian(logmass, point) ≈ ForwardDiff.hessian(reference, point) rtol=1e-7

    # Bins may straddle either edge, including one-sided truncation.
    for d in (truncated(Normal(), 0., 2.), truncated(Normal(); lower=0.),
              truncated(Normal(); upper=2.))
        actual = exp(ScientificFitting._bin_logmass(d, -1., 3., bins))
        @test actual ≈ cdf(d, 3.) - cdf(d, -1.) rtol=1e-12
    end
    @test ScientificFitting._bin_logmass(factory(point), -1., 0., bins) == -Inf
    @test ScientificFitting._bin_logmass(factory(point), 10., 11., bins) == -Inf
end

@testset "PDF-only bin integration retains normalization derivatives" begin
    edges = [0., 0.2, 0.6, 1.]
    bins = ScientificFitting.DistributionHistogram(edges, [4., 8., 13.], 25., :auto, 1e-9)
    factory(p) = MixtureModel([PDFOnlySlope(p[1]), Uniform()], [p[2], 1-p[2]])
    cost(p) = ScientificFitting._distribution_cost(factory(p), bins)
    # Exact exponential bin probabilities, mixed with a uniform background.
    expected(p) = [25*(p[2]*(exp(p[1]*b)-exp(p[1]*a))/expm1(p[1]) + (1-p[2])*(b-a))
                   for (a, b) in zip(edges[1:end-1], edges[2:end])]
    reference(p) = ScientificFitting._poisson_minus2loglik_terms(bins.counts, expected(p))
    p = [1.2, 0.4]
    @test cost(p) ≈ reference(p) atol=1e-8
    @test ForwardDiff.gradient(cost, p) ≈ ForwardDiff.gradient(reference, p) rtol=1e-6
    @test ForwardDiff.hessian(cost, p) ≈ ForwardDiff.hessian(reference, p) rtol=1e-5
    @test_throws ArgumentError fit_distribution(p -> PDFOnlySlope(p[1]), edges, [4, 8, 13];
        p0=[1.], total_count=25., integration=:cdf)
end
