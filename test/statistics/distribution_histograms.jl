using ScientificFitting
using Distributions
using ForwardDiff
using Test

struct PDFOnlySlope{T} <: ContinuousUnivariateDistribution
    rate::T
end
Distributions.pdf(d::PDFOnlySlope, x::Real) =
    0 <= x <= 1 ? exp(d.rate*x)*d.rate/expm1(d.rate) : zero(x)

struct RoundedCDF{T,E} <: ContinuousUnivariateDistribution
    normal::Normal{T}
    offset::E
end
Distributions.pdf(d::RoundedCDF, x::Real) = pdf(d.normal, x)
Distributions.cdf(d::RoundedCDF, x::Real) = cdf(d.normal, x) + d.offset

struct PDFOnlyOscillation{T} <: ContinuousUnivariateDistribution
    amplitude::T
end
function Distributions.pdf(d::PDFOnlyOscillation, x::Real)
    p = d.amplitude
    0 <= x <= 1 || return zero(p+x)
    normalization = 1 + p*(1-cos(73.))/73 + p^2*sin(113.)/113
    return (1 + p*sin(73x) + p^2*cos(113x)) / normalization
end

@testset "Quadrature preserves boundary and integrand derivatives" begin
    SF = ScientificFitting
    bins = SF.DistributionHistogram([-.5,.5], [1.], 1., :quadgk, 1e-10)
    model(p) = truncated(Normal(p[1], exp(p[2])), p[3], 2.)
    actual(p) = SF._bin_logmass(model(p), -.5, .5, bins)
    reference(p) = logdiffcdf(Normal(p[1], exp(p[2])), .5, p[3]) -
                   logdiffcdf(Normal(p[1], exp(p[2])), 2., p[3])
    p = [.1, -.4, -.2]
    @test actual(p) ≈ reference(p) atol=1e-10
    @test ForwardDiff.gradient(actual, p) ≈ ForwardDiff.gradient(reference, p) atol=1e-9
    @test ForwardDiff.hessian(actual, p) ≈ ForwardDiff.hessian(reference, p) atol=1e-8

    # At amplitude zero the density is constant, but its derivatives oscillate.
    # Adaptive accuracy must include those derivatives, not only the density.
    a, b = .13, .92
    integral(p) = exp(SF._bin_logmass(PDFOnlyOscillation(p[1]), a, b, bins))
    exact(p) = ((b-a) + p[1]*(cos(73a)-cos(73b))/73 +
        p[1]^2*(sin(113b)-sin(113a))/113) /
        (1 + p[1]*(1-cos(73.))/73 + p[1]^2*sin(113.)/113)
    @test integral([0.]) ≈ exact([0.]) atol=1e-10
    @test ForwardDiff.gradient(integral, [0.]) ≈ ForwardDiff.gradient(exact, [0.]) atol=1e-9
    @test ForwardDiff.hessian(integral, [0.]) ≈ ForwardDiff.hessian(exact, [0.]) atol=1e-8
end

@testset "CDF roundoff uses integration, not probability clipping" begin
    SF = ScientificFitting
    bins(mode) = SF.DistributionHistogram([8., 8.1], [1.], 1., mode, 1e-10)
    model(p) = RoundedCDF(Normal(p[1], exp(p[2])), 4eps())
    actual(p) = SF._bin_logmass(model(p), 8., 8.1, bins(:auto))
    reference(p) = logdiffcdf(Normal(p[1], exp(p[2])), 8.1, 8.)
    @test actual([0., 0.]) ≈ reference([0., 0.]) atol=1e-10
    @test ForwardDiff.gradient(actual, [0., 0.]) ≈ ForwardDiff.gradient(reference, [0., 0.]) atol=1e-8
    @test ForwardDiff.hessian(actual, [0., 0.]) ≈ ForwardDiff.hessian(reference, [0., 0.]) atol=1e-7
    @test_throws ArgumentError SF._bin_logmass(model([0., 0.]), 8., 8.1, bins(:cdf))
    for offset in (1e-3, NaN, Inf)
        @test_throws ArgumentError SF._bin_logmass(RoundedCDF(Normal(), offset), 8., 8.1, bins(:auto))
    end
end

@testset "Binned mixtures batch components without losing derivatives" begin
    SF = ScientificFitting
    edges = [-3., -1., -0.2, 0.4, 1., 3.]
    bins = SF.DistributionHistogram(edges, [1., 3., 5., 4., 2.], 15., :cdf, 1e-9)
    model(p) = MixtureModel([
        truncated(MixtureModel([Normal(p[1], exp(p[2])), Normal(p[1], 2exp(p[2]))],
            [0.7, 0.3]), -1., 1.), Normal(2., 3.)], [p[3], 1-p[3]])
    actual(p) = SF._distribution_cost(model(p), bins)
    # Preserve the independently evaluated scalar path, including clipping at
    # both selection edges and bins outside the signal's support.
    function reference(p)
        logs = [log(15.) + SF._bin_logmass(model(p), edges[i], edges[i+1], bins)
                for i in eachindex(bins.counts)]
        SF._poisson_minus2loglik_terms(bins.counts, exp.(logs); log_mu=logs)
    end
    for point in ([0.2, -0.3, 0.4], [0.2, -0.3, 0.])
        @test actual(point) ≈ reference(point) rtol=1e-12
        @test ForwardDiff.gradient(actual, point) ≈ ForwardDiff.gradient(reference, point) atol=1e-9
        @test ForwardDiff.hessian(actual, point) ≈ ForwardDiff.hessian(reference, point) atol=1e-8
    end
    @test ForwardDiff.hessian(p -> actual([p[1], p[2], p[3]^2]), [0.2, -0.3, 0.]) ≈
        ForwardDiff.hessian(p -> reference([p[1], p[2], p[3]^2]), [0.2, -0.3, 0.]) atol=1e-8

    # Dispatch and scratch storage belong to each component batch, not every
    # component in every bin. This is an allocation contract, not a time limit.
    n = 10_000
    many = SF.DistributionHistogram(collect(range(-3., 3.; length=n+1)),
        ones(n), Float64(n), :cdf, 1e-9)
    fixed_model = model([0.2, -0.3, 0.4])
    evaluate() = SF._distribution_cost(fixed_model, many)
    evaluate()
    @test (@allocated evaluate()) < 256n + 100_000
end

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

    # Exponential's lower support is already zero. A redundant truncation there
    # can poison the upstream log-normalizer's AD, despite valid scalar values.
    edges, counts = [0., 0.5, 1., 2.], [20, 12, 10]
    upper(p) = truncated(Exponential(p[1]); upper=2.)
    redundant(p) = truncated(Exponential(p[1]), 0., 2.)
    a = fit_distribution(upper, edges, counts; p0=[1.], total_count=42., bounds=([0.1], [10.]))
    b = fit_distribution(redundant, edges, counts; p0=[1.], total_count=42.,
        bounds=([0.1], [10.]), derivatives=:finite)
    @test a.converged && b.converged
    @test a.params ≈ b.params atol=1e-5
    @test a.param_covariance ≈ b.param_covariance rtol=1e-4
    @test a.stats.cost_min ≈ b.stats.cost_min atol=1e-8
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
