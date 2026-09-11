using ScientificFitting
using Distributions
using ForwardDiff
using LinearAlgebra
using Statistics
using Test

@testset "Distribution factories use the complete event likelihood" begin
    data = [-0.8, -0.3, 0.1, 0.6, 0.9, 1.4]
    make_normal(p) = Normal(p[1], exp(p[2]))
    r = fit_distribution(make_normal, data; p0=[0., 0.], parameter_names=["mu", "log_sigma"])
    mu, variance = mean(data), var(data; corrected=false)
    @test r.converged
    @test r.params ≈ [mu, log(sqrt(variance))] atol=1e-6
    @test r.param_covariance ≈ Diagonal([variance/length(data), 1/(2length(data))]) atol=1e-6
    @test r.stats.cost_min ≈ -2loglikelihood(make_normal(r.params), data)
    @test r.problem.nobs == length(data)
    @test isnan(r.stats.pvalue)
    @test occursin("mu", report_text(r))
    old_cost = r.problem.objective(r.params)
    data[1] = 100.
    @test r.problem.objective(r.params) == old_cost

    counts = [0, 1, 2, 0, 3, 5, 1]
    poisson = fit_distribution(p -> Poisson(exp(p[1])), counts; p0=[0.])
    @test poisson.converged
    @test typeof(poisson) === typeof(r)
    @test exp(poisson.params[1]) ≈ mean(counts) atol=1e-6
    @test poisson.param_covariance[1, 1] ≈ 1/sum(counts) rtol=1e-5

    # Each column is one jointly modelled event, not two independent scalars.
    samples = [0.1 0.8 -0.2 1.2; 1.1 0.5 1.8 1.0]
    cov = [1. 0.4; 0.4 2.]
    joint(p) = MvNormal(p, cov)
    columns = fit_distribution(joint, samples; p0=[0., 0.], obsdim=2)
    rows = fit_distribution(joint, permutedims(samples); p0=[0., 0.], obsdim=1)
    @test columns.params ≈ vec(mean(samples; dims=2)) atol=1e-6
    @test columns.param_covariance ≈ cov/4 atol=1e-6
    @test columns.params ≈ rows.params atol=1e-8
    @test columns.problem.nobs == 4
    @test columns.stats.cost_min ≈ rows.stats.cost_min

    # Stable logpdf matters even when pdf itself underflows to zero.
    @test pdf(Normal(), 40.) == 0
    @test isfinite(ScientificFitting._distribution_cost(Normal(), [40.]))
    @test ScientificFitting._distribution_cost(Uniform(), [-1.]) == Inf
    @test_throws ArgumentError fit_distribution(joint, samples; p0=[0., 0.])
    @test_throws ArgumentError fit_distribution(joint, [1., 2.]; p0=[0., 0.])
    @test_throws ArgumentError fit_distribution(make_normal, samples; p0=[0., 0.], obsdim=2)
    @test_throws ArgumentError fit_distribution(make_normal, [NaN]; p0=[0., 0.])
    @test_throws ArgumentError fit_distribution(make_normal, Float64[]; p0=[0., 0.])
end

@testset "Fixed error objects preserve independence or joint dependence" begin
    x, y = [0., 1., 2., 3.], [0.1, 1.2, 1.8, 3.4]
    line(x, p) = @. p[1]*x + p[2]
    sigma = [0.2, 0.3, 0.2, 0.4]
    d = Normal.(0., sigma)
    r = fit_likelihood_model(line, x, y; error=d, p0=[1., 0.])
    gaussian = fit_model(line, x, y; sigma_y=sigma, p0=[1., 0.])
    @test r.converged
    @test r.params ≈ gaussian.params atol=1e-6
    @test r.param_covariance ≈ gaussian.param_covariance rtol=1e-5
    @test r.problem.objective([1., 0.]) ≈ -2sum(logpdf.(d, y.-x))

    covariance = Matrix(Diagonal(sigma.^2)) .+ 0.01
    correlated = fit_likelihood_model(line, x, y; error=MvNormal(zeros(4), covariance), p0=[1., 0.])
    reference = fit_model(line, x, y; cov_y=covariance, p0=[1., 0.])
    @test correlated.params ≈ reference.params atol=1e-6
    @test correlated.param_covariance ≈ reference.param_covariance rtol=1e-5

    heavy = fit_likelihood_model(line, x, y; error=TDist(4), p0=[1., 0.])
    @test heavy.converged
    @test heavy.problem.objective([1., 0.]) ≈ -2sum(logpdf.(TDist(4), y.-x))
    @test_throws ArgumentError fit_likelihood_model(line, x, y; p0=[1., 0.])
    @test_throws ArgumentError fit_likelihood_model(line, x, y;
        error=Normal(), logprob=(y, mu, p) -> y, p0=[1., 0.])
    @test_throws ArgumentError fit_likelihood_model(line, x, y; error=[Normal()], p0=[1., 0.])
    @test_throws ArgumentError fit_likelihood_model(line, x, y; error=MvNormal(zeros(2), I), p0=[1., 0.])
end
