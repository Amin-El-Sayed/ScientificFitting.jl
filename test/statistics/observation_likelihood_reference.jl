using Test, ScientificFitting, Distributions, LinearAlgebra, Statistics

@testset "User-defined observation likelihoods" begin
    line(x, p) = p[1] .* x .+ p[2]
    x = [-1.0, -0.4, 0.0, 0.5, 1.2, 1.8]
    y = [-1.1, -0.23, 0.31, 0.94, 2.04, 2.91]
    sigma = [0.12, 0.2, 0.14, 0.18, 0.11, 0.25]
    gaussian(y, mu, p) = logpdf.(Normal.(mu, sigma), y)
    # Same normalized likelihood, but independent Gaussian and generic entry points.
    native = fit_model(line, x, y; p0=[1., 0.], sigma_y=sigma, cost=:gaussian_likelihood)
    for mode in (:auto, :finite)
        result = fit_likelihood_model(line, x, y; p0=[1., 0.], logprob=gaussian, derivatives=mode)
        @test result.converged
        @test result.params ≈ native.params atol=2e-6
        @test result.param_covariance ≈ native.param_covariance rtol=2e-5
        @test result.stats.cost_min ≈ native.stats.cost_min atol=1e-8
        @test isnan(result.stats.pvalue)
        scan = profile(result, 1; values=result.params[1] .+ [-1., 0., 1.] .* result.param_stderr[1])
        @test scan.delta_cost ≈ [1., 0., 1.] atol=3e-5
    end

    @testset "Fitted scale includes its normalization" begin
        data = [-0.8, -0.2, 0.1, 0.6, 1.4, 1.9]
        constant(x, p) = fill(p[1], length(x))
        result = fit_likelihood_model(constant, collect(1.:6.), data;
            p0=[0., 1.], logprob=(y, mu, p) -> logpdf.(Normal.(mu, p[2]), y),
            bounds=([-3., 0.01], [3., 5.]))
        variance = var(data; corrected=false)
        @test result.converged
        @test result.params ≈ [mean(data), sqrt(variance)] atol=2e-6
        @test result.param_covariance ≈ Diagonal([variance/6, variance/12]) atol=2e-6
    end

    @testset "Discrete binomial measurements" begin
        trials = [8, 12, 10, 15, 9, 20]
        successes = [2, 5, 4, 7, 3, 9]
        expected = sum(successes) / sum(trials)
        result = fit_likelihood_model((x, p) -> fill(p[1], length(x)), collect(1.:6.), successes;
            logprob=(y, probability, p) -> logpdf.(Binomial.(trials, probability), y),
            p0=[0.3], bounds=([0.01], [0.99]))
        @test result.converged
        @test result.params[1] ≈ expected atol=2e-7
        @test result.param_covariance[1,1] ≈ expected*(1-expected)/sum(trials) rtol=2e-5
        @test isnan(result.stats.chi2)
    end

    @testset "Heavy tails have their own curvature, not Gaussian curvature" begin
        data = [-8., -0.8, -0.2, 0., 0.2, 0.8, 8.]
        nu = 4.
        result = fit_likelihood_model((x, p) -> fill(p[1], length(x)), collect(1.:7.), data;
            logprob=(y, mu, p) -> logpdf.(TDist(nu), y .- mu), p0=[0.1])
        hessian = sum(2*(nu+1)*(nu-r^2)/(nu+r^2)^2 for r in data)
        @test result.converged
        @test abs(result.params[1]) < 1e-6
        @test result.param_covariance[1,1] ≈ 2/hessian rtol=2e-5
        scan = profile(result, 1; values=[-1., 0., 1.])
        @test scan.delta_cost[1] ≈ scan.delta_cost[3] atol=1e-8
        @test scan.delta_cost[1] ≈ -2sum(logpdf.(TDist(nu), data .- 1)) - result.stats.cost_min
    end

    @testset "Invalid probabilities are not repaired" begin
        for invalid in ((y, mu, p) -> [0.], (y, mu, p) -> 0.,
                        (y, mu, p) -> fill(NaN, length(y)), (y, mu, p) -> fill(Inf, length(y)))
            @test_throws ArgumentError fit_likelihood_model(line, x, y; p0=[1., 0.],
                logprob=invalid, fixed_parameters=[1=>1., 2=>0.])
        end
        @test_throws ArgumentError fit_likelihood_model(line, x[1:2], y; p0=[1., 0.], logprob=gaussian)
        @test_throws ArgumentError fit_likelihood_model(line, Float64[], Float64[];
            p0=[1., 0.], logprob=gaussian)
        restricted(y, mu, p) = p[1] > 0 ? gaussian(y, mu, p) : fill(-Inf, length(y))
        result = fit_likelihood_model(line, x, y; p0=[1., 0.], logprob=restricted,
            fixed_parameters=[1=>1., 2=>0.])
        @test result.problem.objective([-1., 0.]) == Inf
        @test isfinite(result.stats.cost_min)
    end
end
