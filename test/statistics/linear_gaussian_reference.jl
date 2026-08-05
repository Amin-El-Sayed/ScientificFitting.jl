using Distributions
using JuFitter
using LinearAlgebra
using SparseArrays
using Test

function _linear_design(x)
    return hcat(Float64.(x), ones(length(x)))
end

function _weighted_linear_reference(x, y, sigma_y)
    X = _linear_design(x)
    w = 1.0 ./ abs2.(sigma_y)
    fisher = X' * (w .* X)
    rhs = X' * (w .* y)
    params = Symmetric(fisher) \ rhs
    cov = Symmetric(fisher) \ Matrix{Float64}(I, size(fisher))
    residuals = y .- X * params
    chi2 = sum(abs2, residuals ./ sigma_y)
    return (; params, cov, residuals, chi2)
end

function _generalized_linear_reference(x, y, cov_y)
    X = _linear_design(x)
    F = cholesky(Symmetric(Matrix(cov_y)))
    whitened_X = F.L \ X
    whitened_y = F.L \ y
    fisher = whitened_X' * whitened_X
    rhs = whitened_X' * whitened_y
    params = Symmetric(fisher) \ rhs
    cov = Symmetric(fisher) \ Matrix{Float64}(I, size(fisher))
    residuals = y .- X * params
    z = F.L \ residuals
    chi2 = sum(abs2, z)
    logdet = 2.0 * sum(log, diag(F.L))
    return (; params, cov, residuals, chi2, logdet)
end

@testset "Linear Gaussian statistical references" begin
    @testset "Diagonal weighted least squares matches analytic solution" begin
        x = collect(range(-2.0, 3.0; length=17))
        sigma_y = @. 0.11 + 0.015 * (x + 2.0)
        y = @. 1.4 * x - 0.35 + 0.07 * sin(2.3 * x)
        model(x, p) = @. p[1] * x + p[2]

        ref = _weighted_linear_reference(x, y, sigma_y)
        result = fit_model(model, x, y; p0=[0.0, 0.0], sigma_y=sigma_y, scale_covariance=:never)
        ndf = length(x) - 2
        expected_minus2loglik = length(x) * log(2pi) + 2.0 * sum(log, sigma_y) + ref.chi2

        @test result.backend == :lsqfit
        @test result.stats.cost == :chi2
        @test isapprox(result.params, ref.params; atol=2e-9, rtol=2e-9)
        @test isapprox(result.param_covariance, ref.cov; atol=2e-11, rtol=2e-11)
        @test isapprox(result.param_stderr, sqrt.(diag(ref.cov)); atol=2e-11, rtol=2e-11)
        @test isapprox(result.residuals, ref.residuals; atol=2e-9, rtol=2e-9)
        @test isapprox(result.weighted_residuals, ref.residuals ./ sigma_y; atol=2e-9, rtol=2e-9)
        @test isapprox(result.stats.chi2, ref.chi2; atol=2e-11, rtol=2e-11)
        @test isapprox(result.stats.cost_min, ref.chi2; atol=2e-11, rtol=2e-11)
        @test isapprox(result.stats.minus2loglik_min, expected_minus2loglik; atol=2e-11, rtol=2e-11)
        @test !isapprox(result.stats.cost_min, result.stats.minus2loglik_min; atol=1e-8, rtol=1e-8)
        @test isapprox(result.stats.aic, expected_minus2loglik + 2 * length(result.params); atol=2e-11, rtol=2e-11)
        @test isapprox(
            result.stats.bic,
            expected_minus2loglik + log(length(x)) * length(result.params);
            atol=2e-11,
            rtol=2e-11,
        )
        @test result.stats.ndf == ndf
        @test isapprox(result.stats.chi2_ndf, ref.chi2 / ndf; atol=2e-11, rtol=2e-11)
        @test isapprox(result.stats.pvalue, ccdf(Chisq(ndf), ref.chi2); atol=2e-11, rtol=2e-11)
    end

    @testset "Full covariance Gaussian -2log(L) matches generalized least squares" begin
        x = collect(range(-1.5, 2.5; length=12))
        sigma = @. 0.08 + 0.01 * (x + 1.5)
        rho = 0.42
        cov_y = [sigma[i] * sigma[j] * rho^abs(i - j) for i in eachindex(x), j in eachindex(x)]
        y = @. -0.8 * x + 0.9 + 0.03 * cos(1.7 * x)
        model(x, p) = @. p[1] * x + p[2]

        ref = _generalized_linear_reference(x, y, cov_y)
        result = fit_model(model, x, y; p0=[0.0, 0.0], cov_y=cov_y, cost=:gaussian_likelihood, scale_covariance=:never)
        ndf = length(x) - 2
        expected_minus2loglik = length(x) * log(2pi) + ref.logdet + ref.chi2

        @test result.backend == :optimization
        @test result.stats.cost == :gaussian_likelihood
        @test isapprox(result.params, ref.params; atol=2e-8, rtol=2e-8)
        @test isapprox(result.param_covariance, ref.cov; atol=2e-7, rtol=2e-7)
        @test isapprox(result.stats.chi2, ref.chi2; atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.minus2loglik_min, expected_minus2loglik; atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.cost_min, expected_minus2loglik; atol=2e-8, rtol=2e-8)
        @test result.stats.ndf == ndf
        @test isapprox(result.stats.pvalue, ccdf(Chisq(ndf), ref.chi2); atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.aic, expected_minus2loglik + 2 * length(result.params); atol=2e-8, rtol=2e-8)
        @test isapprox(
            result.stats.bic,
            expected_minus2loglik + log(length(x)) * length(result.params);
            atol=2e-8,
            rtol=2e-8,
        )
    end

    @testset "Sparse covariance keeps the static whitening path sparse" begin
        n = 40
        x = collect(range(-2.0, 3.0; length=n))
        covariance = spdiagm(
            -1 => fill(0.006, n - 1),
            0 => (@. 0.055 + 0.002 * (x + 2.0)),
            1 => fill(0.006, n - 1),
        )
        y = @. 0.82 * x - 0.35 + 0.015 * cos(1.6 * x)
        model(x, p) = @. p[1] * x + p[2]
        jacobian(x, p) = hcat(x, ones(length(x)))

        ref = _generalized_linear_reference(x, y, covariance)
        result = fit_model(
            model,
            x,
            y;
            p0=[0.0, 0.0],
            cov_y=covariance,
            jacobian=jacobian,
            scale_covariance=:never,
        )
        cache = JuFitter._prepare_fit_cache(result.problem)

        @test result.backend == :lsqfit
        @test result.problem.cov_y isa SparseMatrixCSC
        @test cache.covariance.factor isa SparseArrays.CHOLMOD.Factor
        @test isapprox(result.params, ref.params; atol=2e-8, rtol=2e-8)
        @test isapprox(result.param_covariance, ref.cov; atol=2e-7, rtol=2e-7)
        @test isapprox(result.stats.chi2, ref.chi2; atol=2e-8, rtol=2e-8)
        @test isapprox(JuFitter._yerror_for_plot(result.problem), sqrt.(diag(covariance)); atol=1e-14)

        @test_throws ArgumentError fit_model(
            model,
            x,
            y;
            p0=[0.0, 0.0],
            cov_y=covariance,
            bounds=([-2.0, -2.0], [2.0, 2.0]),
        )
    end
end
