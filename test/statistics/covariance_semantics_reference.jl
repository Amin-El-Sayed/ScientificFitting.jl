using Distributions
using JuFitter
using LinearAlgebra
using Test

function _linear_reference_from_covariance(x, y, covariance)
    X = hcat(Float64.(x), ones(length(x)))
    F = cholesky(Symmetric(Matrix(covariance)))
    whitened_X = F.L \ X
    whitened_y = F.L \ y
    fisher = whitened_X' * whitened_X
    rhs = whitened_X' * whitened_y
    params = Symmetric(fisher) \ rhs
    cov = Symmetric(fisher) \ Matrix{Float64}(I, size(fisher))
    residuals = y .- X * params
    chi2 = sum(abs2, F.L \ residuals)
    return (; params, cov, residuals, chi2)
end

@testset "Covariance and constraint semantics" begin
    @testset "Covariance scaling policy matches residual variance estimate" begin
        x = collect(range(-1.0, 2.0; length=18))
        y = @. 0.9 * x + 0.2 + 0.05 * cos(2.1 * x)
        model(x, p) = @. p[1] * x + p[2]

        unit_cov = Matrix{Float64}(I, length(x), length(x))
        ref = _linear_reference_from_covariance(x, y, unit_cov)
        ndf = length(x) - 2

        auto_scaled = fit_model(model, x, y; p0=[0.0, 0.0])
        never_scaled = fit_model(model, x, y; p0=[0.0, 0.0], scale_covariance=:never)

        @test isapprox(auto_scaled.params, ref.params; atol=2e-9, rtol=2e-9)
        @test isapprox(never_scaled.param_covariance, ref.cov; atol=2e-11, rtol=2e-11)
        @test isapprox(auto_scaled.param_covariance, ref.cov .* (ref.chi2 / ndf); atol=2e-11, rtol=2e-11)
        @test auto_scaled.stats.ndf == ndf
        @test isapprox(auto_scaled.stats.chi2_ndf, ref.chi2 / ndf; atol=2e-11, rtol=2e-11)
    end

    @testset "Static error components add as covariance terms" begin
        x = collect(range(0.2, 2.4; length=14))
        y = @. 1.1 * x - 0.3 + 0.02 * sin(3.0 * x)
        model(x, p) = @. p[1] * x + p[2]
        sigma_abs = 0.04
        sigma_extra = @. 0.015 + 0.002 * x
        rel = 0.012
        total_variance = @. sigma_abs^2 + (rel * abs(y))^2 + sigma_extra^2

        ref = _linear_reference_from_covariance(x, y, Diagonal(total_variance))
        result = fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            error_components=[
                (name=:absolute, target=:y, mode=:absolute, values=sigma_abs),
                (name=:relative, target=:y, mode=:relative, values=rel),
                (name=:extra, target=:y, mode=:covariance, values=sigma_extra),
            ],
            scale_covariance=:never,
        )
        ndf = length(x) - 2

        @test result.backend == :optimization
        @test result.stats.cost == :chi2
        @test isapprox(result.params, ref.params; atol=2e-8, rtol=2e-8)
        @test isapprox(result.param_covariance, ref.cov; atol=2e-7, rtol=2e-7)
        @test isapprox(result.stats.chi2, ref.chi2; atol=2e-8, rtol=2e-8)
        @test result.stats.ndf == ndf
        @test isapprox(result.stats.pvalue, ccdf(Chisq(ndf), ref.chi2); atol=2e-8, rtol=2e-8)
    end

    @testset "Correlated parameter constraints match analytic posterior" begin
        x = collect(range(-1.5, 1.5; length=13))
        sigma_y = fill(0.08, length(x))
        y = @. 1.25 * x - 0.2 + 0.03 * x^2
        model(x, p) = @. p[1] * x + p[2]

        X = hcat(x, ones(length(x)))
        w = 1.0 ./ abs2.(sigma_y)
        constraint_mean = [1.2, -0.1]
        constraint_cov = [0.06^2 0.45 * 0.06 * 0.09; 0.45 * 0.06 * 0.09 0.09^2]
        constraint_precision = Symmetric(constraint_cov) \ Matrix{Float64}(I, 2, 2)
        fisher = X' * (w .* X) + constraint_precision
        rhs = X' * (w .* y) + constraint_precision * constraint_mean
        params = Symmetric(fisher) \ rhs
        cov = Symmetric(fisher) \ Matrix{Float64}(I, 2, 2)
        residuals = y .- X * params
        data_chi2 = sum(abs2, residuals ./ sigma_y)
        constraint_delta = params .- constraint_mean
        constraint_chi2 = dot(constraint_delta, constraint_precision * constraint_delta)
        expected_chi2 = data_chi2 + constraint_chi2

        result = fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=sigma_y,
            parameter_constraints=(indices=[1, 2], mean=constraint_mean, covariance=constraint_cov),
            scale_covariance=:never,
        )

        @test result.backend == :optimization
        @test isapprox(result.params, params; atol=2e-8, rtol=2e-8)
        @test isapprox(result.param_covariance, cov; atol=2e-7, rtol=2e-7)
        @test isapprox(result.stats.chi2, expected_chi2; atol=2e-8, rtol=2e-8)
        @test result.stats.ndf == length(x)
    end
end
