using JuFitter
using LinearAlgebra
using Test

function _ar1_whitening(n::Int, sigma::Float64, rho::Float64)
    innovation_sigma = sigma * sqrt(1 - rho^2)
    function whiten!(out, residual)
        out[1] = residual[1] / sigma
        @inbounds for i in 2:length(residual)
            out[i] = (residual[i] - rho * residual[i - 1]) / innovation_sigma
        end
        return nothing
    end
    logdet_covariance = 2n * log(sigma) + (n - 1) * log1p(-rho^2)
    return WhiteningOperator(whiten!; logdet_covariance, marginal_sigma=sigma)
end

@testset "Structured whitening reference" begin
    n = 60
    sigma = 0.20
    rho = 0.65
    x = collect(range(-1.0, 2.0; length=n))
    y = @. 1.4 * x - 0.3 + 0.04 * sin(2.3 * x)
    covariance = [sigma^2 * rho^abs(i - j) for i in 1:n, j in 1:n]
    whitening = _ar1_whitening(n, sigma, rho)

    model(x, p) = @. p[1] * x + p[2]
    jacobian(x, p) = hcat(x, ones(length(x)))
    function model!(out, x, p)
        @. out = p[1] * x + p[2]
        return nothing
    end
    function jacobian!(J, x, p)
        J[:, 1] .= x
        J[:, 2] .= 1
        return nothing
    end

    dense = fit_model(
        model,
        x,
        y;
        p0=[1.0, 0.0],
        cov_y=covariance,
        jacobian,
        scale_covariance=:never,
    )
    structured = fit_model(
        model,
        x,
        y;
        p0=[1.0, 0.0],
        whitening,
        jacobian,
        scale_covariance=:never,
    )
    mutating = fit_model(
        model!,
        x,
        y;
        p0=[1.0, 0.0],
        whitening,
        jacobian=jacobian!,
        inplace=true,
        scale_covariance=:never,
    )

    bounded_dense = fit_model(
        model,
        x,
        y;
        p0=[1.0, 0.0],
        cov_y=covariance,
        bounds=([-3.0, -3.0], [3.0, 3.0]),
        cost=:gaussian_nll,
        scale_covariance=:never,
    )
    bounded_structured = fit_model(
        model,
        x,
        y;
        p0=[1.0, 0.0],
        whitening,
        bounds=([-3.0, -3.0], [3.0, 3.0]),
        cost=:gaussian_nll,
        scale_covariance=:never,
    )
    scan = profile(structured, 1; npoints=7, nsigma=1.5, adaptive=false)

    @test structured.backend == :lsqfit
    @test mutating.backend == :lsqfit
    @test bounded_structured.backend == :optimization
    @test bounded_structured.converged
    @test structured.problem.whitening === whitening
    @test structured.problem.cov_y === nothing
    @test isapprox(structured.params, dense.params; rtol=1e-10, atol=1e-11)
    @test isapprox(mutating.params, dense.params; rtol=1e-10, atol=1e-11)
    @test isapprox(structured.param_covariance, dense.param_covariance; rtol=1e-9, atol=1e-11)
    @test isapprox(mutating.param_covariance, dense.param_covariance; rtol=1e-9, atol=1e-11)
    @test isapprox(structured.weighted_residuals, dense.weighted_residuals; rtol=1e-10, atol=1e-11)
    @test isapprox(structured.stats.chi2, dense.stats.chi2; rtol=1e-10, atol=1e-11)
    @test isapprox(structured.stats.nll_min, dense.stats.nll_min; rtol=1e-10, atol=1e-10)
    @test isapprox(structured.stats.aic, dense.stats.aic; rtol=1e-10, atol=1e-10)
    @test isapprox(structured.stats.bic, dense.stats.bic; rtol=1e-10, atol=1e-10)
    @test isapprox(bounded_structured.params, bounded_dense.params; rtol=1e-9, atol=1e-10)
    @test isapprox(bounded_structured.param_covariance, bounded_dense.param_covariance; rtol=1e-8, atol=1e-10)
    @test isapprox(bounded_structured.stats.nll_min, bounded_dense.stats.nll_min; rtol=1e-10, atol=1e-10)
    @test scan isa ProfileResult
    @test all(isfinite, scan.delta_cost)

    @test_throws ArgumentError WhiteningOperator((out, residual) -> nothing; logdet_covariance=Inf)
    @test_throws ArgumentError WhiteningOperator(
        (out, residual) -> copyto!(out, residual);
        logdet_covariance=0.0,
        marginal_sigma=0.0,
    )
    @test_throws ArgumentError WhiteningOperator(
        (out, residual) -> copyto!(out, residual);
        logdet_covariance=0.0,
        marginal_sigma=[1.0, Inf],
    )
    @test_throws ArgumentError fit_model(model, x, y; p0=[1.0, 0.0], whitening=:invalid)
    @test_throws ArgumentError fit_model(
        model,
        x,
        y;
        p0=[1.0, 0.0],
        whitening=WhiteningOperator(residual -> residual; logdet_covariance=0.0),
    )
    @test_throws ArgumentError fit_model(
        model,
        x,
        y;
        p0=[1.0, 0.0],
        whitening=WhiteningOperator((out, residual) -> (out[1] = residual[1]); logdet_covariance=0.0),
    )
    @test_throws ArgumentError fit_model(
        model,
        x,
        y;
        p0=[1.0, 0.0],
        whitening=WhiteningOperator((out, residual) -> fill!(out, Inf); logdet_covariance=0.0),
    )
    @test_throws ArgumentError fit_model(
        model,
        x,
        y;
        p0=[1.0, 0.0],
        whitening=WhiteningOperator((out, residual) -> fill!(out, 0); logdet_covariance=0.0),
    )
    @test_throws ArgumentError whitening(zeros(n - 1), ones(n))
    vector_only!(out::Vector{Float64}, residual::Vector{Float64}) = copyto!(out, residual)
    @test_throws ArgumentError fit_model(
        model,
        x,
        y;
        p0=[1.0, 0.0],
        whitening=WhiteningOperator(vector_only!; logdet_covariance=0.0),
    )
    @test_throws ArgumentError fit_model(
        model,
        x,
        y;
        p0=[1.0, 0.0],
        whitening=WhiteningOperator(
            (out, residual) -> copyto!(out, residual);
            logdet_covariance=0.0,
            marginal_sigma=ones(n - 1),
        ),
    )
    @test_throws ArgumentError fit_model(model, x, y; p0=[1.0, 0.0], whitening, sigma_y=fill(1.0, n))
    @test_throws ArgumentError fit_model(model, x, y; p0=[1.0, 0.0], whitening, cov_y=I(n))
    @test_throws ArgumentError fit_model(model, x, y; p0=[1.0, 0.0], whitening, sigma_x=fill(0.1, n))
    @test_throws ArgumentError fit_model(
        model,
        x,
        y;
        p0=[1.0, 0.0],
        whitening,
        error_components=(name=:extra, target=:y, mode=:absolute, values=0.1),
    )
end
