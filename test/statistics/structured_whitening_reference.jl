using ScientificFitting
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
    factor = cholesky(Symmetric(covariance))
    design = hcat(x, ones(n))
    gls_params = (factor.L \ design) \ (factor.L \ y)

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
        cost=:gaussian_likelihood,
        scale_covariance=:never,
    )
    bounded_structured = fit_model(
        model,
        x,
        y;
        p0=[1.0, 0.0],
        whitening,
        bounds=([-3.0, -3.0], [3.0, 3.0]),
        cost=:gaussian_likelihood,
        scale_covariance=:never,
    )
    scan = profile(structured, 1; npoints=7, nsigma=1.5, adaptive=false)

    @test structured.backend == :lsqfit
    @test mutating.backend == :lsqfit
    @test bounded_structured.backend == :optimization
    @test bounded_structured.converged
    @test structured.problem.whitening === whitening
    @test structured.problem.cov_y === nothing
    for result in (dense, structured, mutating)
        @test result.converged
        # Iterative fits need solver tolerance; whitening at the same point does not.
        @test isapprox(result.params, gls_params; rtol=1e-8, atol=1e-10)
        @test isapprox(
            result.weighted_residuals,
            factor.L \ (y - design * result.params);
            rtol=1e-10,
            atol=1e-11,
        )
    end
    @test isapprox(structured.param_covariance, dense.param_covariance; rtol=1e-9, atol=1e-11)
    @test isapprox(mutating.param_covariance, dense.param_covariance; rtol=1e-9, atol=1e-11)
    @test isapprox(structured.stats.chi2, dense.stats.chi2; rtol=1e-10, atol=1e-11)
    @test isapprox(structured.stats.minus2loglik_min, dense.stats.minus2loglik_min; rtol=1e-10, atol=1e-10)
    @test isapprox(structured.stats.aic, dense.stats.aic; rtol=1e-10, atol=1e-10)
    @test isapprox(structured.stats.bic, dense.stats.bic; rtol=1e-10, atol=1e-10)
    @test isapprox(bounded_structured.params, bounded_dense.params; rtol=1e-9, atol=1e-10)
    @test isapprox(bounded_structured.param_covariance, bounded_dense.param_covariance; rtol=1e-8, atol=1e-10)
    @test isapprox(
        bounded_structured.stats.minus2loglik_min,
        bounded_dense.stats.minus2loglik_min;
        rtol=1e-10,
        atol=1e-10,
    )
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
    @test_throws DimensionMismatch fit_model(
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

@testset "Exponential whitening reference" begin
    gaps = [0.30, 0.05, 0.75, 0.80, 0.10, 0.75, 0.85, 0.02, 0.88, 1.30, 0.20]
    x = cumsum([0.0; gaps])
    sigma = 0.4
    correlation_length = 0.8
    covariance = [sigma^2 * exp(-abs(xi - xj) / correlation_length) for xi in x, xj in x]
    factor = cholesky(Symmetric(covariance))
    whitening = exponential_whitening(x, sigma, correlation_length)

    @test isapprox(whitening.logdet_covariance, logdet(factor); rtol=1e-12)
    @test whitening.marginal_sigma == sigma
    for residual in (sin.(1:12), collect(range(-1.0, 2.0; length=12)), [(-1.0)^i * i / 12 for i in 1:12])
        whitened = whitening(similar(residual), residual)
        @test isapprox(sum(abs2, whitened), dot(residual, factor \ residual); rtol=1e-12)
    end

    model(x, p) = @. p[1] * x + p[2]
    n = 40
    fit_sigma = 0.2

    fit_x = cumsum(0.05 .+ 0.1 .* (1 .+ sin.(1:n)))
    fit_y = @. 1.4 * fit_x - 0.3 + 0.05 * sin(7.1 * fit_x)
    fit_length = 0.3
    fit_covariance = [fit_sigma^2 * exp(-abs(xi - xj) / fit_length) for xi in fit_x, xj in fit_x]
    structured = fit_model(
        model,
        fit_x,
        fit_y;
        p0=[1.0, 0.0],
        whitening=exponential_whitening(fit_x, fit_sigma, fit_length),
    )
    dense = fit_model(model, fit_x, fit_y; p0=[1.0, 0.0], cov_y=fit_covariance)
    @test structured.converged
    @test dense.converged
    @test isapprox(structured.params, dense.params; rtol=1e-8)
    @test isapprox(structured.param_stderr, dense.param_stderr; rtol=1e-8)
    @test isapprox(structured.stats.chi2, dense.stats.chi2; rtol=1e-8)

    # exp(-1 / 1e-9) underflows to zero, so this covariance is exactly sigma^2 * I.
    unit_x = collect(1.0:n)
    unit_y = @. 0.7 * unit_x + 2.0 + 0.1 * cos(1.9 * unit_x)
    limit = fit_model(
        model,
        unit_x,
        unit_y;
        p0=[1.0, 0.0],
        whitening=exponential_whitening(unit_x, fit_sigma, 1e-9),
    )
    independent = fit_model(model, unit_x, unit_y; p0=[1.0, 0.0], sigma_y=fill(fit_sigma, n))
    @test isapprox(limit.params, independent.params; rtol=1e-6)
    @test isapprox(limit.param_stderr, independent.param_stderr; rtol=1e-6)
    @test isapprox(limit.stats.chi2, independent.stats.chi2; rtol=1e-6)

    @test_throws DomainError exponential_whitening(x, 0.0, correlation_length)
    @test_throws DomainError exponential_whitening(x, sigma, -1.0)
    @test_throws DomainError exponential_whitening([0.0, 1.0, 1.0], sigma, correlation_length)
    @test_throws DomainError exponential_whitening([0.0, 2.0, 1.0], sigma, correlation_length)
    @test_throws ArgumentError exponential_whitening([0.0, NaN, 1.0], sigma, correlation_length)
    @test_throws ArgumentError exponential_whitening(Float64[], sigma, correlation_length)
    @test_throws DimensionMismatch whitening(zeros(11), zeros(11))
end
