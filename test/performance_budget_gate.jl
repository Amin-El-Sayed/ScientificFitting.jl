using ScientificFitting
using LinearAlgebra
using Statistics
using Test

const PERFORMANCE_BUDGET_SCALE = parse(Float64, get(ENV, "SCIENTIFICFITTING_PERFORMANCE_BUDGET_SCALE", "1.0"))

function _median_elapsed(f; samples::Int=5)
    f() # compile and warm caches before measuring steady-state behavior
    times = Float64[]
    for _ in 1:samples
        GC.gc()
        push!(times, @elapsed f())
    end
    return median(times)
end

function _linear_problem(n::Int)
    x = collect(range(0.0, 10.0; length=n))
    model(x, p) = @. p[1] * x + p[2]
    jacobian(x, p) = hcat(x, ones(length(x)))
    sigma_y = fill(0.2, n)
    y = model(x, [2.0, 1.0]) .+ sigma_y .* sin.(1.7 .* x)
    return model, jacobian, x, y, sigma_y
end

function _linear_model!(out, x, p)
    @. out = p[1] * x + p[2]
    return nothing
end

function _linear_jacobian!(J, x, p)
    J[:, 1] .= x
    J[:, 2] .= 1
    return nothing
end

function _dense_covariance_problem(n::Int)
    model, _, x, _, _ = _linear_problem(n)
    sigma = @. 0.1 + 0.01 * abs(sin(x))
    cov = [sigma[i] * sigma[j] * 0.25^abs(i - j) for i in eachindex(x), j in eachindex(x)]
    y = model(x, [2.0, 1.0]) .+ 0.02 .* sin.(x)
    return model, x, y, cov
end

function _structured_covariance_problem(n::Int)
    model, jacobian, x, _, _ = _linear_problem(n)
    sigma = 0.12
    rho = 0.72
    innovation_sigma = sigma * sqrt(1 - rho^2)
    function whiten!(out, residual)
        out[1] = residual[1] / sigma
        @inbounds for i in 2:length(residual)
            out[i] = (residual[i] - rho * residual[i - 1]) / innovation_sigma
        end
        return nothing
    end
    logdet_covariance = 2n * log(sigma) + (n - 1) * log1p(-rho^2)
    whitening = WhiteningOperator(whiten!; logdet_covariance)
    y = model(x, [2.0, 1.0]) .+ 0.03 .* sin.(0.9 .* x)
    return model, jacobian, x, y, whitening
end

@testset "Performance budget gate" begin
    @test Base.get_extension(ScientificFitting, :ScientificFittingCairoMakieExt) === nothing

    model, jacobian, x, y, sigma_y = _linear_problem(10_000)
    fast_result = fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=sigma_y, jacobian=jacobian)
    @test fast_result.backend == :lsqfit

    fast_time = _median_elapsed() do
        fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=sigma_y, jacobian=jacobian)
    end
    @test fast_time < 0.25 * PERFORMANCE_BUDGET_SCALE

    inplace_result = fit_model(
        _linear_model!,
        x,
        y;
        p0=[1.0, 0.0],
        sigma_y=sigma_y,
        jacobian=_linear_jacobian!,
        inplace=true,
    )
    @test inplace_result.backend == :lsqfit

    inplace_time = _median_elapsed() do
        fit_model(
            _linear_model!,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=sigma_y,
            jacobian=_linear_jacobian!,
            inplace=true,
        )
    end
    @test inplace_time < 0.25 * PERFORMANCE_BUDGET_SCALE

    fast_bytes = @allocated fit_model(
        model,
        x,
        y;
        p0=[1.0, 0.0],
        sigma_y=sigma_y,
        jacobian=jacobian,
    )
    inplace_bytes = @allocated fit_model(
        _linear_model!,
        x,
        y;
        p0=[1.0, 0.0],
        sigma_y=sigma_y,
        jacobian=_linear_jacobian!,
        inplace=true,
    )
    @test inplace_bytes < fast_bytes

    noop_result = fit_model(
        model,
        x,
        y;
        p0=[1.0, 0.0],
        sigma_y=sigma_y,
        jacobian=jacobian,
        bounds=([-Inf, -Inf], [Inf, Inf]),
    )
    @test noop_result.backend == :lsqfit

    noop_time = _median_elapsed() do
        fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=sigma_y,
            jacobian=jacobian,
            bounds=([-Inf, -Inf], [Inf, Inf]),
        )
    end
    @test noop_time < 0.25 * PERFORMANCE_BUDGET_SCALE

    dense_model, dense_x, dense_y, dense_cov = _dense_covariance_problem(300)
    dense_result = fit_model(
        dense_model,
        dense_x,
        dense_y;
        p0=[1.0, 0.0],
        cov_y=dense_cov,
        scale_covariance=:never,
    )
    @test dense_result.backend == :lsqfit

    dense_time = _median_elapsed() do
        fit_model(
            dense_model,
            dense_x,
            dense_y;
            p0=[1.0, 0.0],
            cov_y=dense_cov,
            scale_covariance=:never,
        )
    end
    @test dense_time < 0.75 * PERFORMANCE_BUDGET_SCALE

    structured_model, structured_jacobian, structured_x, structured_y, whitening =
        _structured_covariance_problem(50_000)
    structured_result = fit_model(
        structured_model,
        structured_x,
        structured_y;
        p0=[1.0, 0.0],
        whitening,
        jacobian=structured_jacobian,
        scale_covariance=:never,
    )
    @test structured_result.backend == :lsqfit
    @test structured_result.problem.cov_y === nothing

    structured_time = _median_elapsed() do
        fit_model(
            structured_model,
            structured_x,
            structured_y;
            p0=[1.0, 0.0],
            whitening,
            jacobian=structured_jacobian,
            scale_covariance=:never,
        )
    end
    @test structured_time < 0.75 * PERFORMANCE_BUDGET_SCALE

    small_model, small_jacobian, small_x, small_y, small_whitening =
        _structured_covariance_problem(10_000)
    small_bytes = @allocated fit_model(
        small_model,
        small_x,
        small_y;
        p0=[1.0, 0.0],
        whitening=small_whitening,
        jacobian=small_jacobian,
        scale_covariance=:never,
    )
    large_bytes = @allocated fit_model(
        structured_model,
        structured_x,
        structured_y;
        p0=[1.0, 0.0],
        whitening,
        jacobian=structured_jacobian,
        scale_covariance=:never,
    )
    @test large_bytes < 8 * small_bytes
end
