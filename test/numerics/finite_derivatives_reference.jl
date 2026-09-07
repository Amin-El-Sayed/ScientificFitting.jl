using Test
using ScientificFitting
using LinearAlgebra
using SparseArrays

# Float64-only callbacks enforce the same contract as NumPy, without requiring
# Python for the Julia core suite. The Python suite exercises the actual bridge.
finite_linear(x, p::Vector{Float64}) = @. p[1] * x + p[2]

@testset "Foreign-model differentiation" begin
    x = collect(range(-1.0, 1.0; length=15))
    y = finite_linear(x, [1.7, 0.4]) .+ 0.03 .* sin.(1:15)
    sigma = fill(0.12, length(x))
    design = hcat(x, ones(length(x)))
    expected = design \ y
    covariance = inv(design' * design) * 0.12^2

    @testset "Reusable typed callback boundary" begin
        first = ScientificFitting._TypedCallback{Vector{Float64}}(finite_linear)
        second = ScientificFitting._TypedCallback{Vector{Float64}}((x, p) -> p[1] .* x .+ p[2])
        @test typeof(first) == typeof(second)
        @test (@inferred first(x, [1.7, 0.4])) == second(x, [1.7, 0.4])
        for callback in (first, second)
            result = fit_model(callback, x, y; p0=[1.2, 0.2], sigma_y=sigma, derivatives=:finite)
            @test result.converged
            @test result.params ≈ expected atol=2e-6
            @test result.param_covariance ≈ covariance rtol=2e-6
        end
        wrong_type = ScientificFitting._TypedCallback{Vector{Float64}}((x, p) -> 1.0)
        @test_throws TypeError wrong_type(x, [1.0, 0.0])
        failing = ScientificFitting._TypedCallback{Float64}(p -> throw(ArgumentError("callback failed")))
        @test_throws ArgumentError failing([1.0])
    end

    @testset "Bounds, covariance, predictions, and refits" begin
        result = fit_model(finite_linear, x, y; p0=[1.2, 0.2], sigma_y=sigma,
            bounds=([0.0, -1.0], [3.0, 1.0]), derivatives=:finite)
        @test result.converged
        @test result.options.tol == 1e-6
        @test (@inferred ScientificFitting._optimization_ad(result.problem)) isa ScientificFitting.AutoFiniteDiff
        @test result.params ≈ expected atol=2e-6
        @test result.param_covariance ≈ covariance rtol=2e-6
        prediction = predict(result, [-0.5, 0.0, 0.5]; uncertainty=true)
        @test prediction.mean ≈ finite_linear([-0.5, 0.0, 0.5], expected) atol=2e-6
        @test prediction.sigma.^2 ≈ diag(hcat([-0.5, 0.0, 0.5], ones(3)) * covariance *
            hcat([-0.5, 0.0, 0.5], ones(3))') rtol=2e-6
        @test predict(result) ≈ result.model_y
        values = result.params[1] .+ result.param_stderr[1] .* [-1.0, 0.0, 1.0]
        prof = profile(result, 1; values=values)
        @test prof.delta_cost ≈ [1.0, 0.0, 1.0] atol=2e-5
        @test ScientificFitting._refit_with_fixed(result, [FixedParameter(1, expected[1])]).problem.derivatives == :finite
        @test_throws ArgumentError predict(result, [NaN])
        @test_throws ArgumentError FitProblem(finite_linear, x, y; p0=[1.0, 0.0], derivatives=:invalid)
    end

    @testset "Parameter-dependent dense XY covariance" begin
        native(x, p) = @. p[1] * x + p[2]
        cov_x = 0.02^2 .* [0.4^abs(i-j) for i in eachindex(x), j in eachindex(x)]
        options = (; p0=[1.2, 0.2], sigma_y=sigma, cov_x=cov_x,
            bounds=([0.1, -1.0], [3.0, 1.0]))
        reference = fit_model(native, x, y; options...)
        @test (@inferred ScientificFitting._optimization_ad(reference.problem)) isa ScientificFitting.AutoForwardDiff
        result = fit_model(finite_linear, x, y; options..., derivatives=:finite)
        @test result.converged
        @test result.params ≈ reference.params atol=2e-5
        @test result.param_covariance ≈ reference.param_covariance rtol=4e-3
        @test result.stats.cost_min ≈ reference.stats.cost_min atol=1e-7
    end

    @testset "Nonlinear parameter constraints" begin
        eq(p::Vector{Float64}) = [p[1]^2 + p[2]^2 - 1.0]
        objective(p::Vector{Float64}) = sum(abs2, p .- [2.0, 0.0])
        result = fit_custom(objective; p0=[0.8, 0.3], nobs=10,
            constraints=(eq=eq,), derivatives=:finite, tol=1e-10)
        @test result.converged
        @test result.params ≈ [1.0, 0.0] atol=3e-5
        @test abs(only(eq(result.params))) < 1e-6
        @test result.param_covariance ≈ Matrix{Float64}(I, 2, 2) rtol=1e-5
    end

    @testset "Sparse covariance with finite derivatives" begin
        n = length(x)
        cy = spdiagm(-1 => fill(0.002, n-1), 0 => sigma.^2, 1 => fill(0.002, n-1))
        weighted_design = cy \ design
        expected_cov = inv(design' * weighted_design)
        expected_params = expected_cov * weighted_design' * y
        result = fit_model(finite_linear, x, y; p0=[1.2, 0.2], cov_y=cy,
            bounds=([0.0, -1.0], [3.0, 1.0]), derivatives=:finite)
        @test result.converged
        @test result.problem.cov_y isa SparseMatrixCSC
        @test result.params ≈ expected_params atol=2e-6
        @test result.param_covariance ≈ expected_cov rtol=2e-5
        scan = profile(result, 1; values=result.params[1] .+ result.param_stderr[1] .* [-1, 0, 1])
        @test scan.delta_cost ≈ [1, 0, 1] atol=2e-5

        # Parameter-dependent sparse x covariance uses the full NLL, too.
        cx = 0.02^2 * spdiagm(0 => ones(n), 1 => fill(0.2, n-1), -1 => fill(0.2, n-1))
        native(x, p) = @. p[1] * x + p[2]
        dense = fit_model(native, x, y; p0=[1.2, 0.2], cov_y=Matrix(cy), cov_x=Matrix(cx))
        dynamic = fit_model(finite_linear, x, y; p0=[1.2, 0.2], cov_y=cy, cov_x=cx,
            derivatives=:finite)
        @test dynamic.params ≈ dense.params atol=2e-5
        @test dynamic.param_covariance ≈ dense.param_covariance rtol=4e-3
        @test dynamic.stats.cost_min ≈ dense.stats.cost_min atol=1e-7

        # Warm validation first; a dense copy here would allocate 128 MiB.
        component = ErrorComponent(:readout, :y, :covariance, spdiagm(0 => ones(4096)))
        ScientificFitting._normalize_error_components(component, 4096)
        bytes = @allocated ScientificFitting._normalize_error_components(component, 4096)
        @test bytes < 1_000_000
        @test_throws ArgumentError fit_model(native, x, y; p0=[1.0, 0.0], cov_y=cy,
            bounds=([0.0, -1.0], [3.0, 1.0])) # CHOLMOD is still not an AD backend.
    end

    @testset "Poisson Hessian and likelihood profiles" begin
        counts = [9, 13, 10, 12, 8, 11]
        constant(x, p::Vector{Float64}) = fill(p[1], length(x))
        result = fit_poisson_model(constant, collect(1:6), counts; p0=[9.0],
            bounds=([0.1], [30.0]), derivatives=:finite)
        @test result.params[1] ≈ sum(counts)/6 atol=1e-5
        @test result.param_covariance[1, 1] ≈ sum(counts)/36 rtol=2e-5
        prof = profile(result, 1; values=[9.0, 10.5, 12.0])
        @test all(isfinite, prof.delta_cost)
        @test prof.delta_cost[2] ≈ 0.0 atol=1e-7
    end

    @testset "In-place models and analytic overrides" begin
        model!(out::Vector{Float64}, x, p::Vector{Float64}) = (out .= finite_linear(x, p))
        result = fit_model(model!, x, y; p0=[1.0, 0.0], sigma_y=sigma,
            inplace=true, derivatives=:finite, bounds=([0.0, -1.0], [3.0, 1.0]))
        @test result.params ≈ expected atol=2e-6
        calls = Ref(0)
        analytic(x, p) = (calls[] += 1; hcat(x, ones(length(x))))
        result = fit_model(finite_linear, x, y; p0=[1.0, 0.0], sigma_y=sigma,
            jacobian=analytic, derivatives=:finite)
        before = calls[]
        predict(result; uncertainty=true)
        @test calls[] == before + 1
        explicit = fit_model(finite_linear, x, y; p0=[1.0, 0.0], sigma_y=sigma,
            derivatives=:finite, tol=1e-9)
        @test explicit.options.tol == 1e-9
    end
end
