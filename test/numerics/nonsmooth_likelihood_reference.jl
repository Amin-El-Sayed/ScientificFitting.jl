using ScientificFitting
using Distributions
using LinearAlgebra
using Statistics
using Test

@testset "Derivative-free likelihoods and explicit uncertainty policy" begin
    @testset "Laplace errors recover the median without invented curvature" begin
        y = [-1.2, -0.1, 0.2, 0.4, 0.8, 1.3, 5.0]
        location(x, p) = fill(p[1], length(x))
        logprob(y, mu, p) = -abs.(y .- mu) .- log(2.0)
        result = fit_likelihood_model(location, collect(eachindex(y)), y;
            logprob, p0=[0.1], optimizer=:nelder_mead, tol=1e-10)
        @test result.converged
        @test result.params[1] ≈ median(y) atol=1e-7
        @test result.options.optimizer == :nelder_mead
        @test result.options.parameter_covariance == :none
        @test all(isnan, result.param_stderr)
        @test all(isnan, result.param_covariance)
        @test ismissing(result.iterations) # NLopt reports evaluations, not iterations.
        @test :covariance_not_computed in getproperty.(diagnose(result).findings, :code)
        @test !(:invalid_local_covariance in getproperty.(diagnose(result).findings, :code))
        @test occursin("NaN", report_text(result))
        grid = [-0.2, 0.0, 0.4, 0.7, 1.0]
        scan = profile(result, 1; values=grid, on_failure=:throw)
        reference = [2*sum(abs, y .- value) - 2*sum(abs, y .- median(y)) for value in grid]
        @test scan.delta_cost ≈ reference atol=1e-7
    end

    @testset "Finite support remains impossible outside its boundary" begin
        y = [0.3, 0.5, 0.9, 1.2, 2.1]
        logprob(y, mu, p) = [value >= center ? -(value-center) : -Inf for (value, center) in zip(y, mu)]
        result = fit_likelihood_model((x, p) -> fill(p[1], length(x)), collect(eachindex(y)), y;
            logprob, p0=[0.0], bounds=([-1.0], [1.0]), optimizer=:nelder_mead, tol=1e-10)
        @test result.converged
        @test result.params[1] <= minimum(y)
        @test result.params[1] ≈ minimum(y) atol=1e-7
        @test isnan(result.param_stderr[1])
        scan = profile(result, 1; values=[0.0, 0.2, 0.3, 0.4])
        @test scan.delta_cost[1:3] ≈ 2*length(y).*(0.3 .- [0.0, 0.2, 0.3]) atol=1e-6
        @test isinf(scan.delta_cost[4])
        @test_throws ArgumentError fit_likelihood_model((x, p) -> fill(p[1], length(x)), [1, 2], [0.3, 0.5];
            logprob, p0=[0.8], optimizer=:nelder_mead)
    end

    @testset "Bounds, fixed parameters, priors and nuisance refits" begin
        calls = Vector{Vector{Float64}}()
        cost = function (p)
            # A dual input would fail this assertion, including after minimization.
            @assert p isa Vector{Float64}
            push!(calls, copy(p))
            return 2*abs(p[1]-0.7) + (p[2]-p[1])^2 + p[3]^2
        end
        result = fit_custom(cost; p0=[0.2, 0.3, 0.0], nobs=10,
            bounds=([0., 0., 0.], [2., 2., 2.]), fixed_parameters=[FixedParameter(3, 0.4)],
            parameter_priors=[ParameterPrior(2, 0.7, 1.0)], optimizer=:nelder_mead, tol=1e-10)
        @test result.converged
        @test result.params ≈ [0.7, 0.7, 0.4] atol=2e-5
        @test isnan(result.param_stderr[1]) && result.param_stderr[3] == 0
        @test all(iszero, result.param_covariance[3, :])
        scan = profile(result, 1; values=[0.5, 0.7, 0.9], on_failure=:throw)
        @test scan.delta_cost ≈ [0.42, 0.0, 0.42] atol=2e-6
        @test length(calls) > 5
        @test all(p -> 0 <= p[1] <= 2 && 0 <= p[2] <= 2 && p[3] == 0.4, calls)
        stopped = fit_custom(p -> sum(abs2, p .- [1., 2.]); p0=[0., 0.],
            nobs=10, optimizer=:nelder_mead, maxiters=1)
        @test !stopped.converged
    end

    @testset "Covariance is independently selectable for smooth objectives" begin
        cost(p) = (p[1]-1.5)^2/0.04
        with_errors = fit_custom(cost; p0=[0.2], nobs=10, optimizer=:nelder_mead,
            parameter_covariance=:hessian, tol=1e-10)
        @test with_errors.converged
        @test with_errors.param_stderr ≈ [0.2] atol=1e-9
        no_errors = fit_custom(cost; p0=[0.2], nobs=10, parameter_covariance=:none)
        @test no_errors.converged && isnan(no_errors.param_stderr[1])
        @test no_errors.options.optimizer == :lbfgs
        @test_throws ArgumentError fit_custom(cost; p0=[0.2], nobs=10, optimizer=:unknown)
        @test_throws ArgumentError fit_custom(cost; p0=[0.2], nobs=10, parameter_covariance=:unknown)
        @test_throws ArgumentError fit_custom(cost; p0=[0.2], nobs=10,
            constraints=(ineq=p -> [p[1]-1],), optimizer=:nelder_mead)
    end

    @testset "Simplex contraction, not equal costs or a zero-valued optimum" begin
        for center in (0.0, 1.5), tolerance in (1e-6, 1e-8)
            # Symmetric vertices 1.4 and 1.6 have equal costs, but neither is 1.5.
            result = fit_custom(p -> (p[1]-center)^2/0.04; p0=[0.2], nobs=10,
                optimizer=:nelder_mead, tol=tolerance, maxiters=200)
            @test result.converged
            @test result.params[1] ≈ center atol=2*tolerance
        end
    end
end
