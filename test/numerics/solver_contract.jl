using ScientificFitting
using LinearAlgebra
using Test

# A third-party adapter only sees the standard objective and free coordinates.
struct RecordingFitSolver{S} <: AbstractFitSolver
    inner::S
    dimensions::Vector{Int}
end
ScientificFitting.solver_capabilities(s::RecordingFitSolver) = solver_capabilities(s.inner)
function ScientificFitting.solve_fit(s::RecordingFitSolver, problem; kwargs...)
    @test problem isa ScientificFitting.Optimization.OptimizationProblem
    push!(s.dimensions, length(problem.u0))
    return solve_fit(s.inner, problem; kwargs...)
end

struct InvalidFitSolver <: AbstractFitSolver
    bad::Symbol
end
ScientificFitting.solver_capabilities(::InvalidFitSolver) =
    (; bounds=false, constraints=false, gradient=false, hessian=false)
function ScientificFitting.solve_fit(s::InvalidFitSolver, problem; parameter_indices, kwargs...)
    p = s.bad == :nonfinite ? fill(NaN, length(problem.u0)) : copy(problem.u0)
    indices = s.bad == :coordinates ? reverse(parameter_indices) : parameter_indices
    return FitSolverResult(p; backend=:test, converged=false, message="test",
                           parameter_indices=indices)
end

@testset "Open scalar solver contract" begin
    optim = ScientificFitting.OptimizationOptimJL
    nlopt = ScientificFitting.OptimizationNLopt.NLopt
    objective(p) = (p[1]-2)^2 + (p[2]+1)^2 / 4
    for alg in (optim.BFGS(), optim.LBFGS(), nlopt.LN_NELDERMEAD)
        solver = OptimizationSolver(alg)
        r = fit_custom(objective; p0=[0., 0.], nobs=10, solver,
                       parameter_covariance=:hessian, maxiters=500, tol=1e-9)
        @test r.converged
        @test r.params ≈ [2., -1.] atol=1e-5
        @test r.param_covariance ≈ [1. 0.; 0. 4.] atol=1e-8
        @test r.options.solver === solver
        @test r.solver_result.parameter_indices == [1, 2]
        @test r.solver_result.raw !== nothing
        @test r.message == string(r.solver_result.raw.retcode)
    end

    solver = RecordingFitSolver(OptimizationSolver(optim.BFGS(); store_trace=true), Int[])
    r = fit_custom(p -> (p[1]-2)^2 + (p[2]-p[1])^2 + p[3]^2;
        p0=[0., 0., 0.], nobs=10, fixed_parameters=[FixedParameter(3, 0.5)],
        solver, initial_guesses=[[1., 1., 0.]], multistart=2)
    @test r.params ≈ [2., 2., 0.5] atol=1e-7
    scan = profile(r, 1; values=[1., 2., 3.], on_failure=:throw)
    @test scan.delta_cost ≈ [1., 0., 1.] atol=1e-7
    @test solver.dimensions == [2, 2, 1, 1, 1]
    refit = ScientificFitting._refit_with_fixed(r, [FixedParameter(1, 1.)])
    @test refit.options.solver === solver
    @test refit.solver_result.parameter_indices == [2]

    # Both Gaussian and likelihood routes honor the same adapter on refits.
    x, y = [0., 1., 2., 3.], [1.1, 2.9, 5.2, 6.8]
    line(x, p) = @. p[1]*x + p[2]
    gaussian = fit_model(line, x, y; p0=[1., 0.], sigma_y=fill(0.2, 4), solver)
    ref = fit_model(line, x, y; p0=[1., 0.], sigma_y=fill(0.2, 4))
    @test gaussian.params ≈ ref.params atol=1e-7
    @test gaussian.param_covariance ≈ ref.param_covariance rtol=1e-6
    @test gaussian.stats.cost_min ≈ ref.stats.cost_min rtol=1e-8
    fixed = ScientificFitting._refit_with_fixed(gaussian, [FixedParameter(1, 2.)])
    @test fixed.options.solver === solver
    @test fixed.solver_result.parameter_indices == [2]

    constrained = fit_custom(p -> sum(abs2, p .- 2); p0=[0.2, 0.3], nobs=10,
        constraints=(ineq=p -> [p[1]+p[2]-2],), solver=OptimizationSolver(optim.IPNewton()))
    @test constrained.converged
    @test constrained.params ≈ [1., 1.] atol=1e-5
    @test sum(constrained.params) <= 2 + 1e-8
    @test solver_capabilities(OptimizationSolver(optim.IPNewton())).hessian

    @test_throws ArgumentError fit_custom(objective; p0=[0., 0.], nobs=10, solver,
                                          optimizer=:nelder_mead)
    @test_throws ArgumentError fit_model(line, x, y; p0=[0., 0.], solver, backend=:lsqfit)
    @test_throws ArgumentError OptimizationSolver(optim.BFGS(); maxiters=3)
    @test_throws ArgumentError OptimizationSolver(nlopt.Opt(:LN_NELDERMEAD, 2))
    if Base.get_extension(ScientificFitting, :ScientificFittingNativeMinuitExt) === nothing
        @test_throws ArgumentError fit_custom(objective; p0=[0., 0.], nobs=10,
                                               solver=NativeMinuitSolver())
    end
    for bad in (:nonfinite, :coordinates)
        @test_throws ArgumentError fit_custom(objective; p0=[0., 0.], nobs=10,
                                               solver=InvalidFitSolver(bad))
    end
    @test_throws ArgumentError fit_custom(objective; p0=[0., 0.], nobs=10,
        bounds=([-1., -1.], [1., 1.]), solver=InvalidFitSolver(:nonfinite))
    @test_throws ArgumentError fit_custom(objective; p0=[0., 0.], nobs=10,
        constraints=(eq=p -> [p[1]-p[2]],), solver=OptimizationSolver(nlopt.LN_NELDERMEAD))
end

@testset "Likelihood helpers forward solver configuration" begin
    solver = OptimizationSolver(ScientificFitting.OptimizationOptimJL.BFGS())
    opts = (; p0=[0.1], solver, parameter_covariance=:none)
    constant(x, p) = fill(exp(p[1]), length(x))
    normal(x, p) = exp(-(x-p[1])^2/2) / sqrt(2pi)
    results = [
        fit_poisson_model(constant, [0., 1., 2.], [2, 3, 4]; opts...),
        fit_histogram_model((edges, p) -> fill(exp(p[1]), length(edges)-1),
                            [0., 1., 2., 3.], [2, 3, 4]; opts...),
        fit_histogram_density(normal, [-4., 0., 4.], [4, 6]; opts...),
        fit_unbinned_model(normal, [-0.2, 0.1, 0.4]; opts...),
        fit_extended_unbinned_model((x, p) -> exp(p[1]), [0.2, 0.4, 0.7],
                                    (0., 1.); opts...),
        fit_likelihood_model(constant, [0., 1., 2.], [2., 3., 4.];
                             logprob=(y, mu, p) -> -(y.-mu).^2/2, opts...),
        fit_indexed_model(constant, [:a, :b, :c], [2., 3., 4.]; opts...),
        fit_multi_model([constant, constant], [[0., 1.], [0., 1.]],
                        [[2., 3.], [3., 4.]]; opts...),
    ]
    for r in results
        @test r.converged
        @test r.options.solver === solver
        @test r.options.parameter_covariance == :none
        @test all(isnan, r.param_stderr)
    end
end
