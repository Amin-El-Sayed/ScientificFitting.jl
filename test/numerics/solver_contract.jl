using ScientificFitting
using LinearAlgebra
using Test

# A third-party adapter only sees the standard objective and free coordinates.
struct RecordingFitSolver{S} <: AbstractFitSolver
    inner::S
    dimensions::Vector{Int}
    tolerances::Vector{Float64}
end

# Model an adapter that reports success too early, independently of any package.
struct PrematureFitSolver <: AbstractFitSolver
    recovery::Symbol
    budgets::Vector{Int}
    tolerances::Vector{Float64}
end
ScientificFitting.solver_capabilities(::PrematureFitSolver) =
    (; bounds=true, constraints=false, gradient=true, hessian=false)
ScientificFitting._remaining_solver_budget(::PrematureFitSolver, answer, budget) =
    max(0, budget - answer.raw.calls)
function ScientificFitting.solve_fit(s::PrematureFitSolver, problem; maxiters, tol,
                                      parameter_indices, kwargs...)
    push!(s.budgets, maxiters)
    push!(s.tolerances, tol)
    length(s.budgets) > 1 && s.recovery == :error && error("test restart failure")
    q = s.recovery == :success && length(s.budgets) > 1 ? fill(2., length(problem.u0)) : copy(problem.u0)
    return FitSolverResult(q; backend=:test, converged=true, message="reported success",
                           parameter_indices, raw=(calls=3,))
end

@testset "Independent stationarity and one budgeted restart" begin
    for gaussian in (false, true), recover in (false, true)
        solver = PrematureFitSolver(recover ? :success : :stalled, Int[], Float64[])
        r = if gaussian
            fit_model((x,p) -> fill(p[1], length(x)), [1.,2.,3.], [2.,2.,2.];
                p0=[0.], sigma_y=ones(3), solver, maxiters=20, tol=1e-8)
        else
            fit_custom(p -> (p[1]-2)^2; p0=[0.], nobs=10, solver,
                maxiters=20, tol=1e-8)
        end
        @test r.converged == recover
        @test solver.budgets == [20, 17]
        @test solver.tolerances == [1e-8, 1e-8]
        @test r.options.maxiters == 20
        @test r.solver_result.converged # Keep the original backend claim inspectable.
        code = recover ? :optimizer_restarted : :not_stationary
        @test any(f -> f.code == code, r.diagnostics.findings)
        if recover
            @test r.params == [2.]
        else
            @test occursin("Fresh-curvature EDM", r.message)
        end
    end

    solver = PrematureFitSolver(:success, Int[], Float64[])
    r = fit_custom(p -> (p[1]-2)^2; p0=[0.], nobs=10, solver, maxiters=3)
    @test !r.converged
    @test solver.budgets == [3] # Never manufacture a fresh full budget.

    # The gradient need not vanish at a constrained optimum.
    solver = PrematureFitSolver(:stalled, Int[], Float64[])
    bound = fit_custom(p -> (p[1]-2)^2; p0=[1.], nobs=10,
        bounds=([0.], [1.]), solver)
    @test bound.converged
    @test solver.budgets == [1000]
    @test any(f -> f.code == :active_bounds, bound.diagnostics.findings)

    # Rescaling physical coordinates must not hide a displaced solution.
    for unit in (1e-6, 1., 1e6)
        solver = PrematureFitSolver(:stalled, Int[], Float64[])
        r = fit_custom(p -> (p[1]/unit-2)^2; p0=[0.], nobs=10, solver)
        @test !r.converged
        @test any(f -> f.code == :not_stationary, r.diagnostics.findings)
    end

    solver = PrematureFitSolver(:error, Int[], Float64[])
    stopped = fit_custom(p -> (p[1]-2)^2; p0=[0.], nobs=10, solver)
    @test !stopped.converged
    @test stopped.params == [0.]
    @test any(f -> f.code == :optimizer_restart_failed &&
              occursin("test restart failure", f.evidence), stopped.diagnostics.findings)

    # NativeMinuit's acceptance limit is 10 times its nominal EDM goal.
    options = FitOptions(solver=NativeMinuitSolver(), tol=1e-3)
    for (edm, rejected) in ((2.2e-6, false), (2.2e-5, true))
        finding = ScientificFitting._stationarity_finding(stopped.problem, options,
            [2+sqrt(edm)], ones(1,1), edm, q -> (q[1]-2)^2)
        @test (finding !== nothing) == rejected
    end

    solver = PrematureFitSolver(:stalled, Int[], Float64[])
    base = fit_custom(p -> p[1]^2 + (p[2]-p[1]-2)^2; p0=[0.,2.], nobs=10, solver)
    @test base.converged
    scan = profile(base, 1; values=[-1., 0., 1.])
    @test isinf(scan.delta_cost[1]) && isinf(scan.delta_cost[3])
    @test scan.delta_cost[2] == 0
    @test_throws ErrorException profile(base, 1; values=[-1.,0.,1.], on_failure=:throw)
end
ScientificFitting.solver_capabilities(s::RecordingFitSolver) = solver_capabilities(s.inner)
ScientificFitting.default_fit_tolerance(::RecordingFitSolver, ::Symbol) = 3e-9
function ScientificFitting.solve_fit(s::RecordingFitSolver, problem; kwargs...)
    @test problem isa ScientificFitting.Optimization.OptimizationProblem
    push!(s.dimensions, length(problem.u0))
    push!(s.tolerances, kwargs[:tol])
    return solve_fit(s.inner, problem; kwargs...)
end

struct InvalidFitSolver <: AbstractFitSolver
    bad::Symbol
end
ScientificFitting.solver_capabilities(::InvalidFitSolver) =
    (; bounds=false, constraints=false, gradient=false, hessian=false)
ScientificFitting.default_fit_tolerance(s::InvalidFitSolver, ::Symbol) =
    s.bad == :tolerance ? 0.0 : 1e-10
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

    solver = RecordingFitSolver(OptimizationSolver(optim.BFGS(); store_trace=true), Int[], Float64[])
    r = fit_custom(p -> (p[1]-2)^2 + (p[2]-p[1])^2 + p[3]^2;
        p0=[0., 0., 0.], nobs=10, fixed_parameters=[FixedParameter(3, 0.5)],
        solver, initial_guesses=[[1., 1., 0.]], multistart=2)
    @test r.params ≈ [2., 2., 0.5] atol=1e-7
    @test r.options.tol == 3e-9
    scan = profile(r, 1; values=[1., 2., 3.], on_failure=:throw)
    @test scan.delta_cost ≈ [1., 0., 1.] atol=1e-7
    @test solver.dimensions == [2, 2, 1, 1, 1]
    @test solver.tolerances == fill(3e-9, 5)
    refit = ScientificFitting._refit_with_fixed(r, [FixedParameter(1, 1.)])
    @test refit.options.solver === solver
    @test refit.options.tol == r.options.tol
    @test refit.solver_result.parameter_indices == [2]
    explicit = fit_custom(objective; p0=[0., 0.], nobs=10, solver, tol=1e-8)
    @test explicit.options.tol == last(solver.tolerances) == 1e-8

    # Both Gaussian and likelihood routes honor the same adapter on refits.
    x, y = [0., 1., 2., 3.], [1.1, 2.9, 5.2, 6.8]
    line(x, p) = @. p[1]*x + p[2]
    gaussian = fit_model(line, x, y; p0=[1., 0.], sigma_y=fill(0.2, 4), solver)
    ref = fit_model(line, x, y; p0=[1., 0.], sigma_y=fill(0.2, 4))
    @test gaussian.params ≈ ref.params atol=1e-7
    @test gaussian.options.tol == 3e-9
    @test ref.options.tol == 1e-10
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
    for bad in (:nonfinite, :coordinates, :tolerance)
        @test_throws ArgumentError fit_custom(objective; p0=[0., 0.], nobs=10,
                                               solver=InvalidFitSolver(bad))
    end
    @test_throws ArgumentError fit_custom(objective; p0=[0., 0.], nobs=10,
        bounds=([-1., -1.], [1., 1.]), solver=InvalidFitSolver(:nonfinite))
    @test_throws ArgumentError fit_custom(objective; p0=[0., 0.], nobs=10,
        constraints=(eq=p -> [p[1]-p[2]],), solver=OptimizationSolver(nlopt.LN_NELDERMEAD))
end

@testset "Invalid initial derivatives never enter the solver" begin
    optim = ScientificFitting.OptimizationOptimJL
    solver = RecordingFitSolver(OptimizationSolver(optim.BFGS()), Int[], Float64[])
    # Values are an ordinary quadratic, but sqrt(0) gives an undefined AD term.
    bad_ad(p) = (p[1]-2)^2 + sqrt(zero(p[1]))
    err = try
        fit_custom(bad_ad; p0=[0.], nobs=10, solver)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("gradient contains NaN or Inf", sprint(showerror, err))
    @test occursin("derivatives=:finite", sprint(showerror, err))
    @test isempty(solver.dimensions)
    finite = fit_custom(bad_ad; p0=[0.], nobs=10, solver, derivatives=:finite)
    @test finite.converged
    @test finite.params ≈ [2.] atol=1e-6
    @test finite.param_covariance ≈ ones(1, 1) rtol=1e-5

    # A finite first derivative is not sufficient for a Hessian-based solver.
    second = RecordingFitSolver(OptimizationSolver(optim.IPNewton()), Int[], Float64[])
    err = try
        fit_custom(p -> p[1]^1.5; p0=[0.], nobs=10, solver=second)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("Hessian contains NaN or Inf", sprint(showerror, err))
    @test isempty(second.dimensions)

    # Derivative-free optimization must not acquire an implicit AD requirement.
    direct = fit_custom(bad_ad; p0=[0.], nobs=10, optimizer=:nelder_mead)
    @test direct.converged
    @test direct.params ≈ [2.] atol=1e-6
end

@testset "Likelihood helpers forward solver configuration" begin
    solver = RecordingFitSolver(OptimizationSolver(ScientificFitting.OptimizationOptimJL.BFGS()),
                                Int[], Float64[])
    opts = (; p0=[0.1], solver, tol=nothing, parameter_covariance=:none)
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
        @test r.options.tol == 3e-9
        @test r.options.parameter_covariance == :none
        @test all(isnan, r.param_stderr)
    end
    @test solver.tolerances == fill(3e-9, length(results))
end
