"""
    AbstractFitSolver

Optional scalar-minimization adapter. Implement [`solver_capabilities`](@ref)
and [`solve_fit`](@ref); optionally specialize [`default_fit_tolerance`](@ref).
The statistical model stays in ScientificFitting.
"""
abstract type AbstractFitSolver end

"""
    OptimizationSolver(algorithm; kwargs...)

Use an Optimization.jl algorithm with its native solve keywords. Load the
corresponding Optimization solver package first. Set `maxiters` and `tol` on
the fit, not here. Other keywords (for example `g_tol` or `store_trace`) are
preserved in profile refits. Unsupported bounds/constraints are rejected using
SciML's algorithm traits. This does not turn a scalar loss into least squares.
For NLopt, pass an algorithm enum (for example `NLopt.LN_NELDERMEAD`), not a
dimension-bound `NLopt.Opt`: profile refits change the number of free parameters.
"""
struct OptimizationSolver{A, K} <: AbstractFitSolver
    algorithm::A
    kwargs::K
end

function OptimizationSolver(algorithm; kwargs...)
    algorithm isa OptimizationNLopt.NLopt.Opt && throw(ArgumentError(
        "pass an NLopt algorithm enum, not an Opt object; profile refits change the free dimension",
    ))
    any(k -> k in (:maxiters, :abstol, :reltol), keys(kwargs)) && throw(ArgumentError(
        "set maxiters and tol on the fit; reserve OptimizationSolver keywords for native options",
    ))
    return OptimizationSolver(algorithm, (; kwargs...))
end

"""
    NativeMinuitSolver(; steps=nothing, kwargs...)

Select MIGRAD after `import NativeMinuit` (Julia 1.11+); importing only the module
avoids a name conflict with its exported `profile`. `steps` is a positive
scalar or a vector in the original complete parameter order; it specifies
initial numerical step sizes, not statistical errors or priors. Native Minuit
constructor options such as `strategy=2` and `check_gradient=false` are passed
through. Bounds, fixed parameters, derivatives, and `errordef=1` belong to the
fit and cannot be overridden here, including through native `fix_*`, `limit_*`
or `error_*` aliases. Set numerical initial steps with `steps` instead.

`maxiters` is the native function-call budget, not an iteration count; `tol`
is MIGRAD's EDM tolerance parameter, defaulting to the native `0.1` (target
EDM `0.002 * tol` on our `errordef=1` cost scale). An explicit `tol` is never
rescaled or relaxed. One MIGRAD pass is run per ScientificFitting
multistart candidate. The native object is retained in `result.solver_result.raw`;
its coordinates follow `result.solver_result.parameter_indices`. Symmetric
errors in the ScientificFitting result still follow its covariance policy.
"""
struct NativeMinuitSolver{S, K} <: AbstractFitSolver
    steps::S
    kwargs::K

    function NativeMinuitSolver(steps, kwargs::NamedTuple)
        reserved = (:error, :errors, :name, :names, :limits, :fixed, :up, :errordef,
                    :grad, :tol, :maxfcn)
        # Native parameter aliases take precedence over the complete controls.
        any(k -> k in reserved || occursin(r"^(error|fix|limit)_", String(k)), keys(kwargs)) &&
            throw(ArgumentError(
                "parameter controls, error scale, derivatives, and stopping limits belong to the fit; " *
                "native error_*, fix_* and limit_* aliases are not allowed (use steps for numerical step sizes)",
            ))
        if steps !== nothing
            steps isa Real || steps isa AbstractVector{<:Real} || throw(ArgumentError(
                "steps must be a positive scalar or a vector in full parameter order",
            ))
            all(v -> isfinite(v) && v > 0, steps isa Real ? (steps,) : steps) ||
                throw(ArgumentError("initial parameter steps must be finite and positive"))
        end
        stored = steps isa AbstractVector ? collect(Float64, steps) : steps
        return new{typeof(stored), typeof(kwargs)}(stored, kwargs)
    end
end

NativeMinuitSolver(; steps=nothing, kwargs...) = NativeMinuitSolver(steps, (; kwargs...))

"""
    default_fit_tolerance(solver, derivatives::Symbol) -> Real

Stopping tolerance used when a fit omits `tol` or passes `tol=nothing`.
`derivatives` is `:auto` or `:finite`; `solver=nothing` selects the legacy path.
LsqFit/Optimization use `1e-10`, or `1e-6` for noisier finite differences.
NativeMinuit uses its native EDM tolerance `0.1` in either derivative mode.

Third-party solvers may specialize this method to match their stopping rule.
Return a finite positive value. The fit stores the resolved tolerance in
`result.options.tol` and reuses it in multistart and profile refits. Explicit
tolerances bypass this method. These are numerical criteria, not statistical
error guarantees; equal values need not mean equal accuracy across solvers.
"""
default_fit_tolerance(solver, derivatives::Symbol) =
    _validate_derivatives(derivatives) == :finite ? 1e-6 : 1e-10
function default_fit_tolerance(::NativeMinuitSolver, derivatives::Symbol)
    _validate_derivatives(derivatives)
    return 0.1
end

"""
    solver_capabilities(solver) -> NamedTuple

Declare `bounds`, `constraints`, `gradient`, and `hessian` Boolean flags.
The last two indicate required derivatives. Missing adapter implementations
raise `MethodError`, rather than silently assuming capabilities.
"""
function solver_capabilities end

function solver_capabilities(solver::OptimizationSolver)
    alg = solver.algorithm
    api = Optimization.SciMLBase
    return (; bounds=api.allowsbounds(alg), constraints=api.allowsconstraints(alg),
            gradient=api.requiresgradient(alg), hessian=api.requireshessian(alg))
end
solver_capabilities(::NativeMinuitSolver) =
    (; bounds=true, constraints=false, gradient=true, hessian=false)

"""
    FitSolverResult(params; backend, converged, message, iterations=missing,
                    parameter_indices=eachindex(params), raw=nothing)

Solver output in free coordinates, separate from statistical inference.
`parameter_indices` maps these coordinates to the full scientific parameter
vector. `raw` preserves native status/covariance/trace information; mutating it
does not update the already constructed ScientificFitting result.
"""
struct FitSolverResult{R}
    params::Vector{Float64}
    backend::Symbol
    converged::Bool
    iterations::Union{Int, Missing}
    message::String
    parameter_indices::Vector{Int}
    raw::R
end

function FitSolverResult(params; backend::Symbol, converged::Bool, message,
                         iterations=missing, parameter_indices=eachindex(params), raw=nothing)
    return FitSolverResult(collect(Float64, params), backend, converged, iterations,
                          string(message), collect(Int, parameter_indices), raw)
end

"""
    solve_fit(solver, problem::OptimizationProblem; maxiters, tol,
              parameter_indices, parameter_count, parameter_names) -> FitSolverResult

Public solver-extension boundary. `problem.f(q, problem.p)` is the complete
validated cost; `q` contains only free parameters. Bounds and nonlinear
constraints are already mapped to that order. Do not add priors, rescale the
objective, modify inputs, or compute ScientificFitting's statistical summaries.
`parameter_count` is the original full dimension (needed for parameter-specific
solver settings). Retain unavailable iteration counts as `missing` and failed
termination as `converged=false`. The core checks capabilities before dispatch.
`parameter_names` are unique labels in free-coordinate order, suitable for
named native solver operations. Unnamed problems use `p1`, `p2`, etc.
`tol` is already a finite positive number, resolved by
[`default_fit_tolerance`](@ref) unless the user supplied it explicitly.
"""
function solve_fit end

function solve_fit(solver::NativeMinuitSolver, problem; kwargs...)
    throw(ArgumentError("NativeMinuitSolver requires loading NativeMinuit with `using NativeMinuit`"))
end

function solve_fit(solver::OptimizationSolver, problem; maxiters, tol,
                   parameter_indices, parameter_count, parameter_names)
    # NLopt exposes an evaluation count, but OptimizationStats.iterations is zero.
    nlopt = solver.algorithm isa OptimizationNLopt.NLopt.Algorithm
    sol = if nlopt
        alg = OptimizationNLopt.NLopt.Opt(solver.algorithm, length(problem.u0))
        # Use parameter stopping: equal costs need not mean a contracted simplex.
        alg.xtol_rel = tol
        alg.xtol_abs = tol
        solve(problem, alg; maxiters, solver.kwargs...)
    else
        solve(problem, solver.algorithm; maxiters, abstol=tol, reltol=tol, solver.kwargs...)
    end
    iterations = !nlopt && hasproperty(sol, :stats) && hasproperty(sol.stats, :iterations) ?
                 Int(sol.stats.iterations) : missing
    return FitSolverResult(sol.u; backend=:optimization,
        converged=Optimization.SciMLBase.successful_retcode(sol.retcode),
        iterations, message=sol.retcode, parameter_indices, raw=sol)
end

"""Resolve legacy shortcuts without changing the default least-squares path."""
function _scalar_solver(problem, options)
    options.solver !== nothing && return options.solver
    if options.optimizer == :nelder_mead
        return OptimizationSolver(OptimizationNLopt.NLopt.LN_NELDERMEAD)
    end
    alg = has_constraints(problem.constraints) ? OptimizationOptimJL.IPNewton() :
          OptimizationOptimJL.LBFGS()
    return OptimizationSolver(alg)
end

"""Build one scalar solver problem for both statistical problem families."""
function _minimize_scalar(problem, options, objective, cache)
    solver = _scalar_solver(problem, options)
    caps = solver_capabilities(solver)
    bounds = _free_bounds(problem)
    constraints = _free_constraints(problem.constraints, problem)
    has_cons = has_constraints(constraints)
    bounds !== nothing && !caps.bounds && throw(ArgumentError("solver does not support parameter bounds"))
    has_cons && !caps.constraints && throw(ArgumentError("solver does not support nonlinear constraints; they cannot be dropped"))
    free = _free_indices(problem)
    q0 = _free_p0(problem)
    context = problem isa LikelihoodFitProblem ? "initial likelihood cost" : "initial cost"
    isfinite(objective(q0, cache)) || throw(ArgumentError(
        "$context must be finite; choose a starting point inside the model's support",
    ))
    _check_initial_derivatives(problem, q -> objective(q, cache), q0, caps)
    lb, ub = bounds === nothing ? (nothing, nothing) : bounds
    cons!, lcons, ucons = has_cons ? _build_constraint_system(constraints, problem) :
                         (nothing, nothing, nothing)
    optf = if caps.gradient || caps.hessian
        OptimizationFunction(objective, _optimization_ad(problem; second_order=caps.hessian); cons=cons!)
    else
        OptimizationFunction(objective; cons=cons!)
    end
    optprob = OptimizationProblem(optf, q0, cache; lb, ub, lcons, ucons)
    names = hasproperty(problem, :parameter_names) ? problem.parameter_names : nothing
    free_names = names === nothing ? ["p$i" for i in free] : names[free]
    answer = solve_fit(solver, optprob; maxiters=options.maxiters, tol=options.tol,
                       parameter_indices=free, parameter_count=length(problem.p0),
                       parameter_names=free_names)
    answer isa FitSolverResult || throw(ArgumentError("solve_fit must return a FitSolverResult"))
    length(answer.params) == length(free) && answer.parameter_indices == free ||
        throw(ArgumentError("solver returned inconsistent free parameter coordinates"))
    all(isfinite, answer.params) || throw(ArgumentError("solver returned non-finite parameters"))
    return answer
end
