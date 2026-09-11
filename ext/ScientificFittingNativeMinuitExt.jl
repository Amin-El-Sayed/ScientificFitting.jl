module ScientificFittingNativeMinuitExt

using ScientificFitting
using NativeMinuit: Minuit, migrad!
import Optimization
import DifferentiationInterface as DI

function ScientificFitting.solve_fit(solver::NativeMinuitSolver,
                                    problem::Optimization.OptimizationProblem;
                                    maxiters, tol, parameter_indices, parameter_count, parameter_names)
    objective = q -> problem.f(q, problem.p)
    ad = problem.f.adtype
    prepared = DI.prepare_gradient(objective, ad, problem.u0)
    grad = q -> DI.gradient(objective, prepared, ad, q)
    limits = problem.lb === nothing ? nothing : [
        (isfinite(lo) ? lo : nothing, isfinite(hi) ? hi : nothing)
        for (lo, hi) in zip(problem.lb, problem.ub)
    ]
    steps = solver.steps
    if steps isa AbstractVector
        length(steps) == parameter_count || throw(ArgumentError(
            "NativeMinuit steps must have one entry per full model parameter",
        ))
        steps = steps[parameter_indices]
    end
    # ScientificFitting supplies -2logL/chi-square, never Minuit's half-scale NLL.
    native = Minuit(objective, copy(problem.u0); grad, limits, error=steps,
                    names=parameter_names, errordef=1.0,
                    solver.kwargs...)
    migrad!(native; maxfcn=maxiters, tol, iterate=1)
    fm = native.fmin.internal
    message = "MIGRAD: valid=$(native.valid), call_limit=$(fm.reached_call_limit), " *
              "above_max_edm=$(fm.above_max_edm), edm=$(native.edm)"
    return FitSolverResult(native.values; backend=:native_minuit,
                           converged=native.valid, message, parameter_indices, raw=native)
end

end
