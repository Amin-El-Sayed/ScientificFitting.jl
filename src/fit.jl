function _fit_with_lsqfit(problem::FitProblem, options::FitOptions)
    free_p0 = _free_p0(problem)
    weighted_model = (x, q) -> problem.model(x, _expand_free_parameters(problem, q))
    weighted_y = problem.y
    weighted_jacobian = problem.jacobian === nothing ? nothing : (x, q) -> _free_jacobian_from_full(problem, problem.jacobian(x, _expand_free_parameters(problem, q)))

    if problem.sigma_y !== nothing
        sigma = problem.sigma_y
        weighted_model = (x, q) -> problem.model(x, _expand_free_parameters(problem, q)) ./ sigma
        weighted_y = problem.y ./ sigma
        if problem.jacobian !== nothing
            weighted_jacobian = (x, q) -> _free_jacobian_from_full(problem, problem.jacobian(x, _expand_free_parameters(problem, q))) ./ reshape(sigma, :, 1)
        end
    elseif problem.cov_y !== nothing
        F = _stable_cholesky(problem.cov_y)
        weighted_model = (x, q) -> _whiten_with_factor(
            F,
            problem.model(x, _expand_free_parameters(problem, q)),
        )
        weighted_y = _whiten_with_factor(F, problem.y)
        if problem.jacobian !== nothing
            weighted_jacobian = (x, q) -> _whiten_with_factor(
                F,
                _free_jacobian_from_full(
                    problem,
                    problem.jacobian(x, _expand_free_parameters(problem, q)),
                ),
            )
        end
    end

    fit_result = if weighted_jacobian === nothing
        LsqFit.curve_fit(weighted_model, problem.x, weighted_y, free_p0)
    else
        LsqFit.curve_fit(weighted_model, weighted_jacobian, problem.x, weighted_y, free_p0)
    end
    params = _expand_free_parameters(problem, LsqFit.coef(fit_result))

    converged = true
    iterations = hasproperty(fit_result, :iterations) ? getproperty(fit_result, :iterations) : options.maxiters
    message = "Converged with LsqFit"

    # LsqFit stores the weighted model Jacobian. JuFitter stores residual
    # Jacobians, hence the sign flip at result construction.
    return params, converged, iterations, message, Matrix{Float64}(fit_result.jacobian)
end

function _fit_with_optimization(problem::FitProblem, options::FitOptions)
    if _static_effective_covariance_available(problem)
        cov = _effective_covariance(problem, problem.p0)
        if cov isa SparseMatrixCSC
            throw(ArgumentError(
                "sparse covariance currently supports the unbounded least-squares backend; " *
                "use a dense covariance for constrained or Gaussian-NLL fits until " *
                "matrix-free whitening operators are added",
            ))
        end
    end

    cache = _prepare_fit_cache(problem)
    objective = (q, cache) -> _cost_value(cache, _expand_free_parameters(cache.problem, q), options.cost)

    lb = nothing
    ub = nothing
    free_bounds = _free_bounds(problem)
    if free_bounds !== nothing
        lb, ub = free_bounds
    end

    free_constraints = _free_constraints(problem.constraints, problem)
    has_cons = has_constraints(free_constraints)
    if has_cons
        cons!, lcons, ucons = _build_constraint_system(free_constraints, problem)
        ad = DifferentiationInterface.SecondOrder(Optimization.AutoForwardDiff(), Optimization.AutoForwardDiff())
        optf = OptimizationFunction(objective, ad; cons=cons!)
        optprob = OptimizationProblem(optf, _free_p0(problem), cache; lb=lb, ub=ub, lcons=lcons, ucons=ucons)
        sol = solve(
            optprob,
            OptimizationOptimJL.IPNewton();
            maxiters=options.maxiters,
            abstol=options.tol,
            reltol=options.tol,
        )
    else
        optf = OptimizationFunction(objective, Optimization.AutoForwardDiff())
        optprob = OptimizationProblem(optf, _free_p0(problem), cache; lb=lb, ub=ub)
        sol = solve(
            optprob,
            OptimizationOptimJL.LBFGS();
            maxiters=options.maxiters,
            abstol=options.tol,
            reltol=options.tol,
        )
    end

    params = _expand_free_parameters(problem, sol.u)
    retcode_text = string(sol.retcode)
    converged = occursin("Success", retcode_text) || occursin("Default", retcode_text)
    iterations = hasproperty(sol, :stats) && hasproperty(sol.stats, :iterations) ? sol.stats.iterations : options.maxiters
    message = string(sol.retcode)

    return params, converged, iterations, message, nothing
end

function _build_fit_result(
    problem::FitProblem,
    options::FitOptions,
    backend::Symbol,
    params::Vector{Float64},
    converged::Bool,
    iterations::Int,
    message::String,
    backend_jacobian=nothing,
)
    cache = _prepare_fit_cache(problem)
    yhat = _model_values(problem, params)
    residuals = problem.y .- yhat
    weighted_residuals = _weighted_residual(cache, params)
    chi2 = _chi2_cost(cache, params)
    cost_min = Float64(_cost_value(cache, params, options.cost))
    nll_min = Float64(_gaussian_nll(cache, params))

    nconstraint_obs = sum((length(c.indices) for c in problem.parameter_constraints); init=0)
    nobs = length(problem.y) + length(problem.parameter_priors) + nconstraint_obs
    npar = length(_free_indices(problem))
    ndf = nobs - npar
    chi2_ndf = ndf > 0 ? chi2 / ndf : NaN
    pvalue = ndf > 0 ? ccdf(Chisq(ndf), chi2) : NaN
    aic = nll_min + 2.0 * npar
    bic = nll_min + log(nobs) * npar

    Jw = if backend_jacobian === nothing
        _weighted_jacobian(cache, params)
    else
        _full_jacobian_from_free(problem, -backend_jacobian, length(weighted_residuals))
    end
    free_idx = _free_indices(problem)
    cov = if isempty(free_idx)
        _embed_free_covariance(problem, zeros(Float64, 0, 0))
    elseif _resolve_cost(problem, options.cost) == :gaussian_nll
        _covariance_from_cost_hessian(cache, params, options.cost)
    else
        scale = _should_scale_covariance(problem, options.scale_covariance)
        free_cov = _covariance_from_weighted_jacobian(Jw[:, free_idx], chi2, ndf, scale)
        _embed_free_covariance(problem, free_cov)
    end
    stderr = sqrt.(clamp.(diag(cov), 0.0, Inf))
    corr = _correlation_from_covariance(cov)

    stats = FitStatistics(_resolve_cost(problem, options.cost), cost_min, nll_min, chi2, chi2_ndf, ndf, pvalue, aic, bic)
    hessian = if _resolve_cost(problem, options.cost) == :gaussian_nll && !isempty(free_idx)
        ForwardDiff.hessian(q -> _cost_value(cache, _expand_free_parameters(problem, q), options.cost), params[free_idx])
    else
        nothing
    end
    diagnostics = _fit_diagnostics(problem, params, cov, converged, ndf; hessian=hessian, gof=chi2)

    return FitResult(
        problem,
        options,
        backend,
        converged,
        iterations,
        message,
        params,
        stderr,
        cov,
        corr,
        yhat,
        residuals,
        weighted_residuals,
        Jw,
        stats,
        diagnostics,
    )
end

"""
    fit(problem::FitProblem; backend=:auto, cost=:auto, maxiters=500, tol=1e-10, ci_level=0.6827, scale_covariance=:auto, initial_guesses=nothing, multistart=1)

Run a fit and return a `FitResult`.

`cost` accepts:
- `:auto`: choose `:gaussian_nll` for parameter-dependent covariance, otherwise `:chi2`
- `:chi2`: weighted least-squares cost
- `:gaussian_nll`: full Gaussian negative log-likelihood, including log-det terms

`scale_covariance` controls post-fit parameter covariance scaling:
- `:auto`: scale only when no data uncertainties were provided
- `:never`: trust the supplied uncertainties
- `:always`: multiply by `chi2/ndf`

`initial_guesses` and `multistart` can be used for difficult nonlinear fits.
"""
function fit(
    problem::FitProblem;
    backend::Symbol=:auto,
    cost::Symbol=:auto,
    maxiters::Int=500,
    tol::Real=1e-10,
    ci_level::Real=0.6827,
    scale_covariance=:auto,
    initial_guesses=nothing,
    multistart::Int=1,
)
    multistart > 0 || throw(ArgumentError("multistart must be >= 1"))
    options = FitOptions(
        backend=backend,
        cost=_resolve_cost(problem, cost),
        maxiters=maxiters,
        tol=Float64(tol),
        ci_level=Float64(ci_level),
        scale_covariance=_normalize_scale_covariance(scale_covariance),
        multistart=multistart,
    )

    candidates = _initial_candidates(problem, initial_guesses, multistart)
    best_result = nothing
    best_cost = Inf
    last_error = nothing

    for candidate in candidates
        candidate_problem = _with_p0(problem, candidate)
        try
            result = if isempty(_free_indices(candidate_problem))
                params = _expand_free_parameters(candidate_problem, Float64[])
                _build_fit_result(candidate_problem, options, :fixed, params, true, 0, "All parameters fixed", nothing)
            else
                chosen_backend = _solve_backend(candidate_problem, backend, options.cost)
                params, converged, iterations, message, backend_jacobian = if chosen_backend == :lsqfit
                    _fit_with_lsqfit(candidate_problem, options)
                elseif chosen_backend == :optimization
                    _fit_with_optimization(candidate_problem, options)
                else
                    throw(ArgumentError("unsupported backend: $chosen_backend (use :auto, :lsqfit, or :optimization)"))
                end
                _build_fit_result(candidate_problem, options, chosen_backend, params, converged, iterations, message, backend_jacobian)
            end

            current_cost = result.stats.cost_min
            if result.converged && isfinite(current_cost) && current_cost < best_cost
                best_result = result
                best_cost = current_cost
            elseif best_result === nothing && isfinite(current_cost)
                best_result = result
                best_cost = current_cost
            end
        catch err
            last_error = err
        end
    end

    if best_result === nothing
        last_error === nothing || throw(last_error)
        throw(ErrorException("fit failed for all initial guesses"))
    end

    return best_result
end

"""
    fit_model(model, x, y; kwargs...)

Convenience wrapper that constructs `FitProblem` then calls `fit`.

Useful parameter-control kwargs:
- `parameter_priors=(index=i, mean=mu, sigma=sigma)`
- `parameter_priors=(index=i, mean=mu, sigma_minus=sminus, sigma_plus=splus)`
- `parameter_constraints=(indices=[i, j], mean=[mu_i, mu_j], covariance=cov)`
- `fixed_parameters=(index=i, value=value, sigma=sigma)`
- `fixed_parameters=(index=i, value=value, sigma_minus=sminus, sigma_plus=splus)`

For fits with x uncertainty, `x_derivative=(x, p) -> dy_dx` supplies a
vectorized model derivative with respect to x. This avoids the default
point-by-point AD path and is the preferred route for large datasets.
"""
function fit_model(
    model,
    x::AbstractVector,
    y::AbstractVector;
    p0::AbstractVector,
    sigma_y=nothing,
    sigma_x=nothing,
    cov_y=nothing,
    cov_x=nothing,
    error_components=nothing,
    bounds=nothing,
    constraints=nothing,
    parameter_priors=nothing,
    parameter_constraints=nothing,
    fixed_parameters=nothing,
    jacobian=nothing,
    x_derivative=nothing,
    backend::Symbol=:auto,
    cost::Symbol=:auto,
    maxiters::Int=500,
    tol::Real=1e-10,
    ci_level::Real=0.6827,
    scale_covariance=:auto,
    initial_guesses=nothing,
    multistart::Int=1,
)
    problem = FitProblem(
        model,
        x,
        y;
        p0=p0,
        sigma_y=sigma_y,
        sigma_x=sigma_x,
        cov_y=cov_y,
        cov_x=cov_x,
        error_components=error_components,
        bounds=bounds,
        constraints=constraints,
        parameter_priors=parameter_priors,
        parameter_constraints=parameter_constraints,
        fixed_parameters=fixed_parameters,
        jacobian=jacobian,
        x_derivative=x_derivative,
    )

    return fit(
        problem;
        backend=backend,
        cost=cost,
        maxiters=maxiters,
        tol=tol,
        ci_level=ci_level,
        scale_covariance=scale_covariance,
        initial_guesses=initial_guesses,
        multistart=multistart,
    )
end
