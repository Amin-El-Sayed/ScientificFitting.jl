function _apply_lsqfit_weight!(values::AbstractVector, sigma, factor, operator, scratch)
    if sigma !== nothing
        @inbounds for i in eachindex(values, sigma)
            values[i] /= sigma[i]
        end
    elseif factor !== nothing
        copyto!(values, _whiten_with_factor(factor, values))
    elseif operator !== nothing
        _whiten_with_operator!(scratch, operator, values)
        copyto!(values, scratch)
    end
    return values
end

function _apply_lsqfit_weight!(values::AbstractMatrix, sigma, factor, operator, scratch)
    if sigma !== nothing
        @inbounds for j in axes(values, 2), i in axes(values, 1)
            values[i, j] /= sigma[i]
        end
    elseif factor !== nothing
        copyto!(values, _whiten_with_factor(factor, values))
    elseif operator !== nothing
        _whiten_with_operator!(scratch, operator, values)
        copyto!(values, scratch)
    end
    return values
end

function _fit_with_lsqfit(problem::FitProblem, options::FitOptions)
    free_p0 = _free_p0(problem)
    sigma = problem.sigma_y
    factor = problem.cov_y === nothing ? nothing : _stable_cholesky(problem.cov_y)
    operator = problem.whitening
    weighted_y = sigma !== nothing ? problem.y ./ sigma :
                 factor !== nothing ? _whiten_with_factor(factor, problem.y) :
                 operator !== nothing ? _whiten_with_operator(operator, problem.y) : problem.y
    model_whitening_scratch = operator === nothing ? nothing : similar(problem.y)

    fit_result = if problem.model isa _InPlaceModel
        weighted_model! = function (out, x, q)
            problem.model(out, x, _expand_free_parameters(problem, q))
            _apply_lsqfit_weight!(out, sigma, factor, operator, model_whitening_scratch)
            return nothing
        end

        if problem.jacobian === nothing
            LsqFit.curve_fit(weighted_model!, problem.x, weighted_y, free_p0; inplace=true)
        else
            free_idx = _free_indices(problem)
            full_jacobian = length(free_idx) == length(problem.p0) ? nothing :
                            Matrix{Float64}(undef, length(problem.x), length(problem.p0))
            jacobian_whitening_scratch = operator === nothing ? nothing :
                                         Matrix{Float64}(undef, length(problem.x), length(free_idx))
            weighted_jacobian! = function (out, x, q)
                params = _expand_free_parameters(problem, q)
                if full_jacobian === nothing
                    problem.jacobian(out, x, params)
                else
                    problem.jacobian(full_jacobian, x, params)
                    copyto!(out, view(full_jacobian, :, free_idx))
                end
                _apply_lsqfit_weight!(out, sigma, factor, operator, jacobian_whitening_scratch)
                return nothing
            end
            LsqFit.curve_fit(
                weighted_model!,
                weighted_jacobian!,
                problem.x,
                weighted_y,
                free_p0;
                inplace=true,
            )
        end
    else
        weighted_model = (x, q) -> begin
            values = problem.model(x, _expand_free_parameters(problem, q))
            sigma !== nothing && return values ./ sigma
            factor !== nothing && return _whiten_with_factor(factor, values)
            operator !== nothing && return _whiten_with_operator(operator, values)
            return values
        end
        weighted_jacobian = problem.jacobian === nothing ? nothing : (x, q) -> begin
            full = problem.jacobian(x, _expand_free_parameters(problem, q))
            values = _free_jacobian_from_full(problem, full)
            sigma !== nothing && return values ./ reshape(sigma, :, 1)
            factor !== nothing && return _whiten_with_factor(factor, values)
            operator !== nothing && return _whiten_with_operator(operator, values)
            return values
        end

        if weighted_jacobian === nothing
            LsqFit.curve_fit(weighted_model, problem.x, weighted_y, free_p0)
        else
            LsqFit.curve_fit(weighted_model, weighted_jacobian, problem.x, weighted_y, free_p0)
        end
    end
    params = _expand_free_parameters(problem, LsqFit.coef(fit_result))

    converged = true
    iterations = hasproperty(fit_result, :iterations) ?
                 Int(getproperty(fit_result, :iterations)) : missing
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
                "use a dense covariance or an AD-compatible WhiteningOperator for " *
                "constrained or Gaussian-likelihood fits",
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
    iterations = hasproperty(sol, :stats) && hasproperty(sol.stats, :iterations) ?
                 Int(sol.stats.iterations) : missing
    message = string(sol.retcode)

    return params, converged, iterations, message, nothing
end

function _build_fit_result(
    problem::FitProblem,
    options::FitOptions,
    backend::Symbol,
    params::Vector{Float64},
    converged::Bool,
    iterations::Union{Int, Missing},
    message::String,
    backend_jacobian=nothing,
)
    cache = _prepare_fit_cache(problem)
    yhat = _model_values(problem, params)
    residuals = problem.y .- yhat
    weighted_residuals = _weighted_residual(cache, params)
    chi2 = _chi2_cost(cache, params)
    cost_min = Float64(_cost_value(cache, params, options.cost))
    minus2loglik_min = Float64(_gaussian_minus2loglik(cache, params))

    nconstraint_obs = sum((length(c.indices) for c in problem.parameter_constraints); init=0)
    nobs = length(problem.y) + length(problem.parameter_priors) + nconstraint_obs
    npar = length(_free_indices(problem))
    ndf = nobs - npar
    chi2_ndf = ndf > 0 ? chi2 / ndf : NaN
    pvalue = ndf > 0 ? ccdf(Chisq(ndf), chi2) : NaN
    aic = minus2loglik_min + 2.0 * npar
    bic = minus2loglik_min + log(nobs) * npar

    Jw = if backend_jacobian === nothing
        _weighted_jacobian(cache, params)
    else
        _full_jacobian_from_free(problem, -backend_jacobian, length(weighted_residuals))
    end
    free_idx = _free_indices(problem)
    cov = if isempty(free_idx)
        _embed_free_covariance(problem, zeros(Float64, 0, 0))
    elseif _resolve_cost(problem, options.cost) == :gaussian_likelihood
        _covariance_from_cost_hessian(cache, params, options.cost)
    else
        scale = _should_scale_covariance(problem, options.scale_covariance)
        free_cov = _covariance_from_weighted_jacobian(Jw[:, free_idx], chi2, ndf, scale)
        _embed_free_covariance(problem, free_cov)
    end
    stderr = _standard_errors_from_covariance(cov)
    corr = _correlation_from_covariance(cov)

    stats = FitStatistics(
        _resolve_cost(problem, options.cost),
        cost_min,
        minus2loglik_min,
        chi2,
        chi2_ndf,
        ndf,
        pvalue,
        aic,
        bic,
    )
    hessian = if _resolve_cost(problem, options.cost) == :gaussian_likelihood && !isempty(free_idx)
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
    fit(problem::FitProblem; backend=:auto, cost=:auto, maxiters=500,
        tol=1e-10, scale_covariance=:auto, initial_guesses=nothing,
        multistart=1) -> FitResult

Fit a validated Gaussian `FitProblem`.

Keyword contracts:
- `backend`: `:auto`, `:lsqfit`, or `:optimization`. `:auto` uses LsqFit only
  when static chi-square least squares represents the complete problem.
- `cost`: `:auto`, `:chi2`, or `:gaussian_likelihood`. `:auto` selects the
  normalized Gaussian `-2 log(L)` cost for parameter-dependent covariance.
- `maxiters`, `tol`: positive solver limits used for every candidate.
- `scale_covariance`: `:auto`, `:never`, or `:always`. `:auto` estimates a
  residual scale only when no observation uncertainty was supplied.
- `initial_guesses`: additional complete parameter vectors in `p0` order.
- `multistart`: number of deterministic candidates, including `problem.p0`.

The converged finite candidate with the lowest cost is returned. If no candidate
converges but one remains finite, it is returned with `converged == false`; use
`diagnostic_dashboard(result)` before interpreting it. If every candidate
fails, the last model, validation, or solver error is rethrown.

An explicit incompatible `backend=:lsqfit` request raises `ArgumentError`
instead of dropping bounds, parameter terms, constraints, active error
components, or parameter-dependent covariance.
"""
function fit(
    problem::FitProblem;
    backend::Symbol=:auto,
    cost::Symbol=:auto,
    maxiters::Int=500,
    tol::Real=1e-10,
    scale_covariance=:auto,
    initial_guesses=nothing,
    multistart::Int=1,
)
    options = FitOptions(
        backend=backend,
        cost=_resolve_cost(problem, cost),
        maxiters=maxiters,
        tol=Float64(tol),
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
    fit_model(model, x, y; p0, kwargs...) -> FitResult

Fit scalar observations `y` measured at `x` with a model satisfying
`model(x, p) -> yhat`. `x`, `y`, and `p0` are copied to finite `Float64`
storage; `length(x)` must equal `length(y)` and the model must return one finite
prediction per observation.

Observation uncertainty is supplied by exactly the applicable combination of
`sigma_y`, `cov_y`, `sigma_x`, `cov_x`, named `error_components`, or one
complete static `WhiteningOperator`. A whitening operator is mutually exclusive
with every other observation-uncertainty keyword. With no uncertainty,
unweighted least squares is used.

For fits with x uncertainty, `x_derivative=(x, p) -> dy_dx` supplies a
vectorized model derivative with respect to x. This avoids the default
point-by-point AD path and is the preferred route for large datasets.

For large static correlated datasets, `whitening=WhiteningOperator(...)`
supplies the complete covariance through a matrix-free operation. It cannot be
combined with other observation-uncertainty keywords. The operator must accept
generic `AbstractVector` inputs and AD element types when the general optimizer
is used.

Set `inplace=true` for `model!(out, x, p)`. The unbounded least-squares backend
uses LsqFit's native in-place model interface; generic optimizer paths preserve
the same contract with a type-correct output buffer. An optional analytic
Jacobian must then use `jacobian!(J, x, p)`.

Parameter control uses `bounds`, `constraints`, `parameter_priors`,
`parameter_constraints`, and `fixed_parameters`. Solver keywords are forwarded
to `fit(::FitProblem)` with defaults `backend=:auto`, `cost=:auto`,
`maxiters=500`, `tol=1e-10`, `scale_covariance=:auto`, and `multistart=1`.

Returns a `FitResult`. Invalid dimensions, non-finite values, non-positive
standard deviations, contradictory uncertainty inputs, invalid covariance,
bounds, or parameter controls raise `ArgumentError` before a result is
constructed.

# Example

```julia
model(x, p) = @. p[1] * x + p[2]
result = fit_model(model, [0.0, 1.0, 2.0], [0.1, 1.2, 1.9];
    p0=[1.0, 0.0], sigma_y=fill(0.2, 3))
```
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
    whitening=nothing,
    error_components=nothing,
    bounds=nothing,
    constraints=nothing,
    parameter_priors=nothing,
    parameter_constraints=nothing,
    fixed_parameters=nothing,
    jacobian=nothing,
    x_derivative=nothing,
    inplace::Bool=false,
    backend::Symbol=:auto,
    cost::Symbol=:auto,
    maxiters::Int=500,
    tol::Real=1e-10,
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
        whitening=whitening,
        error_components=error_components,
        bounds=bounds,
        constraints=constraints,
        parameter_priors=parameter_priors,
        parameter_constraints=parameter_constraints,
        fixed_parameters=fixed_parameters,
        jacobian=jacobian,
        x_derivative=x_derivative,
        inplace=inplace,
    )

    return fit(
        problem;
        backend=backend,
        cost=cost,
        maxiters=maxiters,
        tol=tol,
        scale_covariance=scale_covariance,
        initial_guesses=initial_guesses,
        multistart=multistart,
    )
end
