"""
    LikelihoodFitProblem(objective, gof, p0; nobs, cost_name, kwargs...)

Internal problem representation for likelihood and custom-objective fits.
`objective(p)` is minimized directly, while optional `gof(p)` supplies a
chi-square-like goodness-of-fit statistic for reduced statistics and p-values.
Most users create this through `fit_poisson_model`, `fit_histogram_model`,
`fit_unbinned_model`, `fit_extended_unbinned_model`, or `fit_custom`.
"""
struct LikelihoodFitProblem{TF, TG}
    objective::TF
    gof::TG
    p0::Vector{Float64}
    bounds::Union{Nothing, Tuple{Vector{Float64}, Vector{Float64}}}
    constraints::ConstraintSpec
    parameter_priors::Vector{ParameterPrior}
    parameter_constraints::Vector{ParameterConstraint}
    fixed_parameters::Vector{FixedParameter}
    nobs::Int
    cost_name::Symbol
    parameter_names::Union{Nothing, Vector{String}}
end

"""
    LikelihoodFitResult

Result of a likelihood or custom-objective fit. It mirrors the parameter,
covariance, statistics, and diagnostics fields of `FitResult`, but stores a
`LikelihoodFitProblem` instead of x-y residual data. Plotting support depends on
the specific likelihood workflow because not every objective has a natural
curve representation.
"""
struct LikelihoodFitResult
    problem::LikelihoodFitProblem
    options::FitOptions
    backend::Symbol
    converged::Bool
    iterations::Int
    message::String
    params::Vector{Float64}
    param_stderr::Vector{Float64}
    param_covariance::Matrix{Float64}
    param_correlation::Matrix{Float64}
    stats::FitStatistics
    diagnostics::FitDiagnostics
end

function LikelihoodFitProblem(
    objective,
    gof,
    p0::AbstractVector;
    bounds=nothing,
    constraints=nothing,
    parameter_priors=nothing,
    parameter_constraints=nothing,
    fixed_parameters=nothing,
    nobs::Integer,
    cost_name::Symbol,
    parameter_names=nothing,
)
    p0_vec = _float_vector(p0)
    length(p0_vec) > 0 || throw(ArgumentError("p0 must contain at least one parameter"))
    nobs > 0 || throw(ArgumentError("nobs must be positive"))
    names = parameter_names === nothing ? nothing : collect(String, parameter_names)
    names === nothing || length(names) == length(p0_vec) || throw(ArgumentError("parameter_names length must match p0"))

    return LikelihoodFitProblem(
        objective,
        gof,
        p0_vec,
        _normalize_bounds(bounds, length(p0_vec)),
        _normalize_constraints(constraints),
        _normalize_parameter_priors(parameter_priors, length(p0_vec)),
        _normalize_parameter_constraints(parameter_constraints, length(p0_vec)),
        _normalize_fixed_parameters(fixed_parameters, length(p0_vec)),
        Int(nobs),
        cost_name,
        names,
    )
end

function _with_p0(problem::LikelihoodFitProblem, p0::AbstractVector)
    return LikelihoodFitProblem(
        problem.objective,
        problem.gof,
        p0;
        bounds=problem.bounds,
        constraints=problem.constraints,
        parameter_priors=problem.parameter_priors,
        parameter_constraints=problem.parameter_constraints,
        fixed_parameters=problem.fixed_parameters,
        nobs=problem.nobs,
        cost_name=problem.cost_name,
        parameter_names=problem.parameter_names,
    )
end

function _likelihood_cost(problem::LikelihoodFitProblem, p::AbstractVector)
    return problem.objective(p) + _prior_nll(problem, p) + _parameter_constraint_nll(problem, p)
end

struct LikelihoodEvaluationCache{TP, TC}
    problem::TP
    parameter_constraints::TC
end

function _prepare_likelihood_cache(problem::LikelihoodFitProblem)
    return LikelihoodEvaluationCache(problem, _prepare_parameter_constraints(problem))
end

function _likelihood_cost(cache::LikelihoodEvaluationCache, p::AbstractVector)
    problem = cache.problem
    return problem.objective(p) + _prior_nll(problem, p) + _parameter_constraint_nll(cache.parameter_constraints, p)
end

function _likelihood_gof(problem::LikelihoodFitProblem, p::AbstractVector)
    return problem.gof === nothing ? NaN : Float64(problem.gof(p))
end

function _fit_likelihood_problem(problem::LikelihoodFitProblem, options::FitOptions)
    cache = _prepare_likelihood_cache(problem)
    objective = (q, cache) -> _likelihood_cost(cache, _expand_free_parameters(cache.problem, q))
    free_constraints = _free_constraints(problem.constraints, problem)
    free_bounds = _free_bounds(problem)
    lb = free_bounds === nothing ? nothing : free_bounds[1]
    ub = free_bounds === nothing ? nothing : free_bounds[2]

    if has_constraints(free_constraints)
        cons!, lcons, ucons = _build_constraint_system(free_constraints, problem)
        ad = DifferentiationInterface.SecondOrder(Optimization.AutoForwardDiff(), Optimization.AutoForwardDiff())
        optf = OptimizationFunction(objective, ad; cons=cons!)
        optprob = OptimizationProblem(optf, _free_p0(problem), cache; lb=lb, ub=ub, lcons=lcons, ucons=ucons)
        sol = solve(optprob, OptimizationOptimJL.IPNewton(); maxiters=options.maxiters, abstol=options.tol, reltol=options.tol)
    else
        optf = OptimizationFunction(objective, Optimization.AutoForwardDiff())
        optprob = OptimizationProblem(optf, _free_p0(problem), cache; lb=lb, ub=ub)
        sol = solve(optprob, OptimizationOptimJL.LBFGS(); maxiters=options.maxiters, abstol=options.tol, reltol=options.tol)
    end

    params = _expand_free_parameters(problem, sol.u)
    retcode_text = string(sol.retcode)
    converged = occursin("Success", retcode_text) || occursin("Default", retcode_text)
    iterations = hasproperty(sol, :stats) && hasproperty(sol.stats, :iterations) ? sol.stats.iterations : options.maxiters
    return params, converged, iterations, retcode_text
end

function _likelihood_covariance(problem::LikelihoodFitProblem, params::Vector{Float64})
    return _likelihood_covariance(_prepare_likelihood_cache(problem), params)
end

function _likelihood_covariance(cache::LikelihoodEvaluationCache, params::Vector{Float64})
    problem = cache.problem
    free_idx = _free_indices(problem)
    if isempty(free_idx)
        return _embed_free_covariance(problem, zeros(Float64, 0, 0))
    end

    q = params[free_idx]
    H = ForwardDiff.hessian(qq -> _likelihood_cost(cache, _expand_free_parameters(problem, qq)), q)
    free_cov = 2.0 .* _stable_symmetric_inverse(H)
    return _embed_free_covariance(problem, free_cov)
end

function _build_likelihood_result(
    problem::LikelihoodFitProblem,
    options::FitOptions,
    params::Vector{Float64},
    converged::Bool,
    iterations::Int,
    message::String,
)
    cache = _prepare_likelihood_cache(problem)
    cost_min = Float64(_likelihood_cost(cache, params))
    gof = _likelihood_gof(problem, params)
    nconstraint_obs = sum((length(c.indices) for c in problem.parameter_constraints); init=0)
    nobs = problem.nobs + length(problem.parameter_priors) + nconstraint_obs
    npar = length(_free_indices(problem))
    ndf = nobs - npar
    gof_ndf = isfinite(gof) && ndf > 0 ? gof / ndf : NaN
    pvalue = isfinite(gof) && ndf > 0 ? ccdf(Chisq(ndf), gof) : NaN
    aic = cost_min + 2.0 * npar
    bic = cost_min + log(nobs) * npar
    cov = _likelihood_covariance(cache, params)
    stderr = sqrt.(clamp.(diag(cov), 0.0, Inf))
    corr = _correlation_from_covariance(cov)
    stats = FitStatistics(problem.cost_name, cost_min, cost_min, gof, gof_ndf, ndf, pvalue, aic, bic)
    free_idx = _free_indices(problem)
    hessian = isempty(free_idx) ? nothing : ForwardDiff.hessian(q -> _likelihood_cost(cache, _expand_free_parameters(problem, q)), params[free_idx])
    diagnostics = _fit_diagnostics(problem, params, cov, converged, ndf; hessian=hessian, gof=gof)

    return LikelihoodFitResult(problem, options, :optimization, converged, iterations, message, params, stderr, cov, corr, stats, diagnostics)
end

function fit(
    problem::LikelihoodFitProblem;
    maxiters::Int=1000,
    tol::Real=1e-10,
    ci_level::Real=0.6827,
    initial_guesses=nothing,
    multistart::Int=1,
    kwargs...,
)
    multistart > 0 || throw(ArgumentError("multistart must be >= 1"))
    options = FitOptions(
        backend=:optimization,
        cost=problem.cost_name,
        maxiters=maxiters,
        tol=Float64(tol),
        ci_level=Float64(ci_level),
        scale_covariance=:never,
        multistart=multistart,
    )

    candidates = _initial_candidates(problem, initial_guesses, multistart)
    best_result = nothing
    best_cost = Inf
    last_error = nothing

    for candidate in candidates
        candidate_problem = _with_p0(problem, candidate)
        try
            params, converged, iterations, message = if isempty(_free_indices(candidate_problem))
                (_expand_free_parameters(candidate_problem, Float64[]), true, 0, "All parameters fixed")
            else
                _fit_likelihood_problem(candidate_problem, options)
            end
            result = _build_likelihood_result(candidate_problem, options, params, converged, iterations, message)
            if result.converged && isfinite(result.stats.cost_min) && result.stats.cost_min < best_cost
                best_result = result
                best_cost = result.stats.cost_min
            elseif best_result === nothing && isfinite(result.stats.cost_min)
                best_result = result
                best_cost = result.stats.cost_min
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
    fit_custom(objective; p0, gof=nothing, nobs, cost_name=:custom, kwargs...)

Fit a user-defined scalar objective. `objective(p)` is minimized directly.
If `gof(p)` is supplied it is used for reduced goodness-of-fit statistics and
p-values; otherwise these fields are `NaN`.
"""
function fit_custom(
    objective;
    p0::AbstractVector,
    gof=nothing,
    nobs::Integer,
    cost_name::Symbol=:custom,
    bounds=nothing,
    constraints=nothing,
    parameter_priors=nothing,
    parameter_constraints=nothing,
    fixed_parameters=nothing,
    parameter_names=nothing,
    maxiters::Int=1000,
    tol::Real=1e-10,
    ci_level::Real=0.6827,
    initial_guesses=nothing,
    multistart::Int=1,
)
    problem = LikelihoodFitProblem(
        objective,
        gof,
        p0;
        bounds=bounds,
        constraints=constraints,
        parameter_priors=parameter_priors,
        parameter_constraints=parameter_constraints,
        fixed_parameters=fixed_parameters,
        nobs=nobs,
        cost_name=cost_name,
        parameter_names=parameter_names,
    )
    return fit(problem; maxiters=maxiters, tol=tol, ci_level=ci_level, initial_guesses=initial_guesses, multistart=multistart)
end

function _assert_finite_observations(name::AbstractString, values::AbstractVector)
    _assert_finite_vector(name, values)
    return values
end

function _assert_count_observations(name::AbstractString, values::AbstractVector)
    _assert_finite_vector(name, values)
    all(>=(0.0), values) || throw(ArgumentError("$name must be non-negative"))
    all(isinteger, values) || throw(ArgumentError("$name must contain integer-valued counts"))
    return values
end

function _positive_expectation(mu, n::Int)
    values = collect(mu)
    length(values) == n || throw(ArgumentError("model expectation length must match observations"))
    all(isfinite, values) || throw(ArgumentError("model expectation contains non-finite values"))
    all(values .> 0) || throw(ArgumentError("model expectation must be strictly positive"))
    return values
end

function _poisson_nll_terms(counts::Vector{Float64}, mu::AbstractVector)
    total = zero(eltype(mu))
    @inbounds for i in eachindex(counts)
        n = counts[i]
        total += 2.0 * (mu[i] - n * log(mu[i]) + loggamma(n + 1.0))
    end
    return total
end

function _poisson_deviance(counts::Vector{Float64}, mu::AbstractVector)
    total = zero(eltype(mu))
    @inbounds for i in eachindex(counts)
        n = counts[i]
        if n == 0
            total += 2.0 * mu[i]
        else
            total += 2.0 * (mu[i] - n + n * log(n / mu[i]))
        end
    end
    return total
end

"""
    fit_poisson_model(model, x, counts; p0, kwargs...)

Fit count data with a Poisson likelihood. `model(x, p)` must return the
strictly positive expected counts for each observation. `counts` must contain
finite non-negative integer-valued observations.
"""
function fit_poisson_model(
    model,
    x::AbstractVector,
    counts::AbstractVector;
    p0::AbstractVector,
    bounds=nothing,
    constraints=nothing,
    parameter_priors=nothing,
    parameter_constraints=nothing,
    fixed_parameters=nothing,
    parameter_names=nothing,
    maxiters::Int=1000,
    tol::Real=1e-10,
    ci_level::Real=0.6827,
    initial_guesses=nothing,
    multistart::Int=1,
)
    x_vec = _float_vector(x)
    counts_vec = _float_vector(counts)
    length(x_vec) == length(counts_vec) || throw(ArgumentError("x and counts must have equal length"))
    _assert_finite_observations("x", x_vec)
    _assert_count_observations("counts", counts_vec)

    objective = p -> _poisson_nll_terms(counts_vec, _positive_expectation(model(x_vec, p), length(counts_vec)))
    gof = p -> _poisson_deviance(counts_vec, _positive_expectation(model(x_vec, p), length(counts_vec)))
    problem = LikelihoodFitProblem(
        objective,
        gof,
        p0;
        bounds=bounds,
        constraints=constraints,
        parameter_priors=parameter_priors,
        parameter_constraints=parameter_constraints,
        fixed_parameters=fixed_parameters,
        nobs=length(counts_vec),
        cost_name=:poisson_nll,
        parameter_names=parameter_names,
    )
    return fit(problem; maxiters=maxiters, tol=tol, ci_level=ci_level, initial_guesses=initial_guesses, multistart=multistart)
end

"""
    fit_histogram_model(expected_counts, edges, counts; p0, kwargs...)

Fit binned counts. `expected_counts(edges, p)` must return one positive expected
count per bin. Histogram edges must be finite and strictly increasing; `counts`
must contain finite non-negative integer-valued observations.
"""
function fit_histogram_model(
    expected_counts,
    edges::AbstractVector,
    counts::AbstractVector;
    p0::AbstractVector,
    bounds=nothing,
    constraints=nothing,
    parameter_priors=nothing,
    parameter_constraints=nothing,
    fixed_parameters=nothing,
    parameter_names=nothing,
    maxiters::Int=1000,
    tol::Real=1e-10,
    ci_level::Real=0.6827,
    initial_guesses=nothing,
    multistart::Int=1,
)
    edges_vec = _float_vector(edges)
    counts_vec = _float_vector(counts)
    length(edges_vec) == length(counts_vec) + 1 || throw(ArgumentError("edges length must be count length + 1"))
    _assert_finite_observations("histogram edges", edges_vec)
    _assert_count_observations("counts", counts_vec)
    any(diff(edges_vec) .<= 0) && throw(ArgumentError("histogram edges must be strictly increasing"))

    objective = p -> _poisson_nll_terms(counts_vec, _positive_expectation(expected_counts(edges_vec, p), length(counts_vec)))
    gof = p -> _poisson_deviance(counts_vec, _positive_expectation(expected_counts(edges_vec, p), length(counts_vec)))
    problem = LikelihoodFitProblem(
        objective,
        gof,
        p0;
        bounds=bounds,
        constraints=constraints,
        parameter_priors=parameter_priors,
        parameter_constraints=parameter_constraints,
        fixed_parameters=fixed_parameters,
        nobs=length(counts_vec),
        cost_name=:histogram_poisson_nll,
        parameter_names=parameter_names,
    )
    return fit(problem; maxiters=maxiters, tol=tol, ci_level=ci_level, initial_guesses=initial_guesses, multistart=multistart)
end

"""
    fit_histogram_density(pdf, edges, counts; p0, total_count=sum(counts), kwargs...)

Fit binned counts from a probability density. Expected bin counts are computed
with adaptive Gauss-Kronrod quadrature. `total_count` and `rtol` must be finite
and positive.
"""
function fit_histogram_density(
    pdf,
    edges::AbstractVector,
    counts::AbstractVector;
    p0::AbstractVector,
    total_count::Real=sum(counts),
    rtol::Real=1e-8,
    bounds=nothing,
    constraints=nothing,
    parameter_priors=nothing,
    parameter_constraints=nothing,
    fixed_parameters=nothing,
    parameter_names=nothing,
    maxiters::Int=1000,
    tol::Real=1e-10,
    ci_level::Real=0.6827,
    initial_guesses=nothing,
    multistart::Int=1,
)
    total = Float64(total_count)
    isfinite(total) && total > 0 || throw(ArgumentError("total_count must be finite and > 0"))
    isfinite(rtol) && rtol > 0 || throw(ArgumentError("rtol must be finite and > 0"))

    expected_counts = function (edge_values, p)
        mu = Vector{eltype(p)}(undef, length(edge_values) - 1)
        @inbounds for i in eachindex(mu)
            integral, _ = quadgk(x -> pdf(x, p), edge_values[i], edge_values[i + 1]; rtol=rtol)
            mu[i] = total * integral
        end
        return mu
    end
    return fit_histogram_model(
        expected_counts,
        edges,
        counts;
        p0=p0,
        bounds=bounds,
        constraints=constraints,
        parameter_priors=parameter_priors,
        parameter_constraints=parameter_constraints,
        fixed_parameters=fixed_parameters,
        parameter_names=parameter_names,
        maxiters=maxiters,
        tol=tol,
        ci_level=ci_level,
        initial_guesses=initial_guesses,
        multistart=multistart,
    )
end

"""
    fit_unbinned_model(pdf, data; p0, kwargs...)

Fit independent unbinned observations with a normalized positive density
`pdf(x, p)`. Observations must be finite.
"""
function fit_unbinned_model(
    pdf,
    data::AbstractVector;
    p0::AbstractVector,
    bounds=nothing,
    constraints=nothing,
    parameter_priors=nothing,
    parameter_constraints=nothing,
    fixed_parameters=nothing,
    parameter_names=nothing,
    maxiters::Int=1000,
    tol::Real=1e-10,
    ci_level::Real=0.6827,
    initial_guesses=nothing,
    multistart::Int=1,
)
    data_vec = _float_vector(data)
    _assert_finite_observations("unbinned data", data_vec)
    objective = function (p)
        total = zero(eltype(p))
        @inbounds for x in data_vec
            density = pdf(x, p)
            density > 0 || throw(ArgumentError("pdf must be strictly positive at all data points"))
            total += -2.0 * log(density)
        end
        return total
    end
    problem = LikelihoodFitProblem(
        objective,
        nothing,
        p0;
        bounds=bounds,
        constraints=constraints,
        parameter_priors=parameter_priors,
        parameter_constraints=parameter_constraints,
        fixed_parameters=fixed_parameters,
        nobs=length(data_vec),
        cost_name=:unbinned_nll,
        parameter_names=parameter_names,
    )
    return fit(problem; maxiters=maxiters, tol=tol, ci_level=ci_level, initial_guesses=initial_guesses, multistart=multistart)
end

"""
    fit_extended_unbinned_model(rate, data, domain; p0, kwargs...)

Fit an inhomogeneous Poisson point process. `rate(x, p)` is the event intensity,
not a normalized density. `domain=(a, b)` defines the integration range for the
expected total event count. The domain endpoints must be finite, and all
observations must lie inside the domain.
"""
function fit_extended_unbinned_model(
    rate,
    data::AbstractVector,
    domain::Tuple{<:Real, <:Real};
    p0::AbstractVector,
    rtol::Real=1e-8,
    bounds=nothing,
    constraints=nothing,
    parameter_priors=nothing,
    parameter_constraints=nothing,
    fixed_parameters=nothing,
    parameter_names=nothing,
    maxiters::Int=1000,
    tol::Real=1e-10,
    ci_level::Real=0.6827,
    initial_guesses=nothing,
    multistart::Int=1,
)
    data_vec = _float_vector(data)
    a, b = Float64(domain[1]), Float64(domain[2])
    isfinite(a) && isfinite(b) || throw(ArgumentError("domain endpoints must be finite"))
    a < b || throw(ArgumentError("domain must satisfy domain[1] < domain[2]"))
    isfinite(rtol) && rtol > 0 || throw(ArgumentError("rtol must be finite and > 0"))
    _assert_finite_observations("unbinned data", data_vec)
    all(x -> a <= x <= b, data_vec) || throw(ArgumentError("extended unbinned data must lie inside the domain"))

    objective = function (p)
        expected, _ = quadgk(x -> rate(x, p), a, b; rtol=rtol)
        expected > 0 || throw(ArgumentError("integrated rate must be positive"))
        total = 2.0 * expected
        @inbounds for x in data_vec
            lambda = rate(x, p)
            lambda > 0 || throw(ArgumentError("rate must be strictly positive at all data points"))
            total -= 2.0 * log(lambda)
        end
        return total
    end

    problem = LikelihoodFitProblem(
        objective,
        nothing,
        p0;
        bounds=bounds,
        constraints=constraints,
        parameter_priors=parameter_priors,
        parameter_constraints=parameter_constraints,
        fixed_parameters=fixed_parameters,
        nobs=length(data_vec),
        cost_name=:extended_unbinned_nll,
        parameter_names=parameter_names,
    )
    return fit(problem; maxiters=maxiters, tol=tol, ci_level=ci_level, initial_guesses=initial_guesses, multistart=multistart)
end

function _gaussian_chi2_from_residual(residual::AbstractVector, sigma_y, cov_y)
    if sigma_y !== nothing && cov_y !== nothing
        throw(ArgumentError("use either sigma_y or cov_y, not both"))
    elseif sigma_y !== nothing
        sigma = _float_vector(sigma_y)
        length(sigma) == length(residual) || throw(ArgumentError("sigma_y length must match observations"))
        _assert_positive_sigma("sigma_y", sigma)
        return sum(abs2, residual ./ sigma)
    elseif cov_y !== nothing
        cov = _float_matrix(cov_y)
        size(cov) == (length(residual), length(residual)) || throw(ArgumentError("cov_y must be n x n"))
        _assert_covariance_matrix("cov_y", cov)
        z = _stable_cholesky(cov).L \ residual
        return sum(abs2, z)
    end
    return sum(abs2, residual)
end

function _normalize_indexed_uncertainty(y_vec::Vector{Float64}, sigma_y, cov_y)
    if sigma_y !== nothing && cov_y !== nothing
        throw(ArgumentError("use either sigma_y or cov_y, not both"))
    elseif sigma_y !== nothing
        sigma = _float_vector(sigma_y)
        length(sigma) == length(y_vec) || throw(ArgumentError("sigma_y length must match observations"))
        _assert_positive_sigma("sigma_y", sigma)
        return sigma, nothing
    elseif cov_y !== nothing
        cov = _float_matrix(cov_y)
        size(cov) == (length(y_vec), length(y_vec)) || throw(ArgumentError("cov_y must be n x n"))
        _assert_covariance_matrix("cov_y", cov)
        return nothing, cov
    end
    return nothing, nothing
end

function _normalize_multi_sigma_sets(sigma_y, y_sets::AbstractVector)
    ndatasets = length(y_sets)
    sigma_y === nothing && return [nothing for _ in 1:ndatasets]

    length(sigma_y) == ndatasets || throw(ArgumentError("sigma_y length must match models"))
    sigma_sets = Vector{Union{Nothing, Vector{Float64}}}(undef, ndatasets)
    for i in 1:ndatasets
        if sigma_y[i] === nothing
            sigma_sets[i] = nothing
        else
            sigma = _float_vector(sigma_y[i])
            length(sigma) == length(y_sets[i]) || throw(ArgumentError("sigma_y[$i] length mismatch"))
            _assert_positive_sigma("sigma_y[$i]", sigma)
            sigma_sets[i] = sigma
        end
    end
    return sigma_sets
end

"""
    fit_indexed_model(model, indices, y; p0, sigma_y=nothing, cov_y=nothing, kwargs...)

Fit observations addressed by arbitrary indices. `model(indices, p)` must return
one model value per index. This is useful when the independent variable is not a
numeric 1D x-axis. Optional `sigma_y` entries must be finite and positive;
optional `cov_y` must be a finite symmetric positive-definite covariance matrix.
"""
function fit_indexed_model(
    model,
    indices,
    y::AbstractVector;
    p0::AbstractVector,
    sigma_y=nothing,
    cov_y=nothing,
    bounds=nothing,
    constraints=nothing,
    parameter_priors=nothing,
    parameter_constraints=nothing,
    fixed_parameters=nothing,
    parameter_names=nothing,
    maxiters::Int=1000,
    tol::Real=1e-10,
    ci_level::Real=0.6827,
    initial_guesses=nothing,
    multistart::Int=1,
)
    y_vec = _float_vector(y)
    _assert_finite_observations("y", y_vec)
    length(indices) == length(y_vec) || throw(ArgumentError("indices and y must have equal length"))
    sigma_vec, cov_mat = _normalize_indexed_uncertainty(y_vec, sigma_y, cov_y)

    objective = function (p)
        yhat = model(indices, p)
        length(yhat) == length(y_vec) || throw(ArgumentError("model output length must match y"))
        residual = y_vec .- yhat
        return _gaussian_chi2_from_residual(residual, sigma_vec, cov_mat)
    end

    problem = LikelihoodFitProblem(
        objective,
        objective,
        p0;
        bounds=bounds,
        constraints=constraints,
        parameter_priors=parameter_priors,
        parameter_constraints=parameter_constraints,
        fixed_parameters=fixed_parameters,
        nobs=length(y_vec),
        cost_name=:indexed_chi2,
        parameter_names=parameter_names,
    )
    return fit(problem; maxiters=maxiters, tol=tol, ci_level=ci_level, initial_guesses=initial_guesses, multistart=multistart)
end

"""
    fit_multi_model(models, xs, ys; p0, sigma_y, kwargs...)

Fit multiple datasets simultaneously with one shared global parameter vector.
By default each `models[i](xs[i], p)` receives the full parameter vector `p`.
With `parameter_map`, model `i` receives `p[parameter_map[i]]` instead.
Per-dataset `sigma_y` entries are treated as physical standard deviations and
must be finite and positive.
"""
function fit_multi_model(
    models::AbstractVector,
    xs::AbstractVector,
    ys::AbstractVector;
    p0::AbstractVector,
    sigma_y=nothing,
    bounds=nothing,
    constraints=nothing,
    parameter_priors=nothing,
    parameter_constraints=nothing,
    fixed_parameters=nothing,
    parameter_names=nothing,
    parameter_map=nothing,
    maxiters::Int=1000,
    tol::Real=1e-10,
    ci_level::Real=0.6827,
    initial_guesses=nothing,
    multistart::Int=1,
)
    ndatasets = length(models)
    ndatasets > 0 || throw(ArgumentError("at least one dataset is required"))
    length(xs) == ndatasets || throw(ArgumentError("xs length must match models"))
    length(ys) == ndatasets || throw(ArgumentError("ys length must match models"))
    maps = parameter_map === nothing ? [nothing for _ in 1:ndatasets] : parameter_map
    length(maps) == ndatasets || throw(ArgumentError("parameter_map length must match models"))
    for i in 1:ndatasets
        maps[i] === nothing && continue
        all(j -> 1 <= j <= length(p0), maps[i]) || throw(ArgumentError("parameter_map[$i] contains an out-of-range parameter index"))
    end

    x_sets = [_float_vector(x) for x in xs]
    y_sets = [_float_vector(y) for y in ys]
    for i in 1:ndatasets
        _assert_finite_observations("x dataset $i", x_sets[i])
        _assert_finite_observations("y dataset $i", y_sets[i])
        length(x_sets[i]) == length(y_sets[i]) || throw(ArgumentError("dataset $i has mismatched x/y lengths"))
    end
    sigma_sets = _normalize_multi_sigma_sets(sigma_y, y_sets)

    objective = function (p)
        total = zero(eltype(p))
        for i in 1:ndatasets
            local_p = maps[i] === nothing ? p : p[maps[i]]
            yhat = models[i](x_sets[i], local_p)
            length(yhat) == length(y_sets[i]) || throw(ArgumentError("model $i output length mismatch"))
            residual = y_sets[i] .- yhat
            if sigma_sets[i] === nothing
                total += sum(abs2, residual)
            else
                total += sum(abs2, residual ./ sigma_sets[i])
            end
        end
        return total
    end

    nobs = sum(length, y_sets)
    problem = LikelihoodFitProblem(
        objective,
        objective,
        p0;
        bounds=bounds,
        constraints=constraints,
        parameter_priors=parameter_priors,
        parameter_constraints=parameter_constraints,
        fixed_parameters=fixed_parameters,
        nobs=nobs,
        cost_name=:multi_chi2,
        parameter_names=parameter_names,
    )
    return fit(problem; maxiters=maxiters, tol=tol, ci_level=ci_level, initial_guesses=initial_guesses, multistart=multistart)
end
