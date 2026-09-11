"""
    LikelihoodFitProblem(objective, gof, p0; nobs, cost_name, kwargs...)

Low-level public problem representation for likelihood and custom-objective
fits. `objective(p)` is minimized directly, while optional `gof(p)` supplies a
chi-square-like data goodness-of-fit statistic for reduced statistics and
p-values. ScientificFitting adds the chi-square contributions from Gaussian parameter
terms to that statistic, matching their treatment as auxiliary observations in
the degrees of freedom.
Nonlinear constraint callbacks receive the same complete `p` vector, including
fixed parameters; ScientificFitting handles the reduced optimizer coordinates
internally.
Starting values must be finite. Fixed values are validated against declared
bounds while the problem is constructed.
For covariance, profile thresholds, and information criteria to have their
documented interpretation, `objective` must use the `-2 log(L)` scale.
Most users should prefer `fit_poisson_model`, `fit_histogram_model`,
`fit_unbinned_model`, `fit_extended_unbinned_model`, or `fit_custom`; construct
this type directly when the objective must be stored, inspected, or refitted.
`derivatives=:finite` selects numerical derivatives for foreign/Float64-only
objectives and constraint callbacks, including the covariance Hessian and all
profile refits. The default `:auto` uses ForwardDiff. Both modes require a
smooth objective near the evaluation point. Fitting defaults to `tol=1e-6`
in finite mode and `1e-10` otherwise; explicit tolerances are preserved.
"""
struct LikelihoodFitProblem{TF, TG, DM}
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
    derivatives::Symbol
end

_derivative_mode(::LikelihoodFitProblem{TF, TG, DM}) where {TF, TG, DM} = DM

"""
    LikelihoodFitResult

Result of a likelihood or custom-objective fit. It mirrors the parameter,
covariance, statistics, and diagnostics fields of `FitResult`, but stores a
`LikelihoodFitProblem` instead of x-y residual data. Plotting support depends on
the specific likelihood workflow because not every objective has a natural
curve representation.

`iterations` is `missing` when the selected optimizer does not expose an
iteration count.
`solver_result` retains the [`FitSolverResult`](@ref) and its native details;
it is `nothing` when every parameter is fixed. Native solver coordinates follow
`solver_result.parameter_indices`, not necessarily the full model vector.
"""
struct LikelihoodFitResult
    problem::LikelihoodFitProblem
    options::FitOptions
    backend::Symbol
    converged::Bool
    iterations::Union{Int, Missing}
    message::String
    params::Vector{Float64}
    param_stderr::Vector{Float64}
    param_covariance::Matrix{Float64}
    param_correlation::Matrix{Float64}
    stats::FitStatistics
    diagnostics::FitDiagnostics
    # Do not specialize reports/plots on every native solution and model type.
    solver_result::Union{Nothing, FitSolverResult}
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
    derivatives::Symbol=:auto,
)
    _validate_derivatives(derivatives)
    p0_vec = _float_vector(p0)
    length(p0_vec) > 0 || throw(ArgumentError("p0 must contain at least one parameter"))
    _assert_finite_vector("p0", p0_vec)
    nobs > 0 || throw(ArgumentError("nobs must be positive"))
    names = parameter_names === nothing ? nothing : collect(String, parameter_names)
    names === nothing || length(names) == length(p0_vec) || throw(ArgumentError("parameter_names length must match p0"))
    names === nothing || (all(name -> !isempty(name), names) && allunique(names)) ||
        throw(ArgumentError("parameter_names must be nonempty and unique"))

    normalized_bounds = _normalize_bounds(bounds, length(p0_vec))
    normalized_fixed = _normalize_fixed_parameters(fixed_parameters, length(p0_vec))
    _assert_fixed_parameters_within_bounds(normalized_fixed, normalized_bounds)

    return LikelihoodFitProblem{typeof(objective), typeof(gof), derivatives}(
        objective,
        gof,
        p0_vec,
        normalized_bounds,
        _normalize_constraints(constraints),
        _normalize_parameter_priors(parameter_priors, length(p0_vec)),
        _normalize_parameter_constraints(parameter_constraints, length(p0_vec)),
        normalized_fixed,
        Int(nobs),
        cost_name,
        names,
        derivatives,
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
        derivatives=problem.derivatives,
    )
end

function _likelihood_cost(problem::LikelihoodFitProblem, p::AbstractVector)
    return problem.objective(p) + _prior_minus2loglik(problem, p) + _parameter_constraint_minus2loglik(problem, p)
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
    return problem.objective(p) +
           _prior_minus2loglik(problem, p) +
           _parameter_constraint_minus2loglik(cache.parameter_constraints, p)
end

function _likelihood_gof(cache::LikelihoodEvaluationCache, p::AbstractVector)
    problem = cache.problem
    problem.gof === nothing && return NaN
    return Float64(
        problem.gof(p) +
        _prior_chi2(problem, p) +
        _parameter_constraint_chi2(cache.parameter_constraints, p),
    )
end

function _fit_likelihood_problem(problem::LikelihoodFitProblem, options::FitOptions)
    cache = _prepare_likelihood_cache(problem)
    # Retain the objective type in its AD tag, just as for Gaussian fits.
    objective = (q, cache) -> _likelihood_cost(cache, _expand_free_parameters(problem, q))
    return _minimize_scalar(problem, options, objective, cache)
end

"""Compute local curvature once, or retain explicit missing free-parameter errors."""
function _likelihood_geometry(cache::LikelihoodEvaluationCache, params::Vector{Float64}, method::Symbol)
    problem = cache.problem
    free_idx = _free_indices(problem)
    if isempty(free_idx) || method == :none
        return _embed_free_covariance(problem, fill(NaN, length(free_idx), length(free_idx))), nothing
    end

    q = params[free_idx]
    H = _derivative_hessian(problem, qq -> _likelihood_cost(cache, _expand_free_parameters(problem, qq)), q)
    free_cov = 2.0 .* _stable_symmetric_inverse(H)
    return _embed_free_covariance(problem, free_cov), H
end

function _build_likelihood_result(
    problem::LikelihoodFitProblem,
    options::FitOptions,
    params::Vector{Float64},
    converged::Bool,
    iterations::Union{Int, Missing},
    message::String,
    solver_result=nothing,
)
    cache = _prepare_likelihood_cache(problem)
    cost_min = Float64(_likelihood_cost(cache, params))
    gof = _likelihood_gof(cache, params)
    nconstraint_obs = sum((length(c.indices) for c in problem.parameter_constraints); init=0)
    nobs = problem.nobs + length(problem.parameter_priors) + nconstraint_obs
    npar = length(_free_indices(problem))
    ndf = nobs - npar
    gof_ndf = isfinite(gof) && ndf > 0 ? gof / ndf : NaN
    pvalue = isfinite(gof) && ndf > 0 ? ccdf(Chisq(ndf), gof) : NaN
    aic = cost_min + 2.0 * npar
    bic = cost_min + log(nobs) * npar
    cov, hessian = _likelihood_geometry(cache, params, options.parameter_covariance)
    stderr = _standard_errors_from_covariance(cov)
    corr = _correlation_from_covariance(cov)
    stats = FitStatistics(problem.cost_name, cost_min, cost_min, gof, gof_ndf, ndf, pvalue, aic, bic)
    diagnostics = _fit_diagnostics(problem, params, cov, converged, ndf;
        hessian, gof, covariance_computed=options.parameter_covariance != :none)

    backend = solver_result === nothing ? :optimization : solver_result.backend
    return LikelihoodFitResult(problem, options, backend, converged, iterations, message,
                               params, stderr, cov, corr, stats, diagnostics, solver_result)
end

"""
    fit(problem::LikelihoodFitProblem; maxiters=1000, tol=1e-10,
        initial_guesses=nothing, multistart=1, optimizer=:auto,
        parameter_covariance=:auto, solver=nothing) -> LikelihoodFitResult

Minimize a validated likelihood-scale or custom objective problem. Bounds,
fixed parameters, Gaussian parameter terms, nonlinear constraints, observation
count, and cost name are already stored in `problem`.

`initial_guesses` may contain additional complete parameter vectors. `multistart`
is the total candidate budget including `p0`; its default of one uses only `p0`.
Remaining slots use deterministic bounded candidates when finite bounds are
available. ScientificFitting returns the converged candidate with the lowest finite cost;
if only non-converged finite candidates remain, the best one is returned with
`converged == false`. If no candidate produces a finite result, the last
objective, validation, or solver error is raised.

`optimizer=:auto` selects LBFGS, or IPNewton with nonlinear constraints.
Alternatively pass `solver=OptimizationSolver(algorithm)` or
`solver=NativeMinuitSolver()` while leaving `optimizer=:auto`. Solver objects
and their native settings are preserved in profile refits. This selection is
independent of the requested local covariance calculation.
`:nelder_mead` uses NLopt's bounded derivative-free simplex method for non-smooth
or support-limited costs. It retains bounds, fixed parameters and Gaussian
parameter terms, but rejects nonlinear constraints; those require `:ipnewton`.
Start at finite cost inside the distribution's support. All methods are local
optimizers of continuous parameters, not guarantees of a global minimum.

`parameter_covariance=:auto` chooses `:none` with Nelder-Mead and `:hessian`
otherwise. Explicit `:hessian` computes `2 * inv(H)` for the complete cost;
use it only for locally smooth likelihoods on the documented `-2 log(L)` scale.
`:none` skips curvature and stores NaN for free-parameter errors (fixed errors
remain zero). Profile refits preserve both options. Supply explicit profile
ranges when local errors are unavailable; profile thresholds still require
statistical justification for non-regular models.

`maxiters` and `tol` must be positive. With Nelder-Mead, `maxiters` limits
objective evaluations, `tol` sets NLopt's absolute/relative parameter tolerances,
and `iterations` is `missing`, not an invented count. Scale parameters accordingly;
`tol` is not a statistical error. Function-value stopping is disabled because
equal costs need not mean a contracted simplex. A budget-limited solve returns
`converged == false`.
"""
function fit(
    problem::LikelihoodFitProblem;
    maxiters::Int=1000,
    tol::Real=_default_fit_tolerance(problem.derivatives),
    initial_guesses=nothing,
    multistart::Int=1,
    optimizer::Symbol=:auto,
    parameter_covariance::Symbol=:auto,
    solver=nothing,
)
    solver !== nothing && optimizer != :auto && throw(ArgumentError(
        "choose solver or optimizer, not both; leave optimizer=:auto with an explicit solver",
    ))
    constrained = has_constraints(problem.constraints)
    optimizer = solver === nothing && optimizer == :auto ? (constrained ? :ipnewton : :lbfgs) : optimizer
    caps = solver isa AbstractFitSolver ? solver_capabilities(solver) : nothing
    derivative_free = caps === nothing ? optimizer == :nelder_mead : !(caps.gradient || caps.hessian)
    parameter_covariance = parameter_covariance == :auto ?
        (derivative_free ? :none : :hessian) : parameter_covariance
    options = FitOptions(
        backend=:optimization,
        cost=problem.cost_name,
        maxiters=maxiters,
        tol=Float64(tol),
        scale_covariance=:never,
        multistart=multistart,
        optimizer=optimizer,
        parameter_covariance=parameter_covariance,
        solver=solver,
    )
    solver === nothing && constrained && optimizer != :ipnewton && throw(ArgumentError(
        "nonlinear constraints require optimizer=:auto or :ipnewton; they cannot be dropped",
    ))
    solver === nothing && !constrained && optimizer == :ipnewton && throw(ArgumentError(
        "optimizer=:ipnewton requires nonlinear constraints; use :auto, :lbfgs, or :nelder_mead",
    ))

    candidates = _initial_candidates(problem, initial_guesses, multistart)
    best_result = nothing
    last_error = nothing

    for candidate in candidates
        candidate_problem = _with_p0(problem, candidate)
        try
            result = if isempty(_free_indices(candidate_problem))
                params = _expand_free_parameters(candidate_problem, Float64[])
                _build_likelihood_result(candidate_problem, options, params, true, 0, "All parameters fixed")
            else
                answer = _fit_likelihood_problem(candidate_problem, options)
                params = _expand_free_parameters(candidate_problem, answer.params)
                _build_likelihood_result(candidate_problem, options, params, answer.converged,
                                         answer.iterations, answer.message, answer)
            end
            _prefer_fit(result, best_result) && (best_result = result)
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
    fit_custom(objective; p0, nobs, gof=nothing, cost_name=:custom,
               kwargs...) -> LikelihoodFitResult

Fit a user-defined scalar objective. `objective(p)` is minimized directly.
If `gof(p)` is supplied it is used for reduced goodness-of-fit statistics and
p-values; otherwise these fields are `NaN`. `nobs` controls degrees of freedom
and the BIC sample-size term, and must represent the number of statistically
independent observations used by the objective.

When requested, ScientificFitting computes local covariance as `2 * inv(H)`, where `H` is the objective
Hessian. AIC, BIC, and that covariance therefore have their standard
interpretation only when `objective` is a normalized `-2log(L)` cost (Gaussian
chi-square is on the same scale). For an arbitrary loss, the optimizer result
remains usable but those inferential fields are only arithmetic summaries.

Common keywords are `bounds`, `constraints`, `parameter_priors`,
`parameter_constraints`, `fixed_parameters`, `parameter_names`, `maxiters`,
`tol`, `initial_guesses`, `multistart`, and `derivatives=:auto` (or `:finite`).
`optimizer` and `parameter_covariance` independently control minimization and
local errors; see [`fit(::LikelihoodFitProblem)`](@ref) for choices and limits.
Parameter callbacks receive the
complete vector in `p0` order. `nobs` must be positive; invalid parameter
controls or a non-finite objective fail with an error rather than producing a
reportable result.

# Example

```julia
cost(p) = ((p[1] - 2.0) / 0.3)^2
result = fit_custom(cost; p0=[0.0], nobs=1, cost_name=:calibration_chi2)
```
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
    derivatives::Symbol=:auto,
    maxiters::Int=1000,
    tol::Real=_default_fit_tolerance(derivatives),
    initial_guesses=nothing,
    multistart::Int=1,
    optimizer::Symbol=:auto,
    parameter_covariance::Symbol=:auto,
    solver=nothing,
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
        derivatives=derivatives,
    )
    return fit(problem; maxiters, tol, initial_guesses, multistart, optimizer, parameter_covariance, solver)
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

"""
    fit_likelihood_model(model, x, y; logprob=nothing, error=nothing, p0, kwargs...)
        -> LikelihoodFitResult

Fit independent observations with a user-defined continuous or discrete
observation distribution. `model(x, p)` returns one prediction per observation;
`logprob(y, prediction, p)` returns a vector of normalized log densities or
log probability masses. The two callbacks are evaluated once per objective
evaluation, allowing vectorized Julia or Python models. Capture per-observation
scales, trial counts, or other known inputs in the `logprob` closure.

Alternatively, `error=distribution` models the additive residual `y - model(x,p)`
directly with a fixed Distributions.jl distribution. A univariate object applies
independently to every residual; a vector of objects gives point-specific errors.
A multivariate object describes the complete residual vector jointly, including
dependence between measurements. Its dimension must equal the number of values.
Use exactly one of `error` and `logprob`; no Gaussian approximation is introduced.
For prediction-dependent distributions such as Poisson counts, keep `logprob`;
these describe observations, not additive continuous errors.

The core minimizes `-2 * sum(logprob(...))`. Include all normalization terms,
especially those depending on fitted parameters. A log density may be positive;
`-Inf` denotes zero probability and yields infinite cost, without clipping.
NaN, positive infinity, wrong output dimensions, or non-finite observations
raise `ArgumentError`. Begin at finite cost inside the distribution's support.
The package cannot verify that an arbitrary callback is normalized.

The default gradient/Hessian backend requires a smooth objective near evaluated
parameters, even for discrete *observations*. Choose `optimizer=:nelder_mead`
for non-smooth or moving-support likelihoods; it does not compute Hessian errors
unless `parameter_covariance=:hessian` is explicitly requested. This does not
support discrete fitted parameters. Use `fit_custom` for dependent observations with a joint
likelihood rather than multiplying their marginal probabilities.

Parameter controls, derivatives, and solver options follow `fit_custom`.
Requested local covariance and profiles use the same complete likelihood. There is no
universal chi-square goodness statistic: `chi2`, `chi2_ndf`, and `pvalue` are
NaN unless a justified `gof(p)` is supplied. Profile coverage is asymptotic,
not automatically guaranteed for every distribution or parameter boundary.

# Example

```julia
using Distributions
model(x, p) = p[1] .* x .+ p[2]
sigma = [0.2, 0.3, 0.2, 0.4]
# Student-t errors allow heavier tails; sigma is a scale, not a standard deviation.
logprob(y, mu, p) = logpdf.(TDist(4), (y .- mu) ./ sigma) .- log.(sigma)
result = fit_likelihood_model(model, [0., 1., 2., 3.], [0.1, 1.2, 1.9, 3.4];
    logprob=logprob, p0=[1., 0.])
```
"""
function fit_likelihood_model(
    model, x::AbstractVector, y::AbstractVector;
    logprob=nothing, error=nothing, p0::AbstractVector,
    cost_name::Symbol=:observation_likelihood, kwargs...,
)
    x_vec, y_vec = _float_vector(x), _float_vector(y)
    length(x_vec) == length(y_vec) || throw(ArgumentError("x and y must have equal length"))
    _assert_finite_observations("x", x_vec)
    _assert_finite_observations("y", y_vec)
    (logprob === nothing) != (error === nothing) || throw(ArgumentError(
        "provide exactly one of logprob or error",
    ))
    error_model = error === nothing ? nothing : _prepare_error_distribution(error, length(y_vec))
    objective = function (p)
        prediction = model(x_vec, p)
        prediction isa AbstractVector && length(prediction) == length(y_vec) ||
            throw(ArgumentError("model must return one prediction per observation"))
        all(isfinite, prediction) || throw(ArgumentError("model predictions must be finite"))
        error_model === nothing || return _minus2logprob(_error_loglikelihood(error_model, y_vec .- prediction))
        terms = logprob(y_vec, prediction, p)
        terms isa AbstractVector && length(terms) == length(y_vec) ||
            throw(ArgumentError("logprob must return one log probability per observation"))
        # Zero support stays impossible; never replace it with an arbitrary floor.
        all(v -> v isa Real && (isfinite(v) || v == -Inf), terms) ||
            throw(ArgumentError("logprob values must be finite or -Inf (zero probability)"))
        return -2 * sum(terms)
    end
    return fit_custom(objective; p0=p0, nobs=length(y_vec), cost_name=cost_name, kwargs...)
end

function _positive_expectation(mu, n::Int)
    values = collect(mu)
    length(values) == n || throw(ArgumentError("model expectation length must match observations"))
    all(isfinite, values) || throw(ArgumentError("model expectation contains non-finite values"))
    all(values .> 0) || throw(ArgumentError("model expectation must be strictly positive"))
    return values
end

function _poisson_minus2loglik_terms(counts::Vector{Float64}, mu::AbstractVector)
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
    fit_poisson_model(model, x, counts; p0, kwargs...) -> LikelihoodFitResult

Fit count data with a Poisson likelihood. `model(x, p)` must return the
strictly positive expected counts for each observation. `counts` must contain
finite non-negative integer-valued observations; `x` must be finite and have
the same length.

The minimized objective is normalized Poisson `-2 log(L)`. `stats.chi2` stores
the Poisson deviance and its p-value uses the asymptotic chi-square reference.
Common parameter-control and solver keywords are listed under `fit_custom`.
Invalid counts, dimensions, expectations, or parameter controls raise
`ArgumentError`. Returns `LikelihoodFitResult`.
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
    derivatives::Symbol=:auto,
    maxiters::Int=1000,
    tol::Real=_default_fit_tolerance(derivatives),
    initial_guesses=nothing,
    multistart::Int=1,
    optimizer::Symbol=:auto,
    parameter_covariance::Symbol=:auto,
    solver=nothing,
)
    x_vec = _float_vector(x)
    counts_vec = _float_vector(counts)
    length(x_vec) == length(counts_vec) || throw(ArgumentError("x and counts must have equal length"))
    _assert_finite_observations("x", x_vec)
    _assert_count_observations("counts", counts_vec)

    objective = p -> _poisson_minus2loglik_terms(counts_vec, _positive_expectation(model(x_vec, p), length(counts_vec)))
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
        cost_name=:poisson_likelihood,
        parameter_names=parameter_names,
        derivatives=derivatives,
    )
    return fit(problem; maxiters, tol, initial_guesses, multistart, optimizer, parameter_covariance, solver)
end

"""
    fit_histogram_model(expected_counts, edges, counts; p0, kwargs...)
        -> LikelihoodFitResult

Fit binned counts. `expected_counts(edges, p)` must return one positive expected
count per bin. Histogram edges must be finite and strictly increasing; `counts`
must contain finite non-negative integer-valued observations.

The objective is normalized Poisson `-2 log(L)` and `stats.chi2` is Poisson
deviance. The caller is responsible for integrating any continuous density over
the bins; use `fit_histogram_density` when ScientificFitting should perform that
integration. Common parameter-control and solver keywords are listed under
`fit_custom`. Invalid edges, counts, expectations, or controls raise
`ArgumentError`. Returns `LikelihoodFitResult`.
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
    derivatives::Symbol=:auto,
    maxiters::Int=1000,
    tol::Real=_default_fit_tolerance(derivatives),
    initial_guesses=nothing,
    multistart::Int=1,
    optimizer::Symbol=:auto,
    parameter_covariance::Symbol=:auto,
    solver=nothing,
)
    edges_vec = _float_vector(edges)
    counts_vec = _float_vector(counts)
    length(edges_vec) == length(counts_vec) + 1 || throw(ArgumentError("edges length must be count length + 1"))
    _assert_finite_observations("histogram edges", edges_vec)
    _assert_count_observations("counts", counts_vec)
    any(diff(edges_vec) .<= 0) && throw(ArgumentError("histogram edges must be strictly increasing"))

    objective = p -> _poisson_minus2loglik_terms(
        counts_vec,
        _positive_expectation(expected_counts(edges_vec, p), length(counts_vec)),
    )
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
        cost_name=:histogram_poisson_likelihood,
        parameter_names=parameter_names,
        derivatives=derivatives,
    )
    return fit(problem; maxiters, tol, initial_guesses, multistart, optimizer, parameter_covariance, solver)
end

"""Validate the vectorized density contract before reduction or buffer copying."""
function _density_batch(pdf, x, p)
    values = pdf(x, p)
    values isa AbstractVector && axes(values) == axes(x) ||
        throw(ArgumentError("vectorized density must return one value per input point"))
    return values
end

"""Sum the same event log likelihood for scalar and batched model evaluation."""
function _density_logcost(pdf, data, p, vectorized::Bool)
    values = vectorized ? _density_batch(pdf, data, p) : (pdf(x, p) for x in data)
    total = zero(eltype(p))
    @inbounds for density in values
        density isa Real && isfinite(density) && density > 0 ||
            throw(ArgumentError("density must be finite and strictly positive at all data points"))
        total -= 2 * log(density)
    end
    return total
end

"""Use QuadGK's native batching and adaptive error control, retaining parameter duals."""
function _density_integrand(pdf, p, vectorized::Bool)
    vectorized || return x -> pdf(x, p)
    evaluate!(out, x) = copyto!(out, _density_batch(pdf, x, p))
    return BatchIntegrand{eltype(p), Float64}(evaluate!)
end

"""
    fit_histogram_density(pdf, edges, counts; p0, total_count=sum(counts),
                          rtol=1e-8, vectorized=false, kwargs...) -> LikelihoodFitResult

Fit binned counts from a probability density. Expected bin counts are computed
with adaptive Gauss-Kronrod quadrature. `pdf(x, p)` must be normalized on its
intended physical domain; ScientificFitting does not renormalize it over the supplied
bins. `total_count` and `rtol` must be finite and positive.

With `vectorized=true`, `pdf(xs, p)` receives a vector of quadrature nodes and
must return one value per node. QuadGK batches these evaluations while retaining
adaptive error control for each bin; the scalar callback remains the default.

Each expectation is `total_count * integral(pdf, edge[i], edge[i+1])`.
Consequently, bins outside the supplied range still carry probability unless
the density is normalized on that range. Common parameter-control and solver
keywords are listed under `fit_custom`. Quadrature/model failures are propagated;
invalid histogram input raises `ArgumentError`. Returns `LikelihoodFitResult`.
"""
function fit_histogram_density(
    pdf,
    edges::AbstractVector,
    counts::AbstractVector;
    p0::AbstractVector,
    total_count::Real=sum(counts),
    rtol::Real=1e-8,
    vectorized::Bool=false,
    bounds=nothing,
    constraints=nothing,
    parameter_priors=nothing,
    parameter_constraints=nothing,
    fixed_parameters=nothing,
    parameter_names=nothing,
    derivatives::Symbol=:auto,
    maxiters::Int=1000,
    tol::Real=_default_fit_tolerance(derivatives),
    initial_guesses=nothing,
    multistart::Int=1,
    optimizer::Symbol=:auto,
    parameter_covariance::Symbol=:auto,
    solver=nothing,
)
    total = Float64(total_count)
    isfinite(total) && total > 0 || throw(ArgumentError("total_count must be finite and > 0"))
    isfinite(rtol) && rtol > 0 || throw(ArgumentError("rtol must be finite and > 0"))

    expected_counts = function (edge_values, p)
        mu = Vector{eltype(p)}(undef, length(edge_values) - 1)
        # Reuse the batch buffers across bins within this objective evaluation.
        integrand = _density_integrand(pdf, p, vectorized)
        @inbounds for i in eachindex(mu)
            integral, _ = quadgk(integrand, edge_values[i], edge_values[i + 1]; rtol=rtol)
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
        derivatives=derivatives,
        maxiters=maxiters,
        tol=tol,
        initial_guesses=initial_guesses,
        multistart=multistart,
        optimizer=optimizer,
        parameter_covariance=parameter_covariance,
        solver=solver,
    )
end

"""
    fit_unbinned_model(pdf, data; p0, vectorized=false, kwargs...) -> LikelihoodFitResult

Fit independent unbinned observations with a normalized positive density
`pdf(x, p)`. ScientificFitting checks positivity at the observations but cannot infer or
verify the normalization domain. Observations must be finite.
With `vectorized=true`, `pdf(data, p)` evaluates all observations in one call
and must return one density per observation. The scalar callback is the default.

The objective is `-2 * sum(log(pdf(x_i, p)))`. No universal chi-square
goodness-of-fit statistic exists, so `chi2`, `chi2_ndf`, and `pvalue` are `NaN`.
Common parameter-control and solver keywords are listed under `fit_custom`.
Non-positive densities and invalid controls raise `ArgumentError`. Returns
`LikelihoodFitResult`.
"""
function fit_unbinned_model(
    pdf,
    data::AbstractVector;
    p0::AbstractVector,
    vectorized::Bool=false,
    bounds=nothing,
    constraints=nothing,
    parameter_priors=nothing,
    parameter_constraints=nothing,
    fixed_parameters=nothing,
    parameter_names=nothing,
    derivatives::Symbol=:auto,
    maxiters::Int=1000,
    tol::Real=_default_fit_tolerance(derivatives),
    initial_guesses=nothing,
    multistart::Int=1,
    optimizer::Symbol=:auto,
    parameter_covariance::Symbol=:auto,
    solver=nothing,
)
    data_vec = _float_vector(data)
    _assert_finite_observations("unbinned data", data_vec)
    objective = p -> _density_logcost(pdf, data_vec, p, vectorized)
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
        cost_name=:unbinned_likelihood,
        parameter_names=parameter_names,
        derivatives=derivatives,
    )
    return fit(problem; maxiters, tol, initial_guesses, multistart, optimizer, parameter_covariance, solver)
end

"""
    fit_extended_unbinned_model(rate, data, domain; p0, rtol=1e-8,
                                vectorized=false, kwargs...) -> LikelihoodFitResult

Fit an inhomogeneous Poisson point process. `rate(x, p)` is the event intensity,
not a normalized density. `domain=(a, b)` defines the integration range for the
expected total event count. The domain endpoints must be finite, and all
observations must lie inside the domain.
With `vectorized=true`, `rate(xs, p)` returns one intensity per input point:
all observations are evaluated together, and QuadGK batches the quadrature nodes.

The objective is `2 * integral(rate, domain) - 2 * sum(log(rate(x_i, p)))`.
No generic chi-square p-value is reported. `rtol` controls Gauss-Kronrod
integration and must be positive. Common parameter-control and solver keywords
are listed under `fit_custom`. Invalid domains, non-positive rates, and
integration/model failures raise an error. Returns `LikelihoodFitResult`.
"""
function fit_extended_unbinned_model(
    rate,
    data::AbstractVector,
    domain::Tuple{<:Real, <:Real};
    p0::AbstractVector,
    rtol::Real=1e-8,
    vectorized::Bool=false,
    bounds=nothing,
    constraints=nothing,
    parameter_priors=nothing,
    parameter_constraints=nothing,
    fixed_parameters=nothing,
    parameter_names=nothing,
    derivatives::Symbol=:auto,
    maxiters::Int=1000,
    tol::Real=_default_fit_tolerance(derivatives),
    initial_guesses=nothing,
    multistart::Int=1,
    optimizer::Symbol=:auto,
    parameter_covariance::Symbol=:auto,
    solver=nothing,
)
    data_vec = _float_vector(data)
    a, b = Float64(domain[1]), Float64(domain[2])
    isfinite(a) && isfinite(b) || throw(ArgumentError("domain endpoints must be finite"))
    a < b || throw(ArgumentError("domain must satisfy domain[1] < domain[2]"))
    isfinite(rtol) && rtol > 0 || throw(ArgumentError("rtol must be finite and > 0"))
    _assert_finite_observations("unbinned data", data_vec)
    all(x -> a <= x <= b, data_vec) || throw(ArgumentError("extended unbinned data must lie inside the domain"))

    objective = function (p)
        expected, _ = quadgk(_density_integrand(rate, p, vectorized), a, b; rtol=rtol)
        expected > 0 || throw(ArgumentError("integrated rate must be positive"))
        return 2 * expected + _density_logcost(rate, data_vec, p, vectorized)
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
        cost_name=:extended_unbinned_likelihood,
        parameter_names=parameter_names,
        derivatives=derivatives,
    )
    return fit(problem; maxiters, tol, initial_guesses, multistart, optimizer, parameter_covariance, solver)
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
    fit_indexed_model(model, indices, y; p0, sigma_y=nothing, cov_y=nothing,
                      kwargs...) -> LikelihoodFitResult

Fit observations addressed by arbitrary indices. `model(indices, p)` must return
one model value per index. This is useful when the independent variable is not a
numeric 1D x-axis. Optional `sigma_y` entries must be finite and positive;
optional `cov_y` must be a finite symmetric positive-definite covariance matrix.

The minimized objective is chi-square without additive Gaussian normalization
constants. AIC/BIC may compare models fit to the same observations and the same
uncertainty model, but not different uncertainty scales or datasets.

Use either `sigma_y` or `cov_y`, never both. With neither, the fit is
unweighted. Common parameter-control and solver keywords are listed under
`fit_custom`. Invalid dimensions, uncertainty, covariance, model output, or
controls raise `ArgumentError`. Returns `LikelihoodFitResult`.
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
    derivatives::Symbol=:auto,
    maxiters::Int=1000,
    tol::Real=_default_fit_tolerance(derivatives),
    initial_guesses=nothing,
    multistart::Int=1,
    optimizer::Symbol=:auto,
    parameter_covariance::Symbol=:auto,
    solver=nothing,
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
        derivatives=derivatives,
    )
    return fit(problem; maxiters, tol, initial_guesses, multistart, optimizer, parameter_covariance, solver)
end

"""
    fit_multi_model(models, xs, ys; p0, sigma_y=nothing,
                    parameter_map=nothing, kwargs...) -> LikelihoodFitResult

Fit multiple datasets simultaneously with one shared global parameter vector.
By default each `models[i](xs[i], p)` receives the full parameter vector `p`.
With `parameter_map`, model `i` receives `p[parameter_map[i]]` instead.
Per-dataset `sigma_y` entries are treated as physical standard deviations and
must be finite and positive.

The minimized objective is the summed chi-square without additive Gaussian
normalization constants. AIC/BIC may compare models fit to the same datasets
and uncertainty model, but not different uncertainty scales or datasets.

`sigma_y` is either `nothing` or one uncertainty vector (or `nothing`) per
dataset. `parameter_map[i]` contains one-based indices into the global `p0` and
defines the local parameter order passed to model `i`. Common parameter-control
and solver keywords are listed under `fit_custom`. Empty or mismatched dataset
collections, invalid maps or uncertainties, and wrong model-output lengths
raise `ArgumentError`. Returns `LikelihoodFitResult`.
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
    derivatives::Symbol=:auto,
    parameter_map=nothing,
    maxiters::Int=1000,
    tol::Real=_default_fit_tolerance(derivatives),
    initial_guesses=nothing,
    multistart::Int=1,
    optimizer::Symbol=:auto,
    parameter_covariance::Symbol=:auto,
    solver=nothing,
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
        derivatives=derivatives,
    )
    return fit(problem; maxiters, tol, initial_guesses, multistart, optimizer, parameter_covariance, solver)
end
