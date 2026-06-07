has_x_uncertainty(problem::FitProblem) = !(problem.sigma_x === nothing && problem.cov_x === nothing) || any(c -> c.active && c.target == :x, problem.error_components)
has_y_uncertainty(problem::FitProblem) = !(problem.sigma_y === nothing && problem.cov_y === nothing) || any(c -> c.active && c.target == :y, problem.error_components)
has_parameter_priors(problem) = !isempty(problem.parameter_priors)
has_parameter_constraints(problem) = !isempty(problem.parameter_constraints)

abstract type PreparedCovariance end

struct DynamicPreparedCovariance <: PreparedCovariance end
struct NoPreparedCovariance <: PreparedCovariance end

struct DiagonalPreparedCovariance{V<:Vector{Float64}} <: PreparedCovariance
    variances::V
    invsqrt::V
    logdet::Float64
end

struct DensePreparedCovariance{F} <: PreparedCovariance
    factor::F
    logdet::Float64
end

"""
    PreparedParameterConstraint

Internal cached representation of a correlated Gaussian parameter constraint.
The covariance factorization and log determinant are independent of the fitted
parameters, so they must be prepared once instead of recomputed inside every
objective call.
"""
struct PreparedParameterConstraint{F}
    indices::Vector{Int}
    mean::Vector{Float64}
    factor::F
    logdet::Float64
end

struct FitEvaluationCache{C<:PreparedCovariance, P}
    problem::FitProblem
    covariance::C
    parameter_constraints::P
end

# Validation checks need plain values, but numerical paths must keep Duals.
_finite_value(x) = x
_finite_value(x::ForwardDiff.Dual) = _finite_value(ForwardDiff.value(x))

function _has_model_relative_y_uncertainty(problem::FitProblem)
    return any(c -> c.active && c.target == :y && c.mode == :model_relative, problem.error_components)
end

function _static_effective_covariance_available(problem::FitProblem)
    has_x_uncertainty(problem) && return false
    _has_model_relative_y_uncertainty(problem) && return false
    return true
end

function _prepare_covariance(cov, n::Int)
    cov === nothing && return NoPreparedCovariance()

    if cov isa AbstractVector
        variances = collect(Float64, cov)
        length(variances) == n || throw(ArgumentError("covariance vector length must match observations"))
        any(variances .<= 0.0) && throw(ArgumentError("all effective variances must be positive"))
        return DiagonalPreparedCovariance(variances, inv.(sqrt.(variances)), sum(log, variances))
    end

    F = _stable_cholesky(cov)
    logdet = 2.0 * sum(log, diag(F.L))
    return DensePreparedCovariance(F, logdet)
end

function _prepare_parameter_constraint(constraint::ParameterConstraint)
    factor = _stable_cholesky(constraint.covariance)
    logdet = 2.0 * sum(log, diag(factor.L))
    return PreparedParameterConstraint(constraint.indices, constraint.mean, factor, logdet)
end

function _prepare_parameter_constraints(problem)
    isempty(problem.parameter_constraints) && return PreparedParameterConstraint[]
    return [_prepare_parameter_constraint(constraint) for constraint in problem.parameter_constraints]
end

function _prepare_fit_cache(problem::FitProblem)
    prepared_constraints = _prepare_parameter_constraints(problem)
    if !_static_effective_covariance_available(problem)
        return FitEvaluationCache(problem, DynamicPreparedCovariance(), prepared_constraints)
    end
    cov = _effective_covariance(problem, problem.p0)
    return FitEvaluationCache(problem, _prepare_covariance(cov, length(problem.y)), prepared_constraints)
end

function _model_values(problem::FitProblem, p::AbstractVector; x::AbstractVector=problem.x)
    values = problem.model(x, p)
    length(values) == length(x) || throw(ArgumentError("model output length must match x length"))
    collected = collect(values)
    finite_values = _finite_value.(collected)
    all(isfinite, finite_values) || throw(ArgumentError("model output must contain only finite values"))
    return collected
end

function _model_scalar(problem::FitProblem, x::Real, p::AbstractVector)
    y = problem.model([x], p)
    length(y) == 1 || throw(ArgumentError("model([x], p) must return exactly one value"))
    return y[1]
end

function _model_dydx(problem::FitProblem, p::AbstractVector; x::AbstractVector=problem.x)
    if problem.x_derivative !== nothing
        values = problem.x_derivative(x, p)
        length(values) == length(x) || throw(ArgumentError("x_derivative output length must match x length"))
        collected = collect(values)
        finite_values = _finite_value.(collected)
        all(isfinite, finite_values) || throw(ArgumentError("x_derivative output must contain only finite values"))
        return collected
    end

    return [ForwardDiff.derivative(t -> _model_scalar(problem, t, p), xi) for xi in x]
end

function _component_scale(values, n::Int)
    values isa Number && return fill(values, n)
    return values
end

function _component_covariance(component::ErrorComponent, problem::FitProblem, p::AbstractVector, target::Symbol)
    component.active || return nothing
    component.target == target || return nothing

    n = length(problem.y)
    values = component.values
    if component.mode == :covariance
        if values isa AbstractVector
            return values .^ 2
        end
        return values
    elseif component.mode == :absolute
        sigma = _component_scale(values, n)
        return sigma .^ 2
    elseif component.mode == :relative
        scale = _component_scale(values, n)
        reference = target == :y ? abs.(problem.y) : abs.(problem.x)
        return (scale .* reference) .^ 2
    elseif component.mode == :model_relative
        scale = _component_scale(values, n)
        yhat = _model_values(problem, p)
        return (scale .* abs.(yhat)) .^ 2
    end

    throw(ArgumentError("unsupported error component mode: $(component.mode)"))
end

function _component_covariance_sum(problem::FitProblem, p::AbstractVector, target::Symbol)
    total = nothing
    for component in problem.error_components
        cov = _component_covariance(component, problem, p, target)
        total = _combine_covariances(total, cov)
    end
    return total
end

function _base_y_covariance(problem::FitProblem, p::AbstractVector)
    base = nothing
    if problem.cov_y !== nothing
        base = problem.cov_y
    end
    if problem.sigma_y !== nothing
        base = _combine_covariances(base, problem.sigma_y .^ 2)
    end
    return _combine_covariances(base, _component_covariance_sum(problem, p, :y))
end

function _base_x_covariance(problem::FitProblem, p::AbstractVector)
    base = nothing
    if problem.cov_x !== nothing
        base = problem.cov_x
    end
    if problem.sigma_x !== nothing
        base = _combine_covariances(base, problem.sigma_x .^ 2)
    end
    return _combine_covariances(base, _component_covariance_sum(problem, p, :x))
end

function _x_covariance_term(problem::FitProblem, p::AbstractVector)
    has_x_uncertainty(problem) || return nothing

    dydx = _model_dydx(problem, p)
    base_x = _base_x_covariance(problem, p)
    base_x === nothing && return nothing

    if base_x isa AbstractVector
        return (dydx .^ 2) .* base_x
    end

    diag_dydx = Diagonal(dydx)
    return diag_dydx * base_x * diag_dydx
end

function _combine_covariances(base_cov, x_cov_term)
    if base_cov === nothing
        return x_cov_term
    end
    if x_cov_term === nothing
        return base_cov
    end

    if base_cov isa AbstractVector && x_cov_term isa AbstractVector
        return base_cov .+ x_cov_term
    end

    base_mat = base_cov isa AbstractVector ? Diagonal(base_cov) : base_cov
    x_mat = x_cov_term isa AbstractVector ? Diagonal(x_cov_term) : x_cov_term
    return base_mat + x_mat
end

function _effective_covariance(problem::FitProblem, p::AbstractVector)
    base_cov = _base_y_covariance(problem, p)
    x_cov_term = _x_covariance_term(problem, p)
    return _combine_covariances(base_cov, x_cov_term)
end

function _stable_cholesky(cov::AbstractMatrix)
    cov_values = Matrix{Float64}(_finite_value.(cov))
    all(isfinite, cov_values) || throw(ArgumentError("covariance matrix must contain only finite values"))
    isapprox(cov_values, cov_values'; rtol=1e-12, atol=1e-12) ||
        throw(ArgumentError("covariance matrix must be symmetric"))
    cov_sym = Symmetric(Matrix(cov))
    try
        return cholesky(cov_sym)
    catch
        throw(ArgumentError("covariance matrix must be symmetric positive definite"))
    end
end

function _whiten_residual(problem::FitProblem, p::AbstractVector, residual::AbstractVector)
    cov = _effective_covariance(problem, p)
    if cov === nothing
        return collect(residual)
    end

    if cov isa AbstractVector
        cov_values = _finite_value.(cov)
        any(cov_values .<= 0.0) && throw(ArgumentError("all effective variances must be positive"))
        return collect(residual ./ sqrt.(cov))
    end

    F = _stable_cholesky(cov)
    return collect(F.L \ residual)
end

function _whiten_residual(cache::FitEvaluationCache{NoPreparedCovariance}, p::AbstractVector, residual::AbstractVector)
    return collect(residual)
end

function _whiten_residual(cache::FitEvaluationCache{<:DiagonalPreparedCovariance}, p::AbstractVector, residual::AbstractVector)
    invsqrt = cache.covariance.invsqrt
    out = Vector{promote_type(eltype(residual), Float64)}(undef, length(residual))
    @inbounds for i in eachindex(residual)
        out[i] = residual[i] * invsqrt[i]
    end
    return out
end

function _whiten_residual(cache::FitEvaluationCache{<:DensePreparedCovariance}, p::AbstractVector, residual::AbstractVector)
    return collect(cache.covariance.factor.L \ residual)
end

function _whiten_residual(cache::FitEvaluationCache{DynamicPreparedCovariance}, p::AbstractVector, residual::AbstractVector)
    return _whiten_residual(cache.problem, p, residual)
end

function _residual(problem::FitProblem, p::AbstractVector)
    yhat = _model_values(problem, p)
    return problem.y .- yhat
end

function _weighted_residual(problem::FitProblem, p::AbstractVector)
    r = _residual(problem, p)
    rw_data = _whiten_residual(problem, p, r)
    if !has_parameter_priors(problem) && !has_parameter_constraints(problem)
        return rw_data
    end

    T = eltype(rw_data)
    n_constraint_terms = sum((length(c.indices) for c in problem.parameter_constraints); init=0)
    rw_priors = Vector{T}(undef, length(problem.parameter_priors) + n_constraint_terms)
    cursor = 1
    @inbounds for (i, prior) in enumerate(problem.parameter_priors)
        sigma = _asymmetric_sigma(p[prior.index], prior.mean, prior.sigma_minus, prior.sigma_plus)
        rw_priors[cursor] = (p[prior.index] - prior.mean) / sigma
        cursor += 1
    end
    @inbounds for constraint in problem.parameter_constraints
        delta = p[constraint.indices] .- constraint.mean
        z = _stable_cholesky(constraint.covariance).L \ delta
        rw_priors[cursor:(cursor + length(z) - 1)] .= z
        cursor += length(z)
    end
    return vcat(rw_data, rw_priors)
end

function _weighted_residual(cache::FitEvaluationCache, p::AbstractVector)
    problem = cache.problem
    r = _residual(problem, p)
    rw_data = _whiten_residual(cache, p, r)
    if !has_parameter_priors(problem) && isempty(cache.parameter_constraints)
        return rw_data
    end

    T = eltype(rw_data)
    n_constraint_terms = sum((length(c.indices) for c in cache.parameter_constraints); init=0)
    rw_priors = Vector{T}(undef, length(problem.parameter_priors) + n_constraint_terms)
    cursor = 1
    @inbounds for prior in problem.parameter_priors
        sigma = _asymmetric_sigma(p[prior.index], prior.mean, prior.sigma_minus, prior.sigma_plus)
        rw_priors[cursor] = (p[prior.index] - prior.mean) / sigma
        cursor += 1
    end
    @inbounds for constraint in cache.parameter_constraints
        delta = p[constraint.indices] .- constraint.mean
        z = constraint.factor.L \ delta
        rw_priors[cursor:(cursor + length(z) - 1)] .= z
        cursor += length(z)
    end
    return vcat(rw_data, rw_priors)
end

function _chi2(problem::FitProblem, p::AbstractVector)
    rw = _weighted_residual(problem, p)
    return sum(abs2, rw)
end

function _chi2(cache::FitEvaluationCache, p::AbstractVector)
    rw = _weighted_residual(cache, p)
    return sum(abs2, rw)
end

function _parameter_jacobian(problem::FitProblem, p::AbstractVector; x::AbstractVector=problem.x)
    if problem.jacobian !== nothing
        J = problem.jacobian(x, p)
        size(J, 1) == length(x) || throw(ArgumentError("jacobian row count must match x length"))
        size(J, 2) == length(p) || throw(ArgumentError("jacobian column count must match parameter count"))
        return Matrix{Float64}(J)
    end

    jac = ForwardDiff.jacobian(pp -> _model_values(problem, pp; x=x), p)
    return Matrix{Float64}(jac)
end

function _weighted_jacobian(problem::FitProblem, p::AbstractVector)
    jac = ForwardDiff.jacobian(pp -> _weighted_residual(problem, pp), p)
    return Matrix{Float64}(jac)
end

function _weighted_jacobian(cache::FitEvaluationCache, p::AbstractVector)
    jac = ForwardDiff.jacobian(pp -> _weighted_residual(cache, pp), p)
    return Matrix{Float64}(jac)
end

function _free_weighted_jacobian(cache::FitEvaluationCache, params::AbstractVector)
    problem = cache.problem
    free_idx = _free_indices(problem)
    jac = ForwardDiff.jacobian(q -> _weighted_residual(cache, _expand_free_parameters(problem, q)), params[free_idx])
    return Matrix{Float64}(jac)
end

function _covariance_from_weighted_jacobian(
    Jw::AbstractMatrix,
    chi2::Float64,
    ndf::Int,
    scale_covariance::Bool,
)
    fisher = Symmetric(Jw' * Jw)
    cov = _stable_symmetric_inverse(fisher)

    if scale_covariance && ndf > 0
        cov .*= (chi2 / ndf)
    end

    return cov
end

function _normalize_scale_covariance(scale_covariance)
    if scale_covariance === true
        return :always
    elseif scale_covariance === false
        return :never
    elseif scale_covariance in (:auto, :always, :never)
        return scale_covariance
    end
    throw(ArgumentError("scale_covariance must be :auto, :always, :never, true, or false"))
end

function _should_scale_covariance(problem::FitProblem, policy::Symbol)
    if policy == :always
        return true
    elseif policy == :never
        return false
    elseif policy == :auto
        return !has_y_uncertainty(problem) && !has_x_uncertainty(problem)
    end
    throw(ArgumentError("unsupported covariance scaling policy: $policy"))
end

function _covariance_from_cost_hessian(problem::FitProblem, p::AbstractVector, cost::Symbol)
    free_idx = _free_indices(problem)
    q = p[free_idx]
    H = ForwardDiff.hessian(qq -> _cost_value(problem, _expand_free_parameters(problem, qq), cost), q)
    cov = 2.0 .* _stable_symmetric_inverse(H)
    return _embed_free_covariance(problem, cov)
end

function _covariance_from_cost_hessian(cache::FitEvaluationCache, p::AbstractVector, cost::Symbol)
    problem = cache.problem
    free_idx = _free_indices(problem)
    q = p[free_idx]
    H = ForwardDiff.hessian(qq -> _cost_value(cache, _expand_free_parameters(problem, qq), cost), q)
    cov = 2.0 .* _stable_symmetric_inverse(H)
    return _embed_free_covariance(problem, cov)
end

function _correlation_from_covariance(cov::AbstractMatrix)
    sigma = sqrt.(clamp.(diag(cov), 0.0, Inf))
    corr = zeros(Float64, size(cov))
    for i in axes(cov, 1), j in axes(cov, 2)
        denom = sigma[i] * sigma[j]
        corr[i, j] = denom > 0 ? cov[i, j] / denom : 0.0
    end
    return corr
end

function _solve_backend(problem::FitProblem, backend::Symbol, cost::Symbol)
    if backend != :auto
        return backend
    end

    if _resolve_cost(problem, cost) != :chi2
        return :optimization
    end

    if has_parameter_priors(problem) || has_parameter_constraints(problem)
        # Priors add parameter-space penalty terms, so we use the generic objective path.
        return :optimization
    end

    if has_constraints(problem.constraints)
        return :optimization
    end

    if problem.bounds !== nothing
        return :optimization
    end

    if any(c -> c.active, problem.error_components)
        return :optimization
    end

    if _has_parameter_dependent_covariance(problem)
        # Effective-variance terms depend on parameters and are not static weights.
        return :optimization
    end

    return :lsqfit
end

function _constraint_vectors(spec::ConstraintSpec, p::AbstractVector)
    if spec.eq === nothing
        eq_vals = Vector{eltype(p)}()
    else
        raw = spec.eq(p)
        eq_vals = raw isa Number ? [raw] : collect(raw)
    end

    if spec.ineq === nothing
        ineq_vals = Vector{eltype(p)}()
    else
        raw = spec.ineq(p)
        ineq_vals = raw isa Number ? [raw] : collect(raw)
    end

    return eq_vals, ineq_vals
end

function _build_constraint_system(spec::ConstraintSpec, problem::FitProblem)
    eq0, ineq0 = _constraint_vectors(spec, _free_p0(problem))
    neq = length(eq0)
    nineq = length(ineq0)

    lcons = vcat(zeros(Float64, neq), fill(-Inf, nineq))
    ucons = vcat(zeros(Float64, neq), zeros(Float64, nineq))

    function cons!(res, p, data)
        eq_vals, ineq_vals = _constraint_vectors(spec, p)
        if neq > 0
            res[1:neq] .= eq_vals
        end
        if nineq > 0
            res[(neq + 1):(neq + nineq)] .= ineq_vals
        end
        return nothing
    end

    return cons!, lcons, ucons
end

function _diag_error(cov)
    cov === nothing && return nothing
    if cov isa AbstractVector
        return sqrt.(clamp.(cov, 0.0, Inf))
    end
    return sqrt.(clamp.(diag(Matrix(cov)), 0.0, Inf))
end

function _yerror_for_plot(problem::FitProblem, p::AbstractVector=problem.p0)
    return _diag_error(_base_y_covariance(problem, p))
end

function _xerror_for_plot(problem::FitProblem, p::AbstractVector=problem.p0)
    return _diag_error(_base_x_covariance(problem, p))
end
