const LOG2PI = log(2.0 * pi)

function _resolve_cost(problem::FitProblem, cost::Symbol)
    if cost == :auto
        return _has_parameter_dependent_covariance(problem) ? :gaussian_nll : :chi2
    elseif cost in (:chi2, :least_squares)
        return :chi2
    elseif cost in (:gaussian_nll, :nll)
        return :gaussian_nll
    end

    throw(ArgumentError("unsupported cost: $cost (use :auto, :chi2, or :gaussian_nll)"))
end

function _asymmetric_sigma(value, center, sigma_minus, sigma_plus)
    return value < center ? sigma_minus : sigma_plus
end

function _weighted_data_residual(problem::FitProblem, p::AbstractVector)
    return _whiten_residual(problem, p, _residual(problem, p))
end

function _prior_chi2(problem, p::AbstractVector)
    isempty(problem.parameter_priors) && return zero(eltype(p))

    total = zero(eltype(p))
    @inbounds for prior in problem.parameter_priors
        sigma = _asymmetric_sigma(p[prior.index], prior.mean, prior.sigma_minus, prior.sigma_plus)
        z = (p[prior.index] - prior.mean) / sigma
        total += abs2(z)
    end
    return total
end

function _parameter_constraint_chi2(problem, p::AbstractVector)
    isempty(problem.parameter_constraints) && return zero(eltype(p))

    return _parameter_constraint_chi2(_prepare_parameter_constraints(problem), p)
end

function _parameter_constraint_chi2(prepared_constraints::AbstractVector{<:PreparedParameterConstraint}, p::AbstractVector)
    isempty(prepared_constraints) && return zero(eltype(p))

    total = zero(eltype(p))
    @inbounds for constraint in prepared_constraints
        delta = p[constraint.indices] .- constraint.mean
        z = constraint.factor.L \ delta
        total += sum(abs2, z)
    end
    return total
end

function _prior_nll(problem, p::AbstractVector)
    isempty(problem.parameter_priors) && return zero(eltype(p))

    total = zero(eltype(p))
    @inbounds for prior in problem.parameter_priors
        sigma = _asymmetric_sigma(p[prior.index], prior.mean, prior.sigma_minus, prior.sigma_plus)
        z = (p[prior.index] - prior.mean) / sigma
        total += LOG2PI + log(abs2(sigma)) + abs2(z)
    end
    return total
end

function _parameter_constraint_nll(problem, p::AbstractVector)
    isempty(problem.parameter_constraints) && return zero(eltype(p))

    return _parameter_constraint_nll(_prepare_parameter_constraints(problem), p)
end

function _parameter_constraint_nll(prepared_constraints::AbstractVector{<:PreparedParameterConstraint}, p::AbstractVector)
    isempty(prepared_constraints) && return zero(eltype(p))

    total = zero(eltype(p))
    @inbounds for constraint in prepared_constraints
        k = length(constraint.indices)
        delta = p[constraint.indices] .- constraint.mean
        z = constraint.factor.L \ delta
        total += k * LOG2PI + constraint.logdet + sum(abs2, z)
    end
    return total
end

function _data_chi2(problem::FitProblem, p::AbstractVector)
    rw = _weighted_data_residual(problem, p)
    return sum(abs2, rw)
end

function _data_chi2(cache::FitEvaluationCache, p::AbstractVector)
    residual = _residual(cache.problem, p)
    rw = _whiten_residual(cache, p, residual)
    return sum(abs2, rw)
end

function _chi2_cost(problem::FitProblem, p::AbstractVector)
    return _data_chi2(problem, p) + _prior_chi2(problem, p) + _parameter_constraint_chi2(problem, p)
end

function _chi2_cost(cache::FitEvaluationCache, p::AbstractVector)
    problem = cache.problem
    return _data_chi2(cache, p) + _prior_chi2(problem, p) + _parameter_constraint_chi2(cache.parameter_constraints, p)
end

function _covariance_logdet(cov, n::Int)
    if cov === nothing
        return zero(Float64)
    elseif cov isa AbstractVector
        cov_values = _finite_value.(cov)
        any(cov_values .<= 0.0) && throw(ArgumentError("all effective variances must be positive"))
        return sum(log, cov)
    end

    F = _stable_cholesky(cov)
    return _cholesky_logdet(F)
end

function _gaussian_data_nll(problem::FitProblem, p::AbstractVector)
    n = length(problem.y)
    if problem.whitening !== nothing
        return n * LOG2PI + problem.whitening.logdet_covariance + _data_chi2(problem, p)
    end
    cov = _effective_covariance(problem, p)
    return n * LOG2PI + _covariance_logdet(cov, n) + _data_chi2(problem, p)
end

function _gaussian_data_nll(cache::FitEvaluationCache{NoPreparedCovariance}, p::AbstractVector)
    return length(cache.problem.y) * LOG2PI + _data_chi2(cache, p)
end

function _gaussian_data_nll(cache::FitEvaluationCache{<:Union{DiagonalPreparedCovariance, DensePreparedCovariance, OperatorPreparedCovariance}}, p::AbstractVector)
    return length(cache.problem.y) * LOG2PI + cache.covariance.logdet + _data_chi2(cache, p)
end

function _gaussian_data_nll(cache::FitEvaluationCache{DynamicPreparedCovariance}, p::AbstractVector)
    return _gaussian_data_nll(cache.problem, p)
end

function _gaussian_nll(problem::FitProblem, p::AbstractVector)
    return _gaussian_data_nll(problem, p) + _prior_nll(problem, p) + _parameter_constraint_nll(problem, p)
end

function _gaussian_nll(cache::FitEvaluationCache, p::AbstractVector)
    problem = cache.problem
    return _gaussian_data_nll(cache, p) + _prior_nll(problem, p) + _parameter_constraint_nll(cache.parameter_constraints, p)
end

function _cost_value(problem::FitProblem, p::AbstractVector, cost::Symbol)
    resolved = _resolve_cost(problem, cost)
    if resolved == :chi2
        return _chi2_cost(problem, p)
    elseif resolved == :gaussian_nll
        return _gaussian_nll(problem, p)
    end
    throw(ArgumentError("unsupported resolved cost: $resolved"))
end

function _cost_value(cache::FitEvaluationCache, p::AbstractVector, cost::Symbol)
    resolved = _resolve_cost(cache.problem, cost)
    if resolved == :chi2
        return _chi2_cost(cache, p)
    elseif resolved == :gaussian_nll
        return _gaussian_nll(cache, p)
    end
    throw(ArgumentError("unsupported resolved cost: $resolved"))
end

function _cost_uses_parameter_dependent_covariance(problem::FitProblem, cost::Symbol)
    return _resolve_cost(problem, cost) == :gaussian_nll && _has_parameter_dependent_covariance(problem)
end

function _has_parameter_dependent_covariance(problem::FitProblem)
    has_x_uncertainty(problem) && return true
    return any(c -> c.active && c.target == :y && c.mode == :model_relative, problem.error_components)
end
