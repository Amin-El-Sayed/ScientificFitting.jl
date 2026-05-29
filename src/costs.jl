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

    total = zero(eltype(p))
    @inbounds for constraint in problem.parameter_constraints
        delta = p[constraint.indices] .- constraint.mean
        z = _stable_cholesky(constraint.covariance).L \ delta
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

    total = zero(eltype(p))
    @inbounds for constraint in problem.parameter_constraints
        k = length(constraint.indices)
        delta = p[constraint.indices] .- constraint.mean
        F = _stable_cholesky(constraint.covariance)
        z = F.L \ delta
        logdet = 2.0 * sum(log, diag(F.L))
        total += k * LOG2PI + logdet + sum(abs2, z)
    end
    return total
end

function _data_chi2(problem::FitProblem, p::AbstractVector)
    rw = _weighted_data_residual(problem, p)
    return sum(abs2, rw)
end

function _chi2_cost(problem::FitProblem, p::AbstractVector)
    return _data_chi2(problem, p) + _prior_chi2(problem, p) + _parameter_constraint_chi2(problem, p)
end

function _covariance_logdet(cov, n::Int)
    if cov === nothing
        return zero(Float64)
    elseif cov isa AbstractVector
        cov_values = ForwardDiff.value.(cov)
        any(cov_values .<= 0.0) && throw(ArgumentError("all effective variances must be positive"))
        return sum(log, cov)
    end

    F = _stable_cholesky(cov)
    return 2.0 * sum(log, diag(F.L))
end

function _gaussian_data_nll(problem::FitProblem, p::AbstractVector)
    n = length(problem.y)
    cov = _effective_covariance(problem, p)
    return n * LOG2PI + _covariance_logdet(cov, n) + _data_chi2(problem, p)
end

function _gaussian_nll(problem::FitProblem, p::AbstractVector)
    return _gaussian_data_nll(problem, p) + _prior_nll(problem, p) + _parameter_constraint_nll(problem, p)
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

function _cost_uses_parameter_dependent_covariance(problem::FitProblem, cost::Symbol)
    return _resolve_cost(problem, cost) == :gaussian_nll && _has_parameter_dependent_covariance(problem)
end

function _has_parameter_dependent_covariance(problem::FitProblem)
    has_x_uncertainty(problem) && return true
    return any(c -> c.active && c.target == :y && c.mode == :model_relative, problem.error_components)
end
