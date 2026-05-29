function _safe_condition_number(matrix)
    isempty(matrix) && return NaN
    values = svdvals(Matrix{Float64}(matrix))
    isempty(values) && return NaN
    smallest = minimum(values)
    return smallest > 0 ? maximum(values) / smallest : Inf
end

function _stable_symmetric_inverse(matrix; rtol::Real=1e-12)
    sym = Symmetric(Matrix{Float64}(matrix))
    try
        return Matrix{Float64}(inv(cholesky(sym)))
    catch
        eigen_decomp = eigen(sym)
        values = eigen_decomp.values
        vectors = eigen_decomp.vectors
        scale = maximum(abs, values; init=0.0)
        cutoff = max(scale * rtol, eps(Float64))
        inv_values = [abs(v) > cutoff ? inv(v) : 0.0 for v in values]
        return Matrix{Float64}(vectors * Diagonal(inv_values) * vectors')
    end
end

function _active_bound_indices(bounds, params::AbstractVector; atol::Real=1e-8)
    bounds === nothing && return Int[]
    lower, upper = bounds
    active = Int[]
    for i in eachindex(params)
        if isfinite(lower[i]) && abs(params[i] - lower[i]) <= atol * max(1.0, abs(lower[i]))
            push!(active, i)
        elseif isfinite(upper[i]) && abs(params[i] - upper[i]) <= atol * max(1.0, abs(upper[i]))
            push!(active, i)
        end
    end
    return active
end

function _diagnostic_warnings(converged::Bool, ndf::Int, cov_cond::Float64, hess_cond::Float64, active_bounds::Vector{Int})
    warnings = String[]
    converged || push!(warnings, "optimizer did not report convergence")
    ndf <= 0 && push!(warnings, "non-positive degrees of freedom; p-values and reduced statistics are not meaningful")
    isfinite(cov_cond) && cov_cond > 1e12 && push!(warnings, "parameter covariance is ill-conditioned")
    isfinite(hess_cond) && hess_cond > 1e12 && push!(warnings, "cost Hessian is ill-conditioned")
    !isempty(active_bounds) && push!(warnings, "one or more parameters are at active bounds; local errors and p-values may be unreliable")
    return warnings
end

function _fit_diagnostics(problem, params::AbstractVector, cov::AbstractMatrix, converged::Bool, ndf::Int; hessian=nothing)
    cov_cond = _safe_condition_number(cov)
    hess_cond = hessian === nothing ? NaN : _safe_condition_number(hessian)
    active_bounds = _active_bound_indices(problem.bounds, params)
    warnings = _diagnostic_warnings(converged, ndf, cov_cond, hess_cond, active_bounds)
    return FitDiagnostics(warnings, cov_cond, hess_cond, active_bounds)
end
