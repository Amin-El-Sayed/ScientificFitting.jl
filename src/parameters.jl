function _fixed_index_set(problem)
    isempty(problem.fixed_parameters) && return Set{Int}()
    return Set(fp.index for fp in problem.fixed_parameters)
end

function _free_indices(problem)
    isempty(problem.fixed_parameters) && return collect(eachindex(problem.p0))
    fixed = _fixed_index_set(problem)
    return [i for i in eachindex(problem.p0) if !(i in fixed)]
end

function _fixed_lookup(problem)
    isempty(problem.fixed_parameters) && return Dict{Int, FixedParameter}()
    return Dict(fp.index => fp for fp in problem.fixed_parameters)
end

function _free_p0(problem)
    isempty(problem.fixed_parameters) && return copy(problem.p0)
    return problem.p0[_free_indices(problem)]
end

function _expand_free_parameters(problem, free_params::AbstractVector)
    if isempty(problem.fixed_parameters)
        length(free_params) == length(problem.p0) || throw(ArgumentError("free parameter vector has wrong length"))
        return collect(free_params)
    end

    free_idx = _free_indices(problem)
    length(free_params) == length(free_idx) || throw(ArgumentError("free parameter vector has wrong length"))

    T = promote_type(eltype(problem.p0), eltype(free_params))
    full = Vector{T}(undef, length(problem.p0))
    full .= problem.p0
    fixed = _fixed_lookup(problem)
    for (idx, fp) in fixed
        full[idx] = fp.value
    end
    full[free_idx] .= free_params
    return full
end

function _free_bounds(problem)
    problem.bounds === nothing && return nothing
    free_idx = _free_indices(problem)
    lower, upper = problem.bounds
    return (lower[free_idx], upper[free_idx])
end

function _free_constraints(spec::ConstraintSpec, problem)
    has_constraints(spec) || return spec

    ineq = spec.ineq === nothing ? nothing : q -> spec.ineq(_expand_free_parameters(problem, q))
    eq = spec.eq === nothing ? nothing : q -> spec.eq(_expand_free_parameters(problem, q))
    return ConstraintSpec(; ineq=ineq, eq=eq)
end

function _embed_free_covariance(problem, free_cov::AbstractMatrix)
    n = length(problem.p0)
    free_idx = _free_indices(problem)
    size(free_cov) == (length(free_idx), length(free_idx)) || throw(ArgumentError("free covariance has wrong shape"))

    cov = zeros(Float64, n, n)
    cov[free_idx, free_idx] .= free_cov

    @inbounds for fp in problem.fixed_parameters
        # A single local covariance cannot represent asymmetry; use the larger side conservatively.
        cov[fp.index, fp.index] = max(fp.sigma_minus, fp.sigma_plus)^2
    end

    return cov
end

function _free_jacobian_from_full(problem::FitProblem, full_jacobian::AbstractMatrix)
    return Matrix{Float64}(full_jacobian[:, _free_indices(problem)])
end

function _full_jacobian_from_free(problem, free_jacobian::AbstractMatrix, nrows::Int)
    free_idx = _free_indices(problem)
    size(free_jacobian) == (nrows, length(free_idx)) || throw(ArgumentError("free jacobian has wrong shape"))
    jacobian = zeros(Float64, nrows, length(problem.p0))
    jacobian[:, free_idx] .= free_jacobian
    return jacobian
end

function _free_weighted_jacobian(problem::FitProblem, params::AbstractVector)
    free_idx = _free_indices(problem)
    jac = ForwardDiff.jacobian(q -> _weighted_residual(problem, _expand_free_parameters(problem, q)), params[free_idx])
    return Matrix{Float64}(jac)
end

function _with_p0(problem::FitProblem, p0::AbstractVector)
    return FitProblem(
        problem.model,
        problem.x,
        problem.y;
        p0=p0,
        sigma_y=problem.sigma_y,
        sigma_x=problem.sigma_x,
        cov_y=problem.cov_y,
        cov_x=problem.cov_x,
        error_components=problem.error_components,
        bounds=problem.bounds,
        constraints=problem.constraints,
        parameter_priors=problem.parameter_priors,
        parameter_constraints=problem.parameter_constraints,
        fixed_parameters=problem.fixed_parameters,
        jacobian=problem.jacobian,
    )
end

function _clip_to_bounds(p::Vector{Float64}, bounds)
    bounds === nothing && return p
    lower, upper = bounds
    return clamp.(p, lower, upper)
end

function _initial_candidates(problem, initial_guesses, multistart::Int)
    candidates = Vector{Vector{Float64}}()
    push!(candidates, copy(problem.p0))

    if initial_guesses !== nothing
        raw = initial_guesses isa AbstractVector && !isempty(initial_guesses) && first(initial_guesses) isa Number ? [initial_guesses] : initial_guesses
        for guess in raw
            push!(candidates, _float_vector(guess))
        end
    end

    if multistart > 1
        free_idx = _free_indices(problem)
        if problem.bounds !== nothing
            lower, upper = problem.bounds
            finite = isfinite.(lower) .& isfinite.(upper)
            center = copy(problem.p0)
            center[finite] .= (lower[finite] .+ upper[finite]) ./ 2
            push!(candidates, center)
            for idx in free_idx
                finite[idx] || continue
                low = copy(center)
                high = copy(center)
                low[idx] = lower[idx] + 0.25 * (upper[idx] - lower[idx])
                high[idx] = lower[idx] + 0.75 * (upper[idx] - lower[idx])
                push!(candidates, low)
                push!(candidates, high)
                length(candidates) >= multistart && break
            end
        else
            for scale in (0.5, 2.0, -1.0)
                candidate = copy(problem.p0)
                candidate[free_idx] .= scale .* candidate[free_idx]
                push!(candidates, candidate)
                length(candidates) >= multistart && break
            end
        end
    end

    fixed = _fixed_lookup(problem)
    for candidate in candidates
        length(candidate) == length(problem.p0) || throw(ArgumentError("each initial guess must match parameter count"))
        for (idx, fp) in fixed
            candidate[idx] = fp.value
        end
    end

    unique_candidates = Vector{Vector{Float64}}()
    seen = Set{Tuple{Vararg{Float64}}}()
    for candidate in candidates
        clipped = _clip_to_bounds(candidate, problem.bounds)
        key = Tuple(clipped)
        if !(key in seen)
            push!(seen, key)
            push!(unique_candidates, clipped)
        end
        length(unique_candidates) >= multistart && break
    end

    return unique_candidates
end
