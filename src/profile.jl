struct ProfileResult
    parameter_index::Int
    values::Vector{Float64}
    cost_values::Vector{Float64}
    delta_cost::Vector{Float64}
    threshold::Float64
    best_value::Float64
end

struct ContourResult
    parameter_indices::Tuple{Int, Int}
    x_values::Vector{Float64}
    y_values::Vector{Float64}
    cost_values::Matrix{Float64}
    delta_cost::Matrix{Float64}
    levels::Vector{Float64}
end

struct ProfileInterval
    parameter_index::Int
    lower::Float64
    upper::Float64
    uncertainty_minus::Float64
    uncertainty_plus::Float64
    threshold::Float64
    profile_result::ProfileResult
end

function _merge_fixed_parameters(existing::Vector{FixedParameter}, added::Vector{FixedParameter})
    merged = Dict(fp.index => fp for fp in existing)
    for fp in added
        merged[fp.index] = fp
    end
    return [merged[i] for i in sort!(collect(keys(merged)))]
end

function _refit_with_fixed(result::FitResult, fixed::Vector{FixedParameter})
    problem = result.problem
    p0 = copy(result.params)
    for fp in fixed
        p0[fp.index] = fp.value
    end

    refit_problem = FitProblem(
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
        fixed_parameters=_merge_fixed_parameters(problem.fixed_parameters, fixed),
        jacobian=problem.jacobian,
    )

    return fit(
        refit_problem;
        backend=:auto,
        cost=result.options.cost,
        maxiters=result.options.maxiters,
        tol=result.options.tol,
        ci_level=result.options.ci_level,
        scale_covariance=result.options.scale_covariance,
    )
end

function _refit_with_fixed(result::LikelihoodFitResult, fixed::Vector{FixedParameter})
    problem = result.problem
    p0 = copy(result.params)
    for fp in fixed
        p0[fp.index] = fp.value
    end

    refit_problem = LikelihoodFitProblem(
        problem.objective,
        problem.gof,
        p0;
        bounds=problem.bounds,
        constraints=problem.constraints,
        parameter_priors=problem.parameter_priors,
        parameter_constraints=problem.parameter_constraints,
        fixed_parameters=_merge_fixed_parameters(problem.fixed_parameters, fixed),
        nobs=problem.nobs,
        cost_name=problem.cost_name,
        parameter_names=problem.parameter_names,
    )

    return fit(
        refit_problem;
        maxiters=result.options.maxiters,
        tol=result.options.tol,
        ci_level=result.options.ci_level,
    )
end

function _default_profile_grid(result, index::Int; npoints::Int, nsigma::Real)
    center = result.params[index]
    sigma = result.param_stderr[index]
    if !isfinite(sigma) || sigma <= 0
        sigma = max(abs(center), 1.0) * 0.1
    end
    return collect(range(center - nsigma * sigma, center + nsigma * sigma; length=npoints))
end

"""
    profile(result, index; values=nothing, npoints=61, nsigma=3, threshold=1.0)

Profile the fitted cost function in one parameter by fixing that parameter to
grid values and re-minimizing all remaining free parameters.
"""
function profile(
    result,
    index::Int;
    values=nothing,
    npoints::Int=61,
    nsigma::Real=3,
    threshold::Real=1.0,
)
    1 <= index <= length(result.params) || throw(ArgumentError("profile index out of range"))
    grid = values === nothing ? _default_profile_grid(result, index; npoints=npoints, nsigma=nsigma) : collect(Float64, values)
    costs = Vector{Float64}(undef, length(grid))

    for (i, value) in enumerate(grid)
        profiled = _refit_with_fixed(result, [FixedParameter(index, value)])
        costs[i] = profiled.stats.cost_min
    end

    delta = costs .- result.stats.cost_min
    return ProfileResult(index, grid, costs, delta, Float64(threshold), result.params[index])
end

function _linear_crossing(x1, y1, x2, y2, threshold)
    y2 == y1 && return (x1 + x2) / 2
    t = (threshold - y1) / (y2 - y1)
    return x1 + t * (x2 - x1)
end

function _profile_crossings(profile_result::ProfileResult)
    values = profile_result.values
    delta = profile_result.delta_cost
    threshold = profile_result.threshold
    center = argmin(delta)

    lower = NaN
    for i in (center - 1):-1:1
        if delta[i] >= threshold && delta[i + 1] <= threshold
            lower = _linear_crossing(values[i], delta[i], values[i + 1], delta[i + 1], threshold)
            break
        end
    end

    upper = NaN
    for i in center:(length(values) - 1)
        if delta[i] <= threshold && delta[i + 1] >= threshold
            upper = _linear_crossing(values[i], delta[i], values[i + 1], delta[i + 1], threshold)
            break
        end
    end

    return lower, upper
end

"""
    profile_interval(result, index; threshold=1.0, npoints=121, nsigma=5)

Compute a profile-based asymmetric interval by finding the profile-cost
crossings at `delta_cost = threshold`.
"""
function profile_interval(
    result,
    index::Int;
    threshold::Real=1.0,
    npoints::Int=121,
    nsigma::Real=5,
    values=nothing,
)
    prof = profile(result, index; values=values, npoints=npoints, nsigma=nsigma, threshold=threshold)
    lower, upper = _profile_crossings(prof)
    center = result.params[index]
    minus = isfinite(lower) ? center - lower : NaN
    plus = isfinite(upper) ? upper - center : NaN
    return ProfileInterval(index, lower, upper, minus, plus, Float64(threshold), prof)
end

function _default_contour_grid(result::FitResult, index::Int; npoints::Int, nsigma::Real)
    return _default_profile_grid(result, index; npoints=npoints, nsigma=nsigma)
end

"""
    contour(result, i, j; xvalues=nothing, yvalues=nothing, npoints=31, nsigma=3, levels=[2.30, 6.18])

Compute a two-parameter profile-likelihood contour grid. At each grid point,
parameters `i` and `j` are fixed and all remaining free parameters are
re-minimized.
"""
function contour(
    result,
    i::Int,
    j::Int;
    xvalues=nothing,
    yvalues=nothing,
    npoints::Int=31,
    nsigma::Real=3,
    levels::AbstractVector=[2.30, 6.18],
)
    i != j || throw(ArgumentError("contour requires two distinct parameter indices"))
    1 <= i <= length(result.params) || throw(ArgumentError("first contour index out of range"))
    1 <= j <= length(result.params) || throw(ArgumentError("second contour index out of range"))

    xs = xvalues === nothing ? _default_contour_grid(result, i; npoints=npoints, nsigma=nsigma) : collect(Float64, xvalues)
    ys = yvalues === nothing ? _default_contour_grid(result, j; npoints=npoints, nsigma=nsigma) : collect(Float64, yvalues)
    costs = Matrix{Float64}(undef, length(xs), length(ys))

    for ix in eachindex(xs), iy in eachindex(ys)
        profiled = _refit_with_fixed(result, [FixedParameter(i, xs[ix]), FixedParameter(j, ys[iy])])
        costs[ix, iy] = profiled.stats.cost_min
    end

    delta = costs .- result.stats.cost_min
    return ContourResult((i, j), xs, ys, costs, delta, collect(Float64, levels))
end
