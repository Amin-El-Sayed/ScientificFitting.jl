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

function _profile_refit_cost(result, fixed::Vector{FixedParameter}; on_failure::Symbol)
    try
        profiled = _refit_with_fixed(result, fixed)
        return Float64(profiled.stats.cost_min)
    catch err
        on_failure == :throw && rethrow(err)
        on_failure == :inf && return Inf
        throw(ArgumentError("on_failure must be :inf or :throw"))
    end
end

"""
    profile(result, index; values=nothing, npoints=61, nsigma=3, threshold=1.0, on_failure=:inf)

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
    on_failure::Symbol=:inf,
)
    1 <= index <= length(result.params) || throw(ArgumentError("profile index out of range"))
    grid = values === nothing ? _default_profile_grid(result, index; npoints=npoints, nsigma=nsigma) : collect(Float64, values)
    costs = Vector{Float64}(undef, length(grid))

    for (i, value) in enumerate(grid)
        costs[i] = _profile_refit_cost(result, [FixedParameter(index, value)]; on_failure=on_failure)
    end

    delta = costs .- result.stats.cost_min
    return ProfileResult(index, grid, costs, delta, Float64(threshold), result.params[index])
end

function _profile_refit_failure_findings(profile_result::ProfileResult)
    failed = count(!isfinite, profile_result.cost_values)
    failed == 0 && return DiagnosticFinding[]
    return DiagnosticFinding[
        _finding(
            :warning,
            :profile_refit_failed,
            "Some profile refits failed",
            "$failed of $(length(profile_result.cost_values)) profile grid point(s) produced non-finite costs.",
            "Inspect bounds, constraints, starting values, and scan range. Treat intervals across failed regions as unreliable.",
        ),
    ]
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

function _profile_threshold_bracket_findings(profile_result::ProfileResult)
    lower, upper = _profile_crossings(profile_result)
    findings = DiagnosticFinding[]
    if !isfinite(lower) || !isfinite(upper)
        side = !isfinite(lower) && !isfinite(upper) ? "both sides" : (!isfinite(lower) ? "lower side" : "upper side")
        push!(
            findings,
            _finding(
                :warning,
                :profile_threshold_not_bracketed,
                "Profile interval is not fully bracketed",
                "The threshold $(profile_result.threshold) was not crossed on the $side of the scan.",
                "Extend the profile range with a larger nsigma or explicit values before quoting a profile interval.",
            ),
        )
    end
    return findings
end

function _profile_parabolicity_findings(profile_result::ProfileResult, local_sigma; tolerance::Real)
    local_sigma === nothing && return DiagnosticFinding[]
    sigma = Float64(local_sigma)
    sigma > 0 || throw(ArgumentError("local_sigma must be positive"))
    tolerance >= 0 || throw(ArgumentError("tolerance must be non-negative"))

    local_delta = @. abs2((profile_result.values - profile_result.best_value) / sigma)
    relevant = (profile_result.delta_cost .<= max(4 * profile_result.threshold, profile_result.threshold + 3))
    any(relevant) || return DiagnosticFinding[]
    deviation = maximum(abs.(profile_result.delta_cost[relevant] .- local_delta[relevant]); init=0.0)

    if deviation > tolerance
        return DiagnosticFinding[
            _finding(
                :warning,
                :profile_not_parabolic,
                "Profile is not well described by the local parabola",
                "max |profile - local parabola| = $(_fmt_scientific(deviation)) near the minimum.",
                "Use profile intervals instead of symmetric local errors, and inspect bounds, correlations, scaling, or model nonlinearity.",
            ),
        ]
    end
    return DiagnosticFinding[]
end

"""
    diagnose(profile_result::ProfileResult; local_sigma=nothing, tolerance=0.25)

Diagnose an already computed one-parameter profile. With `local_sigma`, the
actual profile is compared to the local covariance parabola. This is the
machine-readable counterpart of overlaying both curves in `plot_profile`.
"""
function diagnose(profile_result::ProfileResult; local_sigma=nothing, tolerance::Real=0.25)
    findings = DiagnosticFinding[]
    append!(findings, _profile_refit_failure_findings(profile_result))
    append!(findings, _profile_threshold_bracket_findings(profile_result))
    append!(findings, _profile_parabolicity_findings(profile_result, local_sigma; tolerance=tolerance))
    findings = _sort_findings(_deduplicate_findings(findings))
    return DiagnosticReport(findings, _diagnostic_summary(findings))
end

function _default_contour_grid(result::FitResult, index::Int; npoints::Int, nsigma::Real)
    return _default_profile_grid(result, index; npoints=npoints, nsigma=nsigma)
end

"""
    contour(result, i, j; xvalues=nothing, yvalues=nothing, npoints=31, nsigma=3, levels=[2.30, 6.18], on_failure=:inf)

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
    on_failure::Symbol=:inf,
)
    i != j || throw(ArgumentError("contour requires two distinct parameter indices"))
    1 <= i <= length(result.params) || throw(ArgumentError("first contour index out of range"))
    1 <= j <= length(result.params) || throw(ArgumentError("second contour index out of range"))

    xs = xvalues === nothing ? _default_contour_grid(result, i; npoints=npoints, nsigma=nsigma) : collect(Float64, xvalues)
    ys = yvalues === nothing ? _default_contour_grid(result, j; npoints=npoints, nsigma=nsigma) : collect(Float64, yvalues)
    costs = Matrix{Float64}(undef, length(xs), length(ys))

    for ix in eachindex(xs), iy in eachindex(ys)
        costs[ix, iy] = _profile_refit_cost(result, [FixedParameter(i, xs[ix]), FixedParameter(j, ys[iy])]; on_failure=on_failure)
    end

    delta = costs .- result.stats.cost_min
    return ContourResult((i, j), xs, ys, costs, delta, collect(Float64, levels))
end

function _contour_refit_failure_findings(contour_result::ContourResult)
    failed = count(!isfinite, contour_result.cost_values)
    failed == 0 && return DiagnosticFinding[]
    return DiagnosticFinding[
        _finding(
            :warning,
            :contour_refit_failed,
            "Some contour refits failed",
            "$failed of $(length(contour_result.cost_values)) contour grid point(s) produced non-finite costs.",
            "Inspect bounds, constraints, starting values, and scan range. Do not interpret contour topology through failed grid regions.",
        ),
    ]
end

function _contour_level_bracket_findings(contour_result::ContourResult)
    findings = DiagnosticFinding[]
    finite_delta = contour_result.delta_cost[isfinite.(contour_result.delta_cost)]
    isempty(finite_delta) && return findings
    min_delta = minimum(finite_delta)
    max_delta = maximum(finite_delta)
    missing_levels = [level for level in contour_result.levels if !(min_delta <= level <= max_delta)]
    if !isempty(missing_levels)
        push!(
            findings,
            _finding(
                :warning,
                :contour_levels_not_bracketed,
                "Contour scan does not cover all requested levels",
                "Requested level(s) $(join(_fmt_scientific.(missing_levels), ", ")) are outside the scanned delta-cost range [$(_fmt_scientific(min_delta)), $(_fmt_scientific(max_delta))].",
                "Increase npoints/nsigma or pass wider xvalues/yvalues before interpreting missing contour levels.",
            ),
        )
    end
    return findings
end

function _contour_center(contour_result::ContourResult, local_center)
    if local_center !== nothing
        raw = collect(local_center)
        length(raw) == 2 || throw(ArgumentError("local_center must contain exactly two values"))
        return (Float64(raw[1]), Float64(raw[2]))
    end

    finite_delta = replace(contour_result.delta_cost, NaN => Inf)
    idx = argmin(vec(finite_delta))
    ix, iy = Tuple(CartesianIndices(contour_result.delta_cost)[idx])
    return (contour_result.x_values[ix], contour_result.y_values[iy])
end

function _contour_local_covariance(contour_result::ContourResult, local_covariance)
    local_covariance === nothing && return nothing
    cov = Matrix{Float64}(local_covariance)
    if size(cov) == (2, 2)
        return cov
    end
    i, j = contour_result.parameter_indices
    size(cov, 1) >= max(i, j) && size(cov, 2) >= max(i, j) ||
        throw(ArgumentError("local_covariance must be 2x2 or the full parameter covariance matrix"))
    return cov[[i, j], [i, j]]
end

function _contour_ellipticity_findings(contour_result::ContourResult, local_covariance, local_center; tolerance::Real)
    cov = _contour_local_covariance(contour_result, local_covariance)
    cov === nothing && return DiagnosticFinding[]
    tolerance >= 0 || throw(ArgumentError("tolerance must be non-negative"))

    center = _contour_center(contour_result, local_center)
    precision = Symmetric(cov) \ Matrix{Float64}(I, 2, 2)
    local_delta = Matrix{Float64}(undef, length(contour_result.x_values), length(contour_result.y_values))
    for ix in eachindex(contour_result.x_values), iy in eachindex(contour_result.y_values)
        delta = [contour_result.x_values[ix] - center[1], contour_result.y_values[iy] - center[2]]
        local_delta[ix, iy] = dot(delta, precision * delta)
    end

    first_level = isempty(contour_result.levels) ? 2.30 : minimum(contour_result.levels)
    relevant = contour_result.delta_cost .<= max(4 * first_level, first_level + 3)
    any(relevant) || return DiagnosticFinding[]
    deviation = maximum(abs.(contour_result.delta_cost[relevant] .- local_delta[relevant]); init=0.0)

    if deviation > tolerance
        return DiagnosticFinding[
            _finding(
                :warning,
                :contour_not_elliptic,
                "Contour is not well described by the local covariance ellipse",
                "max |profile contour - local ellipse| = $(_fmt_scientific(deviation)) near the minimum.",
                "Use profile contours for uncertainty interpretation, and inspect parameter correlations, bounds, scaling, or model nonlinearity.",
            ),
        ]
    end
    return DiagnosticFinding[]
end

"""
    diagnose(contour_result::ContourResult; local_covariance=nothing, local_center=nothing, tolerance=0.5)

Diagnose an already computed two-parameter contour grid. With
`local_covariance`, the actual profiled contour surface is compared to the local
covariance ellipse used by symmetric Gaussian error propagation.
"""
function diagnose(contour_result::ContourResult; local_covariance=nothing, local_center=nothing, tolerance::Real=0.5)
    findings = DiagnosticFinding[]
    append!(findings, _contour_refit_failure_findings(contour_result))
    append!(findings, _contour_level_bracket_findings(contour_result))
    append!(
        findings,
        _contour_ellipticity_findings(
            contour_result,
            local_covariance,
            local_center;
            tolerance=tolerance,
        ),
    )
    findings = _sort_findings(_deduplicate_findings(findings))
    return DiagnosticReport(findings, _diagnostic_summary(findings))
end
