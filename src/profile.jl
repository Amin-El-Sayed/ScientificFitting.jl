"""
    ProfileResult

One-dimensional profile scan of the fitted cost function. `values` are the
fixed parameter values, `cost_values` are the refitted objective values, and
`delta_cost` is measured relative to the original fit minimum. The `threshold`
field records the interval threshold requested by `profile`.
"""
struct ProfileResult
    parameter_index::Int
    values::Vector{Float64}
    cost_values::Vector{Float64}
    delta_cost::Vector{Float64}
    threshold::Float64
    best_value::Float64
end

"""
    ContourResult

Two-parameter profile-contour scan. `x_values` and `y_values` define the scan
grid for `parameter_indices`; `delta_cost` stores the refitted cost increase
relative to the best fit, and `levels` stores the requested contour thresholds.
"""
struct ContourResult
    parameter_indices::Tuple{Int, Int}
    x_values::Vector{Float64}
    y_values::Vector{Float64}
    cost_values::Matrix{Float64}
    delta_cost::Matrix{Float64}
    levels::Vector{Float64}
end

"""
    ProfileInterval

Asymmetric interval extracted from a `ProfileResult`. `lower` and `upper` are
the threshold crossings; `uncertainty_minus` and `uncertainty_plus` are measured
relative to the best-fit parameter value. Non-finite fields indicate that the
scan did not bracket the requested threshold.
"""
struct ProfileInterval
    parameter_index::Int
    lower::Float64
    upper::Float64
    uncertainty_minus::Float64
    uncertainty_plus::Float64
    threshold::Float64
    profile_result::ProfileResult
end

"""
    ProfileMatrixResult

Makie-free diagnostic overview for several fitted parameters. It stores the
one-parameter profiles, lower-triangle pairwise contours, per-panel diagnostic
reports, the selected local fit geometry, and a combined diagnostic report.
Plotting is intentionally separate: `plot_profile_matrix(matrix_result)`
renders this object without repeating any profile or contour refits.
`panel_status` maps diagonal `(i, i)` and lower-triangle `(i, j)` panels to
`:ok`, `:review`, or `:stop`.
"""
struct ProfileMatrixResult
    parameters::Vector{Int}
    parameter_names::Vector{String}
    best_values::Vector{Float64}
    local_stderr::Vector{Float64}
    local_covariance::Matrix{Float64}
    local_correlation::Matrix{Float64}
    profiles::Dict{Int, ProfileResult}
    contours::Dict{Tuple{Int, Int}, ContourResult}
    profile_diagnostics::Dict{Int, DiagnosticReport}
    contour_diagnostics::Dict{Tuple{Int, Int}, DiagnosticReport}
    panel_status::Dict{Tuple{Int, Int}, Symbol}
    report::DiagnosticReport
end

"""
    ProfileMatrixPanelTriage

One row in the Makie-free triage summary returned by
`profile_matrix_triage(...)`. `indices` identify the profile or contour panel,
`parameter_names` are the reader-facing parameter labels, `status` is `:ok`,
`:review`, or `:stop`, `finding_codes` are stable machine-readable diagnostic
codes, and `next_action` is the first recommended action for that panel.
"""
struct ProfileMatrixPanelTriage
    indices::Tuple{Int, Int}
    parameter_names::Tuple{String, String}
    status::Symbol
    severity_counts::Dict{Symbol, Int}
    finding_codes::Vector{Symbol}
    next_action::String
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
        whitening=problem.whitening,
        error_components=problem.error_components,
        bounds=problem.bounds,
        constraints=problem.constraints,
        parameter_priors=problem.parameter_priors,
        parameter_constraints=problem.parameter_constraints,
        fixed_parameters=_merge_fixed_parameters(problem.fixed_parameters, fixed),
        jacobian=problem.jacobian,
        x_derivative=problem.x_derivative,
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

function _validate_adaptive_controls(max_refinements::Int, max_points::Int)
    max_refinements >= 0 || throw(ArgumentError("max_refinements must be non-negative"))
    max_points >= 3 || throw(ArgumentError("max_points must be at least 3"))
    return nothing
end

function _validate_profile_controls(npoints::Int, nsigma::Real, threshold::Real; default_grid::Bool)
    if default_grid
        npoints >= 3 || throw(ArgumentError("npoints must be at least 3"))
        isfinite(nsigma) && nsigma > 0 || throw(ArgumentError("nsigma must be finite and > 0"))
    end
    isfinite(threshold) && threshold > 0 || throw(ArgumentError("threshold must be finite and > 0"))
    return nothing
end

function _validate_contour_controls(npoints::Int, nsigma::Real; default_x_grid::Bool, default_y_grid::Bool)
    if default_x_grid || default_y_grid
        npoints >= 2 || throw(ArgumentError("npoints must be at least 2"))
        isfinite(nsigma) && nsigma > 0 || throw(ArgumentError("nsigma must be finite and > 0"))
    end
    return nothing
end

function _unique_sorted(values::AbstractVector{<:Real})
    cleaned = sort!(collect(Float64, values))
    isempty(cleaned) && return cleaned
    out = Float64[cleaned[1]]
    for value in Iterators.drop(cleaned, 1)
        if value != last(out)
            push!(out, value)
        end
    end
    return out
end

function _validated_scan_values(values, name::AbstractString; min_points::Int)
    cleaned = _unique_sorted(collect(Float64, values))
    length(cleaned) >= min_points || throw(ArgumentError("$name must contain at least $min_points distinct values"))
    all(isfinite, cleaned) || throw(ArgumentError("$name must contain only finite values"))
    return cleaned
end

function _profile_refinement_candidates(values, delta, thresholds)
    candidates = Float64[]
    for i in 1:(length(values) - 1)
        d1 = delta[i]
        d2 = delta[i + 1]
        isfinite(d1) && isfinite(d2) || continue
        for threshold in thresholds
            if min(d1, d2) <= threshold <= max(d1, d2)
                push!(candidates, (values[i] + values[i + 1]) / 2)
                break
            end
        end
    end
    return candidates
end

function _profile_from_grid(result, index::Int, grid::Vector{Float64}, threshold::Float64, on_failure::Symbol)
    costs = Vector{Float64}(undef, length(grid))
    for (i, value) in enumerate(grid)
        costs[i] = _profile_refit_cost(result, [FixedParameter(index, value)]; on_failure=on_failure)
    end
    delta = costs .- result.stats.cost_min
    return ProfileResult(index, grid, costs, delta, threshold, result.params[index])
end

function _adaptive_profile(
    result,
    index::Int,
    grid::Vector{Float64},
    threshold::Float64,
    on_failure::Symbol;
    max_refinements::Int,
    max_points::Int,
)
    _validate_adaptive_controls(max_refinements, max_points)
    grid = _unique_sorted(grid)
    prof = _profile_from_grid(result, index, grid, threshold, on_failure)
    thresholds = Float64[threshold]

    for _ in 1:max_refinements
        length(prof.values) >= max_points && break
        candidates = _profile_refinement_candidates(prof.values, prof.delta_cost, thresholds)
        isempty(candidates) && break
        remaining = max_points - length(prof.values)
        grid = _unique_sorted(vcat(prof.values, candidates[1:min(end, remaining)]))
        length(grid) == length(prof.values) && break
        prof = _profile_from_grid(result, index, grid, threshold, on_failure)
    end

    return prof
end

"""
    profile(result, index; values=nothing, npoints=61, nsigma=3, threshold=1.0, adaptive=false, on_failure=:inf)

Profile the fitted cost function in one parameter by fixing that parameter to
grid values and re-minimizing all remaining free parameters.

With `adaptive=true`, JuFitter refines grid intervals that bracket the requested
profile threshold. This improves interval extraction without forcing a dense
grid over the full scan range.
"""
function profile(
    result,
    index::Int;
    values=nothing,
    npoints::Int=61,
    nsigma::Real=3,
    threshold::Real=1.0,
    adaptive::Bool=false,
    max_refinements::Int=3,
    max_points::Int=241,
    on_failure::Symbol=:inf,
)
    1 <= index <= length(result.params) || throw(ArgumentError("profile index out of range"))
    _validate_profile_controls(npoints, nsigma, threshold; default_grid=values === nothing)
    grid = values === nothing ?
        _validated_scan_values(_default_profile_grid(result, index; npoints=npoints, nsigma=nsigma), "profile values"; min_points=3) :
        _validated_scan_values(values, "profile values"; min_points=3)
    if adaptive
        return _adaptive_profile(
            result,
            index,
            grid,
            Float64(threshold),
            on_failure;
            max_refinements=max_refinements,
            max_points=max_points,
        )
    end
    return _profile_from_grid(result, index, _unique_sorted(grid), Float64(threshold), on_failure)
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
    adaptive::Bool=true,
    max_refinements::Int=3,
    max_points::Int=241,
)
    prof = profile(
        result,
        index;
        values=values,
        npoints=npoints,
        nsigma=nsigma,
        threshold=threshold,
        adaptive=adaptive,
        max_refinements=max_refinements,
        max_points=max_points,
    )
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

function _profile_matrix_parameter_names(result, parameters::Vector{Int}, parameter_names)
    if parameter_names !== nothing
        names = collect(String, parameter_names)
        length(names) == length(parameters) || throw(ArgumentError("parameter_names length must match parameters"))
        return names
    end
    if hasproperty(result.problem, :parameter_names) && result.problem.parameter_names !== nothing
        return result.problem.parameter_names[parameters]
    end
    return ["p$(index)" for index in parameters]
end

function _profile_matrix_parameters(result, parameters)
    selected = parameters === nothing ? collect(eachindex(result.params)) : collect(Int, parameters)
    isempty(selected) && throw(ArgumentError("parameters must contain at least one parameter index"))
    all(index -> 1 <= index <= length(result.params), selected) ||
        throw(ArgumentError("parameters contains an out-of-range parameter index"))
    length(unique(selected)) == length(selected) ||
        throw(ArgumentError("parameters must be unique"))
    return selected
end

function _combine_profile_matrix_findings(profile_diagnostics, contour_diagnostics)
    findings = DiagnosticFinding[]
    for report in values(profile_diagnostics)
        append!(findings, report.findings)
    end
    for report in values(contour_diagnostics)
        append!(findings, report.findings)
    end
    findings = _sort_findings(_deduplicate_findings(findings))
    return DiagnosticReport(findings, _diagnostic_summary(findings))
end

function _profile_matrix_panel_status(profile_diagnostics, contour_diagnostics)
    status = Dict{Tuple{Int, Int}, Symbol}()
    for (index, report) in profile_diagnostics
        status[(index, index)] = _diagnostic_status(report.findings)
    end
    for (indices, report) in contour_diagnostics
        status[indices] = _diagnostic_status(report.findings)
    end
    return status
end

"""
    profile_matrix(result; parameters=nothing, parameter_names=nothing, ...)

Compute a multi-parameter profile/contour diagnostic matrix without loading
Makie. Diagonal entries are one-parameter profile scans; lower-triangle entries
are two-parameter profile contours. Each panel is diagnosed against the local
covariance approximation when local errors or covariance entries are finite.

Use this when a fit has several correlated or nonlinear parameters and you
need a quick, machine-readable answer to: "Are local symmetric covariance
errors enough, or do I need profile/contour intervals?"
"""
function profile_matrix(
    result;
    parameters=nothing,
    parameter_names=nothing,
    npoints_profile::Int=61,
    npoints_contour::Int=25,
    nsigma::Real=3,
    profile_threshold::Real=1.0,
    contour_levels::AbstractVector=[2.30, 6.18],
    adaptive::Bool=false,
    max_refinements::Int=2,
    max_points::Int=1200,
    profile_tolerance::Real=0.25,
    contour_tolerance::Real=0.5,
)
    selected = _profile_matrix_parameters(result, parameters)
    names = _profile_matrix_parameter_names(result, selected, parameter_names)
    profiles = Dict{Int, ProfileResult}()
    contours = Dict{Tuple{Int, Int}, ContourResult}()
    profile_diagnostics = Dict{Int, DiagnosticReport}()
    contour_diagnostics = Dict{Tuple{Int, Int}, DiagnosticReport}()

    for index in selected
        prof = profile(
            result,
            index;
            npoints=npoints_profile,
            nsigma=nsigma,
            threshold=profile_threshold,
            adaptive=adaptive,
            max_refinements=max_refinements,
            max_points=max_points,
        )
        profiles[index] = prof
        sigma = result.param_stderr[index]
        profile_diagnostics[index] = isfinite(sigma) && sigma > 0 ?
            diagnose(prof; local_sigma=sigma, tolerance=profile_tolerance) :
            diagnose(prof)
    end

    for row in 2:length(selected), col in 1:(row - 1)
        xindex = selected[col]
        yindex = selected[row]
        key = (xindex, yindex)
        cont = contour(
            result,
            xindex,
            yindex;
            npoints=npoints_contour,
            nsigma=nsigma,
            levels=contour_levels,
            adaptive=adaptive,
            max_refinements=max_refinements,
            max_points=max_points,
        )
        contours[key] = cont
        local_cov = result.param_covariance[[xindex, yindex], [xindex, yindex]]
        contour_diagnostics[key] = diagnose(
            cont;
            local_covariance=local_cov,
            local_center=result.params[[xindex, yindex]],
            tolerance=contour_tolerance,
        )
    end

    return ProfileMatrixResult(
        selected,
        names,
        Float64.(result.params[selected]),
        Float64.(result.param_stderr[selected]),
        Matrix{Float64}(result.param_covariance[selected, selected]),
        Matrix{Float64}(result.param_correlation[selected, selected]),
        profiles,
        contours,
        profile_diagnostics,
        contour_diagnostics,
        _profile_matrix_panel_status(profile_diagnostics, contour_diagnostics),
        _combine_profile_matrix_findings(profile_diagnostics, contour_diagnostics),
    )
end

diagnose(matrix_result::ProfileMatrixResult) = matrix_result.report

function _profile_matrix_panel_report(matrix_result::ProfileMatrixResult, indices::Tuple{Int, Int})
    if indices[1] == indices[2]
        return matrix_result.profile_diagnostics[indices[1]]
    end
    return matrix_result.contour_diagnostics[indices]
end

function _profile_matrix_parameter_name(matrix_result::ProfileMatrixResult, index::Int)
    position = findfirst(==(index), matrix_result.parameters)
    position === nothing && throw(ArgumentError("parameter index $index is not part of the profile matrix"))
    return matrix_result.parameter_names[position]
end

function _profile_matrix_panel_order(matrix_result::ProfileMatrixResult)
    panels = Tuple{Int, Int}[]
    selected = matrix_result.parameters
    for index in selected
        push!(panels, (index, index))
    end
    for row in 2:length(selected), col in 1:(row - 1)
        push!(panels, (selected[col], selected[row]))
    end
    return panels
end

_profile_matrix_status_rank(status::Symbol) =
    status == :stop ? 1 : status == :review ? 2 : status == :ok ? 3 : 4

"""
    profile_matrix_triage(matrix_result; include_ok=false)

Return a sorted, Makie-free list of profile-matrix panels that need attention.
By default only `:review` and `:stop` panels are returned. Set
`include_ok=true` to include successful panels as well, for example when a
notebook wants to print a complete audit table.

The returned `ProfileMatrixPanelTriage` rows are ordered by severity first and
then by the natural matrix panel order. This is the programmatic counterpart of
the status labels in `plot_profile_matrix(...)`.
"""
function profile_matrix_triage(matrix_result::ProfileMatrixResult; include_ok::Bool=false)
    ordered_rows = Tuple{Int, ProfileMatrixPanelTriage}[]
    for (order, indices) in enumerate(_profile_matrix_panel_order(matrix_result))
        status = matrix_result.panel_status[indices]
        include_ok || status != :ok || continue
        report = _profile_matrix_panel_report(matrix_result, indices)
        names = (
            _profile_matrix_parameter_name(matrix_result, indices[1]),
            _profile_matrix_parameter_name(matrix_result, indices[2]),
        )
        next_action = isempty(report.findings) ?
            "No action required by the current profile/contour checks." :
            first(report.findings).recommendation
        row =
            ProfileMatrixPanelTriage(
                indices,
                names,
                status,
                _severity_counts(report.findings),
                [finding.code for finding in report.findings],
                next_action,
            )
        push!(ordered_rows, (order, row))
    end
    sort!(ordered_rows; by=entry -> (_profile_matrix_status_rank(entry[2].status), entry[1]))
    return [row for (_, row) in ordered_rows]
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

function _contour_refit_cost(result, i::Int, j::Int, xvalue::Float64, yvalue::Float64, on_failure::Symbol)
    return _profile_refit_cost(result, [FixedParameter(i, xvalue), FixedParameter(j, yvalue)]; on_failure=on_failure)
end

function _contour_from_grid(result, i::Int, j::Int, xs::Vector{Float64}, ys::Vector{Float64}, levels::Vector{Float64}, on_failure::Symbol)
    cache = Dict{Tuple{Float64, Float64}, Float64}()
    return _contour_from_grid!(cache, result, i, j, xs, ys, levels, on_failure)
end

function _contour_from_grid!(cache, result, i::Int, j::Int, xs::Vector{Float64}, ys::Vector{Float64}, levels::Vector{Float64}, on_failure::Symbol)
    costs = Matrix{Float64}(undef, length(xs), length(ys))
    for ix in eachindex(xs), iy in eachindex(ys)
        key = (xs[ix], ys[iy])
        costs[ix, iy] = get!(cache, key) do
            _contour_refit_cost(result, i, j, xs[ix], ys[iy], on_failure)
        end
    end
    delta = costs .- result.stats.cost_min
    return ContourResult((i, j), xs, ys, costs, delta, levels)
end

function _contour_refinement_candidates(contour_result::ContourResult)
    xs = contour_result.x_values
    ys = contour_result.y_values
    delta = contour_result.delta_cost
    levels = contour_result.levels
    x_candidates = Float64[]
    y_candidates = Float64[]

    for ix in 1:(length(xs) - 1), iy in 1:(length(ys) - 1)
        corners = (delta[ix, iy], delta[ix + 1, iy], delta[ix, iy + 1], delta[ix + 1, iy + 1])
        all(isfinite, corners) || continue
        lo = minimum(corners)
        hi = maximum(corners)
        if any(level -> lo <= level <= hi, levels)
            push!(x_candidates, (xs[ix] + xs[ix + 1]) / 2)
            push!(y_candidates, (ys[iy] + ys[iy + 1]) / 2)
        end
    end

    return _unique_sorted(x_candidates), _unique_sorted(y_candidates)
end

function _adaptive_contour(
    result,
    i::Int,
    j::Int,
    xs::Vector{Float64},
    ys::Vector{Float64},
    levels::Vector{Float64},
    on_failure::Symbol;
    max_refinements::Int,
    max_points::Int,
)
    _validate_adaptive_controls(max_refinements, max_points)
    xs = _unique_sorted(xs)
    ys = _unique_sorted(ys)
    cache = Dict{Tuple{Float64, Float64}, Float64}()
    cont = _contour_from_grid!(cache, result, i, j, xs, ys, levels, on_failure)

    for _ in 1:max_refinements
        length(cont.x_values) * length(cont.y_values) >= max_points && break
        x_candidates, y_candidates = _contour_refinement_candidates(cont)
        isempty(x_candidates) && isempty(y_candidates) && break

        candidate_xs = copy(cont.x_values)
        candidate_ys = copy(cont.y_values)
        for xvalue in x_candidates
            xvalue in candidate_xs && continue
            (length(candidate_xs) + 1) * length(candidate_ys) <= max_points || break
            push!(candidate_xs, xvalue)
        end
        for yvalue in y_candidates
            yvalue in candidate_ys && continue
            length(candidate_xs) * (length(candidate_ys) + 1) <= max_points || break
            push!(candidate_ys, yvalue)
        end
        candidate_xs = _unique_sorted(candidate_xs)
        candidate_ys = _unique_sorted(candidate_ys)

        if length(candidate_xs) == length(cont.x_values) && length(candidate_ys) == length(cont.y_values)
            break
        end

        cont = _contour_from_grid!(cache, result, i, j, candidate_xs, candidate_ys, levels, on_failure)
    end

    return cont
end

"""
    contour(result, i, j; xvalues=nothing, yvalues=nothing, npoints=31, nsigma=3, levels=[2.30, 6.18], adaptive=false, on_failure=:inf)

Compute a two-parameter profile-likelihood contour grid. At each grid point,
parameters `i` and `j` are fixed and all remaining free parameters are
re-minimized.

With `adaptive=true`, JuFitter refines grid cells whose corner values bracket a
requested contour level. This concentrates expensive refits near meaningful
contour geometry instead of spreading them uniformly across the full rectangle.
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
    adaptive::Bool=false,
    max_refinements::Int=2,
    max_points::Int=2601,
    on_failure::Symbol=:inf,
)
    i != j || throw(ArgumentError("contour requires two distinct parameter indices"))
    1 <= i <= length(result.params) || throw(ArgumentError("first contour index out of range"))
    1 <= j <= length(result.params) || throw(ArgumentError("second contour index out of range"))
    _validate_contour_controls(npoints, nsigma; default_x_grid=xvalues === nothing, default_y_grid=yvalues === nothing)

    xs = xvalues === nothing ?
        _validated_scan_values(_default_contour_grid(result, i; npoints=npoints, nsigma=nsigma), "xvalues"; min_points=2) :
        _validated_scan_values(xvalues, "xvalues"; min_points=2)
    ys = yvalues === nothing ?
        _validated_scan_values(_default_contour_grid(result, j; npoints=npoints, nsigma=nsigma), "yvalues"; min_points=2) :
        _validated_scan_values(yvalues, "yvalues"; min_points=2)
    level_values = collect(Float64, levels)
    isempty(level_values) && throw(ArgumentError("contour levels must not be empty"))
    all(isfinite, level_values) || throw(ArgumentError("contour levels must be finite"))
    all(>(0.0), level_values) || throw(ArgumentError("contour levels must be positive delta-cost thresholds"))
    level_values = sort!(unique!(level_values))
    if adaptive
        return _adaptive_contour(
            result,
            i,
            j,
            xs,
            ys,
            level_values,
            on_failure;
            max_refinements=max_refinements,
            max_points=max_points,
        )
    end
    return _contour_from_grid(result, i, j, _unique_sorted(xs), _unique_sorted(ys), level_values, on_failure)
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

function _contour_local_covariance_precision(cov::AbstractMatrix)
    all(isfinite, cov) || return nothing, _finding(
        :warning,
        :contour_local_covariance_unavailable,
        "Local covariance ellipse is unavailable",
        "The supplied local covariance contains non-finite entries.",
        "Read the profile contour directly; do not interpret symmetric local covariance errors for this parameter pair.",
    )

    if !isapprox(cov, cov'; rtol=1e-10, atol=1e-12)
        return nothing, _finding(
            :warning,
            :contour_local_covariance_unavailable,
            "Local covariance ellipse is unavailable",
            "The supplied local covariance is not symmetric.",
            "Check the covariance construction. Use the profiled contour until the local covariance estimate is valid.",
        )
    end
    cov = (cov + cov') ./ 2

    try
        factor = cholesky(Symmetric(cov); check=true)
        return factor \ Matrix{Float64}(I, 2, 2), nothing
    catch
        return nothing, _finding(
            :warning,
            :contour_local_covariance_unavailable,
            "Local covariance ellipse is unavailable",
            "The supplied local covariance is not positive definite.",
            "Inspect parameter degeneracy, active bounds, or scaling. Use profile contours instead of symmetric local covariance errors.",
        )
    end
end

function _contour_ellipticity_findings(contour_result::ContourResult, local_covariance, local_center; tolerance::Real)
    cov = _contour_local_covariance(contour_result, local_covariance)
    cov === nothing && return DiagnosticFinding[]
    tolerance >= 0 || throw(ArgumentError("tolerance must be non-negative"))

    center = _contour_center(contour_result, local_center)
    precision, covariance_finding = _contour_local_covariance_precision(cov)
    covariance_finding === nothing || return DiagnosticFinding[covariance_finding]
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
