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
        return cholesky(sym) \ Matrix{Float64}(I, size(sym))
    catch
        eigen_decomp = eigen(sym)
        values = eigen_decomp.values
        vectors = eigen_decomp.vectors
        scale = maximum(abs, values; init=0.0)
        cutoff = max(scale * rtol, eps(Float64))
        inv_values = [abs(v) > cutoff ? 1.0 / v : 0.0 for v in values]
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

function _diagnostic_warnings(
    converged::Bool,
    ndf::Int,
    cov_cond::Float64,
    hess_cond::Float64,
    active_bounds::Vector{Int};
    gof=nothing,
)
    warnings = String[]
    converged || push!(warnings, "optimizer did not report convergence")
    ndf <= 0 && push!(warnings, "non-positive degrees of freedom; p-values and reduced statistics are not meaningful")
    gof !== nothing && !isfinite(gof) && push!(warnings, "goodness-of-fit statistic is unavailable; p-values are not meaningful")
    !isnan(cov_cond) && cov_cond > 1e12 && push!(warnings, "parameter covariance is ill-conditioned")
    !isnan(hess_cond) && hess_cond > 1e12 && push!(warnings, "cost Hessian is ill-conditioned")
    !isempty(active_bounds) && push!(warnings, "one or more parameters are at active bounds; local errors and p-values may be unreliable")
    return warnings
end

function _finding(severity::Symbol, code::Symbol, title, evidence, recommendation)
    severity in (:info, :warning, :critical) ||
        throw(ArgumentError("diagnostic severity must be :info, :warning, or :critical"))
    return DiagnosticFinding(severity, code, String(title), String(evidence), String(recommendation))
end

function _basic_diagnostic_findings(
    converged::Bool,
    ndf::Int,
    cov_cond::Float64,
    hess_cond::Float64,
    active_bounds::Vector{Int};
    gof=nothing,
)
    findings = DiagnosticFinding[]

    if !converged
        push!(
            findings,
            _finding(
                :critical,
                :optimizer_not_converged,
                "Optimizer did not converge",
                "The optimizer did not report a successful convergence code.",
                "Do not trust parameter errors yet. Try better initial values, parameter scaling, bounds, multistart, or a simpler model.",
            ),
        )
    end

    if ndf <= 0
        push!(
            findings,
            _finding(
                :warning,
                :nonpositive_ndf,
                "No meaningful degrees of freedom",
                "ndf = $ndf, so reduced chi-square and p-values are not statistically meaningful.",
                "Use more observations, reduce the number of free parameters, or treat this as interpolation/calibration instead of a goodness-of-fit test.",
            ),
        )
    end

    if gof !== nothing && !isfinite(gof)
        push!(
            findings,
            _finding(
                :warning,
                :gof_unavailable,
                "Goodness-of-fit statistic is unavailable",
                "The fit does not provide a finite chi-square-like goodness-of-fit statistic.",
                "Use residual diagnostics, profiles, simulation, or a likelihood-specific goodness-of-fit test instead of interpreting p-values.",
            ),
        )
    end

    if !isnan(cov_cond) && cov_cond > 1e12
        severity = cov_cond > 1e16 ? :critical : :warning
        push!(
            findings,
            _finding(
                severity,
                :ill_conditioned_covariance,
                "Parameter covariance is ill-conditioned",
                "condition(covariance) = $(_fmt_scientific(cov_cond)).",
                "Check parameter degeneracies, units/scaling, too-flexible models, and whether profile intervals are more reliable than local covariance errors.",
            ),
        )
    end

    if !isnan(hess_cond) && hess_cond > 1e12
        severity = hess_cond > 1e16 ? :critical : :warning
        push!(
            findings,
            _finding(
                severity,
                :ill_conditioned_hessian,
                "Cost Hessian is ill-conditioned",
                "condition(Hessian) = $(_fmt_scientific(hess_cond)).",
                "Rescale parameters, inspect correlations and contours, and avoid interpreting symmetric local errors as final uncertainty.",
            ),
        )
    end

    if !isempty(active_bounds)
        push!(
            findings,
            _finding(
                :warning,
                :active_bounds,
                "One or more parameters are at bounds",
                "Active parameter indices: $(join(active_bounds, ", ")).",
                "Check whether the bound is physical. If yes, use profile intervals; if not, widen bounds or improve the model/start values.",
            ),
        )
    end

    return findings
end

function _fit_diagnostics(problem, params::AbstractVector, cov::AbstractMatrix, converged::Bool, ndf::Int; hessian=nothing, gof=nothing)
    cov_cond = _safe_condition_number(cov)
    hess_cond = hessian === nothing ? NaN : _safe_condition_number(hessian)
    active_bounds = _active_bound_indices(problem.bounds, params)
    warnings = _diagnostic_warnings(converged, ndf, cov_cond, hess_cond, active_bounds; gof=gof)
    findings = _basic_diagnostic_findings(converged, ndf, cov_cond, hess_cond, active_bounds; gof=gof)
    return FitDiagnostics(warnings, cov_cond, hess_cond, active_bounds, findings)
end

function _fmt_scientific(x::Real)
    !isfinite(x) && return string(x)
    return string(round(Float64(x); sigdigits=4))
end

function _data_pull_values(result::FitResult)
    n = length(result.problem.y)
    return result.weighted_residuals[1:n]
end

function _run_count(values::AbstractVector)
    isempty(values) && return 0
    signs = sign.(values)
    filtered = [s for s in signs if s != 0]
    isempty(filtered) && return 0
    runs = 1
    for i in 2:length(filtered)
        filtered[i] != filtered[i - 1] && (runs += 1)
    end
    return runs
end

function _lag1_autocorrelation(values::AbstractVector)
    length(values) < 3 && return NaN
    centered = values .- mean(values)
    denom = sum(abs2, centered)
    denom <= 0 && return NaN
    return sum(centered[1:(end - 1)] .* centered[2:end]) / denom
end

function _parameter_correlation_findings(result)
    corr = result.param_correlation
    n = size(corr, 1)
    n < 2 && return DiagnosticFinding[]

    best_abs_corr = 0.0
    best_pair = (0, 0)
    for i in 1:(n - 1), j in (i + 1):n
        value = abs(corr[i, j])
        if isfinite(value) && value > best_abs_corr
            best_abs_corr = value
            best_pair = (i, j)
        end
    end

    best_pair[1] == 0 && return DiagnosticFinding[]
    if best_abs_corr > 0.995
        severity = :critical
    elseif best_abs_corr > 0.95
        severity = :warning
    else
        return DiagnosticFinding[]
    end

    return DiagnosticFinding[
        _finding(
            severity,
            :strong_parameter_correlation,
            "Parameters are strongly correlated",
            "max |correlation| = $(_fmt_scientific(best_abs_corr)) between parameters $(best_pair[1]) and $(best_pair[2]).",
            "Inspect a contour/profile plot. Re-center or rescale the independent variable, reparameterize the model, or add data that breaks the degeneracy.",
        ),
    ]
end

function _xy_diagnostic_findings(result::FitResult)
    pulls = _data_pull_values(result)
    findings = DiagnosticFinding[]
    isempty(pulls) && return findings

    max_pull_index = argmax(abs.(pulls))
    max_abs_pull = abs(pulls[max_pull_index])
    max_pull_x = result.problem.x[max_pull_index]
    pull_evidence = "max |pull| = $(_fmt_scientific(max_abs_pull)) at point $max_pull_index (x = $(_fmt_scientific(max_pull_x)))."
    if isfinite(max_abs_pull) && max_abs_pull > 5
        push!(
            findings,
            _finding(
                :critical,
                :extreme_pull,
                "Extreme pull detected",
                pull_evidence,
                "Inspect the corresponding data point, uncertainty, units, and possible outlier handling before trusting the fit.",
            ),
        )
    elseif isfinite(max_abs_pull) && max_abs_pull > 3
        push!(
            findings,
            _finding(
                :warning,
                :large_pull,
                "Large pull detected",
                pull_evidence,
                "Inspect residuals near the largest pull. One point may dominate the result or the uncertainty model may be too optimistic.",
            ),
        )
    end

    if isfinite(result.stats.chi2_ndf)
        if result.stats.chi2_ndf > 5
            push!(
                findings,
                _finding(
                    :critical,
                    :very_large_reduced_chi2,
                    "Fit is very unlikely under the stated uncertainties",
                    "chi2/ndf = $(_fmt_scientific(result.stats.chi2_ndf)).",
                    "Look for missing physics, underestimated uncertainties, wrong correlations, outliers, or a failed optimizer.",
                ),
            )
        elseif result.stats.chi2_ndf > 2
            push!(
                findings,
                _finding(
                    :warning,
                    :large_reduced_chi2,
                    "Reduced chi-square is high",
                    "chi2/ndf = $(_fmt_scientific(result.stats.chi2_ndf)).",
                    "Check residual structure and uncertainty estimates. If residuals are structured, improve the model before tuning errors.",
                ),
            )
        elseif result.stats.chi2_ndf < 0.2
            push!(
                findings,
                _finding(
                    :warning,
                    :very_small_reduced_chi2,
                    "Data are too good for the assigned uncertainties",
                    "chi2/ndf = $(_fmt_scientific(result.stats.chi2_ndf)).",
                    "Uncertainties may be overestimated, correlations may be ignored, or the data may have been smoothed/averaged.",
                ),
            )
        end
    end

    if isfinite(result.stats.pvalue)
        if result.stats.pvalue < 1e-3
            push!(
                findings,
                _finding(
                    :critical,
                    :tiny_pvalue,
                    "Tiny chi-square probability",
                    "P(chi2) = $(_fmt_scientific(result.stats.pvalue)).",
                    "Under the stated assumptions this fit is statistically implausible. Inspect residuals and the uncertainty model.",
                ),
            )
        elseif result.stats.pvalue < 0.01
            push!(
                findings,
                _finding(
                    :warning,
                    :small_pvalue,
                    "Small chi-square probability",
                    "P(chi2) = $(_fmt_scientific(result.stats.pvalue)).",
                    "Treat the result as suspicious unless you can explain the residual pattern or uncertainty model.",
                ),
            )
        elseif result.stats.pvalue > 0.99
            push!(
                findings,
                _finding(
                    :warning,
                    :huge_pvalue,
                    "Chi-square probability is suspiciously high",
                    "P(chi2) = $(_fmt_scientific(result.stats.pvalue)).",
                    "The uncertainties may be too large, correlations may be ignored, or the data may not be independent.",
                ),
            )
        end
    end

    runs = _run_count(pulls)
    expected_runs = (length(pulls) + 1) / 2
    if length(pulls) >= 12 && runs < 0.45 * expected_runs
        push!(
            findings,
            _finding(
                :warning,
                :structured_residual_signs,
                "Residual signs look structured",
                "Observed $runs sign runs; roughly $(_fmt_scientific(expected_runs)) would be typical for structureless residuals.",
                "Look for missing curvature, drift, hysteresis, time dependence, or an incorrect independent variable transformation.",
            ),
        )
    end

    rho1 = _lag1_autocorrelation(pulls)
    if isfinite(rho1) && abs(rho1) > 0.45 && length(pulls) >= 12
        push!(
            findings,
            _finding(
                :warning,
                :autocorrelated_pulls,
                "Neighboring pulls are correlated",
                "lag-1 correlation = $(_fmt_scientific(rho1)).",
                "Use a covariance model, inspect acquisition order/time dependence, or fit a model with the missing systematic component.",
            ),
        )
    end

    return findings
end

function _deduplicate_findings(findings::Vector{DiagnosticFinding})
    seen = Set{Symbol}()
    out = DiagnosticFinding[]
    for finding in findings
        finding.code in seen && continue
        push!(seen, finding.code)
        push!(out, finding)
    end
    return out
end

function _severity_rank(severity::Symbol)
    severity == :critical && return 3
    severity == :warning && return 2
    return 1
end

function _sort_findings(findings::Vector{DiagnosticFinding})
    return sort(findings; by=f -> (-_severity_rank(f.severity), string(f.code)))
end

function _diagnostic_summary(findings::Vector{DiagnosticFinding})
    critical = count(f -> f.severity == :critical, findings)
    warning = count(f -> f.severity == :warning, findings)
    if critical > 0
        return "$critical critical issue(s), $warning warning(s). Do not treat this fit as publication-ready."
    elseif warning > 0
        return "$warning warning(s). Inspect before trusting uncertainties or conclusions."
    end
    return "No major diagnostic issues detected by the current checks."
end

"""
    diagnose(result)

Return a structured `DiagnosticReport` with actionable findings for a fit
result. The report is designed for quick notebook/lab use: each finding includes
a severity, evidence, and a concrete next step.
"""
function diagnose(result::FitResult)
    findings = DiagnosticFinding[]
    append!(findings, result.diagnostics.findings)
    append!(findings, _parameter_correlation_findings(result))
    append!(findings, _xy_diagnostic_findings(result))
    findings = _sort_findings(_deduplicate_findings(findings))
    return DiagnosticReport(findings, _diagnostic_summary(findings))
end

function diagnose(result)
    hasproperty(result, :diagnostics) || throw(ArgumentError("diagnose expects a fit result with diagnostics"))
    findings = copy(result.diagnostics.findings)
    if hasproperty(result, :param_correlation)
        append!(findings, _parameter_correlation_findings(result))
    end
    findings = _sort_findings(_deduplicate_findings(findings))
    return DiagnosticReport(findings, _diagnostic_summary(findings))
end

function _diagnostic_report_lines(report::DiagnosticReport)
    lines = String["Fit diagnosis", report.summary]
    if isempty(report.findings)
        push!(lines, "No action required by the current diagnostic checks.")
        return lines
    end

    for finding in report.findings
        push!(lines, "")
        push!(lines, "[$(uppercase(String(finding.severity)))] $(finding.title)")
        push!(lines, "  evidence: $(finding.evidence)")
        push!(lines, "  action: $(finding.recommendation)")
    end
    return lines
end

function diagnose_text(report::DiagnosticReport)
    return join(_diagnostic_report_lines(report), "\n")
end

diagnose_text(result) = diagnose_text(diagnose(result))

function Base.show(io::IO, report::DiagnosticReport)
    print(io, diagnose_text(report))
end
