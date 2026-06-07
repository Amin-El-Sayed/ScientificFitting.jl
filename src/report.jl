"""
    ParameterEstimate

One fitted parameter as stored in a `FitReport`. It records the public name,
best-fit value, symmetric display uncertainty, asymmetric lower/upper
uncertainties when available, and whether the parameter was fixed rather than
fitted.
"""
struct ParameterEstimate
    index::Int
    name::String
    value::Float64
    uncertainty::Float64
    uncertainty_minus::Float64
    uncertainty_plus::Float64
    fixed::Bool
end

"""
    FitReport

Serializable summary returned by `fit_report(result)`. It separates parameter
estimates, statistics, covariance/correlation matrices, solver status, and
diagnostics from the full fit object so reports can be printed, tested, or
exported without depending on Makie.
"""
struct FitReport
    parameters::Vector{ParameterEstimate}
    statistics::NamedTuple
    covariance::Matrix{Float64}
    correlation::Matrix{Float64}
    backend::Symbol
    converged::Bool
    iterations::Int
    message::String
    diagnostics::FitDiagnostics
end

function _plain_parameter_name(name)
    raw = name isa LaTeXString ? String(name) : string(name)
    return _strip_math_delims(raw)
end

function _parameter_estimates(
    result,
    parameter_names;
    errors::Symbol=:local,
    profile_threshold::Real=1.0,
    profile_npoints::Int=121,
    profile_nsigma::Real=5,
)
    n = length(result.params)
    names = if parameter_names === nothing
        if hasproperty(result.problem, :parameter_names) && result.problem.parameter_names !== nothing
            result.problem.parameter_names
        else
            ["p$i" for i in 1:n]
        end
    else
        length(parameter_names) == n || throw(ArgumentError("parameter_names length must match parameter count"))
        [_plain_parameter_name(name) for name in parameter_names]
    end

    fixed = _fixed_lookup(result.problem)
    estimates = ParameterEstimate[]
    for i in 1:n
        if haskey(fixed, i)
            fp = fixed[i]
            uncertainty = max(fp.sigma_minus, fp.sigma_plus)
            push!(estimates, ParameterEstimate(i, names[i], result.params[i], uncertainty, fp.sigma_minus, fp.sigma_plus, true))
        else
            uncertainty = result.param_stderr[i]
            if errors == :profile
                interval = profile_interval(
                    result,
                    i;
                    threshold=profile_threshold,
                    npoints=profile_npoints,
                    nsigma=profile_nsigma,
                )
                uncertainty_minus = isfinite(interval.uncertainty_minus) ? interval.uncertainty_minus : uncertainty
                uncertainty_plus = isfinite(interval.uncertainty_plus) ? interval.uncertainty_plus : uncertainty
                uncertainty = max(uncertainty_minus, uncertainty_plus)
                push!(estimates, ParameterEstimate(i, names[i], result.params[i], uncertainty, uncertainty_minus, uncertainty_plus, false))
            elseif errors == :local
                push!(estimates, ParameterEstimate(i, names[i], result.params[i], uncertainty, uncertainty, uncertainty, false))
            else
                throw(ArgumentError("errors must be :local or :profile"))
            end
        end
    end
    return estimates
end

"""
    fit_report(result; parameter_names=nothing, errors=:local)

Return an extractable report object for a fit result. Parameters are available as
`report.parameters[i].value` and `report.parameters[i].uncertainty`.

Use `errors=:profile` to compute profile-based asymmetric uncertainties. This
re-runs fits and can be expensive.
"""
function fit_report(
    result;
    parameter_names=nothing,
    errors::Symbol=:local,
    profile_threshold::Real=1.0,
    profile_npoints::Int=121,
    profile_nsigma::Real=5,
)
    statistics = (
        cost=result.stats.cost,
        cost_min=result.stats.cost_min,
        nll_min=result.stats.nll_min,
        chi2=result.stats.chi2,
        chi2_ndf=result.stats.chi2_ndf,
        ndf=result.stats.ndf,
        pvalue=result.stats.pvalue,
        aic=result.stats.aic,
        bic=result.stats.bic,
    )

    return FitReport(
        _parameter_estimates(
            result,
            parameter_names;
            errors=errors,
            profile_threshold=profile_threshold,
            profile_npoints=profile_npoints,
            profile_nsigma=profile_nsigma,
        ),
        statistics,
        copy(result.param_covariance),
        copy(result.param_correlation),
        result.backend,
        result.converged,
        result.iterations,
        result.message,
        result.diagnostics,
    )
end

function _report_lines(report::FitReport; sigdigits::Int=6)
    lines = String[
        "Fit report",
        "backend = $(report.backend)",
        "converged = $(report.converged)",
        "iterations = $(report.iterations)",
        "message = $(report.message)",
        "",
        "Parameters:",
    ]

    for p in report.parameters
        value = _fmt_value(p.value; sigdigits=sigdigits)
        if p.uncertainty_minus != p.uncertainty_plus
            minus = _fmt_value(p.uncertainty_minus; sigdigits=sigdigits)
            plus = _fmt_value(p.uncertainty_plus; sigdigits=sigdigits)
            suffix = p.fixed ? " (fixed)" : ""
            push!(lines, "  $(p.name) = $value -$minus +$plus$suffix")
        else
            uncertainty = _fmt_value(p.uncertainty; sigdigits=sigdigits)
            suffix = p.fixed ? " (fixed)" : ""
            push!(lines, "  $(p.name) = $value +/- $uncertainty$suffix")
        end
    end

    push!(lines, "")
    push!(lines, "Statistics:")
    push!(lines, "  cost = $(report.statistics.cost)")
    push!(lines, "  cost_min = $(_fmt_value(report.statistics.cost_min; sigdigits=sigdigits))")
    push!(lines, "  nll_min = $(_fmt_value(report.statistics.nll_min; sigdigits=sigdigits))")
    push!(lines, "  chi2 = $(_fmt_value(report.statistics.chi2; sigdigits=sigdigits))")
    push!(lines, "  ndf = $(report.statistics.ndf)")
    push!(lines, "  chi2/ndf = $(_fmt_value(report.statistics.chi2_ndf; sigdigits=sigdigits))")
    push!(lines, "  pvalue = $(_fmt_value(report.statistics.pvalue; sigdigits=sigdigits))")
    push!(lines, "  AIC = $(_fmt_value(report.statistics.aic; sigdigits=sigdigits))")
    push!(lines, "  BIC = $(_fmt_value(report.statistics.bic; sigdigits=sigdigits))")

    if !isempty(report.diagnostics.warnings)
        push!(lines, "")
        push!(lines, "Warnings:")
        append!(lines, ["  $warning" for warning in report.diagnostics.warnings])
    end

    if !isempty(report.diagnostics.findings)
        push!(lines, "")
        push!(lines, "Diagnosis:")
        for finding in report.diagnostics.findings
            push!(lines, "  [$(uppercase(String(finding.severity)))] $(finding.title)")
            push!(lines, "    evidence: $(finding.evidence)")
            push!(lines, "    action: $(finding.recommendation)")
        end
    end

    return lines
end

"""
    report_text(report::FitReport; sigdigits=6)
    report_text(result; parameter_names=nothing, sigdigits=6, kwargs...)

Render a `FitReport` or fit result as plain text. The result method first calls
`fit_report`, so keyword arguments such as `errors=:profile` are forwarded to
the report builder.
"""
function report_text(report::FitReport; sigdigits::Int=6)
    return join(_report_lines(report; sigdigits=sigdigits), "\n")
end

function report_text(result; parameter_names=nothing, sigdigits::Int=6, kwargs...)
    return report_text(fit_report(result; parameter_names=parameter_names, kwargs...); sigdigits=sigdigits)
end

function Base.show(io::IO, report::FitReport)
    print(io, report_text(report))
end
