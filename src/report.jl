"""
    ParameterEstimate

One fitted parameter as stored in a `FitReport`. It records the public name,
best-fit value, symmetric display uncertainty, asymmetric lower/upper
uncertainties when available, and whether the parameter was fixed rather than
fitted. `index` is the position in the fit's parameter vector. For a fixed
parameter, the uncertainties are the user-declared values from
`FixedParameter`, not fit results.
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
    iterations::Union{Int, Missing}
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
        length(parameter_names) == n || throw(DimensionMismatch("parameter_names length must match parameter count"))
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
                # Missing crossings are not evidence for a local symmetric error.
                uncertainty_minus = interval.uncertainty_minus
                uncertainty_plus = interval.uncertainty_plus
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
    fit_report(result; parameter_names=nothing, errors=:local,
               profile_threshold=1.0, profile_npoints=121,
               profile_nsigma=5)

Return an extractable report object for a fit result. Parameters are available as
`report.parameters[i].value` and `report.parameters[i].uncertainty`.
`report.statistics` has `cost`, `cost_min`, `minus2loglik_min`, `chi2`,
`chi2_ndf`, `ndf`, `pvalue`, `aic`, and `bic`.

Use `errors=:profile` to compute profile-based asymmetric uncertainties. This
re-runs fits and can be expensive. The keywords map to `profile_interval`'s
`threshold`, `npoints`, and `nsigma`: `profile_threshold` is the delta-cost
level whose crossings define the interval, on the scale where `1.0` is the
one-sigma (68.3%) cut consistent with `param_stderr` (use `k^2` for a
k-sigma interval); `profile_npoints` sets the scan grid size per parameter;
`profile_nsigma` sets the scan half-width in units of the local standard
error. All three are ignored for `errors=:local`. A side where the scan does
not cross the threshold within its range stays `NaN`; it is never silently
replaced with the local symmetric error.
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
        minus2loglik_min=result.stats.minus2loglik_min,
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
    iterations = ismissing(report.iterations) ? "unavailable" : string(report.iterations)
    lines = String[
        "Fit report",
        "backend = $(report.backend)",
        "converged = $(report.converged)",
        "iterations = $iterations",
        "message = $(report.message)",
        "",
        "Parameters:",
    ]

    for p in report.parameters
        # The uncertainty sets the printed precision: two significant digits
        # of the error, value rounded to the same decimal place.
        suffix = p.fixed ? " (fixed)" : ""
        if p.uncertainty_minus != p.uncertainty_plus
            value, minus, plus = _value_error_strings(
                p.value, p.uncertainty_minus, p.uncertainty_plus; sigdigits=sigdigits)
            push!(lines, "  $(p.name) = $value -$minus +$plus$suffix")
        else
            value, uncertainty = _value_error_strings(
                p.value, p.uncertainty; sigdigits=sigdigits)
            push!(lines, "  $(p.name) = $value +/- $uncertainty$suffix")
        end
    end

    push!(lines, "")
    push!(lines, "Statistics:")
    push!(lines, "  cost = $(report.statistics.cost)")
    push!(lines, "  cost_min = $(_fmt_value(report.statistics.cost_min; sigdigits=sigdigits))")
    push!(lines, "  minus2loglik_min = $(_fmt_value(report.statistics.minus2loglik_min; sigdigits=sigdigits))")
    push!(lines, "  chi2 = $(_fmt_value(report.statistics.chi2; sigdigits=sigdigits))")
    push!(lines, "  ndf = $(report.statistics.ndf)")
    push!(lines, "  chi2/ndf = $(_fmt_value(report.statistics.chi2_ndf; sigdigits=sigdigits))")
    push!(lines, "  pvalue = $(_fmt_value(report.statistics.pvalue; sigdigits=sigdigits))")
    push!(lines, "  AIC = $(_fmt_value(report.statistics.aic; sigdigits=sigdigits))")
    push!(lines, "  BIC = $(_fmt_value(report.statistics.bic; sigdigits=sigdigits))")

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

Parameter lines round to the uncertainty: the error keeps two significant
digits and the value is cut at the same decimal place (`15.56 +/- 0.11`).
The result panels of the plot functions follow the same rule. `sigdigits`
controls the statistics lines, and parameter lines fall back to it when the
uncertainty has no compact fixed-point form (non-finite, zero, below `1e-6`,
or at `1e5` and above).
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
