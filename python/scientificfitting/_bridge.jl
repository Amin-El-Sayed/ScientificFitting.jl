using ScientificFitting, PythonCall, SparseArrays

# Fail before fitting when an installed wheel is paired with the old core.
supports_finite_derivatives = hasfield(LikelihoodFitProblem, :derivatives)

"""One vectorized foreign call; no dual numbers or per-observation Python loops."""
vector_model(f) = (x, p) -> pyconvert(Vector{Float64}, f(x, p))
matrix_model(f) = (x, p) -> pyconvert(Matrix{Float64}, f(x, p))
scalar_cost(f) = p -> pyconvert(Float64, f(p))
vector_constraint(f) = p -> pyconvert(Vector{Float64}, f(p))
scalar_model(f) = (x, p) -> pyconvert(Float64, f(x, p))
mutating_model(f) = (out, x, p) -> (f(out, x, p); nothing)
vector(x) = pyconvert(Vector{Float64}, Py(x))
matrix(x) = pyconvert(Matrix{Float64}, Py(x))

"""Reconstruct canonical CSC without allocating an n-by-n dense intermediary."""
function covariance(value::Py)
    pyisinstance(value, pybuiltins.dict) || return matrix(value)
    n, m = pyconvert(Tuple{Int, Int}, value["shape"])
    return SparseMatrixCSC(n, m, pyconvert(Vector{Int}, value["indptr"]) .+ 1,
        pyconvert(Vector{Int}, value["indices"]) .+ 1, vector(value["data"]))
end

"""Copy scalar, vector, or covariance metadata once, not during fit evaluation."""
function uncertainty_values(value::Py)
    pyisinstance(value, pybuiltins.dict) && return covariance(value)
    ndim = pyhasattr(value, "ndim") ? pyconvert(Int, value.ndim) : 0
    return ndim == 0 ? pyconvert(Float64, value) : ndim == 1 ? vector(value) : matrix(value)
end

"""Convert the supported Python keyword boundary once, outside numerical loops."""
function fit_keywords(options)
    result = Dict{Symbol, Any}(:derivatives => :finite)
    values = pyconvert(Dict{String, Py}, Py(options))
    inplace = haskey(values, "inplace") && pyconvert(Bool, values["inplace"])
    for (key, value) in values
        name = Symbol(key)
        pyis(value, pybuiltins.None) && continue
        result[name] = if name in (:backend, :cost, :scale_covariance, :cost_name)
            Symbol(pyconvert(String, value))
        elseif name in (:cov_x, :cov_y)
            covariance(value)
        elseif name == :whitening
            callback = value[0]
            marginal = pyis(value[2], pybuiltins.None) ? nothing : uncertainty_values(value[2])
            WhiteningOperator((out, residual) -> (callback(out, residual); nothing);
                logdet_covariance=pyconvert(Float64, value[1]), marginal_sigma=marginal)
        elseif name == :error_components
            [ErrorComponent(Symbol(pyconvert(String, row[0])), Symbol(pyconvert(String, row[1])),
                Symbol(pyconvert(String, row[2])), uncertainty_values(row[3]);
                active=pyconvert(Bool, row[4])) for row in value]
        elseif name in (:sigma_x, :sigma_y)
            vector(value)
        elseif name == :bounds
            lower, upper = pyconvert(Tuple{Py, Py}, value)
            (vector(lower), vector(upper))
        elseif name == :constraints
            callbacks = pyconvert(Dict{String, Py}, value)
            ConstraintSpec(;
                eq=haskey(callbacks, "eq") ? vector_constraint(callbacks["eq"]) : nothing,
                ineq=haskey(callbacks, "ineq") ? vector_constraint(callbacks["ineq"]) : nothing,
            )
        elseif name in (:parameter_priors, :fixed_parameters)
            constructor = name == :parameter_priors ? ParameterPrior : FixedParameter
            [constructor(Int(row[1]), row[2:end]...) for row in pyconvert(Vector{Vector{Float64}}, value)]
        elseif name == :parameter_constraints
            [ParameterConstraint(pyconvert(Vector{Int}, row[0]), vector(row[1]), matrix(row[2]))
             for row in value]
        elseif name == :parameter_names
            pyconvert(Vector{String}, value)
        elseif name == :initial_guesses
            pyconvert(Vector{Vector{Float64}}, value)
        elseif name == :jacobian
            inplace ? mutating_model(value) : matrix_model(value)
        elseif name == :x_derivative
            vector_model(value)
        elseif name == :gof
            scalar_cost(value)
        elseif name == :logprob
            (y, mu, p) -> pyconvert(Vector{Float64}, value(y, mu, p))
        elseif name in (:maxiters, :multistart, :nobs)
            pyconvert(Int, value)
        elseif name == :inplace
            pyconvert(Bool, value)
        else
            pyconvert(Float64, value)
        end
    end
    return result
end

"""Dispatch to existing Julia fits; the bridge owns conversion, never statistics."""
function run_fit(kind::String, callback::Py, x, y, p0, options)
    kwargs = fit_keywords(options)
    start = vector(p0)
    kind == "custom" && return fit_custom(scalar_cost(callback); p0=start, kwargs...)
    kind == "unbinned" && return fit_unbinned_model(scalar_model(callback), vector(y); p0=start, kwargs...)
    if kind == "extended_unbinned"
        domain = vector(x)
        length(domain) == 2 || throw(ArgumentError("domain must contain exactly two endpoints"))
        return fit_extended_unbinned_model(scalar_model(callback), vector(y), Tuple(domain); p0=start, kwargs...)
    end
    kind == "histogram_density" && return fit_histogram_density(scalar_model(callback), vector(x), vector(y); p0=start, kwargs...)
    kind in ("gaussian", "poisson", "histogram", "indexed", "likelihood") ||
        throw(ArgumentError("unknown fit family: $kind"))
    model = kind == "gaussian" && get(kwargs, :inplace, false) ? mutating_model(callback) : vector_model(callback)
    fit_function = kind == "gaussian" ? fit_model :
                   kind == "poisson" ? fit_poisson_model :
                   kind == "indexed" ? fit_indexed_model :
                   kind == "likelihood" ? fit_likelihood_model : fit_histogram_model
    return fit_function(model, vector(x), vector(y); p0=start, kwargs...)
end

"""Preserve the single global parameter map while converting dataset arrays once."""
function run_multi(callbacks, xs, ys, sigma, maps, p0, options)
    models = [vector_model(f) for f in Py(callbacks)]
    scales = [pyis(s, pybuiltins.None) ? nothing : vector(s) for s in Py(sigma)]
    return fit_multi_model(models, [vector(x) for x in Py(xs)], [vector(y) for y in Py(ys)];
        p0=vector(p0), sigma_y=scales, parameter_map=pyconvert(Vector{Vector{Int}}, Py(maps)),
        fit_keywords(options)...)
end

"""Convert scalar records only; symbols become strings and missing counts become None."""
scalar_value(value) = value isa Symbol ? String(value) : ismissing(value) ? nothing : value
scalar_fields(record) = pydict(String(name) => scalar_value(getproperty(record, name)) for name in propertynames(record))

"""Stored numerical checks, with parameter names instead of one-based indices."""
function numerical_values(diagnostics, names)
    return pydict(warnings=pylist(diagnostics.warnings),
        covariance_condition=diagnostics.covariance_condition, hessian_condition=diagnostics.hessian_condition,
        active_bounds=pylist(names[diagnostics.active_bounds]),
        findings=pylist(scalar_fields(f) for f in diagnostics.findings))
end

"""Transfer actual core reports; Python never infers findings by parsing text."""
function diagnostic_values(report::DiagnosticReport, max_actions::Int=5)
    dashboard = diagnostic_dashboard(report; max_actions)
    return pydict(findings=pylist(scalar_fields(f) for f in report.findings),
        summary=report.summary, status=String(dashboard.status),
        severity_counts=pydict(String(k) => v for (k, v) in dashboard.severity_counts),
        next_actions=pylist(dashboard.next_actions), text=diagnose_text(report),
        dashboard_text=diagnostic_dashboard_text(dashboard))
end

"""Return snapshots, retaining the Julia fit privately for later refits."""
function result_values(result, names)
    labels = pyconvert(Vector{String}, Py(names))
    return pydict(params=Py(result.params), stderr=Py(result.param_stderr),
        covariance=Py(result.param_covariance), correlation=Py(result.param_correlation),
        converged=result.converged, statistics=scalar_fields(result.stats),
        options=scalar_fields(result.options), backend=String(result.backend),
        iterations=scalar_value(result.iterations), message=result.message,
        numerical_diagnostics=numerical_values(result.diagnostics, labels),
        data=result isa FitResult ? pydict(x=Py(result.problem.x), y=Py(result.problem.y),
            model_y=Py(result.model_y), residuals=Py(result.residuals),
            weighted_residuals=Py(result.weighted_residuals), jacobian=Py(result.jacobian)) : pybuiltins.None)
end

function report_values(report::FitReport, names, sigdigits::Int)
    return pydict(parameters=pylist(pydict(name=p.name, value=p.value, uncertainty=p.uncertainty,
            uncertainty_minus=p.uncertainty_minus, uncertainty_plus=p.uncertainty_plus, fixed=p.fixed)
            for p in report.parameters),
        statistics=scalar_fields(report.statistics), covariance=Py(report.covariance),
        correlation=Py(report.correlation), backend=String(report.backend), converged=report.converged,
        iterations=scalar_value(report.iterations), message=report.message,
        numerical_diagnostics=numerical_values(report.diagnostics, pyconvert(Vector{String}, Py(names))),
        text=report_text(report; sigdigits))
end

function run_report(result, names, errors::String, threshold::Float64, npoints::Int, nsigma::Float64)
    return fit_report(result; parameter_names=pyconvert(Vector{String}, Py(names)),
        errors=Symbol(errors), profile_threshold=threshold, profile_npoints=npoints, profile_nsigma=nsigma)
end

result_diagnose(result, max_actions::Int) = diagnostic_values(diagnose(result), max_actions)
prediction(result, x, uncertainty::Bool) = predict(result, vector(x); uncertainty=uncertainty)

"""Pass scan controls without recomputing profile costs in Python."""
function scan_keywords(options)
    result = Dict{Symbol, Any}()
    for (key, value) in pyconvert(Dict{String, Py}, Py(options))
        name = Symbol(key)
        result[name] = if name in (:values, :xvalues, :yvalues, :levels, :contour_levels)
            vector(value)
        elseif name == :on_failure
            Symbol(pyconvert(String, value))
        else
            pyconvert(Any, value)
        end
    end
    return result
end

run_profile(result, index::Int, options) = profile(result, index; scan_keywords(options)...)
run_contour(result, i::Int, j::Int, options) = contour(result, i, j; scan_keywords(options)...)
run_interval(result, index::Int, options) = profile_interval(result, index; scan_keywords(options)...)
function run_matrix(result, indices, names, options)
    return profile_matrix(result; parameters=pyconvert(Vector{Int}, Py(indices)),
        parameter_names=pyconvert(Vector{String}, Py(names)), scan_keywords(options)...)
end

profile_diagnostics(scan::ProfileResult, sigma::Real) = diagnostic_values(
    isfinite(sigma) && sigma > 0 ? diagnose(scan; local_sigma=sigma) : diagnose(scan))
contour_diagnostics(scan::ContourResult, center, covariance) = diagnostic_values(
    diagnose(scan; local_center=vector(center), local_covariance=matrix(covariance)))

"""Keep core ordering and axis orientation while replacing indices with names."""
function matrix_values(result::ProfileMatrixResult)
    labels = Dict(zip(result.parameters, result.parameter_names))
    triage = profile_matrix_triage(result; include_ok=true)
    return pydict(parameters=pylist(result.parameter_names),
        best_values=Py(result.best_values), local_stderr=Py(result.local_stderr),
        local_covariance=Py(result.local_covariance), local_correlation=Py(result.local_correlation),
        profiles=pylist(pytuple((labels[i], Py(scan), diagnostic_values(result.profile_diagnostics[i])))
            for (i, scan) in result.profiles),
        contours=pylist(pytuple((labels[i], labels[j], Py(scan), diagnostic_values(result.contour_diagnostics[(i, j)])))
            for ((i, j), scan) in result.contours),
        panel_status=pydict(pytuple((labels[i], labels[j])) => String(status)
            for ((i, j), status) in result.panel_status),
        diagnostics=diagnostic_values(result.report),
        triage=pylist(pydict(parameters=pytuple(row.parameter_names), status=String(row.status),
            severity_counts=pydict(String(k) => v for (k, v) in row.severity_counts),
            finding_codes=pylist(String.(row.finding_codes)), next_action=row.next_action) for row in triage))
end

plot_errors(result) = (ScientificFitting._xerror_for_plot(result.problem, result.params),
    ScientificFitting._yerror_for_plot(result.problem, result.params))

diagnostic_data(result::FitResult, kind::String) = ScientificFitting._diagnostic_values(result, Symbol(kind))
