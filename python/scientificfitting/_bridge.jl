using ScientificFitting, PythonCall

# Fail before fitting when an installed wheel is paired with the old core.
supports_finite_derivatives = hasfield(LikelihoodFitProblem, :derivatives)

"""One vectorized foreign call; no dual numbers or per-observation Python loops."""
vector_model(f) = (x, p) -> pyconvert(Vector{Float64}, f(x, p))
matrix_model(f) = (x, p) -> pyconvert(Matrix{Float64}, f(x, p))
scalar_cost(f) = p -> pyconvert(Float64, f(p))
vector_constraint(f) = p -> pyconvert(Vector{Float64}, f(p))
vector(x) = pyconvert(Vector{Float64}, Py(x))
matrix(x) = pyconvert(Matrix{Float64}, Py(x))

"""Convert the supported Python keyword boundary once, outside numerical loops."""
function fit_keywords(options)
    result = Dict{Symbol, Any}(:derivatives => :finite)
    for (key, value) in pyconvert(Dict{String, Py}, Py(options))
        name = Symbol(key)
        pyis(value, pybuiltins.None) && continue
        result[name] = if name in (:backend, :cost, :scale_covariance, :cost_name)
            Symbol(pyconvert(String, value))
        elseif name in (:cov_x, :cov_y)
            matrix(value)
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
        elseif name == :initial_guesses
            pyconvert(Vector{Vector{Float64}}, value)
        elseif name == :jacobian
            matrix_model(value)
        elseif name == :x_derivative
            vector_model(value)
        elseif name in (:maxiters, :multistart, :nobs)
            pyconvert(Int, value)
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
    kind in ("gaussian", "poisson", "histogram") || throw(ArgumentError("unknown fit family: $kind"))
    model = vector_model(callback)
    fit_function = kind == "gaussian" ? fit_model :
                   kind == "poisson" ? fit_poisson_model : fit_histogram_model
    return fit_function(model, vector(x), vector(y); p0=start, kwargs...)
end

"""Return ordinary mappings at the Python boundary, retaining the Julia fit for refits."""
function result_values(result)
    stats = Dict(String(name) => getproperty(result.stats, name) for name in
        (:cost_min, :chi2, :chi2_ndf, :ndf, :pvalue, :aic, :bic))
    return pydict(params=Py(result.params), stderr=Py(result.param_stderr),
        covariance=Py(result.param_covariance), correlation=Py(result.param_correlation),
        converged=result.converged, statistics=stats)
end

result_report(result, names) = report_text(result; parameter_names=pyconvert(Vector{String}, Py(names)))
result_diagnose(result) = diagnostic_dashboard_text(result)
prediction(result, x, uncertainty::Bool) = predict(result, vector(x); uncertainty=uncertainty)

"""Pass scan controls without recomputing profile costs in Python."""
function scan_keywords(options)
    result = Dict{Symbol, Any}()
    for (key, value) in pyconvert(Dict{String, Py}, Py(options))
        name = Symbol(key)
        result[name] = if name in (:values, :xvalues, :yvalues, :levels)
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
plot_errors(result) = (ScientificFitting._xerror_for_plot(result.problem, result.params),
    ScientificFitting._yerror_for_plot(result.problem, result.params))
