"""
    ConfidenceInterval

Confidence interval for one fitted parameter, returned by
[`confidence_interval`](@ref). `parameter_index` is the parameter's 1-based
index in `result.params`. `lower` and `upper` are the interval limits; a
profile limit the scan did not bracket is `NaN`. `level` is the nominal
two-sided coverage probability and `method` the method that produced the
limits (`:hessian`, `:profile`, or `:bootstrap`).
"""
struct ConfidenceInterval
    parameter_index::Int
    lower::Float64
    upper::Float64
    level::Float64
    method::Symbol
end

function Base.show(io::IO, interval::ConfidenceInterval)
    print(
        io,
        "ConfidenceInterval(p", interval.parameter_index,
        ": [", _fmt_value(interval.lower), ", ", _fmt_value(interval.upper), "], ",
        Printf.@sprintf("%.4g", 100 * interval.level), "% via ", interval.method, ")",
    )
end

const _PROFILE_INTERVAL_KEYWORDS = (:values, :npoints, :nsigma, :adaptive, :max_refinements, :max_points, :on_failure)
const _BOOTSTRAP_INTERVAL_KEYWORDS = (:nsamples, :rng)

"""
    confidence_interval(result, parameter; method=:profile, level=0.6827,
                        nsamples=999, rng=Random.default_rng(), kwargs...)

Confidence interval for one free parameter of a `FitResult` or
`LikelihoodFitResult`. Returns a [`ConfidenceInterval`](@ref).

`parameter` is the 1-based index into `result.params`, or a name (`Symbol` or
string) looked up in `result.problem.parameter_names`. Only likelihood fits
store names; for a Gaussian `FitProblem`, pass the index.

`level` is the nominal two-sided coverage probability. The default 0.6827 is
the one-sigma convention of `param_stderr` and of the default
[`profile_interval`](@ref) threshold `delta_cost = 1`; 0.9545 corresponds to
two sigma. Coverage is asymptotic for `:hessian` and `:profile` and
approximate for `:bootstrap`.

`method` selects one of three intervals:

- `:profile` (default): profile-likelihood interval. The limits are the
  crossings of the profiled `delta_cost` with `quantile(Chisq(1), level)`
  (1.0000 for the default level, 4.0 for two sigma). `delta_cost` is on the
  scale of `param_stderr`, including a chi2/ndf covariance scaling applied by
  the fit. The keywords `values`, `npoints=121`, `nsigma=5`, `adaptive=true`,
  `max_refinements=3`, `max_points=241`, and `on_failure=:inf` are passed to
  [`profile`](@ref); the defaults are those of `profile_interval`. A limit
  whose crossing the scan does not bracket is `NaN`: the scan range or a
  parameter bound ends first, or a failed refit interrupts the scan. The
  interval is then one-sided or undefined; `diagnose` on a `profile` scan with
  the same keywords and threshold names the cause, and wider `nsigma` or
  explicit `values` extend the scan.
- `:hessian`: symmetric interval `params[i] ± z * param_stderr[i]` with
  `z = quantile(Normal(), (1 + level) / 2)`. It assumes a quadratic cost near
  the minimum; for a model linear in its parameters with a fixed observation
  covariance it coincides with `:profile`.
- `:bootstrap`: percentile interval, the empirical `(1 - level) / 2` and
  `(1 + level) / 2` quantiles of the parameter over the successful refits of
  [`parametric_bootstrap`](@ref) with the given `nsamples` and `rng`; that
  function defines the resampling model. Pass a seeded `rng` for reproducible
  limits.

Keywords apply only to their method: `nsamples` and `rng` to `:bootstrap`,
the profile keywords to `:profile`. The profile threshold follows from
`level` and cannot be passed separately.

Throws `DomainError` when `level` is not in (0, 1) or `nsamples < 1`, and
`ArgumentError` for an index outside `1:length(result.params)`, an unknown
name, a name lookup on a fit without stored names, a fixed parameter (an
input of the fit, not an estimate), an unknown `method`, a keyword the
selected method does not use, and for `:hessian` a non-finite
`param_stderr[i]` (for example after `parameter_covariance=:none`).
"""
function confidence_interval(
    result::Union{FitResult, LikelihoodFitResult},
    parameter::Union{Integer, Symbol, AbstractString};
    method::Symbol=:profile,
    level::Real=0.6827,
    kwargs...,
)
    0 < level < 1 || throw(DomainError(level, "level must lie strictly between 0 and 1"))
    index = _interval_parameter_index(result, parameter)
    lower, upper = _interval_limits(result, index, Float64(level), method; kwargs...)
    return ConfidenceInterval(index, lower, upper, Float64(level), method)
end

function _interval_parameter_index(result, index::Integer)
    nparams = length(result.params)
    1 <= index <= nparams || throw(ArgumentError("parameter index $index is outside 1:$nparams"))
    any(fp -> fp.index == index, result.problem.fixed_parameters) && throw(ArgumentError(
        "parameter $index is fixed; a fixed parameter is an input of the fit and has no confidence interval",
    ))
    return Int(index)
end

function _interval_parameter_index(result, name::Union{Symbol, AbstractString})
    names = hasproperty(result.problem, :parameter_names) ? result.problem.parameter_names : nothing
    names === nothing && throw(ArgumentError(
        "this fit stores no parameter names (a Gaussian FitProblem never does); pass the 1-based index into result.params",
    ))
    index = findfirst(==(String(name)), names)
    index === nothing && throw(ArgumentError(
        "unknown parameter name $(repr(String(name))); stored names: $(join(names, ", "))",
    ))
    return _interval_parameter_index(result, index)
end

function _interval_limits(result, index::Int, level::Float64, method::Symbol; kwargs...)
    if method == :hessian
        _check_interval_keywords(method, kwargs, ())
        return _hessian_limits(result, index, level)
    elseif method == :profile
        _check_interval_keywords(method, kwargs, _PROFILE_INTERVAL_KEYWORDS)
        return _profile_limits(result, index, level; kwargs...)
    elseif method == :bootstrap
        _check_interval_keywords(method, kwargs, _BOOTSTRAP_INTERVAL_KEYWORDS)
        return _bootstrap_limits(result, index, level; kwargs...)
    end
    throw(ArgumentError("method must be :hessian, :profile, or :bootstrap, got $(repr(method))"))
end

function _check_interval_keywords(method::Symbol, kwargs, supported::Tuple)
    unsupported = [key for key in keys(kwargs) if !(key in supported)]
    isempty(unsupported) && return nothing
    accepted = isempty(supported) ? "none" : join(supported, ", ")
    throw(ArgumentError(
        "method=:$method does not use the keyword(s) $(join(unsupported, ", ")); accepted: $accepted",
    ))
end

function _hessian_limits(result, index::Int, level::Float64)
    sigma = result.param_stderr[index]
    isfinite(sigma) || throw(ArgumentError(
        "param_stderr[$index] = $sigma is not finite; use method=:profile or method=:bootstrap",
    ))
    halfwidth = quantile(Normal(), (1 + level) / 2) * sigma
    center = result.params[index]
    return center - halfwidth, center + halfwidth
end

function _profile_limits(result, index::Int, level::Float64; npoints::Int=121, nsigma::Real=5, adaptive::Bool=true, kwargs...)
    scan = profile(
        result,
        index;
        npoints=npoints,
        nsigma=nsigma,
        adaptive=adaptive,
        kwargs...,
        threshold=quantile(Chisq(1), level),
    )
    interval = profile_interval(scan)
    return interval.lower, interval.upper
end

function _bootstrap_limits(result, index::Int, level::Float64; nsamples::Integer=999, rng::Random.AbstractRNG=Random.default_rng())
    nsamples >= 1 || throw(DomainError(nsamples, "nsamples must be at least 1"))
    samples = parametric_bootstrap(result; nsamples=nsamples, rng=rng)
    # Refit failures correlate with where the resample landed, so surviving
    # samples are a selected subsample and their percentiles are biased.
    samples.n_failed == 0 || throw(ErrorException(
        "$(samples.n_failed) of $(samples.n_requested) bootstrap refits failed; " *
        "a percentile interval over the survivors would be selection-biased. " *
        "Inspect parametric_bootstrap(result).n_failed and make every refit " *
        "feasible, for example with physical bounds."))
    lower, upper = quantile(view(samples.params, :, index), ((1 - level) / 2, (1 + level) / 2))
    return Float64(lower), Float64(upper)
end
