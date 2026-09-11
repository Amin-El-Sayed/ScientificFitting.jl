"""Keep zero probability impossible and reject undefined or infinite densities."""
function _minus2logprob(value)
    value isa Real && (isfinite(value) || value == -Inf) || throw(ArgumentError(
        "log likelihood must be finite or -Inf (zero probability)",
    ))
    return -2 * value
end

"""Snapshot observations once; column-oriented joint events match Distributions.jl."""
function _distribution_data(data::AbstractArray{<:Real}; obsdim)
    !isempty(data) && all(isfinite, data) || throw(ArgumentError(
        "distribution observations must be nonempty and finite",
    ))
    if data isa AbstractVector
        obsdim === nothing || throw(ArgumentError("obsdim is only used for a matrix of multivariate events"))
        return collect(data), length(data)
    elseif data isa AbstractMatrix
        obsdim in (1, 2) || throw(ArgumentError(
            "set obsdim=1 for events in rows or obsdim=2 for events in columns",
        ))
        return obsdim == 2 ? Matrix(data) : permutedims(data), size(data, obsdim)
    end
    throw(ArgumentError("use a vector of scalar observations or a matrix of multivariate events"))
end

function _distribution_shape(d::Distribution, data)
    compatible = d isa UnivariateDistribution && data isa AbstractVector ||
                 d isa MultivariateDistribution && data isa AbstractMatrix && size(data, 1) == length(d)
    compatible || throw(ArgumentError("observation shape does not match the distribution's variate dimension"))
    return nothing
end

"""Use upstream log densities, including any specialized batched implementation."""
function _distribution_cost(d::Distribution, data)
    _distribution_shape(d, data)
    return _minus2logprob(_distribution_loglikelihood(d, data))
end

"""Supply the missing log-density interface without adding methods to upstream types."""
struct PDFLogAdapter{S, D} <: UnivariateDistribution{S}
    distribution::D
end
PDFLogAdapter(d::UnivariateDistribution{S}) where {S} = PDFLogAdapter{S, typeof(d)}(d)

# A fixed zero can be omitted; a zero with gradient or Hessian sensitivity cannot.
_structural_zero(x) = iszero(x)
_structural_zero(x::ForwardDiff.Dual) = _structural_zero(ForwardDiff.value(x)) &&
    all(_structural_zero, ForwardDiff.partials(x))

"""Sum weighted probabilities in log space without differentiating log(zero weight)."""
function _log_weighted_probability(logvalues, weights)
    active = ((v, w) for (v, w) in zip(logvalues, weights) if !_structural_zero(w))
    scale = maximum(first, active; init=-Inf)
    scale == -Inf && return scale
    return scale + log(sum(w * exp(v - scale) for (v, w) in active))
end

"""Retain sensitivity to a zero mixture weight; upstream logpdf omits that component."""
struct MixtureLogAdapter{V, S, D} <: Distribution{V, S}
    distribution::D
end
MixtureLogAdapter(d::MixtureModel{V, S}) where {V, S} = MixtureLogAdapter{V, S, typeof(d)}(d)
Base.length(d::MixtureLogAdapter{Multivariate}) = length(d.distribution)

function _weighted_logpdf(d, x)
    return _log_weighted_probability(map(part -> logpdf(part, x), components(d)), probs(d))
end
Distributions.logpdf(d::MixtureLogAdapter{Univariate}, x::Real) = _weighted_logpdf(d.distribution, x)
# Distributions' array-variate interface owns shape checking and batched calls.
Distributions._logpdf(d::MixtureLogAdapter{Multivariate}, x::AbstractVector{<:Real}) = _weighted_logpdf(d.distribution, x)

function Distributions.logpdf(d::PDFLogAdapter, x::Real)
    density = pdf(d.distribution, x)
    density isa Real && isfinite(density) && density >= 0 ||
        throw(ArgumentError("density must be finite and nonnegative"))
    return log(density)
end

_with_logpdf(d, x) = d
_with_logpdf(d::UnivariateDistribution, x::Real) = applicable(logpdf, d, x) ? d : PDFLogAdapter(d)

_with_logpdf(d::MixtureModel{Univariate}, x::Real) = _mixture_with_logpdf(d, x)
_with_logpdf(d::MixtureModel{Multivariate}, x::AbstractVector) = _mixture_with_logpdf(d, x)
function _mixture_with_logpdf(d, x)
    original = components(d)
    prepared = map(part -> _with_logpdf(part, x), original)
    model = all(a === b for (a, b) in zip(original, prepared)) ? d : MixtureModel(prepared, probs(d))
    # Inspect plain values only to choose the evaluation path. The actual sum
    # retains the original weights (and all their first/second derivatives).
    return any(w -> iszero(_finite_value(w)) && !_structural_zero(w), probs(d)) ? MixtureLogAdapter(model) : model
end

# Distributions 0.25 currently returns two product representations from its public
# constructor (vector vs tuple inputs). Keep that storage detail at this boundary.
_product_parts(d::Distributions.Product) = d.v
_product_parts(d::Distributions.VectorOfUnivariateDistribution) = d.dists
function _with_logpdf(d::Union{Distributions.Product, Distributions.VectorOfUnivariateDistribution}, x::AbstractVector)
    original = _product_parts(d)
    length(original) == length(x) || throw(DimensionMismatch("event dimension must match product distribution"))
    prepared = map((part, i) -> _with_logpdf(part, x[i]), original, eachindex(original))
    all(a === b for (a, b) in zip(original, prepared)) && return d
    return prepared isa Tuple ? product_distribution(prepared...) : product_distribution(prepared)
end

_distribution_logpdf(d, x) = logpdf(_with_logpdf(d, x), x)

_distribution_logpdfs(d, data) = logpdf(d, data)
_distribution_logpdfs(d::MixtureLogAdapter, data) = _mixture_logpdfs(d.distribution, data)
_distribution_logpdfs(d::MixtureModel, data) = isconcretetype(eltype(components(d))) ?
    logpdf(d, data) : _mixture_logpdfs(d, data)

"""Combine p = exp(scale) * total without taking the logarithm of a zero weight."""
function _merge_mixture_terms(scale, total, logs, weight)
    next_scale = max.(scale, logs)
    next_total = map(scale, total, logs, next_scale) do a, s, b, m
        # Outside both supports any finite scaled total represents zero probability.
        m == -Inf ? s + weight : s * exp(a-m) + weight * exp(b-m)
    end
    return next_scale, next_total
end

"""Combine component batches in log space with linear storage and one dispatch per component."""
function _weighted_logbatches(batch::F, parts, weights, n) where {F}
    state = nothing
    for (part, weight) in zip(parts, weights)
        isfinite(weight) && weight >= 0 || throw(ArgumentError(
            "mixture weights must be finite and nonnegative",
        ))
        _structural_zero(weight) && continue
        logs = batch(part)
        state = state === nothing ? (logs, fill(weight, length(logs))) :
                _merge_mixture_terms(state..., logs, weight)
    end
    state === nothing && return fill(-Inf, n)
    scale, total = state
    return scale .+ log.(total)
end

_mixture_logpdfs(d::MixtureModel, data) = _weighted_logbatches(
    part -> _distribution_logpdfs(part, data), components(d), probs(d), size(data, ndims(data)))

_event_batch(data::AbstractVector, indices) = view(data, indices)
_event_batch(data::AbstractMatrix, indices) = view(data, :, indices)

"""Bound event-sized AD scratch storage while retaining native component batches."""
function _batched_mixture_loglikelihood(d, data)
    # Only the reduction is blocked. Model construction, normalization and the
    # likelihood are unchanged; every event contributes, including the last block.
    blocks = Iterators.partition(1:size(data, ndims(data)), 4096)
    return sum(indices -> sum(_mixture_logpdfs(d, _event_batch(data, indices))), blocks)
end

_prepared_loglikelihood(d, data) = loglikelihood(d, data)
# Concrete component types already allow native per-event specialization.
_prepared_loglikelihood(d::MixtureModel, data) = isconcretetype(eltype(components(d))) ?
    loglikelihood(d, data) : _batched_mixture_loglikelihood(d, data)
_prepared_loglikelihood(d::MixtureLogAdapter, data) = _batched_mixture_loglikelihood(d.distribution, data)

function _distribution_loglikelihood(d, data)
    events = data isa AbstractMatrix ? eachcol(data) : data
    # Preparation and component dispatch belong outside the event loop.
    return _prepared_loglikelihood(_with_logpdf(d, first(events)), data)
end

_prepare_error_distribution(d::UnivariateDistribution, n) = deepcopy(d)
function _prepare_error_distribution(d::MultivariateDistribution, n)
    length(d) == n || throw(ArgumentError("joint error distribution dimension must match y"))
    return deepcopy(d)
end
function _prepare_error_distribution(ds::AbstractVector, n)
    length(ds) == n && all(d -> d isa UnivariateDistribution, ds) || throw(ArgumentError(
        "error must contain one univariate distribution per observation",
    ))
    return deepcopy(ds)
end
_prepare_error_distribution(d, n) = throw(ArgumentError(
    "error must be a univariate or joint multivariate distribution, or a vector of univariate distributions",
))

_error_loglikelihood(d::UnivariateDistribution, residuals) = _distribution_loglikelihood(d, residuals)
_error_loglikelihood(d::MultivariateDistribution, residuals) = _distribution_logpdf(d, residuals)
_error_loglikelihood(ds::AbstractVector, residuals) = sum(_distribution_logpdf(d, r) for (d, r) in zip(ds, residuals))

"""Keep the model factory available for reconstruction, not just a closed scalar loss."""
struct DistributionObjective{F, D}
    factory::F
    observations::D
end
(objective::DistributionObjective)(p) = _distribution_cost(objective.factory(p), objective.observations)

"""
    fitted_model(result::LikelihoodFitResult)

Reconstruct the distribution or event-intensity model of a [`fit_distribution`](@ref)
result using its fitted full parameter vector. The returned object has the
upstream package's native interface. For BuildConstructors fits this uses the
private constructor snapshot, even if the caller later updates the original.
No refit is performed. Arbitrary `fit_custom` results have no stored model
factory and raise `ArgumentError` rather than guessing how to reconstruct one.
"""
function fitted_model(result::LikelihoodFitResult)
    objective = result.problem.objective
    objective isa DistributionObjective || throw(ArgumentError(
        "fitted_model requires a result from fit_distribution; custom costs have no stored model factory",
    ))
    return objective.factory(copy(result.params))
end

"""
    fit_distribution(make_distribution, data; p0, obsdim=nothing, kwargs...)
        -> LikelihoodFitResult

Fit independent events using `make_distribution(p)`, which returns a normalized
Distributions.jl distribution. The factory runs once per objective evaluation,
not once per event: expensive parameter-dependent normalization belongs there.
ScientificFitting uses the upstream `loglikelihood`/`logpdf` implementation and
minimizes `-2 log(L)`. Heterogeneous mixtures batch the upstream component log
densities in bounded event blocks, so their temporary AD arrays do not grow with
the sample size. Every event contributes; there is no subsampling. Homogeneous
mixtures keep their specialized native loop. A scaled-sum
adapter retains first/second derivatives when a free mixture weight reaches zero.
A distribution exposing only `pdf` uses `log(pdf)`;
extreme-tail underflow then follows that upstream density implementation.
Univariate PDF-only components also work inside native mixtures and products:
an internal adapter supplies `logpdf` without changing the upstream model.
No parameters are inferred from a distribution instance;
declare which are fitted explicitly in the factory and `p0`.

For univariate continuous or discrete observations, pass a vector. For joint
multivariate events pass a matrix and specify `obsdim=1` (rows are events) or
`obsdim=2` (columns are events). Events are independent, but components of each
event may be correlated. `nobs` counts events, not scalar coordinates. Discrete
observations do not enable discrete fitted parameters. Data are copied once.

Distribution support, normalization and derivatives come from the upstream
package. Zero probability yields infinite cost, without clipping. Start inside
the support; use `derivatives=:finite` if a constructor does not support dual
numbers, or a derivative-free solver for nonsmooth parameter dependence.
Bounds, parameter terms, solvers, multistart and covariance options follow
[`fit_custom`](@ref). A universal goodness-of-fit p-value is not assumed.

After `using DistributionsHEP`, factories returning `ExtendedMixtureModel` are
also supported. These retain the extended Poisson yield term rather than being
treated as just a normalized mixture. Component distributions must describe the
observed domain (use upstream truncation for a selected interval).

With `using BuildConstructors`, an `AbstractConstructor` can replace the factory.
Its names, starts, bounds and fixed/shared state are read through the public
metadata API; conflicting shared descriptors are rejected. An optional named
`p0=(mu=0.2,)` overrides free starts or supplies missing ones. Names, bounds and
fixed state have one source of truth: set them on the constructor, not through
duplicate fit keywords. Numerical metadata uncertainties are not priors or
reported fixed-parameter errors. The fit holds a private constructor snapshot;
`fitted_model(result)` reconstructs the fitted model without changing the original
constructor. `BuildConstructors.parameter_values(result)` returns a named tuple
for an explicit `update!` if desired. No Minuit dependency is needed.

# Example

```julia
using Distributions
result = fit_distribution(p -> Normal(p[1], exp(p[2])), [-0.5, 0.2, 0.8, 1.1];
    p0=[0.0, 0.0], parameter_names=["mean", "log_scale"])
```
"""
function fit_distribution(make_distribution, data::AbstractArray{<:Real};
                          p0::AbstractVector, obsdim=nothing, kwargs...)
    observations, nobs = _distribution_data(data; obsdim)
    objective = DistributionObjective(make_distribution, observations)
    return fit_custom(objective; p0, nobs, cost_name=:distribution_likelihood, kwargs...)
end

"""Keep integration policy and observed bins separate from parameter-dependent models."""
struct DistributionHistogram
    edges::Vector{Float64}
    counts::Vector{Float64}
    total_count::Union{Nothing, Float64}
    integration::Symbol
    rtol::Float64
end

function _quadrature_logmass(d::ContinuousUnivariateDistribution, a, b, rtol)
    # Use ordinary quadrature nodes even when truncation bounds carry duals.
    # The affine map retains both moving-boundary terms in the derivatives.
    t = float(_finite_value(a))
    mass = _likelihood_integral(x -> (b-a)*pdf(d, a+(b-a)*x), zero(t), one(t); rtol)
    mass >= 0 || throw(ArgumentError("bin quadrature produced a negative probability"))
    return iszero(mass) ? -Inf : log(mass)
end
_quadrature_logmass(d, a, b, rtol) = throw(ArgumentError(
    "quadrature requires a continuous univariate distribution; use its CDF for discrete data",
))

_valid_cdf_interval(a, b) = 0 <= _finite_value(a) <= _finite_value(b) <= 1

"""Reintegrate roundoff-sized CDF violations; never clip probabilities or accept invalid CDFs."""
function _cdf_roundoff_logmass(d, a, b, bins, fa, fb)
    x, y = _finite_value(fa), _finite_value(fb)
    tolerance = 8eps(float(one(x+y)))
    if bins.integration == :auto && d isa ContinuousUnivariateDistribution &&
       isfinite(x) && isfinite(y) && -tolerance <= x <= y+tolerance && y <= 1+tolerance
        return _quadrature_logmass(d, a, b, bins.rtol)
    end
    throw(ArgumentError("CDF must be finite, nondecreasing, and between zero and one; " *
        "received ($x, $y). For a continuous density with an inaccurate CDF, use integration=:quadgk."))
end

"""Prefer upstream CDF integrals; retain tail log probabilities and handle exact support edges."""
function _bin_logmass(d::UnivariateDistribution, a, b, bins::DistributionHistogram)
    if bins.integration == :quadgk || (bins.integration == :auto && !applicable(cdf, d, a))
        return _quadrature_logmass(d, a, b, bins.rtol)
    end
    applicable(cdf, d, a) || throw(ArgumentError("distribution has no CDF; use integration=:quadgk"))
    fa, fb = cdf(d, a), cdf(d, b)
    _valid_cdf_interval(fa, fb) || return _cdf_roundoff_logmass(d, a, b, bins, fa, fb)
    # Avoid log(0) duals at exact support endpoints. Normal's specialized
    # logdiffcdf still handles tails where both ordinary CDFs underflow.
    fa == 0 && fb > 0 && return log(fb)
    qa, qb = fb == 1 ? (ccdf(d, a), ccdf(d, b)) : (one(fa), one(fb))
    _valid_cdf_interval(qb, qa) || return _cdf_roundoff_logmass(d, a, b, bins, qb, qa)
    qa > 0 && qb == 0 && return log(qa)
    value = logdiffcdf(d, b, a)
    isfinite(value) && return value
    qa > qb && return log(qa - qb)
    if bins.integration == :auto && d isa ContinuousUnivariateDistribution
        return _quadrature_logmass(d, a, b, bins.rtol)
    end
    fa == fb && return -Inf
    throw(ArgumentError("CDF bin integral is undefined; try integration=:quadgk for a continuous density"))
end

"""Integrate inside a continuous selection window, retaining the upstream normalizer."""
function _bin_logmass(d::Truncated{<:ContinuousUnivariateDistribution}, a, b, bins::DistributionHistogram)
    lo = d.lower === nothing ? a : max(a, d.lower)
    hi = d.upper === nothing ? b : min(b, d.upper)
    lo < hi || return -Inf
    # Truncated CDF subtraction can create log(0) duals at the window edges.
    return _bin_logmass(d.untruncated, lo, hi, bins) - d.logtp
end

function _bin_logmass(d::MixtureModel{Univariate}, a, b, bins::DistributionHistogram)
    # Integrate components before mixing: a PDF-only component may need quadrature.
    return _log_weighted_probability(map(part -> _bin_logmass(part, a, b, bins), components(d)), probs(d))
end

"""Specialize the bin loop once for each concrete upstream distribution type."""
function _bin_logmasses(d::UnivariateDistribution, edges, bins::DistributionHistogram)
    # Repeated edges can arise when clipping a selected window. Such bins have
    # exactly zero probability, not a failed CDF or quadrature calculation.
    return [edges[i] < edges[i+1] ? _bin_logmass(d, edges[i], edges[i+1], bins) : -Inf
            for i in 1:length(edges)-1]
end

function _bin_logmasses(d::Truncated{<:ContinuousUnivariateDistribution}, edges, bins::DistributionHistogram)
    clipped = clamp.(edges, something(d.lower, -Inf), something(d.upper, Inf))
    return _bin_logmasses(d.untruncated, clipped, bins) .- d.logtp
end

function _bin_logmasses(d::MixtureModel{Univariate}, edges, bins::DistributionHistogram)
    return _weighted_logbatches(part -> _bin_logmasses(part, edges, bins),
        components(d), probs(d), length(edges)-1)
end

function _bin_logexpectation(d::UnivariateDistribution, bins::DistributionHistogram)
    bins.total_count === nothing && throw(ArgumentError(
        "a normalized distribution requires total_count (expected events on its full support)",
    ))
    return log(bins.total_count) .+ _bin_logmasses(d, bins.edges, bins)
end
_bin_logexpectation(d, bins::DistributionHistogram) = throw(ArgumentError(
    "one-dimensional histogram fitting requires a univariate distribution or extended mixture",
))

function _distribution_cost(d::Distribution, bins::DistributionHistogram)
    log_mu = _bin_logexpectation(d, bins)
    mu = _nonnegative_expectation(exp.(log_mu), length(bins.counts))
    return _poisson_minus2loglik_terms(bins.counts, mu; log_mu)
end

"""
    fit_distribution(make_distribution, edges, counts;
        p0, total_count=nothing, integration=:auto, rtol=1e-8, kwargs...)

Fit independent Poisson bin counts using a univariate distribution factory.
`edges` must be finite and strictly increasing. Bin `i` is `(edges[i], edges[i+1]]`;
for integer data, half-integer edges avoid any endpoint convention ambiguity.
Counts must be nonnegative integers. The factory runs once per objective call.

For a normalized distribution, `total_count` is required: it is the known expected
event count on the distribution's full support, not an automatically fitted or
conditioned observed total. For a fitted rate use an `ExtendedMixtureModel` from
DistributionsHEP and omit `total_count`; its component yields supply the rates.
There is no implicit renormalization over the histogram window. Truncate the
upstream distribution explicitly when it describes a selected sample instead.

`integration=:auto` uses upstream CDF/log-CDF differences where available and
adaptive QuadGK integration otherwise. Roundoff-sized CDF range/order violations
are reintegrated from the PDF, not clipped; larger violations raise an error.
`:cdf` requires valid CDFs and does not apply this recovery. `:quadgk` forces
quadrature for continuous components. Its estimated error uses relative tolerance
`rtol` in the maximum norm of the value and all carried AD coefficients.
Moving truncation bounds retain their derivatives through a fixed-interval map.
Mixture components are integrated in batches, with dispatch outside the bin loop
and temporary storage linear in the number of bins. Their normalization and
parameter derivatives are retained. A missing dual-number implementation requires the usual
explicit `derivatives=:finite` option, not silent loss of derivatives.

The objective is normalized Poisson `-2log(L)`, with tail bin probabilities kept
in log space. Empty support bins contribute zero for zero observations and
infinite cost otherwise. Goodness-of-fit uses asymptotic Poisson deviance only
when all expected bin counts are strictly positive; otherwise those fields are
`NaN`. Low counts and fitted boundaries can also invalidate that approximation.
Other fit controls and [`fitted_model`](@ref) work as for unbinned distribution fits.
BuildConstructors metadata can be used with either form.
"""
function fit_distribution(make_distribution, edges::AbstractVector, counts::AbstractVector;
                          p0::AbstractVector, total_count=nothing,
                          integration::Symbol=:auto, rtol::Real=1e-8, kwargs...)
    integration in (:auto, :cdf, :quadgk) || throw(ArgumentError("integration must be :auto, :cdf, or :quadgk"))
    isfinite(rtol) && rtol > 0 || throw(ArgumentError("rtol must be finite and positive"))
    total_count === nothing || (total_count isa Real && isfinite(total_count) && total_count > 0) ||
        throw(ArgumentError("total_count must be finite and positive"))
    edge_values, observations = _histogram_data(edges, counts)
    bins = DistributionHistogram(edge_values, observations,
        total_count === nothing ? nothing : Float64(total_count), integration, Float64(rtol))
    objective = DistributionObjective(make_distribution, bins)
    gof = function (p)
        log_mu = _bin_logexpectation(make_distribution(p), bins)
        return _poisson_gof(observations, exp.(log_mu); log_mu)
    end
    return fit_custom(objective; p0, nobs=length(observations),
        cost_name=:histogram_poisson_likelihood, gof, kwargs...)
end
