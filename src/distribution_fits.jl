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
    all(a === b for (a, b) in zip(original, prepared)) && return d
    return MixtureModel(prepared, probs(d))
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

function _distribution_loglikelihood(d, data)
    events = data isa AbstractMatrix ? eachcol(data) : data
    # Prepare once per parameter point, not once per event. Upstream batching stays intact.
    return loglikelihood(_with_logpdf(d, first(events)), data)
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
minimizes `-2 log(L)`. A distribution exposing only `pdf` uses `log(pdf)`;
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
