module ScientificFittingDistributionsHEPExt

using ScientificFitting
using Distributions
using DistributionsHEP: ExtendedMixtureModel, yields, total_yield

function validated_yields(model)
    weights = yields(model)
    total = total_yield(model)
    all(w -> isfinite(w) && w >= 0, weights) && isfinite(total) && total > 0 ||
        throw(ArgumentError("extended mixture yields must be finite and nonnegative, with positive total"))
    return weights, total
end

function ScientificFitting._distribution_cost(model::ExtendedMixtureModel, data)
    parts = components(model)
    weights, total = validated_yields(model)
    for part in parts
        ScientificFitting._distribution_shape(part, data)
    end
    # Reuse Distributions' stable log-sum-exp mixture evaluation. Restore the
    # intensity scale explicitly: lambda(x) = total * normalized_pdf(x).
    n = data isa AbstractMatrix ? size(data, 2) : length(data)
    shape_cost = ScientificFitting._distribution_cost(MixtureModel(model), data)
    return shape_cost + 2 * (total - n * log(total))
end

function ScientificFitting._bin_logexpectation(model::ExtendedMixtureModel, bins::ScientificFitting.DistributionHistogram)
    bins.total_count === nothing || throw(ArgumentError("extended mixtures already specify event yields; omit total_count"))
    parts, (weights, _) = components(model), validated_yields(model)
    all(part -> part isa UnivariateDistribution, parts) ||
        throw(ArgumentError("one-dimensional histograms require univariate mixture components"))
    return ScientificFitting._weighted_logbatches(
        part -> ScientificFitting._bin_logmasses(part, bins.edges, bins),
        parts, weights, length(bins.counts))
end

function ScientificFitting._distribution_cost(model::ExtendedMixtureModel, bins::ScientificFitting.DistributionHistogram)
    log_mu = ScientificFitting._bin_logexpectation(model, bins)
    mu = ScientificFitting._nonnegative_expectation(exp.(log_mu), length(bins.counts))
    return ScientificFitting._poisson_minus2loglik_terms(bins.counts, mu; log_mu)
end

end
