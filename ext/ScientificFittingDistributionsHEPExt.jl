module ScientificFittingDistributionsHEPExt

using ScientificFitting
using Distributions
using DistributionsHEP: ExtendedMixtureModel, yields, total_yield

function ScientificFitting._distribution_cost(model::ExtendedMixtureModel, data)
    parts, weights = components(model), yields(model)
    total = total_yield(model)
    all(w -> isfinite(w) && w >= 0, weights) && isfinite(total) && total > 0 ||
        throw(ArgumentError("extended mixture yields must be finite and nonnegative, with positive total"))
    for part in parts
        ScientificFitting._distribution_shape(part, data)
    end
    # Reuse Distributions' stable log-sum-exp mixture evaluation. Restore the
    # intensity scale explicitly: lambda(x) = total * normalized_pdf(x).
    n = data isa AbstractMatrix ? size(data, 2) : length(data)
    shape_cost = ScientificFitting._distribution_cost(MixtureModel(model), data)
    return shape_cost + 2 * (total - n * log(total))
end

end
