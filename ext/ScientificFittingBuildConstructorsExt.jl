module ScientificFittingBuildConstructorsExt

using ScientificFitting
import BuildConstructors as BC

"""Build once per parameter point, preserving the parameter element type for AD."""
struct ConstructorFactory{C, N}
    constructor::C
end
(f::ConstructorFactory{C, N})(p) where {C, N} = BC.build_model(f.constructor, NamedTuple{N}(Tuple(p)))

"""Read the public metadata contract once, validating shared names before deduplication."""
function constructor_inputs(constructor; p0)
    BC.validate_parameters(constructor)
    metadata = BC.parameter_metadata(constructor)
    names = BC.parameter_names(metadata)
    isempty(names) && throw(ArgumentError("constructor has no named parameters to fit"))
    defaults = BC.parameter_values(metadata)
    starts = if p0 === nothing
        collect(values(defaults))
    elseif p0 isa NamedTuple
        all(k -> k in names, keys(p0)) || throw(ArgumentError("p0 contains unknown constructor parameter names"))
        collect(values(merge(defaults, p0)))
    elseif p0 isa AbstractVector
        collect(p0)
    else
        throw(ArgumentError("p0 must be a named tuple or a vector in constructor parameter order"))
    end
    length(starts) == length(names) || throw(ArgumentError("p0 length must match constructor parameter count"))
    for (name, value) in zip(names, starts)
        value isa Real && isfinite(value) || throw(ArgumentError(
            "parameter $name needs a finite starting value; supply it through p0",
        ))
    end
    fixed_values = BC.fixed_values(metadata)
    fixed = FixedParameter[]
    for (i, name) in enumerate(names)
        if haskey(fixed_values, name)
            starts[i] == fixed_values[name] || throw(ArgumentError(
                "p0 cannot change fixed constructor parameter $name; update the constructor explicitly",
            ))
            push!(fixed, FixedParameter(i, starts[i]))
        end
    end
    bounds = (collect(values(BC.parameter_lower_boundaries(metadata))),
              collect(values(BC.parameter_upper_boundaries(metadata))))
    # All numeric values must come from p, even for originally fixed descriptors.
    # The shared core enforces fixing; release only our private constructor copy.
    snapshot = deepcopy(constructor)
    BC.release!(snapshot, names)
    factory = ConstructorFactory{typeof(snapshot), names}(snapshot)
    return (; factory, p0=Float64.(starts), bounds,
            fixed_parameters=fixed, parameter_names=String.(collect(names)))
end

ScientificFitting.fit_distribution(constructor::BC.AbstractConstructor, data::AbstractArray{<:Real}; kwargs...) =
    fit_constructor(constructor, data; kwargs...)
ScientificFitting.fit_distribution(constructor::BC.AbstractConstructor, edges::AbstractVector, counts::AbstractVector; kwargs...) =
    fit_constructor(constructor, edges, counts; kwargs...)

function fit_constructor(constructor, data...; p0=nothing, kwargs...)
    any(k -> k in (:bounds, :fixed_parameters, :parameter_names), keys(kwargs)) &&
        throw(ArgumentError("set names, bounds and fixed state on the constructor, not as duplicate fit keywords"))
    inputs = constructor_inputs(constructor; p0)
    return fit_distribution(inputs.factory, data...; p0=inputs.p0, bounds=inputs.bounds,
        fixed_parameters=inputs.fixed_parameters, parameter_names=inputs.parameter_names, kwargs...)
end

"""
    BuildConstructors.parameter_values(result::LikelihoodFitResult)

Return fitted values as a named tuple, including fixed values. Requires stored
parameter names. Pass this tuple to the explicit `update!` operation to change
constructor defaults; obtaining it does not mutate any object.
"""
function BC.parameter_values(result::LikelihoodFitResult)
    names = result.problem.parameter_names
    names === nothing && throw(ArgumentError("result has no stored parameter names"))
    return NamedTuple{Tuple(Symbol.(names))}(Tuple(result.params))
end

end
