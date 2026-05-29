struct ConstraintSpec{FI, FE}
    ineq::FI
    eq::FE
end

ConstraintSpec(; ineq=nothing, eq=nothing) = ConstraintSpec(ineq, eq)

has_constraints(spec::ConstraintSpec) = !(spec.ineq === nothing && spec.eq === nothing)

struct ParameterPrior
    index::Int
    mean::Float64
    sigma_minus::Float64
    sigma_plus::Float64
end

ParameterPrior(index::Integer, mean::Real, sigma::Real) =
    ParameterPrior(Int(index), Float64(mean), Float64(sigma), Float64(sigma))

ParameterPrior(index::Integer, mean::Real, sigma_minus::Real, sigma_plus::Real) =
    ParameterPrior(Int(index), Float64(mean), Float64(sigma_minus), Float64(sigma_plus))

struct FixedParameter
    index::Int
    value::Float64
    sigma_minus::Float64
    sigma_plus::Float64
end

FixedParameter(index::Integer, value::Real) =
    FixedParameter(Int(index), Float64(value), 0.0, 0.0)

FixedParameter(index::Integer, value::Real, sigma::Real) =
    FixedParameter(Int(index), Float64(value), Float64(sigma), Float64(sigma))

FixedParameter(index::Integer, value::Real, sigma_minus::Real, sigma_plus::Real) =
    FixedParameter(Int(index), Float64(value), Float64(sigma_minus), Float64(sigma_plus))

struct ParameterConstraint
    indices::Vector{Int}
    mean::Vector{Float64}
    covariance::Matrix{Float64}
end

struct ErrorComponent
    name::Symbol
    target::Symbol
    mode::Symbol
    values::Any
    active::Bool
end

ErrorComponent(name::Symbol, target::Symbol, mode::Symbol, values; active::Bool=true) =
    ErrorComponent(name, target, mode, values, active)

Base.@kwdef struct FitOptions
    backend::Symbol = :auto
    cost::Symbol = :auto
    maxiters::Int = 500
    tol::Float64 = 1e-10
    ci_level::Float64 = 0.6827
    scale_covariance::Symbol = :auto
    multistart::Int = 1
end

const CovarianceLike = Union{Matrix{Float64}, SparseMatrixCSC{Float64, Int}}

struct FitProblem{TF}
    model::TF
    x::Vector{Float64}
    y::Vector{Float64}
    p0::Vector{Float64}
    sigma_y::Union{Nothing, Vector{Float64}}
    sigma_x::Union{Nothing, Vector{Float64}}
    cov_y::Union{Nothing, CovarianceLike}
    cov_x::Union{Nothing, CovarianceLike}
    error_components::Vector{ErrorComponent}
    bounds::Union{Nothing, Tuple{Vector{Float64}, Vector{Float64}}}
    constraints::ConstraintSpec
    parameter_priors::Vector{ParameterPrior}
    parameter_constraints::Vector{ParameterConstraint}
    fixed_parameters::Vector{FixedParameter}
    jacobian::Any
end

struct FitStatistics
    cost::Symbol
    cost_min::Float64
    nll_min::Float64
    chi2::Float64
    chi2_ndf::Float64
    ndf::Int
    pvalue::Float64
    aic::Float64
    bic::Float64
end

struct FitDiagnostics
    warnings::Vector{String}
    covariance_condition::Float64
    hessian_condition::Float64
    active_bounds::Vector{Int}
end

struct FitResult
    problem::FitProblem
    options::FitOptions
    backend::Symbol
    converged::Bool
    iterations::Int
    message::String
    params::Vector{Float64}
    param_stderr::Vector{Float64}
    param_covariance::Matrix{Float64}
    param_correlation::Matrix{Float64}
    model_y::Vector{Float64}
    residuals::Vector{Float64}
    weighted_residuals::Vector{Float64}
    jacobian::Matrix{Float64}
    stats::FitStatistics
    diagnostics::FitDiagnostics
end

function _float_vector(v::AbstractVector)
    return collect(Float64, v)
end

function _float_matrix(m::AbstractMatrix)
    if m isa SparseMatrixCSC
        return sparse(Float64.(m))
    end
    return Matrix{Float64}(m)
end

function _normalize_bounds(bounds, nparams::Int)
    bounds === nothing && return nothing
    if !(bounds isa Tuple) || length(bounds) != 2
        throw(ArgumentError("bounds must be a tuple (lower, upper)"))
    end
    lower = _float_vector(bounds[1])
    upper = _float_vector(bounds[2])
    length(lower) == nparams || throw(ArgumentError("lower bounds length must equal parameter count"))
    length(upper) == nparams || throw(ArgumentError("upper bounds length must equal parameter count"))
    any(lower .> upper) && throw(ArgumentError("each lower bound must be <= corresponding upper bound"))
    return (lower, upper)
end

function _normalize_constraints(constraints)
    constraints === nothing && return ConstraintSpec()
    constraints isa ConstraintSpec && return constraints

    if constraints isa NamedTuple
        ineq = haskey(constraints, :ineq) ? constraints[:ineq] : nothing
        eq = haskey(constraints, :eq) ? constraints[:eq] : nothing
        return ConstraintSpec(; ineq=ineq, eq=eq)
    end

    throw(ArgumentError("constraints must be nothing, ConstraintSpec, or NamedTuple with optional :ineq/:eq"))
end

function _normalize_parameter_priors(parameter_priors, nparams::Int)
    parameter_priors === nothing && return ParameterPrior[]

    raw_priors = parameter_priors isa AbstractVector ? parameter_priors : [parameter_priors]
    priors = ParameterPrior[]

    for raw in raw_priors
        prior = if raw isa ParameterPrior
            raw
        elseif raw isa NamedTuple
            haskey(raw, :index) || throw(ArgumentError("parameter prior NamedTuple must include :index"))
            haskey(raw, :mean) || throw(ArgumentError("parameter prior NamedTuple must include :mean"))
            if haskey(raw, :sigma)
                ParameterPrior(Int(raw.index), Float64(raw.mean), Float64(raw.sigma))
            elseif haskey(raw, :sigma_minus) && haskey(raw, :sigma_plus)
                ParameterPrior(Int(raw.index), Float64(raw.mean), Float64(raw.sigma_minus), Float64(raw.sigma_plus))
            else
                throw(ArgumentError("parameter prior NamedTuple must include :sigma or both :sigma_minus/:sigma_plus"))
            end
        else
            throw(ArgumentError("parameter_priors entries must be ParameterPrior or NamedTuple(index, mean, sigma)"))
        end

        1 <= prior.index <= nparams || throw(ArgumentError("parameter prior index $(prior.index) is out of range 1:$nparams"))
        prior.sigma_minus > 0 || throw(ArgumentError("parameter prior sigma_minus must be > 0"))
        prior.sigma_plus > 0 || throw(ArgumentError("parameter prior sigma_plus must be > 0"))
        push!(priors, prior)
    end

    return priors
end

function _normalize_parameter_constraints(parameter_constraints, nparams::Int)
    parameter_constraints === nothing && return ParameterConstraint[]

    raw_constraints = parameter_constraints isa AbstractVector ? parameter_constraints : [parameter_constraints]
    constraints = ParameterConstraint[]

    for raw in raw_constraints
        constraint = if raw isa ParameterConstraint
            raw
        elseif raw isa NamedTuple
            haskey(raw, :indices) || throw(ArgumentError("parameter constraint NamedTuple must include :indices"))
            haskey(raw, :mean) || throw(ArgumentError("parameter constraint NamedTuple must include :mean"))
            haskey(raw, :covariance) || throw(ArgumentError("parameter constraint NamedTuple must include :covariance"))
            ParameterConstraint(
                collect(Int, raw.indices),
                collect(Float64, raw.mean),
                Matrix{Float64}(raw.covariance),
            )
        else
            throw(ArgumentError("parameter_constraints entries must be ParameterConstraint or NamedTuple(indices, mean, covariance)"))
        end

        k = length(constraint.indices)
        k > 0 || throw(ArgumentError("parameter constraint must contain at least one index"))
        length(constraint.mean) == k || throw(ArgumentError("parameter constraint mean length must match indices"))
        size(constraint.covariance) == (k, k) || throw(ArgumentError("parameter constraint covariance must be k x k"))
        all(1 .<= constraint.indices .<= nparams) || throw(ArgumentError("parameter constraint index out of range 1:$nparams"))
        length(unique(constraint.indices)) == k || throw(ArgumentError("parameter constraint indices must be unique"))
        try
            cholesky(Symmetric(constraint.covariance))
        catch
            throw(ArgumentError("parameter constraint covariance must be symmetric positive definite"))
        end
        push!(constraints, constraint)
    end

    return constraints
end

function _normalize_error_component_values(values)
    if values isa Number
        return Float64(values)
    elseif values isa AbstractVector
        return _float_vector(values)
    elseif values isa AbstractMatrix
        return _float_matrix(values)
    end
    throw(ArgumentError("error component values must be a number, vector, or matrix"))
end

function _normalize_error_components(error_components, nobs::Int)
    error_components === nothing && return ErrorComponent[]

    raw_components = error_components isa AbstractVector ? error_components : [error_components]
    components = ErrorComponent[]

    for raw in raw_components
        component = if raw isa ErrorComponent
            raw
        elseif raw isa NamedTuple
            haskey(raw, :name) || throw(ArgumentError("error component NamedTuple must include :name"))
            haskey(raw, :target) || throw(ArgumentError("error component NamedTuple must include :target"))
            haskey(raw, :mode) || throw(ArgumentError("error component NamedTuple must include :mode"))
            haskey(raw, :values) || throw(ArgumentError("error component NamedTuple must include :values"))
            active = haskey(raw, :active) ? Bool(raw.active) : true
            ErrorComponent(Symbol(raw.name), Symbol(raw.target), Symbol(raw.mode), _normalize_error_component_values(raw.values); active=active)
        else
            throw(ArgumentError("error_components entries must be ErrorComponent or NamedTuple(name, target, mode, values)"))
        end

        component.target in (:y, :x) || throw(ArgumentError("error component target must be :y or :x"))
        component.mode in (:absolute, :relative, :model_relative, :covariance) ||
            throw(ArgumentError("error component mode must be :absolute, :relative, :model_relative, or :covariance"))
        component.target == :x && component.mode == :model_relative &&
            throw(ArgumentError("x error components do not support :model_relative mode"))

        if component.values isa AbstractVector
            length(component.values) == nobs || throw(ArgumentError("error component vector length must match data length"))
        elseif component.values isa AbstractMatrix
            size(component.values) == (nobs, nobs) || throw(ArgumentError("error component covariance must be n x n"))
        end

        component.mode == :covariance && component.values isa Number &&
            throw(ArgumentError("covariance error components require a vector or matrix"))
        component.mode != :covariance && component.values isa AbstractMatrix &&
            throw(ArgumentError("non-covariance error components require a scalar or vector"))

        push!(components, component)
    end

    return components
end

function _normalize_fixed_parameters(fixed_parameters, nparams::Int)
    fixed_parameters === nothing && return FixedParameter[]

    raw_fixed = fixed_parameters isa AbstractVector ? fixed_parameters : [fixed_parameters]
    fixed = FixedParameter[]
    seen = Set{Int}()

    for raw in raw_fixed
        fp = if raw isa FixedParameter
            raw
        elseif raw isa Pair
            FixedParameter(Int(raw.first), Float64(raw.second))
        elseif raw isa NamedTuple
            haskey(raw, :index) || throw(ArgumentError("fixed parameter NamedTuple must include :index"))
            haskey(raw, :value) || throw(ArgumentError("fixed parameter NamedTuple must include :value"))
            if haskey(raw, :sigma)
                FixedParameter(Int(raw.index), Float64(raw.value), Float64(raw.sigma))
            elseif haskey(raw, :sigma_minus) && haskey(raw, :sigma_plus)
                FixedParameter(Int(raw.index), Float64(raw.value), Float64(raw.sigma_minus), Float64(raw.sigma_plus))
            else
                FixedParameter(Int(raw.index), Float64(raw.value))
            end
        else
            throw(ArgumentError("fixed_parameters entries must be FixedParameter, Pair, or NamedTuple(index, value, ...)"))
        end

        1 <= fp.index <= nparams || throw(ArgumentError("fixed parameter index $(fp.index) is out of range 1:$nparams"))
        fp.sigma_minus >= 0 || throw(ArgumentError("fixed parameter sigma_minus must be >= 0"))
        fp.sigma_plus >= 0 || throw(ArgumentError("fixed parameter sigma_plus must be >= 0"))
        fp.index in seen && throw(ArgumentError("fixed parameter index $(fp.index) appears more than once"))
        push!(seen, fp.index)
        push!(fixed, fp)
    end

    return fixed
end

"""
    FitProblem(model, x, y; p0, sigma_y, sigma_x, cov_y, cov_x, error_components, bounds, constraints, parameter_priors, parameter_constraints, fixed_parameters, jacobian)

Build a fit problem for 1D `x` and scalar `y` observations.

`constraints` accepts either `ConstraintSpec` or a NamedTuple with:
- `ineq = p -> vector` interpreted as `ineq(p) <= 0`
- `eq = p -> vector` interpreted as `eq(p) == 0`

`parameter_priors` adds Gaussian penalties in parameter space and accepts a
single NamedTuple or vector of NamedTuples:
- `(index=i, mean=mu, sigma=sigma)`
- `(index=i, mean=mu, sigma_minus=sminus, sigma_plus=splus)`

`parameter_constraints` adds correlated Gaussian constraints:
- `(indices=[i, j], mean=[mu_i, mu_j], covariance=cov)`

`error_components` adds named y/x uncertainty sources:
- `(name=:stat, target=:y, mode=:absolute, values=sigma)`
- `(name=:scale, target=:y, mode=:relative, values=0.02)`
- `(name=:model_scale, target=:y, mode=:model_relative, values=0.02)`
- `(name=:corr, target=:y, mode=:covariance, values=cov)`

`fixed_parameters` removes parameters from the optimizer and accepts:
- `i => value`
- `(index=i, value=value)`
- `(index=i, value=value, sigma=sigma)`
- `(index=i, value=value, sigma_minus=sminus, sigma_plus=splus)`
"""
function FitProblem(
    model,
    x::AbstractVector,
    y::AbstractVector;
    p0::AbstractVector,
    sigma_y=nothing,
    sigma_x=nothing,
    cov_y=nothing,
    cov_x=nothing,
    error_components=nothing,
    bounds=nothing,
    constraints=nothing,
    parameter_priors=nothing,
    parameter_constraints=nothing,
    fixed_parameters=nothing,
    jacobian=nothing,
)
    x_vec = _float_vector(x)
    y_vec = _float_vector(y)
    p0_vec = _float_vector(p0)

    n = length(y_vec)
    length(x_vec) == n || throw(ArgumentError("x and y must have equal length"))
    length(p0_vec) > 0 || throw(ArgumentError("p0 must contain at least one parameter"))

    if sigma_y !== nothing && cov_y !== nothing
        throw(ArgumentError("use either sigma_y or cov_y, not both"))
    end
    if sigma_x !== nothing && cov_x !== nothing
        throw(ArgumentError("use either sigma_x or cov_x, not both"))
    end

    sigma_y_vec = sigma_y === nothing ? nothing : _float_vector(sigma_y)
    sigma_x_vec = sigma_x === nothing ? nothing : _float_vector(sigma_x)
    cov_y_mat = cov_y === nothing ? nothing : _float_matrix(cov_y)
    cov_x_mat = cov_x === nothing ? nothing : _float_matrix(cov_x)
    components = _normalize_error_components(error_components, n)

    sigma_y_vec !== nothing && length(sigma_y_vec) != n && throw(ArgumentError("sigma_y length must match y"))
    sigma_x_vec !== nothing && length(sigma_x_vec) != n && throw(ArgumentError("sigma_x length must match x"))

    if cov_y_mat !== nothing
        size(cov_y_mat) == (n, n) || throw(ArgumentError("cov_y must be n x n"))
    end
    if cov_x_mat !== nothing
        size(cov_x_mat) == (n, n) || throw(ArgumentError("cov_x must be n x n"))
    end

    bnd = _normalize_bounds(bounds, length(p0_vec))
    cons = _normalize_constraints(constraints)
    priors = _normalize_parameter_priors(parameter_priors, length(p0_vec))
    par_constraints = _normalize_parameter_constraints(parameter_constraints, length(p0_vec))
    fixed = _normalize_fixed_parameters(fixed_parameters, length(p0_vec))

    return FitProblem(
        model,
        x_vec,
        y_vec,
        p0_vec,
        sigma_y_vec,
        sigma_x_vec,
        cov_y_mat,
        cov_x_mat,
        components,
        bnd,
        cons,
        priors,
        par_constraints,
        fixed,
        jacobian,
    )
end
