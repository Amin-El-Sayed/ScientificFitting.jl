"""
    ConstraintSpec(; ineq=nothing, eq=nothing)

Container for nonlinear parameter constraints passed to the general
`Optimization.jl` backend. `ineq` and `eq` are user functions evaluated on the
free parameter vector and interpreted as inequality and equality constraints.
Use this only when simple bounds, fixed parameters, or Gaussian parameter
constraints are not expressive enough.
"""
struct ConstraintSpec{FI, FE}
    ineq::FI
    eq::FE
end

ConstraintSpec(; ineq=nothing, eq=nothing) = ConstraintSpec(ineq, eq)

has_constraints(spec::ConstraintSpec) = !(spec.ineq === nothing && spec.eq === nothing)

"""
    ParameterPrior(index, mean, sigma)
    ParameterPrior(index, mean, sigma_minus, sigma_plus)

Gaussian prior term for one fitted parameter. The prior contributes a
chi-square-like penalty centered at `mean`; asymmetric uncertainties use
`sigma_minus` below the mean and `sigma_plus` above the mean.
"""
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

"""
    FixedParameter(index, value[, sigma])
    FixedParameter(index, value, sigma_minus, sigma_plus)

Fix one parameter to `value` during the fit. Optional uncertainties describe
the externally known value for reports and downstream uncertainty accounting;
they do not make the parameter free again.
"""
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

"""
    ParameterConstraint(indices, mean, covariance)

Correlated Gaussian constraint on several parameters. The contribution is
formed from the selected parameter vector, the supplied `mean`, and the
positive-definite covariance matrix. This is the parameter-space analogue of a
correlated measurement.
"""
struct ParameterConstraint
    indices::Vector{Int}
    mean::Vector{Float64}
    covariance::Matrix{Float64}
end

"""
    ErrorComponent(name, target, mode, values; active=true)

Named uncertainty contribution used by component-based covariance models.
`target` identifies the affected variable, `mode` describes how `values` are
interpreted, and `active=false` keeps the component documented but excluded
from the current fit.
"""
struct ErrorComponent
    name::Symbol
    target::Symbol
    mode::Symbol
    values::Any
    active::Bool
end

ErrorComponent(name::Symbol, target::Symbol, mode::Symbol, values; active::Bool=true) =
    ErrorComponent(name, target, mode, values, active)

"""
    FitOptions(; backend=:auto, cost=:auto, maxiters=500, tol=1e-10,
                ci_level=0.6827, scale_covariance=:auto, multistart=1)

Normalized solver and reporting options stored in a `FitResult`. User-facing
fit functions expose these as keyword arguments; constructing `FitOptions`
directly is mainly useful for lower-level workflows and tests. Invalid backend,
iteration, tolerance, confidence-level, covariance-scaling, and multistart
settings fail during construction rather than inside a solver.
"""
Base.@kwdef struct FitOptions
    backend::Symbol = :auto
    cost::Symbol = :auto
    maxiters::Int = 500
    tol::Float64 = 1e-10
    ci_level::Float64 = 0.6827
    scale_covariance::Symbol = :auto
    multistart::Int = 1

    function FitOptions(
        backend::Symbol,
        cost::Symbol,
        maxiters::Integer,
        tol::Real,
        ci_level::Real,
        scale_covariance::Symbol,
        multistart::Integer,
    )
        maxiters_value = Int(maxiters)
        tol_value = Float64(tol)
        ci_level_value = Float64(ci_level)
        multistart_value = Int(multistart)
        backend in (:auto, :lsqfit, :optimization) || throw(ArgumentError(
            "backend must be :auto, :lsqfit, or :optimization",
        ))
        scale_covariance in (:auto, :always, :never) || throw(ArgumentError(
            "scale_covariance must be :auto, :always, or :never",
        ))
        maxiters_value > 0 || throw(ArgumentError("maxiters must be > 0"))
        isfinite(tol_value) && tol_value > 0 || throw(ArgumentError(
            "tol must be finite and > 0",
        ))
        isfinite(ci_level_value) && 0 < ci_level_value < 1 || throw(ArgumentError(
            "ci_level must be finite and strictly between 0 and 1",
        ))
        multistart_value > 0 || throw(ArgumentError("multistart must be >= 1"))
        return new(
            backend,
            cost,
            maxiters_value,
            tol_value,
            ci_level_value,
            scale_covariance,
            multistart_value,
        )
    end
end

const CovarianceLike = Union{Matrix{Float64}, SparseMatrixCSC{Float64, Int}}

"""
    WhiteningOperator(whiten!; logdet_covariance, marginal_sigma=nothing)

Represent a complete static data covariance without materializing it.
`whiten!(out, residual)` must apply a linear operator `W` satisfying
`W'W = inv(C)`, where `C` is the observation covariance. The required
`logdet_covariance` is `log(det(C))`; JuFitter uses it for the normalized
Gaussian likelihood, AIC, and BIC.

`marginal_sigma` is optional scalar or pointwise marginal standard deviation.
It does not change the fit; plotting uses it for data error bars and pointwise
prediction bands. Without it, confidence bands remain available but a
prediction band would not have enough information.

The mutating function must write every element of `out` and support the element
types used by automatic differentiation when the general optimizer is needed.
It must accept `AbstractVector` views because JuFitter applies the same operator
columnwise to analytic Jacobians.
"""
struct WhiteningOperator{F, MS}
    whiten!::F
    logdet_covariance::Float64
    marginal_sigma::MS

    function WhiteningOperator{F, MS}(
        whiten!::F,
        logdet_covariance::Float64,
        marginal_sigma::MS,
    ) where {F, MS}
        isfinite(logdet_covariance) || throw(ArgumentError(
            "logdet_covariance must be finite",
        ))
        marginal_sigma isa Union{Nothing, Float64, Vector{Float64}} || throw(ArgumentError(
            "marginal_sigma must be nothing, a scalar, or a vector",
        ))
        if marginal_sigma isa Float64
            isfinite(marginal_sigma) || throw(ArgumentError(
                "marginal_sigma must be finite",
            ))
            marginal_sigma > 0 || throw(ArgumentError(
                "marginal_sigma must be positive",
            ))
        elseif marginal_sigma isa Vector{Float64}
            all(isfinite, marginal_sigma) || throw(ArgumentError(
                "marginal_sigma must contain only finite values",
            ))
            all(>(0), marginal_sigma) || throw(ArgumentError(
                "marginal_sigma must contain only positive values",
            ))
        end
        return new{F, MS}(whiten!, logdet_covariance, marginal_sigma)
    end
end

function WhiteningOperator(whiten!; logdet_covariance, marginal_sigma=nothing)
    normalized_sigma = if marginal_sigma === nothing
        nothing
    elseif marginal_sigma isa Real
        Float64(marginal_sigma)
    elseif marginal_sigma isa AbstractVector
        collect(Float64, marginal_sigma)
    else
        throw(ArgumentError("marginal_sigma must be nothing, a scalar, or a vector"))
    end
    return WhiteningOperator{typeof(whiten!), typeof(normalized_sigma)}(
        whiten!,
        Float64(logdet_covariance),
        normalized_sigma,
    )
end

function (operator::WhiteningOperator)(out::AbstractVector, residual::AbstractVector)
    length(out) == length(residual) || throw(ArgumentError(
        "whitening input and output must have equal length",
    ))
    operator.whiten!(out, residual)
    return out
end

function _validate_whitening_operator(operator::WhiteningOperator, n::Int)
    input = collect(range(-0.75, 1.25; length=n))
    output = fill(NaN, n)
    applicable(operator.whiten!, output, input) || throw(ArgumentError(
        "whiten! must have the signature whiten!(out, residual)",
    ))
    operator(output, input)
    all(isfinite, output) || throw(ArgumentError(
        "whiten! must write a finite value to every output element",
    ))
    any(x -> !iszero(x), output) || throw(ArgumentError(
        "whiten! must not map a nonzero residual vector to zero",
    ))

    # Analytic Jacobians are whitened columnwise through views. Reject
    # Vector-only methods before solver dispatch instead of failing mid-fit.
    input_view = view(input, :)
    output_view = view(output, :)
    applicable(operator.whiten!, output_view, input_view) || throw(ArgumentError(
        "whiten! must accept AbstractVector views for Jacobian whitening",
    ))
    fill!(output_view, NaN)
    operator(output_view, input_view)
    all(isfinite, output_view) || throw(ArgumentError(
        "whiten! must write a finite value to every output element",
    ))

    marginal_sigma = operator.marginal_sigma
    if marginal_sigma isa AbstractVector && length(marginal_sigma) != n
        throw(ArgumentError("marginal_sigma length must match y"))
    end
    return nothing
end

# These adapters keep one model contract inside JuFitter while allowing the
# least-squares backend to call allocation-free user functions directly.
struct _InPlaceModel{F}
    f!::F
end

function (model::_InPlaceModel)(out::AbstractVector, x::AbstractVector, p::AbstractVector)
    model.f!(out, x, p)
    return out
end

function (model::_InPlaceModel)(x::AbstractVector, p::AbstractVector)
    T = promote_type(eltype(x), eltype(p))
    out = Vector{T}(undef, length(x))
    return model(out, x, p)
end

struct _InPlaceJacobian{F}
    f!::F
end

function (jacobian::_InPlaceJacobian)(out::AbstractMatrix, x::AbstractVector, p::AbstractVector)
    jacobian.f!(out, x, p)
    return out
end

function (jacobian::_InPlaceJacobian)(x::AbstractVector, p::AbstractVector)
    T = promote_type(eltype(x), eltype(p))
    out = Matrix{T}(undef, length(x), length(p))
    return jacobian(out, x, p)
end

function _validate_inplace_output!(f!, output, x, p, name::AbstractString)
    fill!(output, NaN)
    f!(output, x, p)
    all(isfinite, output) || throw(ArgumentError(
        "$name must write a finite value to every output element at p0",
    ))
    return nothing
end

struct FitProblem{TF, TW}
    model::TF
    x::Vector{Float64}
    y::Vector{Float64}
    p0::Vector{Float64}
    sigma_y::Union{Nothing, Vector{Float64}}
    sigma_x::Union{Nothing, Vector{Float64}}
    cov_y::Union{Nothing, CovarianceLike}
    cov_x::Union{Nothing, CovarianceLike}
    whitening::TW
    error_components::Vector{ErrorComponent}
    bounds::Union{Nothing, Tuple{Vector{Float64}, Vector{Float64}}}
    constraints::ConstraintSpec
    parameter_priors::Vector{ParameterPrior}
    parameter_constraints::Vector{ParameterConstraint}
    fixed_parameters::Vector{FixedParameter}
    jacobian::Any
    x_derivative::Any
end

"""
    FitStatistics

Goodness-of-fit and information-criterion summary for a fit. The fields include
the minimized objective, negative log-likelihood convention, chi-square-like
goodness of fit, degrees of freedom, p-value, AIC, and BIC. Fields that are not
meaningful for a given likelihood are set to `NaN`.
"""
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

"""
    DiagnosticFinding

One structured diagnostic issue or note. Each finding has a `severity`
(`:info`, `:warning`, or `:critical`), a stable machine-readable `code`,
reader-facing `title`, concrete `evidence`, and a recommended next action.
"""
struct DiagnosticFinding
    severity::Symbol
    code::Symbol
    title::String
    evidence::String
    recommendation::String
end

"""
    FitDiagnostics

Numerical and statistical diagnostics stored with every fit result. It contains
legacy warning strings, covariance/Hessian condition numbers, active-bound
indices, and structured `DiagnosticFinding`s used by `diagnose` and diagnostic
plots.
"""
struct FitDiagnostics
    warnings::Vector{String}
    covariance_condition::Float64
    hessian_condition::Float64
    active_bounds::Vector{Int}
    findings::Vector{DiagnosticFinding}
end

"""
    DiagnosticReport

Structured diagnostic report returned by `diagnose(result)`. It keeps the full
list of findings plus a short summary suitable for notebook output or a lab log.
Use `diagnose_text(report)` for a plain-text representation.
"""
struct DiagnosticReport
    findings::Vector{DiagnosticFinding}
    summary::String
end

"""
    DiagnosticDashboard

Compact, action-oriented summary returned by `diagnostic_dashboard(...)`. It
groups diagnostic findings into an overall status, severity counts, and
deduplicated next actions for quick interactive use.
"""
struct DiagnosticDashboard
    report::DiagnosticReport
    status::Symbol
    severity_counts::Dict{Symbol, Int}
    next_actions::Vector{String}
end

"""
    FitResult

Result of a Gaussian/least-squares style fit. It stores the normalized
`FitProblem`, solver options and status, fitted parameters, local covariance
and correlation estimates, fitted model values, residuals, Jacobian,
statistics, and diagnostics.

The parameter covariance is a local quadratic approximation. For nonlinear
models, active bounds, weak data, or asymmetric likelihoods, inspect
`profile(...)` or `contour(...)` before treating symmetric errors as final.
"""
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

function _assert_finite_vector(name::AbstractString, values::AbstractVector)
    all(isfinite, values) || throw(ArgumentError("$name must contain only finite values"))
    return nothing
end

function _assert_finite_matrix(name::AbstractString, values::AbstractMatrix)
    all(isfinite, values) || throw(ArgumentError("$name must contain only finite values"))
    return nothing
end

function _assert_finite_matrix(name::AbstractString, values::SparseMatrixCSC)
    all(isfinite, nonzeros(values)) || throw(ArgumentError("$name must contain only finite values"))
    return nothing
end

function _assert_positive_sigma(name::AbstractString, values::AbstractVector)
    _assert_finite_vector(name, values)
    all(>(0.0), values) || throw(ArgumentError("$name entries must be > 0"))
    return nothing
end

function _assert_covariance_matrix(name::AbstractString, cov::AbstractMatrix)
    dense = Matrix{Float64}(cov)
    _assert_finite_matrix(name, dense)
    isapprox(dense, dense'; rtol=1e-12, atol=1e-12) ||
        throw(ArgumentError("$name must be symmetric"))
    try
        cholesky(Symmetric(dense))
    catch
        throw(ArgumentError("$name must be symmetric positive definite"))
    end
    return nothing
end

function _assert_covariance_matrix(name::AbstractString, cov::SparseMatrixCSC{Float64, Int})
    _assert_finite_matrix(name, cov)
    isapprox(cov, cov'; rtol=1e-12, atol=1e-12) ||
        throw(ArgumentError("$name must be symmetric"))
    try
        cholesky(Symmetric(cov))
    catch
        throw(ArgumentError("$name must be symmetric positive definite"))
    end
    return nothing
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
    all(isinf.(lower) .& (lower .< 0.0)) && all(isinf.(upper) .& (upper .> 0.0)) && return nothing
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
        isfinite(prior.mean) || throw(ArgumentError("parameter prior mean must be finite"))
        isfinite(prior.sigma_minus) && prior.sigma_minus > 0 ||
            throw(ArgumentError("parameter prior sigma_minus must be finite and > 0"))
        isfinite(prior.sigma_plus) && prior.sigma_plus > 0 ||
            throw(ArgumentError("parameter prior sigma_plus must be finite and > 0"))
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
        _assert_finite_vector("parameter constraint mean", constraint.mean)
        _assert_covariance_matrix("parameter constraint covariance", constraint.covariance)
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
            _assert_finite_vector("error component $(component.name)", component.values)
        elseif component.values isa AbstractMatrix
            size(component.values) == (nobs, nobs) || throw(ArgumentError("error component covariance must be n x n"))
            _assert_finite_matrix("error component $(component.name)", Matrix(component.values))
        elseif component.values isa Number
            isfinite(component.values) || throw(ArgumentError("error component $(component.name) must be finite"))
        end

        component.mode == :covariance && component.values isa Number &&
            throw(ArgumentError("covariance error components require a vector or matrix"))
        component.mode != :covariance && component.values isa AbstractMatrix &&
            throw(ArgumentError("non-covariance error components require a scalar or vector"))
        if component.mode == :covariance && component.values isa AbstractVector
            all(>(0.0), component.values) || throw(ArgumentError("covariance error component vector entries must be > 0"))
        elseif component.mode != :covariance
            if component.values isa Number
                component.values > 0.0 || throw(ArgumentError("error component sigma entries must be > 0"))
            else
                all(>(0.0), component.values) || throw(ArgumentError("error component sigma entries must be > 0"))
            end
        end

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
        isfinite(fp.value) || throw(ArgumentError("fixed parameter value must be finite"))
        isfinite(fp.sigma_minus) && fp.sigma_minus >= 0 ||
            throw(ArgumentError("fixed parameter sigma_minus must be finite and >= 0"))
        isfinite(fp.sigma_plus) && fp.sigma_plus >= 0 ||
            throw(ArgumentError("fixed parameter sigma_plus must be finite and >= 0"))
        fp.index in seen && throw(ArgumentError("fixed parameter index $(fp.index) appears more than once"))
        push!(seen, fp.index)
        push!(fixed, fp)
    end

    return fixed
end

function _assert_fixed_parameters_within_bounds(fixed::Vector{FixedParameter}, bounds)
    bounds === nothing && return nothing
    lower, upper = bounds
    for fp in fixed
        value = fp.value
        lo = lower[fp.index]
        hi = upper[fp.index]
        if value < lo || value > hi
            throw(ArgumentError("fixed parameter $(fp.index) value must satisfy its bounds"))
        end
    end
    return nothing
end

"""
    FitProblem(model, x, y; p0, sigma_y, sigma_x, cov_y, cov_x, whitening,
               error_components, bounds, constraints, parameter_priors,
               parameter_constraints, fixed_parameters, jacobian,
               x_derivative, inplace=false)

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

`x_derivative(x, p)` optionally supplies the vector derivative `dy/dx` used for
effective x-uncertainty propagation. If omitted, JuFitter differentiates the
model with respect to each x value by automatic differentiation.

`whitening=WhiteningOperator(...)` supplies the complete static observation
covariance through a matrix-free whitening operation. It is mutually exclusive
with y/x uncertainties and active `error_components`; combining covariance
models without an explicit derivation would change the statistical model.

Set `inplace=true` when the model has the signature `model!(out, x, p)`. On the
unbounded least-squares path, JuFitter forwards this contract to LsqFit's native
in-place solver interface. Other solver paths use the same model through a
type-preserving output buffer. If `jacobian` is also supplied, its in-place
signature must be `jacobian!(J, x, p)`. Mutating functions used with bounds,
constraints, or parameter-dependent covariance must accept buffers whose
element type is chosen by automatic differentiation; avoid `Float64`-specific
method signatures.
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
    whitening=nothing,
    error_components=nothing,
    bounds=nothing,
    constraints=nothing,
    parameter_priors=nothing,
    parameter_constraints=nothing,
    fixed_parameters=nothing,
    jacobian=nothing,
    x_derivative=nothing,
    inplace::Bool=false,
)
    x_vec = _float_vector(x)
    y_vec = _float_vector(y)
    p0_vec = _float_vector(p0)

    n = length(y_vec)
    n > 0 || throw(ArgumentError("x and y must contain at least one observation"))
    length(x_vec) == n || throw(ArgumentError("x and y must have equal length"))
    length(p0_vec) > 0 || throw(ArgumentError("p0 must contain at least one parameter"))
    _assert_finite_vector("x", x_vec)
    _assert_finite_vector("y", y_vec)
    _assert_finite_vector("p0", p0_vec)

    if inplace
        model_output = similar(y_vec)
        applicable(model, model_output, x_vec, p0_vec) || throw(ArgumentError(
            "inplace=true requires model!(out, x, p)",
        ))
        _validate_inplace_output!(model, model_output, x_vec, p0_vec, "model!")
        if jacobian !== nothing
            J = Matrix{Float64}(undef, n, length(p0_vec))
            applicable(jacobian, J, x_vec, p0_vec) || throw(ArgumentError(
                "inplace=true requires jacobian!(J, x, p)",
            ))
            _validate_inplace_output!(jacobian, J, x_vec, p0_vec, "jacobian!")
        end
    end

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

    if whitening !== nothing
        whitening isa WhiteningOperator || throw(ArgumentError(
            "whitening must be a WhiteningOperator",
        ))
        if sigma_y !== nothing || cov_y !== nothing || sigma_x !== nothing || cov_x !== nothing ||
           any(component -> component.active, components)
            throw(ArgumentError(
                "whitening describes the complete data covariance and cannot be combined " *
                "with sigma_y, cov_y, sigma_x, cov_x, or active error_components",
            ))
        end
        _validate_whitening_operator(whitening, n)
    end

    sigma_y_vec !== nothing && length(sigma_y_vec) != n && throw(ArgumentError("sigma_y length must match y"))
    sigma_x_vec !== nothing && length(sigma_x_vec) != n && throw(ArgumentError("sigma_x length must match x"))
    sigma_y_vec !== nothing && _assert_positive_sigma("sigma_y", sigma_y_vec)
    sigma_x_vec !== nothing && _assert_positive_sigma("sigma_x", sigma_x_vec)

    if cov_y_mat !== nothing
        size(cov_y_mat) == (n, n) || throw(ArgumentError("cov_y must be n x n"))
        _assert_covariance_matrix("cov_y", cov_y_mat)
    end
    if cov_x_mat !== nothing
        size(cov_x_mat) == (n, n) || throw(ArgumentError("cov_x must be n x n"))
        _assert_covariance_matrix("cov_x", cov_x_mat)
    end

    bnd = _normalize_bounds(bounds, length(p0_vec))
    cons = _normalize_constraints(constraints)
    priors = _normalize_parameter_priors(parameter_priors, length(p0_vec))
    par_constraints = _normalize_parameter_constraints(parameter_constraints, length(p0_vec))
    fixed = _normalize_fixed_parameters(fixed_parameters, length(p0_vec))
    _assert_fixed_parameters_within_bounds(fixed, bnd)

    model_impl = inplace ? _InPlaceModel(model) : model
    jacobian_impl = inplace && jacobian !== nothing ? _InPlaceJacobian(jacobian) : jacobian

    return FitProblem(
        model_impl,
        x_vec,
        y_vec,
        p0_vec,
        sigma_y_vec,
        sigma_x_vec,
        cov_y_mat,
        cov_x_mat,
        whitening,
        components,
        bnd,
        cons,
        priors,
        par_constraints,
        fixed,
        jacobian_impl,
        x_derivative,
    )
end
