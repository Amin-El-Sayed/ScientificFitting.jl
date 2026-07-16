using Test

@testset "JuFitter core test suite" begin
    include("regression/current_api.jl")
    include("statistics/covariance_semantics_reference.jl")
    include("statistics/diagnostics_reference.jl")
    include("statistics/linear_gaussian_reference.jl")
    include("statistics/structured_whitening_reference.jl")
    include("statistics/likelihood_reference.jl")
    include("statistics/profile_contour_reference.jl")
    include("numerics/inplace_model_reference.jl")
    include("numerics/torture_inputs.jl")
end
