using Test

@testset "ScientificFitting statistical reference suite" begin
    include("covariance_semantics_reference.jl")
    include("diagnostics_reference.jl")
    include("linear_gaussian_reference.jl")
    include("likelihood_reference.jl")
    include("profile_contour_reference.jl")
end
