using Test

@testset "JuFitter test suite" begin
    include("regression/current_api.jl")
    include("plots/fitplot.jl")
end
