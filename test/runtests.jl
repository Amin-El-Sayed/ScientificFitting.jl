using Test

@testset "JuFitter test suite" begin
    # Keep one authoritative Makie-free inventory. Pkg.test adds only the
    # optional plotting slice so new core references cannot be omitted here.
    include("core_runtests.jl")
    include("plots/fitplot.jl")
end
