using Test

# Keep one authoritative Makie-free inventory. Pkg.test adds only the
# optional plotting slice so new core references cannot be omitted here.
# SCIENTIFICFITTING_TEST_SLICE="plots" restricts a run to the plotting slice;
# CI uses it to check the Pkg.test sandbox and Makie compatibility without
# repeating the core suite, which its own jobs already run per shard. Unset,
# Pkg.test runs everything.
const _SLICE = strip(get(ENV, "SCIENTIFICFITTING_TEST_SLICE", ""))
isempty(_SLICE) || _SLICE == "plots" ||
    error("SCIENTIFICFITTING_TEST_SLICE must be empty or \"plots\"")

@testset "ScientificFitting test suite" begin
    isempty(_SLICE) && include("core_runtests.jl")
    include("plots/fitplot.jl")
    include("plots/legend_theme.jl")
end
