using Test

# Every entry carries a CI shard tag; runtimes are balanced so the three
# shards finish together. SCIENTIFICFITTING_TEST_SHARD selects one shard
# ("1", "2", or "3"); unset or empty runs the complete suite, which is what
# Pkg.test and local runs use.
const CORE_TEST_FILES = (
    ("regression/current_api.jl", 3),
    ("statistics/scaling_consistency_reference.jl", 3),
    ("statistics/nist_strd_reference.jl", 3),
    ("statistics/covariance_semantics_reference.jl", 1),
    ("statistics/diagnostics_reference.jl", 3),
    ("statistics/linear_gaussian_reference.jl", 3),
    ("statistics/structured_whitening_reference.jl", 2),
    ("statistics/likelihood_reference.jl", 3),
    ("statistics/vectorized_density_reference.jl", 3),
    ("statistics/observation_likelihood_reference.jl", 2),
    ("statistics/distribution_objects.jl", 2),
    ("statistics/distribution_histograms.jl", 1),
    ("statistics/profile_contour_reference.jl", 1),
    ("numerics/inplace_model_reference.jl", 2),
    ("numerics/finite_derivatives_reference.jl", 2),
    ("numerics/nonsmooth_likelihood_reference.jl", 3),
    ("numerics/solver_status_reference.jl", 3),
    ("numerics/solver_contract.jl", 2),
    ("numerics/torture_inputs.jl", 1),
)

const _SHARD = strip(get(ENV, "SCIENTIFICFITTING_TEST_SHARD", ""))
isempty(_SHARD) || _SHARD in ("1", "2", "3") ||
    error("SCIENTIFICFITTING_TEST_SHARD must be empty, \"1\", \"2\", or \"3\"")

@testset "ScientificFitting core test suite$(isempty(_SHARD) ? "" : " (shard $_SHARD)")" begin
    for (file, shard) in CORE_TEST_FILES
        isempty(_SHARD) || string(shard) == _SHARD || continue
        # Show progress before compilation; Julia's test summary appears only at the end.
        println("Core test: ", file)
        flush(stdout)
        @time include(file)
    end
end
