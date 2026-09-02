using Test

const ROOT = abspath(joinpath(@__DIR__, ".."))
const BENCHMARK_RUNNER = joinpath(ROOT, "benchmarks", "runbenchmarks.jl")
const STARTUP_PROBE = joinpath(ROOT, "benchmarks", "startup_probe.jl")
const PERFORMANCE_DOC = joinpath(ROOT, "docs", "src", "performance.md")
const CI_WORKFLOW = joinpath(ROOT, ".github", "workflows", "ci.yml")

const CANONICAL_BENCHMARK_COMMAND =
    "julia --project=benchmarks benchmarks/runbenchmarks.jl"

const REQUIRED_BENCHMARK_CASES = [
    "fit/linear_100",
    "fit/linear_10000",
    "fit/linear_10000_inplace",
    "fit/linear_10000_noop_bounds",
    "fit/nonlinear_1000",
    "fit/full_covariance_500",
    "fit/full_covariance_500_bounded",
    "fit/sparse_covariance_5000",
    "fit/structured_whitening_100000",
    "likelihood/poisson_5000",
    "profile/linear_profile",
    "diagnostics/saturation_profile_matrix",
]

const REQUIRED_SUMMARY_FIELDS = [
    "julia_version",
    "os",
    "cpu_name",
    "cpu_threads",
    "threads",
    "blas_threads",
    "machine",
    "git_commit",
    "created_utc",
    "scientificfitting_version",
    "seconds_per_benchmark",
    "benchmark_case_count",
    "unit_time",
    "unit_memory",
    "median_seconds",
    "memory_bytes",
    "allocations",
]

function file_text(path)
    @test isfile(path)
    return read(path, String)
end

@testset "Benchmark contract gate" begin
    runner = file_text(BENCHMARK_RUNNER)
    startup_probe = file_text(STARTUP_PROBE)
    performance = file_text(PERFORMANCE_DOC)
    workflow = file_text(CI_WORKFLOW)

    @test occursin(CANONICAL_BENCHMARK_COMMAND, performance)
    @test occursin("benchmarks/startup_probe.jl", performance)

    for case in REQUIRED_BENCHMARK_CASES
        group, name = split(case, "/"; limit=2)
        @test occursin("SUITE[\"$group\"][\"$name\"]", runner)
    end

    for field in REQUIRED_SUMMARY_FIELDS
        @test occursin(field, runner)
    end

    @test occursin("--save=", runner)
    @test occursin("--compare=", runner)
    @test occursin("--tolerance=", runner)
    @test occursin("--list", runner)
    @test occursin("--plot", runner)
    @test occursin("--allow-metadata-mismatch", runner)
    @test occursin("STRICT_METADATA_FIELDS", runner)
    @test occursin("Use --allow-metadata-mismatch only for exploratory comparisons", runner)
    @test occursin("missing from the current run", runner)
    @test occursin("missing from the baseline", runner)
    @test occursin("finite positive number", runner)
    @test occursin("finite non-negative number", runner)
    @test occursin("requires a non-empty path", runner)

    @test occursin("loaded_plot_modules", startup_probe)
    @test occursin(":Makie", startup_probe)
    @test occursin(":CairoMakie", startup_probe)
    @test occursin("core_without_makie", startup_probe)

    @test occursin("Benchmark contract gate", workflow)
    @test occursin("test/benchmark_contract_gate.jl", workflow)
    @test occursin("Startup probe gate", workflow)
    @test occursin("test/startup_probe_gate.jl", workflow)
end
