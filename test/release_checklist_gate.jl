using Test

const ROOT = abspath(joinpath(@__DIR__, ".."))
const RELEASE_CHECKLIST = joinpath(ROOT, "RELEASE_CHECKLIST.md")

const REQUIRED_RELEASE_CHECKLIST_STRINGS = [
    "do not push, publish, register, deploy documentation",
    "git status --short --branch",
    "git diff --check",
    "julia --project=. --startup-file=no -e 'include(\"test/core_runtests.jl\")'",
    "julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'",
    "test/torture_runtests.jl",
    "test/docs_gallery_gate.jl",
    "test/docs_public_release_gate.jl",
    "test/docs_api_reference_gate.jl",
    "test/docs_link_gate.jl",
    "test/docs_visual_asset_gate.jl",
    "test/docs_visual_snapshot_gate.jl",
    "julia --project=docs --startup-file=no docs/make.jl",
    "test/docs_html_link_gate.jl",
    "test/docs_output_snapshots.jl",
    "test/plots/fitplot.jl",
    "benchmarks/startup_probe.jl",
    "test/benchmark_contract_gate.jl",
    "test/performance_budget_gate.jl",
    "julia --project=benchmarks benchmarks/runbenchmarks.jl --save",
    "julia --project=benchmarks benchmarks/runbenchmarks.jl --compare",
    "structured covariance or custom whitening",
    "Parameter covariance is a local approximation",
    "profile/contour intervals",
    "Parameter-dependent dense `cov_x` propagation",
    "finite-difference value, gradient, and Hessian",
    "JUFITTER_RUN_PYTHON_INTEROP=1",
    "Amin-El-Sayed <78275938+Amin-El-Sayed@users.noreply.github.com>",
    "manual approval",
]

@testset "Pre-release checklist coverage" begin
    @test isfile(RELEASE_CHECKLIST)
    text = read(RELEASE_CHECKLIST, String)
    @test !isempty(strip(text))

    for required in REQUIRED_RELEASE_CHECKLIST_STRINGS
        @test occursin(required, text)
    end
end
