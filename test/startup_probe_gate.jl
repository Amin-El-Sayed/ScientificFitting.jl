using Test
using TOML

const ROOT = abspath(joinpath(@__DIR__, ".."))
const STARTUP_PROBE = joinpath(ROOT, "benchmarks", "startup_probe.jl")

@testset "Startup probe gate" begin
    @test isfile(STARTUP_PROBE)

    output_path = tempname() * ".toml"
    try
        cmd = `$(Base.julia_cmd()) --project=$ROOT --startup-file=no $STARTUP_PROBE --save=$output_path`
        output = read(cmd, String)

        @test occursin("loaded_plot_modules=", output)
        @test occursin("core_without_makie=true", output)
        @test !occursin("loaded_plot_modules=Makie", output)
        @test !occursin("loaded_plot_modules=CairoMakie", output)
        @test isfile(output_path)

        summary = TOML.parsefile(output_path)
        metadata = summary["metadata"]
        startup = summary["startup"]["core_without_makie"]

        @test metadata["jufitter_version"] == "0.1.0"
        @test metadata["unit_time"] == "seconds"
        @test haskey(metadata, "git_commit")
        @test startup["elapsed_seconds"] > 0
        @test occursin("loaded_plot_modules=", startup["stdout"])
        @test occursin("core_without_makie=true", startup["stdout"])
    finally
        rm(output_path; force=true)
    end
end
