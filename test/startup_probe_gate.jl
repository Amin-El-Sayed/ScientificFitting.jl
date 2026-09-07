using Test
using TOML

const ROOT = abspath(joinpath(@__DIR__, ".."))
const STARTUP_PROBE = joinpath(ROOT, "benchmarks", "startup_probe.jl")

@testset "Nested derivatives after precompilation" begin
    # A warmed suite can hide ForwardDiff tag-order collisions. This must be
    # the first fit in a new process, using a model absent from the workload.
    code = raw"""
    using ScientificFitting, Test, LinearAlgebra
    if only(ARGS) == "gaussian"
        model(x, p) = @. p[1] * x + p[2]
        dydx(x, p) = fill(p[1], length(x))
        x = [0., 1., 2., 3.]
        y = [0.8, 3.1, 4.8, 7.2]
        for xerror in ((sigma_x=fill(0.05, 4),),
                       (cov_x=0.05^2 .* [0.4^abs(i-j) for i in 1:4, j in 1:4],))
            options = (; p0=[1., 0.], sigma_y=fill(0.2, 4), xerror...)
            result = fit_model(model, x, y; options...)
            reference = fit_model(model, x, y; options..., x_derivative=dydx)
            @test result.converged && reference.converged
            @test result.params ≈ reference.params atol=1e-8
            @test result.param_covariance ≈ reference.param_covariance rtol=1e-7
            @test result.stats.cost_min ≈ reference.stats.cost_min atol=1e-9
        end
    else
        # Exact cost: (2a - 4)^2 + (b - 1)^2, including a nested derivative.
        objective(p) = (ScientificFitting.ForwardDiff.derivative(t -> p[1]*t^2, 1.) - 4)^2 +
                       (p[2] - 1)^2
        result = fit_custom(objective; p0=[0.5, 0.], nobs=3)
        @test result.converged
        @test result.params ≈ [2., 1.] atol=1e-8
        @test result.param_covariance ≈ Diagonal([0.25, 1.]) rtol=1e-8
        @test result.stats.cost_min < 1e-14
    end
    """
    for family in ("gaussian", "likelihood")
        @test success(`$(Base.julia_cmd()) --project=$ROOT --startup-file=no -e $code $family`)
    end
end

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
        project = TOML.parsefile(joinpath(ROOT, "Project.toml"))
        metadata = summary["metadata"]
        startup = summary["startup"]["core_without_makie"]

        @test metadata["scientificfitting_version"] == project["version"]
        @test metadata["unit_time"] == "seconds"
        @test haskey(metadata, "git_commit")
        @test startup["elapsed_seconds"] > 0
        @test occursin("loaded_plot_modules=", startup["stdout"])
        @test occursin("core_without_makie=true", startup["stdout"])
    finally
        rm(output_path; force=true)
    end
end
