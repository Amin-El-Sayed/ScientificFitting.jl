using Test

const ROOT = abspath(joinpath(@__DIR__, ".."))
const PYTHON_INTEROP_SCRIPT = joinpath(ROOT, "examples", "python", "fit_from_python.py")

function python3_cmd()
    return get(ENV, "PYTHON", "python3")
end

function command_available(cmd::AbstractString)
    try
        run(pipeline(Cmd([cmd, "--version"]); stdout=devnull, stderr=devnull))
        return true
    catch
        return false
    end
end

function juliacall_available(cmd::AbstractString)
    try
        run(pipeline(Cmd([cmd, "-c", "import juliacall"]); stdout=devnull, stderr=devnull))
        return true
    catch
        return false
    end
end

@testset "Python interoperability" begin
    run_interop = get(ENV, "JUFITTER_RUN_PYTHON_INTEROP", "0") == "1"
    py = python3_cmd()

    @test isfile(PYTHON_INTEROP_SCRIPT)

    if !run_interop
        @test true
        @info(
            "Python interoperability gate is opt-in. Set " *
            "JUFITTER_RUN_PYTHON_INTEROP=1 after installing juliacall."
        )
    else
        @test command_available(py)
        @test juliacall_available(py)
        output = read(Cmd([py, PYTHON_INTEROP_SCRIPT]), String)
        @test occursin("JuFitter from Python", output)
        @test occursin("Fit report", output)
        @test occursin("Fit diagnostic dashboard", output)
        @test occursin("slope", output)
        @test occursin("offset", output)
        @test occursin("plot modules loaded = []", output)
        @test !occursin("plot modules loaded = ['Makie']", output)
        @test !occursin("plot modules loaded = ['CairoMakie']", output)
    end
end
