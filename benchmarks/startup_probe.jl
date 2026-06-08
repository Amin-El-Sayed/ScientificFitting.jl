using Dates
using TOML

const BENCHMARK_ROOT = @__DIR__
const PROJECT_ROOT = dirname(BENCHMARK_ROOT)

function _parse_args(args)
    options = Dict{String, Any}(
        "save" => nothing,
    )
    for arg in args
        if startswith(arg, "--save=")
            value = strip(last(split(arg, "="; limit=2)))
            isempty(value) && throw(ArgumentError("--save requires a non-empty path"))
            options["save"] = value
        else
            throw(ArgumentError("unknown startup probe option: $arg"))
        end
    end
    return options
end

function _git_commit()
    try
        return strip(read(`git -C $PROJECT_ROOT rev-parse --short HEAD`, String))
    catch
        return "unknown"
    end
end

function _project_version()
    project = TOML.parsefile(joinpath(PROJECT_ROOT, "Project.toml"))
    return string(get(project, "version", "unknown"))
end

function _core_probe_code()
    return """
        using JuFitter
        loaded_plot_modules = sort!(
            String[
                string(nameof(mod))
                for mod in values(Base.loaded_modules)
                if nameof(mod) in (:Makie, :CairoMakie)
            ],
        )
        println("loaded_plot_modules=", join(loaded_plot_modules, ","))
        isempty(loaded_plot_modules) || exit(2)
        println("core_without_makie=true")
    """
end

function _run_core_probe()
    cmd = `$(Base.julia_cmd()) --project=$PROJECT_ROOT --startup-file=no -e $(_core_probe_code())`
    output = ""
    elapsed = @elapsed output = read(cmd, String)
    return Dict{String, Any}(
        "elapsed_seconds" => elapsed,
        "stdout" => strip(output),
    )
end

function _summary(result)
    return Dict{String, Any}(
        "metadata" => Dict{String, Any}(
            "created_utc" => string(now(UTC)),
            "jufitter_version" => _project_version(),
            "julia_version" => string(VERSION),
            "machine" => Sys.MACHINE,
            "git_commit" => _git_commit(),
            "unit_time" => "seconds",
        ),
        "startup" => Dict{String, Any}(
            "core_without_makie" => result,
        ),
    )
end

function _write_summary(path::AbstractString, summary)
    directory = dirname(path)
    isempty(directory) || mkpath(directory)
    open(path, "w") do io
        TOML.print(io, summary; sorted=true)
    end
    return path
end

function main(args=ARGS)
    options = _parse_args(args)
    result = _run_core_probe()
    summary = _summary(result)

    println(result["stdout"])
    println("core startup probe elapsed seconds = ", result["elapsed_seconds"])

    if options["save"] !== nothing
        path = abspath(options["save"])
        _write_summary(path, summary)
        println("Saved startup probe summary to ", path)
    end
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
