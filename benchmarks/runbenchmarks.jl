using Pkg

const BENCHMARK_ROOT = @__DIR__
const PROJECT_ROOT = dirname(BENCHMARK_ROOT)

if Base.active_project() == joinpath(BENCHMARK_ROOT, "Project.toml")
    Pkg.develop(PackageSpec(path=PROJECT_ROOT))
    Pkg.instantiate()
end

using BenchmarkTools
using Dates
using JuFitter
using LinearAlgebra
using TOML

const SUITE = BenchmarkGroup()
SUITE["fit"] = BenchmarkGroup()
SUITE["likelihood"] = BenchmarkGroup()
SUITE["profile"] = BenchmarkGroup()

function linear_problem(n::Int)
    x = collect(range(0.0, 10.0; length=n))
    model(x, p) = @. p[1] * x + p[2]
    jacobian(x, p) = hcat(x, ones(length(x)))
    sigma_y = fill(0.2, n)
    y = model(x, [2.0, 1.0]) .+ sigma_y .* sin.(1.7 .* x)
    return model, jacobian, x, y, sigma_y
end

function nonlinear_problem(n::Int)
    x = collect(range(0.0, 4.0; length=n))
    model(x, p) = @. p[1] * exp(-p[2] * x) + p[3]
    sigma_y = fill(0.05, n)
    y = model(x, [2.0, 0.7, 0.2]) .+ sigma_y .* cos.(2.3 .* x)
    return model, x, y, sigma_y
end

function full_covariance_problem(n::Int)
    x = collect(range(0.0, 10.0; length=n))
    model(x, p) = @. p[1] * x + p[2]
    sigma = @. 0.1 + 0.01 * abs(sin(x))
    cov = [sigma[i] * sigma[j] * 0.25^abs(i - j) for i in eachindex(x), j in eachindex(x)]
    y = model(x, [2.0, 1.0]) .+ 0.02 .* sin.(x)
    return model, x, y, cov
end

function poisson_problem(n::Int)
    x = collect(range(0.0, 5.0; length=n))
    model(x, p) = @. exp(p[1] + p[2] * x)
    counts = round.(model(x, [1.2, 0.15]))
    return model, x, counts
end

model100, jac100, x100, y100, sy100 = linear_problem(100)
model10k, jac10k, x10k, y10k, sy10k = linear_problem(10_000)
nlmodel, nlx, nly, nlsy = nonlinear_problem(1_000)
fcmodel, fcx, fcy, fccov = full_covariance_problem(500)
pmodel, px, pcounts = poisson_problem(5_000)

SUITE["fit"]["linear_100"] = @benchmarkable fit_model($model100, $x100, $y100; p0=[1.0, 0.0], sigma_y=$sy100, jacobian=$jac100)
SUITE["fit"]["linear_10000"] = @benchmarkable fit_model($model10k, $x10k, $y10k; p0=[1.0, 0.0], sigma_y=$sy10k, jacobian=$jac10k)
SUITE["fit"]["linear_10000_noop_bounds"] = @benchmarkable fit_model($model10k, $x10k, $y10k; p0=[1.0, 0.0], sigma_y=$sy10k, jacobian=$jac10k, bounds=([-Inf, -Inf], [Inf, Inf]))
SUITE["fit"]["nonlinear_1000"] = @benchmarkable fit_model($nlmodel, $nlx, $nly; p0=[1.5, 0.4, 0.0], sigma_y=$nlsy)
SUITE["fit"]["full_covariance_500"] = @benchmarkable fit_model($fcmodel, $fcx, $fcy; p0=[1.0, 0.0], cov_y=$fccov, scale_covariance=:never)
SUITE["fit"]["full_covariance_500_bounded"] = @benchmarkable fit_model($fcmodel, $fcx, $fcy; p0=[1.0, 0.0], cov_y=$fccov, bounds=([0.0, -5.0], [5.0, 5.0]), scale_covariance=:never, maxiters=200)
SUITE["likelihood"]["poisson_5000"] = @benchmarkable fit_poisson_model($pmodel, $px, $pcounts; p0=[1.0, 0.1], bounds=([-10.0, -10.0], [10.0, 10.0]), maxiters=200)

baseline = fit_model(model100, x100, y100; p0=[1.0, 0.0], sigma_y=sy100, jacobian=jac100)
SUITE["profile"]["linear_profile"] = @benchmarkable profile($baseline, 1; npoints=21, nsigma=2)

function _try_enable_plot_benchmarks!(suite)
    try
        @eval using CairoMakie
    catch err
        @info "Skipping plot benchmarks because CairoMakie is not available in this environment" exception=(err, catch_backtrace())
        return false
    end

    suite["plot"] = BenchmarkGroup()
    plot_file = joinpath(tempdir(), "jufitter_benchmark_fit.png")
    suite["plot"]["fit_png"] = @benchmarkable plot_fit($baseline; filename=$plot_file, format=:png)
    return true
end

function _parse_args(args)
    options = Dict{String, Any}(
        "seconds" => 1.0,
        "save" => nothing,
        "compare" => nothing,
        "tolerance" => 0.25,
        "plot" => false,
        "list" => false,
    )
    for arg in args
        if arg == "--plot"
            options["plot"] = true
        elseif arg == "--list"
            options["list"] = true
        elseif startswith(arg, "--seconds=")
            options["seconds"] = _parse_positive_float_option(arg, "--seconds")
        elseif startswith(arg, "--save=")
            options["save"] = _parse_nonempty_path_option(arg, "--save")
        elseif startswith(arg, "--compare=")
            options["compare"] = _parse_nonempty_path_option(arg, "--compare")
        elseif startswith(arg, "--tolerance=")
            options["tolerance"] = _parse_nonnegative_float_option(arg, "--tolerance")
        else
            throw(ArgumentError("unknown benchmark option: $arg"))
        end
    end
    return options
end

function _option_value(arg::AbstractString, name::AbstractString)
    value = last(split(arg, "="; limit=2))
    isempty(value) && throw(ArgumentError("$name requires a non-empty value"))
    return value
end

function _parse_positive_float_option(arg::AbstractString, name::AbstractString)
    value = parse(Float64, _option_value(arg, name))
    isfinite(value) && value > 0 ||
        throw(ArgumentError("$name must be a finite positive number"))
    return value
end

function _parse_nonnegative_float_option(arg::AbstractString, name::AbstractString)
    value = parse(Float64, _option_value(arg, name))
    isfinite(value) && value >= 0 ||
        throw(ArgumentError("$name must be a finite non-negative number"))
    return value
end

function _parse_nonempty_path_option(arg::AbstractString, name::AbstractString)
    value = strip(_option_value(arg, name))
    isempty(value) && throw(ArgumentError("$name requires a non-empty path"))
    return value
end

function _flatten_benchmarks(group::BenchmarkGroup; prefix::String="")
    pairs = Pair{String, Any}[]
    for key in sort!(collect(keys(group)); by=string)
        name = isempty(prefix) ? string(key) : string(prefix, "/", key)
        child = group[key]
        if child isa BenchmarkGroup
            append!(pairs, _flatten_benchmarks(child; prefix=name))
        else
            push!(pairs, name => child)
        end
    end
    return pairs
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

function _cpu_name()
    try
        info = Sys.cpu_info()
        return isempty(info) ? "unknown" : string(first(info).model)
    catch
        return "unknown"
    end
end

function _summarize_trial(trial)
    estimate = median(trial)
    return Dict{String, Any}(
        "median_seconds" => estimate.time / 1.0e9,
        "memory_bytes" => estimate.memory,
        "allocations" => estimate.allocs,
    )
end

function _summarize_results(results::BenchmarkGroup; seconds::Float64)
    summary = Dict{String, Any}()
    benchmark_names = first.(_flatten_benchmarks(results))
    summary["metadata"] = Dict{String, Any}(
        "created_utc" => string(now(UTC)),
        "jufitter_version" => _project_version(),
        "julia_version" => string(VERSION),
        "os" => string(Sys.KERNEL),
        "cpu_name" => _cpu_name(),
        "cpu_threads" => Sys.CPU_THREADS,
        "threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "machine" => Sys.MACHINE,
        "git_commit" => _git_commit(),
        "seconds_per_benchmark" => seconds,
        "benchmark_case_count" => length(benchmark_names),
        "unit_time" => "seconds",
        "unit_memory" => "bytes",
    )

    benchmarks = Dict{String, Any}()
    for (name, trial) in _flatten_benchmarks(results)
        benchmarks[name] = _summarize_trial(trial)
    end
    summary["benchmarks"] = benchmarks
    return summary
end

function _write_summary(path::AbstractString, summary)
    directory = dirname(path)
    isempty(directory) || mkpath(directory)
    open(path, "w") do io
        TOML.print(io, summary; sorted=true)
    end
    return path
end

function _compare_summary(summary, baseline_path::AbstractString; tolerance::Float64)
    baseline = TOML.parsefile(baseline_path)
    current = summary["benchmarks"]
    reference = baseline["benchmarks"]

    failures = String[]
    missing_current = sort!(collect(setdiff(keys(reference), keys(current))))
    missing_reference = sort!(collect(setdiff(keys(current), keys(reference))))
    for name in missing_current
        push!(failures, "$name is present in the baseline but missing from the current run")
    end
    for name in missing_reference
        push!(failures, "$name is new in the current run but missing from the baseline")
    end

    for name in sort!(collect(intersect(keys(current), keys(reference))))
        old = reference[name]["median_seconds"]
        new = current[name]["median_seconds"]
        if old > 0 && new > old * (1 + tolerance)
            push!(failures, "$name median time regressed from $old to $new seconds")
        end
    end

    if !isempty(failures)
        println("Benchmark regressions relative to ", baseline_path, ":")
        foreach(msg -> println("  - ", msg), failures)
        return false
    end

    println("No benchmark regression above ", round(100tolerance; digits=1), "% relative to ", baseline_path)
    return true
end

function main(args=ARGS)
    options = _parse_args(args)

    if options["plot"]
        _try_enable_plot_benchmarks!(SUITE)
    end

    if options["list"]
        foreach(pair -> println(first(pair)), _flatten_benchmarks(SUITE))
        return nothing
    end

    seconds = options["seconds"]
    println("Julia ", VERSION)
    println("Threads: ", Threads.nthreads())
    println("Benchmarks: ", length(_flatten_benchmarks(SUITE)))
    tune!(SUITE)
    results = run(SUITE; seconds=seconds)
    display(results)

    summary = _summarize_results(results; seconds=seconds)
    if options["save"] !== nothing
        path = abspath(options["save"])
        _write_summary(path, summary)
        println("Saved benchmark summary to ", path)
    end
    if options["compare"] !== nothing
        ok = _compare_summary(summary, options["compare"]; tolerance=options["tolerance"])
        ok || exit(1)
    end
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
