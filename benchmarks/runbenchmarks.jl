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
using SparseArrays
using TOML

const SUITE = BenchmarkGroup()
SUITE["fit"] = BenchmarkGroup()
SUITE["likelihood"] = BenchmarkGroup()
SUITE["profile"] = BenchmarkGroup()
SUITE["diagnostics"] = BenchmarkGroup()

function linear_problem(n::Int)
    x = collect(range(0.0, 10.0; length=n))
    model(x, p) = @. p[1] * x + p[2]
    jacobian(x, p) = hcat(x, ones(length(x)))
    sigma_y = fill(0.2, n)
    y = model(x, [2.0, 1.0]) .+ sigma_y .* sin.(1.7 .* x)
    return model, jacobian, x, y, sigma_y
end

function linear_model!(out, x, p)
    @. out = p[1] * x + p[2]
    return nothing
end

function linear_jacobian!(J, x, p)
    J[:, 1] .= x
    J[:, 2] .= 1
    return nothing
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

function sparse_covariance_problem(n::Int)
    x = collect(range(0.0, 20.0; length=n))
    model(x, p) = @. p[1] * x + p[2]
    jacobian(x, p) = hcat(x, ones(length(x)))
    variance = @. 0.08^2 * (1 + 0.05 * sin(x)^2)
    covariance = spdiagm(
        -1 => fill(0.0012, n - 1),
        0 => variance,
        1 => fill(0.0012, n - 1),
    )
    y = model(x, [1.5, -0.2]) .+ 0.03 .* sin.(0.8 .* x)
    return model, jacobian, x, y, covariance
end

function structured_covariance_problem(n::Int)
    x = collect(range(0.0, 20.0; length=n))
    model(x, p) = @. p[1] * x + p[2]
    jacobian(x, p) = hcat(x, ones(length(x)))
    sigma = 0.10
    rho = 0.70
    innovation_sigma = sigma * sqrt(1 - rho^2)
    function whiten!(out, residual)
        out[1] = residual[1] / sigma
        @inbounds for i in 2:length(residual)
            out[i] = (residual[i] - rho * residual[i - 1]) / innovation_sigma
        end
        return nothing
    end
    logdet_covariance = 2n * log(sigma) + (n - 1) * log1p(-rho^2)
    whitening = WhiteningOperator(whiten!; logdet_covariance)
    y = model(x, [1.5, -0.2]) .+ 0.03 .* sin.(0.8 .* x)
    return model, jacobian, x, y, whitening
end

function poisson_problem(n::Int)
    x = collect(range(0.0, 5.0; length=n))
    model(x, p) = @. exp(p[1] + p[2] * x)
    counts = round.(model(x, [1.2, 0.15]))
    return model, x, counts
end

function saturation_problem()
    x = collect(range(0.15, 2.2; length=18))
    model(t, p) = @. p[1] * (1 - exp(-t / p[2])) + p[3]
    sigma_x = @. 0.010 + 0.004 * x
    sigma_y = @. 0.045 + 0.008 * x
    pattern = [
        0.50, -0.90, 0.30, 1.10, -0.70, 0.80, -1.00, 0.40, 0.90,
        -0.60, 0.70, -0.80, 1.00, -0.40, 0.55, -0.75, 0.65, -0.35,
    ]
    y = model(x, [4.8, 3.4, 0.12]) .+ sigma_y .* pattern
    return model, x, y, sigma_x, sigma_y
end

model100, jac100, x100, y100, sy100 = linear_problem(100)
model10k, jac10k, x10k, y10k, sy10k = linear_problem(10_000)
nlmodel, nlx, nly, nlsy = nonlinear_problem(1_000)
fcmodel, fcx, fcy, fccov = full_covariance_problem(500)
scmodel, scjac, scx, scy, sccov = sparse_covariance_problem(5_000)
wcmodel, wcjac, wcx, wcy, wcwhitening = structured_covariance_problem(100_000)
pmodel, px, pcounts = poisson_problem(5_000)
saturation_model_bench, saturation_x, saturation_y, saturation_sigma_x, saturation_sigma_y =
    saturation_problem()

SUITE["fit"]["linear_100"] = @benchmarkable fit_model($model100, $x100, $y100; p0=[1.0, 0.0], sigma_y=$sy100, jacobian=$jac100)
SUITE["fit"]["linear_10000"] = @benchmarkable fit_model($model10k, $x10k, $y10k; p0=[1.0, 0.0], sigma_y=$sy10k, jacobian=$jac10k)
SUITE["fit"]["linear_10000_inplace"] = @benchmarkable fit_model(
    $linear_model!,
    $x10k,
    $y10k;
    p0=[1.0, 0.0],
    sigma_y=$sy10k,
    jacobian=$linear_jacobian!,
    inplace=true,
)
SUITE["fit"]["linear_10000_noop_bounds"] = @benchmarkable fit_model($model10k, $x10k, $y10k; p0=[1.0, 0.0], sigma_y=$sy10k, jacobian=$jac10k, bounds=([-Inf, -Inf], [Inf, Inf]))
SUITE["fit"]["nonlinear_1000"] = @benchmarkable fit_model($nlmodel, $nlx, $nly; p0=[1.5, 0.4, 0.0], sigma_y=$nlsy)
SUITE["fit"]["full_covariance_500"] = @benchmarkable fit_model($fcmodel, $fcx, $fcy; p0=[1.0, 0.0], cov_y=$fccov, scale_covariance=:never)
SUITE["fit"]["full_covariance_500_bounded"] = @benchmarkable fit_model($fcmodel, $fcx, $fcy; p0=[1.0, 0.0], cov_y=$fccov, bounds=([0.0, -5.0], [5.0, 5.0]), scale_covariance=:never, maxiters=200)
SUITE["fit"]["sparse_covariance_5000"] = @benchmarkable fit_model($scmodel, $scx, $scy; p0=[1.0, 0.0], cov_y=$sccov, jacobian=$scjac, scale_covariance=:never)
SUITE["fit"]["structured_whitening_100000"] = @benchmarkable fit_model($wcmodel, $wcx, $wcy; p0=[1.0, 0.0], whitening=$wcwhitening, jacobian=$wcjac, scale_covariance=:never)
SUITE["likelihood"]["poisson_5000"] = @benchmarkable fit_poisson_model($pmodel, $px, $pcounts; p0=[1.0, 0.1], bounds=([-10.0, -10.0], [10.0, 10.0]), maxiters=200)

baseline = fit_model(model100, x100, y100; p0=[1.0, 0.0], sigma_y=sy100, jacobian=jac100)
SUITE["profile"]["linear_profile"] = @benchmarkable profile($baseline, 1; npoints=21, nsigma=2)

saturation_baseline = fit_model(
    saturation_model_bench,
    saturation_x,
    saturation_y;
    p0=[3.0, 2.0, 0.0],
    sigma_y=saturation_sigma_y,
    sigma_x=saturation_sigma_x,
    bounds=([0.1, 0.1, -0.5], [20.0, 20.0, 1.0]),
    parameter_priors=(index=3, mean=0.10, sigma=0.08),
    initial_guesses=[[3.0, 2.0, 0.0], [8.0, 7.0, 0.1], [2.0, 1.0, 0.2]],
)
SUITE["diagnostics"]["saturation_profile_matrix"] = @benchmarkable profile_matrix(
    $saturation_baseline;
    parameters=[1, 2, 3],
    npoints_profile=21,
    npoints_contour=11,
    nsigma=3,
    adaptive=false,
)

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
        "allow_metadata_mismatch" => false,
    )
    for arg in args
        if arg == "--plot"
            options["plot"] = true
        elseif arg == "--list"
            options["list"] = true
        elseif arg == "--allow-metadata-mismatch"
            options["allow_metadata_mismatch"] = true
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

const STRICT_METADATA_FIELDS = [
    "julia_version",
    "os",
    "cpu_name",
    "machine",
    "threads",
    "blas_threads",
    "unit_time",
    "unit_memory",
]

function _metadata_failures(current_metadata, baseline_metadata)
    failures = String[]
    for field in STRICT_METADATA_FIELDS
        current = get(current_metadata, field, nothing)
        reference = get(baseline_metadata, field, nothing)
        current == reference && continue
        push!(failures, "metadata $field differs: baseline=$reference current=$current")
    end
    return failures
end

function _compare_summary(
    summary,
    baseline_path::AbstractString;
    tolerance::Float64,
    allow_metadata_mismatch::Bool=false,
)
    baseline = TOML.parsefile(baseline_path)
    current = summary["benchmarks"]
    reference = baseline["benchmarks"]

    failures = String[]
    if !allow_metadata_mismatch
        append!(
            failures,
            _metadata_failures(
                summary["metadata"],
                get(baseline, "metadata", Dict{String, Any}()),
            ),
        )
    end

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
        if any(startswith(failure, "metadata ") for failure in failures)
            println("Use --allow-metadata-mismatch only for exploratory comparisons, not release evidence.")
        end
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
        ok = _compare_summary(
            summary,
            options["compare"];
            tolerance=options["tolerance"],
            allow_metadata_mismatch=options["allow_metadata_mismatch"],
        )
        ok || exit(1)
    end
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
