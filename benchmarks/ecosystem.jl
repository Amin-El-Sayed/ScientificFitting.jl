# Optional ecosystem throughput probe; run in the resolved docs environment.
# Controlled data are generated outside timing. This is not a solver ranking.
using ScientificFitting, Distributions, DistributionsHEP, BuildConstructors
using Random, Statistics, LinearAlgebra, Printf
import NativeMinuit

BLAS.set_num_threads(1)

# Both entry points construct exactly the same normalized peak and continuum.
peak_model(p, exposure) = ExtendedMixtureModel(
    [truncated(Normal(p[1], p[2]), -5., 5.), Uniform(-5., 5.)], exposure .* p[3:4])

struct PeakConstructor{M, S, A, B} <: BuildConstructors.AbstractConstructor
    mean::M
    sigma::S
    signal_yield::A
    background_yield::B
    exposure::Float64
end
BuildConstructors.build_model(c::PeakConstructor, pars) = peak_model([
    BuildConstructors.value(c.mean; pars), BuildConstructors.value(c.sigma; pars),
    BuildConstructors.value(c.signal_yield; pars), BuildConstructors.value(c.background_yield; pars)], c.exposure)

function measured(f; repeats=3)
    result = f() # Compile and check before measuring warmed runs.
    @assert result.converged "Do not time an unsuccessful fit"
    timings = map(1:repeats) do _
        GC.gc()
        sample = @timed f()
        (; seconds=sample.time, bytes=sample.bytes)
    end
    (; result, seconds=median(t.seconds for t in timings),
       bytes=median(t.bytes for t in timings))
end

function compare(n)
    rng = MersenneTwister(301)
    model = p -> peak_model(p, n)
    truth = model([0.15, 0.8, 0.65, 0.35])
    data = rand(rng, MixtureModel(truth), n)
    # Scaling the yields leaves the likelihood unchanged and avoids comparing
    # poorly conditioned coordinates with a solver's internal rescaling.
    p0 = [0.4, 1.1, 0.5, 0.5]
    bounds = ([-2., 0.2, 0., 0.], [2., 3., 2., 2.])
    names = ["mean", "sigma", "signal_yield", "background_yield"]
    constructor = PeakConstructor((AdvancedParameter(names[i], p0[i];
        boundaries=(bounds[1][i], bounds[2][i])) for i in 1:4)..., Float64(n))
    steps = 0.05
    settings = (; p0, bounds, parameter_names=names, maxiters=1000)
    native = NativeMinuitSolver(; steps)
    optim = OptimizationSolver(ScientificFitting.OptimizationOptimJL.LBFGS())
    # Same inference/report construction for both objective implementations.
    direct_cost(p) = 2extended_negative_log_likelihood(model(p), data)
    jobs = (
        "factory / MIGRAD" => () -> fit_distribution(model, data;
            settings..., solver=native, tol=1e-3),
        "constructor / MIGRAD" => () -> fit_distribution(constructor, data;
            solver=native, tol=1e-3, maxiters=1000),
        "factory / LBFGS" => () -> fit_distribution(model, data;
            settings..., solver=optim, tol=1e-7),
        "upstream NLL / MIGRAD" => () -> fit_custom(direct_cost;
            settings..., nobs=n, solver=native, tol=1e-3),
    )
    println("\nevents = ", n, "; repeats = 3; includes covariance and FitResult; excludes JIT/data generation")
    outcomes = map(jobs) do (name, f)
        sample = measured(f)
        r = sample.result
        @printf("%-24s %9.3f ms  %9.3f MiB  converged=%s  cost=%.9g\n",
            name, sample.seconds * 1e3, sample.bytes / 2.0^20, r.converged, r.stats.cost_min)
        println("  parameters = ", r.params, "; stderr = ", r.param_stderr)
        flush(stdout)
        sample
    end
    reference = outcomes[1].result
    for (job, outcome) in zip(jobs, outcomes)
        r = outcome.result
        pull = maximum(abs.((r.params - reference.params) ./ reference.param_stderr))
        covariance = norm(r.param_covariance - reference.param_covariance) / norm(reference.param_covariance)
        @printf("%-24s max difference/error=%.5g; delta cost=%.5g; covariance relnorm=%.5g\n",
            job.first, pull, r.stats.cost_min-reference.stats.cost_min, covariance)
        @assert r.converged
        @assert pull < 0.01
        @assert abs(r.stats.cost_min-reference.stats.cost_min) < 2e-4
        @assert covariance < 1e-3
    end
end

println("Julia ", VERSION, "; NativeMinuit ", pkgversion(NativeMinuit),
    "; CPU ", Sys.CPU_NAME, "; Julia threads ", Threads.nthreads(),
    "; BLAS threads ", BLAS.get_num_threads())
println("Distributions ", pkgversion(Distributions), "; DistributionsHEP ",
    pkgversion(DistributionsHEP), "; BuildConstructors ", pkgversion(BuildConstructors))
for n in (isempty(ARGS) ? [10_000, 100_000] : parse.(Int, ARGS))
    compare(n)
end
