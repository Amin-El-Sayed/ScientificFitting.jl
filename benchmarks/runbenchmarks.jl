using BenchmarkTools
using JuFitter
using LinearAlgebra
using Random

const SUITE = BenchmarkGroup()
SUITE["fit"] = BenchmarkGroup()
SUITE["likelihood"] = BenchmarkGroup()
SUITE["plot"] = BenchmarkGroup()
SUITE["profile"] = BenchmarkGroup()

function linear_problem(n::Int)
    x = collect(range(0.0, 10.0; length=n))
    model(x, p) = @. p[1] * x + p[2]
    sigma_y = fill(0.2, n)
    y = model(x, [2.0, 1.0]) .+ sigma_y .* sin.(1.7 .* x)
    return model, x, y, sigma_y
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

model100, x100, y100, sy100 = linear_problem(100)
model10k, x10k, y10k, sy10k = linear_problem(10_000)
nlmodel, nlx, nly, nlsy = nonlinear_problem(1_000)
fcmodel, fcx, fcy, fccov = full_covariance_problem(500)
pmodel, px, pcounts = poisson_problem(5_000)

SUITE["fit"]["linear_100"] = @benchmarkable fit_model($model100, $x100, $y100; p0=[1.0, 0.0], sigma_y=$sy100)
SUITE["fit"]["linear_10000"] = @benchmarkable fit_model($model10k, $x10k, $y10k; p0=[1.0, 0.0], sigma_y=$sy10k)
SUITE["fit"]["linear_10000_noop_bounds"] = @benchmarkable fit_model($model10k, $x10k, $y10k; p0=[1.0, 0.0], sigma_y=$sy10k, bounds=([-Inf, -Inf], [Inf, Inf]))
SUITE["fit"]["nonlinear_1000"] = @benchmarkable fit_model($nlmodel, $nlx, $nly; p0=[1.5, 0.4, 0.0], sigma_y=$nlsy)
SUITE["fit"]["full_covariance_500"] = @benchmarkable fit_model($fcmodel, $fcx, $fcy; p0=[1.0, 0.0], cov_y=$fccov, scale_covariance=:never)
SUITE["fit"]["full_covariance_500_bounded"] = @benchmarkable fit_model($fcmodel, $fcx, $fcy; p0=[1.0, 0.0], cov_y=$fccov, bounds=([0.0, -5.0], [5.0, 5.0]), scale_covariance=:never, maxiters=200)
SUITE["likelihood"]["poisson_5000"] = @benchmarkable fit_poisson_model($pmodel, $px, $pcounts; p0=[1.0, 0.1], bounds=([-10.0, -10.0], [10.0, 10.0]), maxiters=200)

baseline = fit_model(model100, x100, y100; p0=[1.0, 0.0], sigma_y=sy100)
SUITE["plot"]["fit_png"] = @benchmarkable plot_fit($baseline; filename=tempname() * ".png", format=:png)
SUITE["profile"]["linear_profile"] = @benchmarkable profile($baseline, 1; npoints=21, nsigma=2)

if abspath(PROGRAM_FILE) == @__FILE__
    println("Julia ", VERSION)
    println("Threads: ", Threads.nthreads())
    tune!(SUITE)
    results = run(SUITE; seconds=1)
    display(results)
end
