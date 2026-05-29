using BenchmarkTools
using JuFitter
using LinearAlgebra
using Random

const SUITE = BenchmarkGroup()
SUITE["fit"] = BenchmarkGroup()
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

model100, x100, y100, sy100 = linear_problem(100)
model10k, x10k, y10k, sy10k = linear_problem(10_000)
nlmodel, nlx, nly, nlsy = nonlinear_problem(1_000)

SUITE["fit"]["linear_100"] = @benchmarkable fit_model($model100, $x100, $y100; p0=[1.0, 0.0], sigma_y=$sy100)
SUITE["fit"]["linear_10000"] = @benchmarkable fit_model($model10k, $x10k, $y10k; p0=[1.0, 0.0], sigma_y=$sy10k)
SUITE["fit"]["nonlinear_1000"] = @benchmarkable fit_model($nlmodel, $nlx, $nly; p0=[1.5, 0.4, 0.0], sigma_y=$nlsy)

baseline = fit_model(model100, x100, y100; p0=[1.0, 0.0], sigma_y=sy100)
SUITE["plot"]["fit_png"] = @benchmarkable plot_fit($baseline; filename=tempname() * ".png", format=:png)
SUITE["profile"]["linear_profile"] = @benchmarkable profile($baseline, 1; npoints=21, nsigma=2)

if abspath(PROGRAM_FILE) == @__FILE__
    tune!(SUITE)
    results = run(SUITE; seconds=1)
    display(results)
end
