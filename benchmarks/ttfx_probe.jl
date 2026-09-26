using ScientificFitting, Printf
x = collect(1.0:30.0); y = 2 .* x .+ 1 .+ 0.01 .* sin.(x)
mA(x, p) = @. p[1] * x + p[2]
mB(x, p) = @. p[1] * exp(-x / p[2])
mC(x, p) = @. p[1] * x + p[2] + 0 * p[1]
yb = 5 .* exp.(-x ./ 7.0)
# :finite ist jetzt der geteilte Explorationspfad — ohne manuelles Wrappen.
t0 = time(); fit_model(mA, x, y; p0=[1.0, 0.0], derivatives=:finite); tA = time() - t0
t0 = time(); rB = fit_model(mB, x, yb; p0=[4.0, 5.0], derivatives=:finite); tB = time() - t0
t0 = time(); fit_model(mC, x, y; p0=[1.0, 0.0], derivatives=:finite); tC = time() - t0
# :auto (AD) bleibt spezialisiert:
t0 = time(); fit_model(mA, x, y; p0=[1.0, 0.0]); tAuto = time() - t0
# Poisson :finite über zwei Ratenmodelle:
cx = collect(1.0:8.0); counts = round.(Int, 40 .* exp.(-cx ./ 4))
r1(x, p) = @. p[1] * exp(-x / p[2])
r2(x, p) = @. p[1] * exp(-x / p[2]) + 0 * p[1]
t0 = time(); fit_poisson_model(r1, cx, counts; p0=[45.0, 3.5], derivatives=:finite); tP1 = time() - t0
t0 = time(); rp = fit_poisson_model(r2, cx, counts; p0=[45.0, 3.5], derivatives=:finite); tP2 = time() - t0
@printf "finite: A=%.2fs B=%.3fs C=%.3fs | auto A=%.2fs | poisson P1=%.2fs P2=%.3fs | B=%s P=%s\n" tA tB tC tAuto tP1 tP2 round.(rB.params; digits=3) round.(rp.params; digits=2)
