using ScientificFitting, LinearAlgebra, Printf

# AR(1) whitening: matrix-free operator vs dense covariance, growing n.
model(x, p) = @. p[1] * x + p[2]
sigma, rho = 0.3, 0.7
for n in (10_000, 10_000, 100_000, 1_000_000)  # first 10k run is compilation
    x = collect(range(0.0, 10.0; length=n))
    y = 2 .* x .+ 1 .+ sigma .* randn(n)
    op = exponential_whitening(1:n, sigma, -1 / log(rho))
    t0 = time(); r1 = fit_model(model, x, y; p0=[1.0, 0.0], whitening=op, derivatives=:finite); t_op = time() - t0
    t_dense = NaN
    if n <= 10_000
        C = [sigma^2 * rho^abs(i - j) for i in 1:n, j in 1:n]
        t0 = time(); r2 = fit_model(model, x, y; p0=[1.0, 0.0], cov_y=C); t_dense = time() - t0
        @assert isapprox(r1.params, r2.params; rtol=1e-6)
    end
    @printf "n=%-8d operator=%.3fs dense=%.3fs params=%s\n" n t_op t_dense round.(r1.params; digits=4)
end
