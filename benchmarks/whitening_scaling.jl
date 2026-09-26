using ScientificFitting, LinearAlgebra, Printf

# AR(1) whitening: matrix-free operator vs dense covariance, growing n.
function ar1_operator(n, sigma, rho)
    logdet_cov = 2n * log(sigma) + (n - 1) * log(1 - rho^2)
    scale = 1 / (sigma * sqrt(1 - rho^2))
    op = function (out, r)
        out[1] = r[1] / sigma
        @inbounds for i in 2:n
            out[i] = (r[i] - rho * r[i-1]) * scale
        end
        return nothing
    end
    return WhiteningOperator(op; logdet_covariance=logdet_cov,
                             marginal_sigma=fill(sigma, n))
end

model(x, p) = @. p[1] * x + p[2]
sigma, rho = 0.3, 0.7
for n in (10_000, 10_000, 100_000, 1_000_000)  # first 10k run is compilation
    x = collect(range(0.0, 10.0; length=n))
    y = 2 .* x .+ 1 .+ sigma .* randn(n)
    op = ar1_operator(n, sigma, rho)
    t0 = time(); r1 = fit_model(model, x, y; p0=[1.0, 0.0], whitening=op, derivatives=:finite); t_op = time() - t0
    t_dense = NaN
    if n <= 10_000
        C = [sigma^2 * rho^abs(i - j) for i in 1:n, j in 1:n]
        t0 = time(); r2 = fit_model(model, x, y; p0=[1.0, 0.0], cov_y=C); t_dense = time() - t0
        @assert isapprox(r1.params, r2.params; rtol=1e-6)
    end
    @printf "n=%-8d operator=%.3fs dense=%.3fs params=%s\n" n t_op t_dense round.(r1.params; digits=4)
end
