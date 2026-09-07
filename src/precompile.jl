# Cache common numerical and reporting kernels at installation, not at first fit.
# These tiny deterministic workloads never load plotting or call foreign runtimes.
@setup_workload begin
    x = collect(0.0:5.0)
    y = [0.8, 2.9, 4.7, 6.8, 8.7, 10.9]
    model = (x, p) -> @. p[1] * x + p[2]
    counts = [4.0, 6.0, 8.0]
    rate_model = (x, p) -> fill(p[1], length(x))

    @compile_workload begin
        for derivatives in (:auto, :finite)
            callback = derivatives == :finite ? _TypedCallback{Vector{Float64}}(model) : model
            result = fit_model(callback, x, y; p0=[1.0, 0.0], sigma_y=fill(0.2, length(x)), derivatives)
            predict(result, x; uncertainty=true)
            report_text(result)
            diagnostic_dashboard_text(result)

            rate_callback = derivatives == :finite ? _TypedCallback{Vector{Float64}}(rate_model) : rate_model
            poisson = fit_poisson_model(rate_callback, x[1:3], counts;
                                        p0=[5.0], derivatives)
            report_text(poisson)
            diagnostic_dashboard_text(poisson)
        end
    end
end
