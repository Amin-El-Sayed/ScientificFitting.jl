using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using JuFitter
include(joinpath(@__DIR__, "_example_utils.jl"))

# Poisson count fit: model returns expected counts, not y-values.
x = collect(1.0:10.0)
counts = [2, 4, 5, 6, 9, 11, 12, 14, 15, 19]
poisson_model(x, p) = @. exp(p[1] + p[2] * x)
poisson_result = fit_poisson_model(
    poisson_model,
    x,
    counts;
    p0=[0.0, 0.1],
    bounds=([-10.0, -10.0], [10.0, 10.0]),
    parameter_names=["log_scale", "slope"],
)
println("Poisson params = ", poisson_result.params)
println("Poisson deviance/ndf = ", poisson_result.stats.chi2_ndf)

# Histogram fit: expected_counts(edges, p) returns one positive expected count per bin.
edges = collect(0.0:1.0:5.0)
hist_counts = [4, 9, 13, 20, 31]
expected_counts(edges, p) = [p[1] * (edges[i + 1] - edges[i]) * exp(p[2] * (edges[i] + edges[i + 1]) / 2) for i in 1:(length(edges) - 1)]
hist_result = fit_histogram_model(
    expected_counts,
    edges,
    hist_counts;
    p0=[3.0, 0.3],
    bounds=([1e-6, -5.0], [100.0, 5.0]),
)
println("Histogram params = ", hist_result.params)
println("Histogram deviance/ndf = ", hist_result.stats.chi2_ndf)

# Unbinned likelihood fit for a normalized density.
data = [-1.1, -0.2, 0.1, 0.3, 0.9, 1.2]
normal_pdf(x, p) = exp(-0.5 * ((x - p[1]) / p[2])^2) / (p[2] * sqrt(2 * pi))
unbinned_result = fit_unbinned_model(
    normal_pdf,
    data;
    p0=[0.0, 1.0],
    bounds=([-5.0, 0.1], [5.0, 5.0]),
    parameter_names=["mu", "sigma"],
)
println("Unbinned normal params = ", unbinned_result.params)

# Extended unbinned fit: rate(x, p) is an intensity, not a normalized PDF.
rate(x, p) = exp(p[1])
extended_result = fit_extended_unbinned_model(
    rate,
    [0.1, 0.2, 0.8, 0.9],
    (0.0, 1.0);
    p0=[0.0],
    parameter_names=["log_rate"],
)
println("Extended unbinned expected count = ", exp(extended_result.params[1]))

# Indexed fit: observations are addressed by labels rather than a numeric x-axis.
indices = [:a, :b, :a, :c]
y_indexed = [1.0, 2.0, 1.1, 3.0]
indexed_model(indices, p) = [p[idx == :a ? 1 : idx == :b ? 2 : 3] for idx in indices]
indexed_result = fit_indexed_model(
    indexed_model,
    indices,
    y_indexed;
    p0=[0.0, 0.0, 0.0],
    sigma_y=fill(0.1, length(y_indexed)),
    parameter_names=["level_a", "level_b", "level_c"],
)
println("Indexed params = ", indexed_result.params)

# Custom fit: directly minimize a user-defined scalar objective.
custom_result = fit_custom(
    p -> sum(abs2, p .- [1.0, 2.0]);
    p0=[0.0, 0.0],
    nobs=4,
    parameter_names=["a", "b"],
)
println("Custom params = ", custom_result.params)

# MultiFit with parameter mapping: both datasets share global p[1], but use different offsets.
x1 = collect(0.0:1.0:5.0)
x2 = collect(0.0:1.0:5.0)
local_linear(x, p) = @. p[1] * x + p[2]
y1 = local_linear(x1, [2.0, 1.0])
y2 = local_linear(x2, [2.0, -1.0])
multi_result = fit_multi_model(
    [local_linear, local_linear],
    [x1, x2],
    [y1, y2];
    p0=[0.0, 0.0, 0.0],
    parameter_map=[[1, 2], [1, 3]],
    sigma_y=[fill(0.1, length(x1)), fill(0.1, length(x2))],
    parameter_names=["shared_slope", "offset_1", "offset_2"],
)
println("Mapped MultiFit params = ", multi_result.params)

# Profiles and reports also work for likelihood fits.
prof = profile(poisson_result, 2; npoints=9)
plot_profile(prof; filename=example_output("poisson_profile"), format=:png)
println("Saved plot to ", example_output("poisson_profile.png"))
println(report_text(poisson_result))
