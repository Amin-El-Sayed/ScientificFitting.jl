using CairoMakie
using Printf
using ScientificFitting
using SpecialFunctions
include(joinpath(@__DIR__, "..", "_example_utils.jl"))

# Poisson decay: synthetic counts in a 10 s acquisition window taken once per
# minute (t in minutes, one Poisson-like draw per window from ~45 exp(-0.15 t) + 2),
# so the fitted decay constant is in 1/min and the half-life in minutes.
t_obs = collect(0.0:1.0:18.0)
counts = [48, 37, 35, 27, 27, 17, 22, 13, 16, 8, 13, 5, 11, 4, 7, 2, 6, 1, 5]
poisson_model(t, p) = @. p[1] * exp(-p[2] * t) + p[3]
poisson_result = fit_poisson_model(
    poisson_model,
    t_obs,
    counts;
    p0=[40.0, 0.15, 3.0],
    bounds=([1e-6, 1e-6, 1e-6], [200.0, 2.0, 50.0]),
    parameter_names=["initial signal", "decay constant", "background"],
    initial_guesses=[[70.0, 0.30, 2.0], [25.0, 0.08, 5.0]], # always tried in addition to p0
)
print_result_summary("Poisson count fit", poisson_result)
half_life = log(2) / poisson_result.params[2]
sigma_half_life = log(2) * poisson_result.param_stderr[2] / poisson_result.params[2]^2
@printf("Half-life = %.2f +/- %.2f min\n", half_life, sigma_half_life)

# Histogram fit: the bins have unequal widths, so the model integrates the
# peak and the background density over each bin instead of evaluating at
# bin centers. Synthetic counts around a Gaussian peak on a flat background.
edges = [0.0, 0.4, 0.9, 1.5, 2.2, 3.0, 4.0, 5.2, 6.6, 8.2, 10.0]
hist_counts = [0, 3, 9, 24, 47, 69, 51, 24, 8, 4]
function expected_counts(edges, p)
    peak_yield, centroid, width, background_density = p
    return [
        peak_yield * 0.5 * (
            erf((edges[i + 1] - centroid) / (sqrt(2) * width)) -
            erf((edges[i] - centroid) / (sqrt(2) * width))
        ) + background_density * (edges[i + 1] - edges[i])
        for i in 1:(length(edges) - 1)
    ]
end
hist_result = fit_histogram_model(
    expected_counts,
    edges,
    hist_counts;
    p0=[210.0, 3.8, 1.0, 1.0],
    bounds=([1e-6, 0.0, 0.05, 1e-6], [1000.0, 10.0, 5.0, 100.0]),
    parameter_names=["peak yield", "centroid", "width", "background density"],
    initial_guesses=[[300.0, 4.2, 1.5, 0.5], [150.0, 3.2, 0.7, 2.0]],
)
print_result_summary("Histogram Poisson fit", hist_result)

# Unbinned likelihood fit for a normalized density. Six points only
# demonstrate the API; the reported Hessian uncertainties are asymptotic
# and not reliable at n = 6.
data = [-1.1, -0.2, 0.1, 0.3, 0.9, 1.2]
normal_pdf(x, p) = exp(-0.5 * ((x - p[1]) / p[2])^2) / (p[2] * sqrt(2 * pi))
unbinned_result = fit_unbinned_model(
    normal_pdf,
    data;
    p0=[0.0, 1.0],
    bounds=([-5.0, 0.1], [5.0, 5.0]),
    parameter_names=["mu", "sigma"],
)
print_result_summary("Unbinned normal fit", unbinned_result)

# Extended unbinned fit: rate(x, p) is an intensity, not a normalized PDF.
rate(x, p) = exp(p[1])
extended_result = fit_extended_unbinned_model(
    rate,
    [0.1, 0.2, 0.8, 0.9],
    (0.0, 1.0);
    p0=[0.0],
    parameter_names=["log_rate"],
)
print_result_summary("Extended unbinned fit", extended_result)

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
print_result_summary("Indexed fit", indexed_result)

# Custom scalar objective on the -2 log L scale: a chi-square comparing the
# two parameters against two unit-sigma reference values, so nobs=2. Two
# observations and two free parameters give ndf = 0; the report flags that
# goodness-of-fit statistics are not meaningful for this calibration-style fit.
custom_result = fit_custom(
    p -> sum(abs2, p .- [1.0, 2.0]);
    p0=[0.0, 0.0],
    nobs=2,
    parameter_names=["a", "b"],
)
print_result_summary("Custom objective fit", custom_result)

# MultiFit with parameter mapping: both datasets share a slope but have
# different offsets. The data are noise-free model values, so the fit
# reproduces the truth exactly (chi2 = 0); this block only demonstrates
# parameter_map.
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
print_result_summary("Mapped multi-dataset fit", multi_result)

prof = profile(poisson_result, 2; npoints=41, nsigma=3)
plot_profile(
    prof;
    filename=example_output("06_decay_constant_profile.pdf"),
    xlabel="decay constant",
    local_sigma=poisson_result.param_stderr[2],
    delta_max=8,
)
println("Saved profile to ", example_output("06_decay_constant_profile.pdf"))
