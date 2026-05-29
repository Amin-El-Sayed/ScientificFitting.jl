using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using JuFitter
include(joinpath(@__DIR__, "_example_utils.jl"))
using CairoMakie
using LaTeXStrings
using Statistics

t = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0]
alpha_1 = [23.00, 22.75, 22.45, 22.40, 22.20, 22.00, 21.75, 21.55, 21.25, 21.10, 20.90, 20.50, 19.95, 19.95, 19.05, 18.80, 18.25, 17.85, 17.30]
alpha_2 = [22.95, 22.6, 22.40, 22.30, 22.15, 21.9, 21.65, 21.35, 21.15, 21.0, 20.85, 20.35, 19.75, 19.25, 18.85, 18.6, 18.1, 17.6, 17.15]
alpha = (alpha_1 .+ alpha_2) ./ 2

sigma_t = fill(1 / 6, length(t))
sigma_alpha_1 = fill(0.1, length(alpha_1))
sigma_alpha_2 = fill(0.1, length(alpha_2))
sigma_alpha = sqrt.(sigma_alpha_1.^2 .+ sigma_alpha_2.^2)

# Deliberately simple linearization:
# ln((alpha_t - alpha_inf) / (alpha_0 - alpha_inf)) = -k t
# We treat alpha_0 and alpha_inf as fixed reference values for the transformation.
alpha_0_ref = alpha[1]
sigma_alpha_0 = sigma_alpha[1]
alpha_inf_ref = 10.575
sigma_alpha_inf = 0.1414

ratio = (alpha .- alpha_inf_ref) ./ (alpha_0_ref - alpha_inf_ref)
valid = isfinite.(ratio) .& (ratio .> 0.0)
all(valid) || error("Linearized logarithm is not defined for all data points with the chosen alpha_inf.")

t_lin = t[valid]
alpha_lin = alpha[valid]
sigma_t_lin = sigma_t[valid]
sigma_alpha_lin = sigma_alpha[valid]
y_lin = log.(ratio[valid])

# First-order uncertainty propagation for the transformed ordinate.
inv_den = 1.0 / (alpha_0_ref - alpha_inf_ref)
dydalpha = 1.0 ./ (alpha_lin .- alpha_inf_ref)
dydalpha0 = fill(-inv_den, length(alpha_lin))
dydalphainf = -1.0 ./ (alpha_lin .- alpha_inf_ref) .+ inv_den
sigma_y_lin = sqrt.(
    (dydalpha .* sigma_alpha_lin).^2 .+
    (dydalpha0 .* sigma_alpha_0).^2 .+
    (dydalphainf .* sigma_alpha_inf).^2
)

linear_model(x, p) = @. p[1] * x + p[2]
result = fit_model(linear_model, t_lin, y_lin; p0=[-0.01, 0.0], sigma_y=sigma_y_lin)

slope, intercept = result.params
sigma_slope, sigma_intercept = result.param_stderr
k_fit = -slope
sigma_k = sigma_slope

y_fit = result.model_y
ss_res = sum((y_lin .- y_fit) .^ 2)
ss_tot = sum((y_lin .- mean(y_lin)) .^ 2)
r2 = ss_tot > 0 ? 1.0 - ss_res / ss_tot : NaN
r = cor(t_lin, y_lin)

t_grid = collect(range(minimum(t_lin), maximum(t_lin); length=500))
y_grid = linear_model(t_grid, result.params)
J = hcat(t_grid, ones(length(t_grid)))
band_sigma = sqrt.(vec(sum((J * result.param_covariance) .* J; dims=2)))

latex_theme = Theme(
    font="CMU Serif",
    fontsize=22,
    figure_padding=18,
    Axis=(
        titlesize=28,
        xlabelsize=26,
        ylabelsize=26,
        xticklabelsize=18,
        yticklabelsize=18,
        xgridvisible=true,
        ygridvisible=true,
        xminorgridvisible=true,
        yminorgridvisible=true,
    ),
    Legend=(
        labelsize=20,
        framevisible=true,
        backgroundcolor=:white,
        patchsize=(28, 20),
    ),
)

fig = with_theme(latex_theme) do
    Figure(size=(900, 500))
end

colgap!(fig.layout, 8)

ax = Axis(
    fig[1, 1],
    title=L"\text{Linearisiert: }\ln\!\left(\frac{\alpha_t-\alpha_{\infty}}{\alpha_0-\alpha_{\infty}}\right)=-kt",
    xlabel=L"\text{Zeit }t\ \text{in min}",
    ylabel=L"\ln\!\left(\frac{\alpha_t-\alpha_{\infty}}{\alpha_0-\alpha_{\infty}}\right)",
)

band = band!(
    ax,
    t_grid,
    y_grid .- band_sigma,
    y_grid .+ band_sigma;
    color=(:dodgerblue, 0.22),
    label=L"1\sigma\ \text{Band}",
)
line = lines!(ax, t_grid, y_grid; color=:darkblue, linewidth=2.5, label=L"\text{Linearer Fit}")
scatter = scatter!(ax, t_lin, y_lin; color=:black, marker=:circle, markersize=10, label=L"\text{Daten}")
errorbars!(ax, t_lin, y_lin, sigma_t_lin; direction=:x, color=:black, whiskerwidth=10)
errorbars!(ax, t_lin, y_lin, sigma_y_lin; color=:black, whiskerwidth=10)
axislegend(ax; position=:rt, framevisible=true, backgroundcolor=:white)

stats_box = Box(
    fig[1, 2];
    color=:white,
    strokecolor=:black,
    strokewidth=1.0,
)

stats_grid = GridLayout(fig[1, 2], tellwidth=true, tellheight=true)

Label(stats_grid[1, 1], L"\textbf{Fit Summary}"; tellwidth=false, halign=:center, valign=:center, fontsize=22)
#Label(stats_grid[2, 1], L"m = %$(round(slope; sigdigits=5)) \pm %$(round(sigma_slope; sigdigits=4))"; tellwidth=false, halign=:center, valign=:center, fontsize=20)
Label(stats_grid[2, 1], L"b = %$(round(intercept; sigdigits=5)) \pm %$(round(sigma_intercept; sigdigits=4))"; tellwidth=false, halign=:center, valign=:center, fontsize=20)
Label(stats_grid[3, 1], L"k = %$(round(k_fit; sigdigits=5)) \pm %$(round(sigma_k; sigdigits=4))"; tellwidth=false, halign=:center, valign=:center, fontsize=20)
Label(stats_grid[4, 1], L"R = %$(round(r; sigdigits=5))"; tellwidth=false, halign=:center, valign=:center, fontsize=20)
Label(stats_grid[5, 1], L"R^2 = %$(round(r2; sigdigits=5))"; tellwidth=false, halign=:center, valign=:center, fontsize=20)
Label(stats_grid[6, 1], L"\alpha_0 = %$(round(alpha_0_ref; sigdigits=5))"; tellwidth=false, halign=:center, valign=:center, fontsize=20)
Label(stats_grid[7, 1], L"\alpha_{\infty} = %$(round(alpha_inf_ref; sigdigits=5))"; tellwidth=false, halign=:center, valign=:center, fontsize=20)
Label(stats_grid[8, 1], ""; tellwidth=false, fontsize=4)

colsize!(fig.layout, 2, Fixed(280))

outfile = example_output("protolysierte_Reaktion_linearisiert_plot.svg")
save(outfile, fig)

println("slope m = ", slope, " ± ", sigma_slope)
println("intercept b = ", intercept, " ± ", sigma_intercept)
println("k = ", k_fit, " ± ", sigma_k)
println("R = ", r)
println("R^2 = ", r2)
println("Saved plot to ", outfile)
