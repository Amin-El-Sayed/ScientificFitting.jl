using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using JuFitter
include(joinpath(@__DIR__, "_example_utils.jl"))
using CairoMakie
using LaTeXStrings
using Statistics

beta = [0.05, 0.1, 0.15, 0.2]
alpha_1 = [19.9, 23.5, 32.3, 37.75]
alpha_2 = [19.8, 23.4, 32.2, 37.6]
alpha = (alpha_1 .+ alpha_2) ./ 2

sigma_beta = fill(0.005, length(beta))
sigma_alpha_1 = fill(0.1, length(alpha_1))
sigma_alpha_2 = fill(0.1, length(alpha_2))
sigma_alpha = sqrt.(sigma_alpha_1.^2 .+ sigma_alpha_2.^2)

# Deliberately simple "linearized" treatment:
# subtract a fixed offset and fit a straight line through the origin
#   alpha_corr = (alpha* * l) * beta
l_ref = 0.2
sigma_l = 0.005
alpha_offset_ref = 11.95
sigma_alpha_offset = 0.1414

alpha_corr = alpha .- alpha_offset_ref
sigma_alpha_corr = sqrt.(sigma_alpha .^ 2 .+ sigma_alpha_offset^2)

linear_model(x, p) = @. p[1] * x
result = fit_model(linear_model, beta, alpha_corr; p0=[100.0], sigma_y=sigma_alpha_corr)

slope = result.params[1]
sigma_slope = result.param_stderr[1]
alpha_star = slope / l_ref
sigma_alpha_star = alpha_star * sqrt((sigma_slope / slope)^2 + (sigma_l / l_ref)^2)

y_fit = result.model_y
ss_res = sum((alpha_corr .- y_fit) .^ 2)
ss_tot = sum((alpha_corr .- mean(alpha_corr)) .^ 2)
r2 = ss_tot > 0 ? 1.0 - ss_res / ss_tot : NaN
r = cor(beta, alpha_corr)

beta_grid = collect(range(minimum(beta), maximum(beta); length=500))
y_grid = linear_model(beta_grid, result.params)
band_sigma = abs.(beta_grid) .* sigma_slope

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
    Figure(size=(700, 500))
end

colgap!(fig.layout, 8)

ax = Axis(
    fig[1, 1],
    title=L"\text{Linearisiert: }\alpha = \alpha^\ast l\,\beta", #\alpha_{\mathrm{offset}}
    xlabel=L"\text{Konzentration }\beta\ \text{in g/mL}",
    ylabel=L"\alpha\ \text{in }^\circ",
)

band!(
    ax,
    beta_grid,
    y_grid .- band_sigma,
    y_grid .+ band_sigma;
    color=(:dodgerblue, 0.22),
    label=L"1\sigma\ \text{Band}",
)
lines!(ax, beta_grid, y_grid; color=:darkblue, linewidth=2.5, label=L"\text{Linearer Fit}")
scatter!(ax, beta, alpha_corr; color=:black, marker=:circle, markersize=10, label=L"\text{Daten}")
errorbars!(ax, beta, alpha_corr, sigma_beta; direction=:x, color=:black, whiskerwidth=10)
errorbars!(ax, beta, alpha_corr, sigma_alpha_corr; color=:black, whiskerwidth=10)
axislegend(ax; position=:lt, framevisible=true, backgroundcolor=:white)

Box(
    fig[1, 2];
    color=:white,
    strokecolor=:black,
    strokewidth=1.0,
)

stats_grid = GridLayout(fig[1, 2], tellwidth=true, tellheight=true)
Label(stats_grid[1, 1], L"\textbf{Fit Summary}"; tellwidth=false, halign=:center, valign=:top, fontsize=22)
#Label(stats_grid[2, 1], L"m = %$(round(slope; sigdigits=5)) \pm %$(round(sigma_slope; sigdigits=4))"; tellwidth=false, halign=:center, valign=:top, fontsize=20)
Label(stats_grid[2, 1], L"\alpha^\ast = %$(round(alpha_star; sigdigits=5)) \pm %$(round(sigma_alpha_star; sigdigits=4))"; tellwidth=false, halign=:center, valign=:top, fontsize=20)
Label(stats_grid[3, 1], L"l = %$(round(l_ref; sigdigits=5)) \pm %$(round(sigma_l; sigdigits=4))"; tellwidth=false, halign=:center, valign=:top, fontsize=20)
#Label(stats_grid[5, 1], L"\alpha_{\mathrm{offset}} = %$(round(alpha_offset_ref; sigdigits=5))"; tellwidth=false, halign=:center, valign=:top, fontsize=20)
Label(stats_grid[4, 1], L"R = %$(round(r; sigdigits=5))"; tellwidth=false, halign=:center, valign=:top, fontsize=20)
Label(stats_grid[5, 1], L"R^2 = %$(round(r2; sigdigits=5))"; tellwidth=false, halign=:center, valign=:top, fontsize=20)
Label(stats_grid[6, 1], ""; tellwidth=false, fontsize=4)

colsize!(fig.layout, 2, Fixed(200))

outfile = example_output("Drehwinkel_Konzentration.svg")
save(outfile, fig)

println("slope m = ", slope, " ± ", sigma_slope)
println("alpha_star = ", alpha_star, " ± ", sigma_alpha_star)
println("R = ", r)
println("R^2 = ", r2)
println("Saved plot to ", outfile)
