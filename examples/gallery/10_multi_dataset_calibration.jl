using CairoMakie
using JuFitter
using LinearAlgebra
using Printf

const MULTI_OUTPUT_DIR = joinpath(@__DIR__, "..", "output")
const MULTI_DOC_ASSET_DIR = joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "gallery")

linear_channel(x, p) = @. p[1] * x + p[2]

# Controlled three-channel calibration record with realistic heteroscedastic scatter.
x_a = collect(0.0:1.0:10.0)
x_b = collect(0.5:1.0:9.5)
x_c = collect(0.0:1.25:10.0)

sigma_a = @. 0.075 + 0.008 * x_a
sigma_b = @. 0.085 + 0.006 * x_b
sigma_c = @. 0.080 + 0.009 * x_c

pattern_a = 1.65 .* [0.2, -0.7, 0.4, -0.5, 0.8, -0.3, 0.1, -0.8, 0.6, -0.2, 0.5]
pattern_b = 1.65 .* [-0.4, 0.7, -0.5, 0.1, 0.8, -0.6, 0.3, -0.2, 0.6, -0.3]
pattern_c = 1.65 .* [0.3, -0.6, 0.7, -0.4, 0.1, 0.8, -0.5, 0.4, -0.7]

y_a = linear_channel(x_a, [1.82, 0.72]) .+ sigma_a .* pattern_a
y_b = linear_channel(x_b, [1.82, -0.46]) .+ sigma_b .* pattern_b
y_c = linear_channel(x_c, [1.91, 0.12]) .+ sigma_c .* pattern_c

models = [linear_channel, linear_channel, linear_channel]
x_sets = [x_a, x_b, x_c]
y_sets = [y_a, y_b, y_c]
sigma_sets = [sigma_a, sigma_b, sigma_c]

# Hypothesis one: all three channels share one physical gain.
all_shared_result = fit_multi_model(
    models,
    x_sets,
    y_sets;
    p0=[1.85, 0.7, -0.4, 0.1],
    sigma_y=sigma_sets,
    parameter_map=[[1, 2], [1, 3], [1, 4]],
    parameter_names=["shared gain", "offset A", "offset B", "offset C"],
)

# Hypothesis two: A and B share a gain, while channel C has its own gain.
partial_shared_result = fit_multi_model(
    models,
    x_sets,
    y_sets;
    p0=[1.82, 0.7, -0.4, 1.90, 0.1],
    sigma_y=sigma_sets,
    parameter_map=[[1, 2], [1, 3], [4, 5]],
    parameter_names=["gain A/B", "offset A", "offset B", "gain C", "offset C"],
)

gain_gap_gradient = [-1.0, 0.0, 0.0, 1.0, 0.0]
gain_gap = partial_shared_result.params[4] - partial_shared_result.params[1]
sigma_gain_gap = sqrt(dot(
    gain_gap_gradient,
    partial_shared_result.param_covariance * gain_gap_gradient,
))

function multi_theme(dark::Bool)
    foreground = dark ? "#edf2f4" : "#14151a"
    muted = dark ? "#b8c1ca" : "#5b6270"
    grid = dark ? ("#2a313a", 0.85) : ("#e9eef4", 0.95)
    background = dark ? "#111318" : "#ffffff"
    return Theme(
        fontsize=20,
        font="TeX Gyre Heros",
        Figure=(backgroundcolor=background,),
        Axis=(
            backgroundcolor=background,
            xlabelsize=27,
            ylabelsize=27,
            titlesize=26,
            xticklabelsize=20,
            yticklabelsize=20,
            xlabelcolor=foreground,
            ylabelcolor=foreground,
            titlecolor=foreground,
            xticklabelcolor=muted,
            yticklabelcolor=muted,
            xtickcolor=muted,
            ytickcolor=muted,
            xgridcolor=grid,
            ygridcolor=grid,
            topspinevisible=false,
            rightspinevisible=false,
            leftspinecolor=muted,
            bottomspinecolor=muted,
        ),
        Legend=(framevisible=false, labelcolor=foreground, labelsize=19, patchsize=(30, 16)),
    )
end

function confidence_band(result, parameter_indices, x)
    local_covariance = result.param_covariance[parameter_indices, parameter_indices]
    jacobian = hcat(x, ones(length(x)))
    variance = vec(sum((jacobian * local_covariance) .* jacobian; dims=2))
    return sqrt.(clamp.(variance, 0.0, Inf))
end

function dataset_pulls(result, parameter_maps)
    return [
        (y_sets[i] .- linear_channel(x_sets[i], result.params[parameter_maps[i]])) ./ sigma_sets[i]
        for i in eachindex(x_sets)
    ]
end

function add_pull_reference!(axis, xmin, xmax, color_1sigma, color_2sigma, zero_color)
    band!(axis, [xmin, xmax], [-2.0, -2.0], [2.0, 2.0]; color=color_2sigma)
    band!(axis, [xmin, xmax], [-1.0, -1.0], [1.0, 1.0]; color=color_1sigma)
    hlines!(axis, [0.0]; color=zero_color, linewidth=1.4)
end

fmt(x, digits=4) = @sprintf("%.*g", digits, x)

function save_multi_dataset_calibration(filename; dark::Bool=false)
    foreground = dark ? "#edf2f4" : "#14151a"
    muted = dark ? "#b8c1ca" : "#5b6270"
    background = dark ? "#111318" : "#ffffff"
    colors = dark ? ["#66d9ef", "#f4b860", "#f28fad"] : ["#007f9e", "#b45309", "#b83280"]
    bands = dark ?
            [("#66d9ef", 0.15), ("#f4b860", 0.14), ("#f28fad", 0.14)] :
            [("#89d5e0", 0.30), ("#f4b183", 0.28), ("#e7a5ca", 0.28)]
    pull_1sigma = dark ? ("#66d9ef", 0.13) : ("#a8dadc", 0.30)
    pull_2sigma = dark ? ("#66d9ef", 0.06) : ("#a8dadc", 0.14)
    markers = [:circle, :rect, :diamond]
    labels = ["channel A", "channel B", "channel C"]
    partial_maps = [[1, 2], [1, 3], [4, 5]]
    all_shared_maps = [[1, 2], [1, 3], [1, 4]]

    figure = with_theme(multi_theme(dark)) do
        Figure(size=(1680, 1040), backgroundcolor=background)
    end
    fit_axis = Axis(
        figure[1, 1];
        title="Three-channel calibration transfer",
        ylabel="channel response y (V)",
    )
    shared_pull_axis = Axis(
        figure[2, 1];
        title="Pulls: all channels forced to share one gain",
        titlealign=:left,
        titlesize=19,
        ylabel="pull rᵢ",
        ylabelsize=22,
    )
    partial_pull_axis = Axis(
        figure[3, 1];
        title="Pulls: channels A/B share a gain; C is independent",
        titlealign=:left,
        titlesize=19,
        xlabel="reference input x",
        ylabel="pull rᵢ",
        ylabelsize=22,
    )

    x_grid = collect(range(0.0, 10.0; length=300))
    for i in eachindex(x_sets)
        partial_parameters = partial_shared_result.params[partial_maps[i]]
        partial_prediction = linear_channel(x_grid, partial_parameters)
        partial_sigma = confidence_band(partial_shared_result, partial_maps[i], x_grid)
        all_shared_prediction = linear_channel(x_grid, all_shared_result.params[all_shared_maps[i]])

        band!(
            fit_axis,
            x_grid,
            partial_prediction .- partial_sigma,
            partial_prediction .+ partial_sigma;
            color=bands[i],
        )
        lines!(
            fit_axis,
            x_grid,
            all_shared_prediction;
            color=(colors[i], 0.48),
            linestyle=:dash,
            linewidth=2.0,
        )
        lines!(
            fit_axis,
            x_grid,
            partial_prediction;
            color=colors[i],
            linewidth=2.8,
        )
        errorbars!(fit_axis, x_sets[i], y_sets[i], sigma_sets[i]; color=(colors[i], 0.50), whiskerwidth=5)
        scatter!(
            fit_axis,
            x_sets[i],
            y_sets[i];
            color=colors[i],
            marker=markers[i],
            markersize=9,
        )
    end
    hidexdecorations!(fit_axis; grid=false)

    shared_pulls = dataset_pulls(all_shared_result, all_shared_maps)
    partial_pulls = dataset_pulls(partial_shared_result, partial_maps)
    add_pull_reference!(shared_pull_axis, 0.0, 10.0, pull_1sigma, pull_2sigma, (foreground, 0.55))
    add_pull_reference!(partial_pull_axis, 0.0, 10.0, pull_1sigma, pull_2sigma, (foreground, 0.55))
    for i in eachindex(x_sets)
        lines!(shared_pull_axis, x_sets[i], shared_pulls[i]; color=(colors[i], 0.42), linewidth=1.2)
        scatter!(shared_pull_axis, x_sets[i], shared_pulls[i]; color=colors[i], marker=markers[i], markersize=6)
        lines!(partial_pull_axis, x_sets[i], partial_pulls[i]; color=(colors[i], 0.42), linewidth=1.2)
        scatter!(partial_pull_axis, x_sets[i], partial_pulls[i]; color=colors[i], marker=markers[i], markersize=6)
    end
    hidexdecorations!(shared_pull_axis; grid=false)
    linkxaxes!(fit_axis, shared_pull_axis, partial_pull_axis)
    ylims!(shared_pull_axis, -3.4, 3.4)
    ylims!(partial_pull_axis, -3.4, 3.4)

    side = GridLayout()
    figure[1:3, 2] = side
    legend_elements = LegendElement[
        MarkerElement(color=colors[i], marker=markers[i], markersize=12)
        for i in eachindex(labels)
    ]
    append!(
        legend_elements,
        [
            LineElement(color=foreground, linewidth=3),
            LineElement(color=(foreground, 0.55), linestyle=:dash, linewidth=2),
            PolyElement(color=dark ? ("#66d9ef", 0.18) : ("#89d5e0", 0.32)),
        ],
    )
    Legend(
        side[1, 1],
        legend_elements,
        [
            labels...,
            "partial-sharing model",
            "all-shared-gain hypothesis",
            "local 1σ fit band",
        ];
        framevisible=false,
        tellheight=true,
        nbanks=2,
        labelsize=17,
    )
    report_axis = Axis(side[2, 1]; backgroundcolor=:transparent)
    hidedecorations!(report_axis)
    hidespines!(report_axis)
    text!(
        report_axis,
        0,
        1;
        text="partial-sharing result\n" *
             "  gain A/B = $(fmt(partial_shared_result.params[1], 6)) ± $(fmt(partial_shared_result.param_stderr[1], 2))\n" *
             "  gain C = $(fmt(partial_shared_result.params[4], 6)) ± $(fmt(partial_shared_result.param_stderr[4], 2))\n" *
             "  gain C − gain A/B\n" *
             "    $(fmt(gain_gap, 5)) ± $(fmt(sigma_gain_gap, 2))\n" *
             "    $(fmt(gain_gap / sigma_gain_gap, 3))σ from zero\n\n" *
             "all-shared-gain hypothesis\n" *
             "  χ²/ndf = $(fmt(all_shared_result.stats.chi2_ndf, 4))\n" *
             "  P(χ²) = $(fmt(all_shared_result.stats.pvalue, 4))\n" *
             "  AIC = $(fmt(all_shared_result.stats.aic, 5))\n\n" *
             "partial-sharing model\n" *
             "  χ²/ndf = $(fmt(partial_shared_result.stats.chi2_ndf, 4))\n" *
             "  P(χ²) = $(fmt(partial_shared_result.stats.pvalue, 4))\n" *
             "  AIC = $(fmt(partial_shared_result.stats.aic, 5))\n\n" *
             "ΔAIC = $(fmt(all_shared_result.stats.aic - partial_shared_result.stats.aic, 5))\n" *
             "Do not transfer channel C's gain.",
        space=:relative,
        align=(:left, :top),
        color=foreground,
        fontsize=18,
        lineheight=1.08,
    )

    rowsize!(side, 1, Auto())
    rowsize!(side, 2, Relative(1))
    rowsize!(figure.layout, 1, Relative(0.58))
    rowsize!(figure.layout, 2, Relative(0.21))
    rowsize!(figure.layout, 3, Relative(0.21))
    colsize!(figure.layout, 2, Fixed(600))
    save(filename, figure)
end

mkpath(MULTI_OUTPUT_DIR)
mkpath(MULTI_DOC_ASSET_DIR)

for (dark, suffix) in ((false, "light"), (true, "dark"))
    save_multi_dataset_calibration(
        joinpath(MULTI_OUTPUT_DIR, "10_multi_dataset_calibration_$(suffix).png");
        dark=dark,
    )
    save_multi_dataset_calibration(
        joinpath(MULTI_DOC_ASSET_DIR, "multi_dataset_shared_slope_$(suffix).png");
        dark=dark,
    )
end

println("All-shared-gain hypothesis")
println(report_text(all_shared_result))
println()
println("Partial-sharing model")
println(report_text(partial_shared_result))
println("gain C - gain A/B = ", gain_gap, " +/- ", sigma_gain_gap)
