const MULTI_SNAPSHOT_ONLY = get(ENV, "JUFITTER_DOC_SNAPSHOT_ONLY", "0") == "1"
const MULTI_RENDER_PLOTS = !MULTI_SNAPSHOT_ONLY
const MULTI_RENDER_DOC_ASSETS = MULTI_RENDER_PLOTS &&
    get(ENV, "JUFITTER_RENDER_DOC_ASSETS", "0") == "1"

if MULTI_RENDER_PLOTS
    using CairoMakie
end

using Distributions
using JuFitter
using LinearAlgebra
using Printf

const MULTI_OUTPUT_DIR = joinpath(@__DIR__, "..", "output")
const MULTI_DOC_ASSET_DIR = joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "gallery")
const MULTI_EMIT_DOC_OUTPUT_SNAPSHOTS = get(ENV, "JUFITTER_DOC_OUTPUT_SNAPSHOTS", "0") == "1"
const MULTI_RENDER_WIDTH = 1280
const MULTI_PX_PER_UNIT = 2.0

function emit_multi_doc_output_snapshot(body::Function, id::AbstractString)
    MULTI_EMIT_DOC_OUTPUT_SNAPSHOTS || return nothing

    println("=== JUFITTER_DOC_OUTPUT_BEGIN ", id, " ===")
    body()
    println("=== JUFITTER_DOC_OUTPUT_END ", id, " ===")
    return nothing
end

linear_channel(x, p) = @. p[1] * x + p[2]

# Controlled teaching record; the arrays are the complete fit input.
x_a = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
y_a = [
    0.744750, 2.444135, 4.420060, 6.098325, 8.141240, 9.763075,
    11.660295, 13.287080, 15.417610, 17.051490, 19.047875,
]
sigma_a = [
    0.075, 0.083, 0.091, 0.099, 0.107, 0.115,
    0.123, 0.131, 0.139, 0.147, 0.155,
]

x_b = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]
y_b = [
    0.391920, 2.378570, 4.007500, 5.927490, 7.877840,
    9.433180, 11.431380, 13.147100, 15.144640, 16.759710,
]
sigma_b = [0.088, 0.094, 0.100, 0.106, 0.112, 0.118, 0.124, 0.130, 0.136, 0.142]

x_c = [0.0, 1.25, 2.5, 3.75, 5.0, 6.25, 7.5, 8.75, 10.0]
y_c = [
    0.159600, 2.417162, 5.013387, 7.207425, 9.690625,
    12.237350, 14.323312, 16.937275, 19.023650,
]
sigma_c = [0.08000, 0.09125, 0.10250, 0.11375, 0.12500,
           0.13625, 0.14750, 0.15875, 0.17000]

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
delta_chi2 = all_shared_result.stats.chi2 - partial_shared_result.stats.chi2
nested_pvalue = ccdf(Chisq(1), delta_chi2)

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

function save_multi_dataset_calibration(
    filename;
    dark::Union{Nothing, Bool}=nothing,
    style::Symbol=:screen,
    appearance::Symbol=dark === nothing ? :light : (dark ? :dark : :light),
)
    MULTI_RENDER_PLOTS || return nothing

    palette = plot_palette(style; appearance=appearance)
    foreground = palette.stats_color
    colors = collect(palette.series_colors[1:3])
    bands = [(color, max(0.10, 0.70 * palette.band_alpha)) for color in colors]
    pull_1sigma = (palette.band_color, max(0.12, 0.70 * palette.band_alpha))
    pull_2sigma = (palette.band_color, max(0.06, 0.35 * palette.band_alpha))
    markers = [:circle, :rect, :diamond]
    labels = ["channel A", "channel B", "channel C"]
    partial_maps = [[1, 2], [1, 3], [4, 5]]
    all_shared_maps = [[1, 2], [1, 3], [1, 4]]

    figure = with_theme(plot_theme(style; appearance=appearance)) do
        Figure(size=(MULTI_RENDER_WIDTH, 860), backgroundcolor=palette.background_color)
    end
    fit_axis = Axis(
        figure[1, 1];
        title="Three-channel calibration transfer",
        ylabel="channel response y (V)",
    )
    shared_pull_axis = Axis(
        figure[2, 1];
        title="Pulls: shared gain",
        titlealign=:left,
        titlesize=palette.subplot_titlesize,
        titlegap=palette.subplot_titlegap,
        ylabel="pull rᵢ",
    )
    partial_pull_axis = Axis(
        figure[3, 1];
        title="Pulls: partial sharing",
        titlealign=:left,
        titlesize=palette.subplot_titlesize,
        titlegap=palette.subplot_titlegap,
        xlabel="reference input x",
        ylabel="pull rᵢ",
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
            linewidth=max(1.5, palette.fit_linewidth - 0.5),
        )
        lines!(
            fit_axis,
            x_grid,
            partial_prediction;
            color=colors[i],
            linewidth=palette.fit_linewidth,
        )
        errorbars!(
            fit_axis,
            x_sets[i],
            y_sets[i],
            sigma_sets[i];
            color=(colors[i], 0.50),
            whiskerwidth=palette.error_whiskerwidth,
        )
        scatter!(
            fit_axis,
            x_sets[i],
            y_sets[i];
            color=colors[i],
            marker=markers[i],
            markersize=palette.data_markersize + 1.2,
        )
    end
    hidexdecorations!(fit_axis; grid=false)

    shared_pulls = dataset_pulls(all_shared_result, all_shared_maps)
    partial_pulls = dataset_pulls(partial_shared_result, partial_maps)
    pull_linewidth = max(1.4, 0.45 * palette.fit_linewidth)
    pull_markersize = max(5.5, 0.60 * palette.data_markersize)
    add_pull_reference!(shared_pull_axis, 0.0, 10.0, pull_1sigma, pull_2sigma, (foreground, 0.55))
    add_pull_reference!(partial_pull_axis, 0.0, 10.0, pull_1sigma, pull_2sigma, (foreground, 0.55))
    for i in eachindex(x_sets)
        lines!(
            shared_pull_axis,
            x_sets[i],
            shared_pulls[i];
            color=(colors[i], 0.52),
            linewidth=pull_linewidth,
        )
        scatter!(
            shared_pull_axis,
            x_sets[i],
            shared_pulls[i];
            color=colors[i],
            marker=markers[i],
            markersize=pull_markersize,
        )
        lines!(
            partial_pull_axis,
            x_sets[i],
            partial_pulls[i];
            color=(colors[i], 0.52),
            linewidth=pull_linewidth,
        )
        scatter!(
            partial_pull_axis,
            x_sets[i],
            partial_pulls[i];
            color=colors[i],
            marker=markers[i],
            markersize=pull_markersize,
        )
    end
    hidexdecorations!(shared_pull_axis; grid=false)
    linkxaxes!(fit_axis, shared_pull_axis, partial_pull_axis)
    ylims!(shared_pull_axis, -3.4, 3.4)
    ylims!(partial_pull_axis, -3.4, 3.4)

    legend_elements = LegendElement[
        MarkerElement(color=colors[i], marker=markers[i], markersize=12)
        for i in eachindex(labels)
    ]
    append!(
        legend_elements,
        [
            LineElement(color=foreground, linewidth=palette.fit_linewidth),
            LineElement(color=(foreground, 0.55), linestyle=:dash, linewidth=max(1.5, palette.fit_linewidth - 0.5)),
            PolyElement(color=(palette.band_color, max(0.12, palette.band_alpha))),
        ],
    )
    legend_labels = [
        labels...,
        "partial-sharing model",
        "all-shared-gain hypothesis",
        "local 1σ fit band",
    ]
    if style == :article
        Legend(
            figure[1:3, 2],
            legend_elements,
            legend_labels;
            framevisible=false,
            tellheight=false,
            halign=:left,
            valign=:top,
            nbanks=1,
            labelsize=palette.legend_labelsize,
            patchsize=palette.legend_patchsize,
            rowgap=palette.legend_rowgap,
        )
    else
        plot_info_panel!(
            figure[1:3, 2];
            legend_plots=legend_elements,
            legend_labels=legend_labels,
            legend_kwargs=(nbanks=2,),
            title="Partial-sharing result",
            parameter_lines=[
                "gain A/B = $(fmt(partial_shared_result.params[1], 6)) ± $(fmt(partial_shared_result.param_stderr[1], 2))",
                "gain C = $(fmt(partial_shared_result.params[4], 6)) ± $(fmt(partial_shared_result.param_stderr[4], 2))",
                "gain C − gain A/B = $(fmt(gain_gap, 5)) ± $(fmt(sigma_gain_gap, 2))",
                "difference = $(fmt(gain_gap / sigma_gain_gap, 3))σ from zero",
            ],
            statistic_lines=[
                "all-shared: χ²/ndf = $(fmt(all_shared_result.stats.chi2_ndf, 4)), P(χ²) = $(fmt(all_shared_result.stats.pvalue, 4))",
                "all-shared: AIC = $(fmt(all_shared_result.stats.aic, 5))",
                "partial-sharing: χ²/ndf = $(fmt(partial_shared_result.stats.chi2_ndf, 4)), P(χ²) = $(fmt(partial_shared_result.stats.pvalue, 4))",
                "partial-sharing: AIC = $(fmt(partial_shared_result.stats.aic, 5))",
                "nested test: Δχ² = $(fmt(delta_chi2, 5)), p = $(fmt(nested_pvalue, 3))",
                "ΔAIC = $(fmt(all_shared_result.stats.aic - partial_shared_result.stats.aic, 5))",
                "Do not transfer channel C's gain.",
            ],
            theme=style,
            appearance=appearance,
        )
    end

    rowsize!(figure.layout, 1, Relative(0.58))
    rowsize!(figure.layout, 2, Relative(0.21))
    rowsize!(figure.layout, 3, Relative(0.21))
    colgap!(figure.layout, 24)
    save(filename, figure; px_per_unit=MULTI_PX_PER_UNIT)
end

if MULTI_RENDER_PLOTS
    mkpath(MULTI_OUTPUT_DIR)

    for (dark, suffix) in ((false, "light"), (true, "dark"))
        save_multi_dataset_calibration(
            joinpath(MULTI_OUTPUT_DIR, "10_multi_dataset_calibration_$(suffix).png");
            dark=dark,
        )
    end
end

if MULTI_RENDER_DOC_ASSETS
    mkpath(MULTI_DOC_ASSET_DIR)
    for style in (:screen, :article), appearance in (:light, :dark)
        save_multi_dataset_calibration(
            joinpath(MULTI_DOC_ASSET_DIR, "multi_dataset_shared_slope_$(style)_$(appearance).png");
            style=style,
            appearance=appearance,
        )
    end
end

println("All-shared-gain hypothesis")
println(report_text(all_shared_result))
println()
println("Partial-sharing model")
println(report_text(partial_shared_result))
@printf("gain C - gain A/B = %.5f +/- %.5f\n", gain_gap, sigma_gain_gap)
@printf("nested test: delta chi2 = %.5f for 1 dof, p = %.4g\n",
        delta_chi2, nested_pvalue)
println()
println("All-shared diagnostic dashboard")
println(diagnostic_dashboard_text(all_shared_result))
println("Partial-sharing diagnostic dashboard")
println(diagnostic_dashboard_text(partial_shared_result))
emit_multi_doc_output_snapshot("multi_dataset") do
    println("All-shared-gain hypothesis")
    println(report_text(all_shared_result))
    println()
    println("Partial-sharing model")
    println(report_text(partial_shared_result))
    @printf("gain C - gain A/B = %.5f +/- %.5f\n", gain_gap, sigma_gain_gap)
    @printf("nested test: delta chi2 = %.5f for 1 dof, p = %.4g\n",
            delta_chi2, nested_pvalue)
    println()
    println("All-shared diagnostic dashboard")
    println(diagnostic_dashboard_text(all_shared_result))
    println("Partial-sharing diagnostic dashboard")
    println(diagnostic_dashboard_text(partial_shared_result))
end
