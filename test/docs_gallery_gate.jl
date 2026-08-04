using Test

const ROOT = abspath(joinpath(@__DIR__, ".."))
const DOCS_SRC = joinpath(ROOT, "docs", "src")
const GALLERY_SRC = joinpath(DOCS_SRC, "gallery")

const GALLERY_PAGES = [
    "linear_calibration.md",
    "xy_uncertainties.md",
    "full_covariance.md",
    "resonance_decay.md",
    "photoelectric_threshold.md",
    "constraints_profiles.md",
    "poisson_histogram.md",
    "multi_dataset.md",
]

function read_gallery_page(name)
    path = joinpath(GALLERY_SRC, name)
    return path, read(path, String)
end

function has_any(text, patterns)
    return any(pattern -> occursin(pattern, text), patterns)
end

function image_sources(text)
    return [match.captures[1] for match in eachmatch(r"<img[^>]+src=\"([^\"]+)\"", text)]
end

function style_groups(text)
    groups = Dict{String, Set{Tuple{String,String}}}()
    pattern = r"data-jufitter-plot-group=\"([^\"]+)\"[^>]+data-jufitter-plot-style=\"([^\"]+)\"[^>]+src=\"([^\"]+)\""
    for match in eachmatch(pattern, text)
        group, style, src = match.captures
        appearance = occursin("_dark", src) ? "dark" : occursin("_light", src) ? "light" : ""
        get!(groups, group, Set{Tuple{String,String}}())
        push!(groups[group], (style, appearance))
    end
    return groups
end

function resolve_doc_asset(page_path, src)
    return normpath(joinpath(dirname(page_path), src))
end

@testset "Documentation gallery release gate" begin
    required_style_pairs = Set(
        (style, appearance)
        for style in ("lab", "modern", "article")
        for appearance in ("light", "dark")
    )

    for page in GALLERY_PAGES
        path, text = read_gallery_page(page)

        @testset "$page" begin
            @test occursin(r"^# "m, text)
            @test has_any(text, [r"^## .*Question"m, r"^## Scientific Question"m])
            @test has_any(text, [r"^## Data"m, r"^## The Measurement"m])
            @test has_any(text, [r"^## .*Model"m, r"^## .*Cost"m, r"^## Poisson Likelihood"m])
            @test has_any(text, [r"```julia"m, r"^## Complete .*Code"m, r"^## Complete .*Fit"m])
            @test occursin("jufitter-cell-output", text)
            @test occursin("Fit diagnostic dashboard", text)
            @test has_any(text, [r"^## Diagnostics"m, r"^## What To Inspect"m, r"^## Diagnose"m])
            @test has_any(text, [r"^## Interpretation"m, r"^## Decision"m, r"^## Read "m, r"^## Reading "m, r"^## Why Compare"m])
            @test has_any(text, [r"^## What Can Go Wrong"m, r"^## Failure Modes"m, r"^## What To Do Before"m])
            @test occursin("1σ", text) || occursin("1-sigma", text) || occursin("one-sigma", text)
            @test !occursin("P1", text)
            @test !occursin("Praktikum", text)
            @test !occursin("placeholder prose", lowercase(text))
            @test !occursin("being rewritten", text)
            @test !occursin("not all of them are finished", text)

            sources = image_sources(text)
            @test !isempty(sources)
            for src in sources
                @test isfile(resolve_doc_asset(path, src))
            end

            groups = style_groups(text)
            @test !isempty(groups)
            for (_, pairs) in groups
                @test required_style_pairs ⊆ pairs
            end
        end
    end

    overview = read(joinpath(DOCS_SRC, "gallery.md"), String)
    @test !occursin("being rewritten", overview)
    @test !occursin("remaining documentation passes", overview)
    @test !occursin("not all of them are finished", overview)
    @test !occursin("toy-like synthetic examples", overview)

    expected_path = [
        "Linear calibration",
        "XY uncertainties",
        "Full covariance",
        "Damped oscillator",
        "Photoelectric work function",
        "Poisson and histograms",
        "Constraints and profiles",
        "Multi-dataset fit",
    ]
    card_positions = [findfirst(">$(title)</a></h3>", overview) for title in expected_path]
    @test all(!isnothing, card_positions)
    @test issorted(first.(card_positions))

    make_source = read(joinpath(ROOT, "docs", "make.jl"), String)
    navigation_pages = [
        "linear_calibration",
        "xy_uncertainties",
        "full_covariance",
        "resonance_decay",
        "photoelectric_threshold",
        "poisson_histogram",
        "constraints_profiles",
        "multi_dataset",
    ]
    navigation_positions = [
        findfirst("gallery/$(page).md", make_source) for page in navigation_pages
    ]
    @test all(!isnothing, navigation_positions)
    @test issorted(first.(navigation_positions))

    for step in (
        "First fit",
        "Measured x",
        "Shared noise",
        "Model criticism",
        "Derived quantity",
        "Count data",
        "Nonlinear uncertainty",
        "Shared hypotheses",
    )
        @test occursin("<strong>$(step)</strong>", overview)
    end

    full_covariance = read(joinpath(GALLERY_SRC, "full_covariance.md"), String)
    @test occursin("exp(-p[2] * t)", full_covariance)
    @test !occursin("exp(p[2] * t)", full_covariance)
    @test occursin("sigma_stat^2 * (i == j)", full_covariance)
    @test occursin("sigma_corr^2 * exp", full_covariance)
    @test occursin("filename=\"full_covariance_decay.pdf\"", full_covariance)

    photoelectric = read(joinpath(GALLERY_SRC, "photoelectric_threshold.md"), String)
    @test occursin("photoelectric_slope = me - mb", photoelectric)
    @test occursin("(nu_THz - reference_frequency_THz)", photoelectric)
    @test occursin("h_fit = photoelectric_slope * elementary_charge", photoelectric)
    @test occursin("work_function_eV = photoelectric_slope * threshold_THz", photoelectric)
    @test !occursin("h_fit = me *", photoelectric)
    @test !occursin("examples/gallery/09_docs_gallery_suite.jl", photoelectric)

    gallery_generator = read(
        joinpath(ROOT, "examples", "gallery", "09_docs_gallery_suite.jl"),
        String,
    )
    @test occursin(
        "photoelectric_slope = emission_slope - baseline_slope",
        gallery_generator,
    )
    @test !occursin("h_fit = emission_slope *", gallery_generator)
end
