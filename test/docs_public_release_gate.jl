using Test

const ROOT = abspath(joinpath(@__DIR__, ".."))
const DOCS_SRC = joinpath(ROOT, "docs", "src")
const DOCS_MAKE = joinpath(ROOT, "docs", "make.jl")

const PUBLIC_DOC_PAGES = [
    "index.md",
    "install.md",
    "quickstart.md",
    "how_jufitter_works.md",
    "gallery.md",
    "gallery/linear_calibration.md",
    "gallery/xy_uncertainties.md",
    "gallery/full_covariance.md",
    "gallery/resonance_decay.md",
    "gallery/photoelectric_threshold.md",
    "gallery/constraints_profiles.md",
    "gallery/poisson_histogram.md",
    "gallery/multi_dataset.md",
    "fitting_for_practitioners.md",
    "plotting_design.md",
    "statistical_foundations.md",
    "api.md",
    "overview.md",
    "backend_design.md",
    "performance.md",
    "maintenance.md",
]

const PUBLIC_TEXT_FILES = vcat(joinpath.(DOCS_SRC, PUBLIC_DOC_PAGES), [joinpath(ROOT, "README.md")])

const FORBIDDEN_PUBLIC_PATTERNS = Pair{String, Regex}[
    "AI/LLM disclosure text" => r"(?i)\b(as an ai|chatgpt|large language model|ai[- ]?generated|ai slop)\b",
    "placeholder marker" => r"(?i)\b(todo|fixme|lorem ipsum|placeholder prose|being rewritten|not all of them are finished)\b",
    "private local path" => r"(?i)(/Users/|Documents/Projekte|private P1|P1-Praktikum|Praktikum)",
    "private author handle in public prose" => r"(?i)\bAmin_El_Sayed\b",
    "course-internal wording" => r"(?i)\b(course[- ]internal|lab-course-internal|private dataset)\b",
    "stale public API identifier" => r"\b(profile_curve|contour_grid)\b",
]

function public_file_text(path)
    isfile(path) || error("public documentation file missing: $(relpath(path, ROOT))")
    return read(path, String)
end

function documenter_navigation_pages()
    text = read(DOCS_MAKE, String)
    pages = String[]
    for match in eachmatch(r"\"([^\"]+\.md)\"", text)
        push!(pages, match.captures[1])
    end
    return sort(unique(pages))
end

@testset "Public documentation release hygiene" begin
    @testset "Documenter navigation coverage" begin
        @test documenter_navigation_pages() == sort(PUBLIC_DOC_PAGES)
    end

    @testset "Configured public files exist" begin
        for path in PUBLIC_TEXT_FILES
            @test isfile(path)
        end
    end

    for path in PUBLIC_TEXT_FILES
        text = public_file_text(path)
        rel = relpath(path, ROOT)
        @testset "$rel" begin
            @test !isempty(strip(text))
            for (label, pattern) in FORBIDDEN_PUBLIC_PATTERNS
                @testset "$label" begin
                    @test !occursin(pattern, text)
                end
            end
        end
    end
end
