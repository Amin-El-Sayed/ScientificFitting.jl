using Test

const ROOT = abspath(joinpath(@__DIR__, ".."))
const DOCS_SRC = joinpath(ROOT, "docs", "src")
const DOCS_MAKE = joinpath(ROOT, "docs", "make.jl")
const MAINTAINERS = joinpath(ROOT, "MAINTAINERS.md")

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
]

const PUBLIC_TEXT_FILES = vcat(joinpath.(DOCS_SRC, PUBLIC_DOC_PAGES), [joinpath(ROOT, "README.md")])

const FORBIDDEN_PUBLIC_PATTERNS = Pair{String, Regex}[
    "AI/LLM disclosure text" => r"(?i)\b(as an ai|chatgpt|large language model|ai[- ]?generated|ai slop)\b",
    "placeholder marker" => r"(?i)\b(todo|fixme|lorem ipsum|placeholder prose|being rewritten|not all of them are finished|work in progress|coming soon|to be written|to be added)\b",
    "draft/tutorial residue" => r"(?i)\b(draft-only|toy example|left as an exercise|why this example matters|synthetic perfect-data)\b",
    "private local path" => r"(?i)(/Users/|Documents/Projekte|private P1|P1-Praktikum|Praktikum)",
    "private author handle in public prose" => r"(?i)\bAmin_El_Sayed\b",
    "course-internal wording" => r"(?i)\b(course[- ]internal|lab-course-internal|private dataset)\b",
    "stale public API identifier" => r"\b(profile_curve|contour_grid)\b",
    "stale likelihood-scale identifier" =>
        r"\b(nll_min|gaussian_nll|poisson_nll|histogram_poisson_nll|unbinned_nll|extended_unbinned_nll)\b",
    "ungrouped TeX operator subscript" => r"_(\\min|\\max)\b",
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

function docs_source_markdown_pages()
    pages = String[]
    for (directory, _, filenames) in walkdir(DOCS_SRC)
        for filename in filenames
            endswith(filename, ".md") || continue
            push!(pages, relpath(joinpath(directory, filename), DOCS_SRC))
        end
    end
    return sort(pages)
end

function documenter_make_text()
    return read(DOCS_MAKE, String)
end

function install_page_text()
    return public_file_text(joinpath(DOCS_SRC, "install.md"))
end

function quickstart_page_text()
    return public_file_text(joinpath(DOCS_SRC, "quickstart.md"))
end

function home_page_text()
    return public_file_text(joinpath(DOCS_SRC, "index.md"))
end

function gallery_page_text()
    return public_file_text(joinpath(DOCS_SRC, "gallery.md"))
end

function architecture_page_text()
    return public_file_text(joinpath(DOCS_SRC, "how_jufitter_works.md"))
end

function practitioner_page_text()
    return public_file_text(joinpath(DOCS_SRC, "fitting_for_practitioners.md"))
end

function plotting_page_text()
    return public_file_text(joinpath(DOCS_SRC, "plotting_design.md"))
end

function statistical_foundations_page_text()
    return public_file_text(joinpath(DOCS_SRC, "statistical_foundations.md"))
end

function reference_map_page_text()
    return public_file_text(joinpath(DOCS_SRC, "overview.md"))
end

function markdown_image_alt_texts(text::AbstractString)
    return [match.captures[1] for match in eachmatch(r"!\[([^\]]*)\]\([^)]+\)", text)]
end

function html_image_tags(text::AbstractString)
    return [match.captures[1] for match in eachmatch(r"<img\b([^>]*)>", text)]
end

function html_image_alt_text(tag::AbstractString)
    match_result = match(r"\balt=\"([^\"]*)\"", tag)
    return isnothing(match_result) ? nothing : match_result.captures[1]
end

@testset "Public documentation release hygiene" begin
    @testset "Documenter navigation coverage" begin
        @test documenter_navigation_pages() == sort(PUBLIC_DOC_PAGES)
        @test docs_source_markdown_pages() == sort(PUBLIC_DOC_PAGES)
        @test !occursin(r"(?m)^\s{8}\"Engineering Notes\"\s*=>", documenter_make_text())
        @test occursin(r"(?m)^\s{12}\"Technical Notes\"\s*=>", documenter_make_text())
        @test !occursin("Maintenance Notes", documenter_make_text())
        @test !isfile(joinpath(DOCS_SRC, "maintenance.md"))
        @test isfile(MAINTAINERS)
        maintainers = read(MAINTAINERS, String)
        @test occursin("## Numerical Invariants", maintainers)
        @test occursin("## Plotting Invariants", maintainers)
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
            @testset "image alt text" begin
                for alt in markdown_image_alt_texts(text)
                    @test !isempty(strip(alt))
                end
                for tag in html_image_tags(text)
                    alt = html_image_alt_text(tag)
                    @test !isnothing(alt)
                    @test !isempty(strip(alt))
                end
            end
        end
    end

    @testset "Python interoperability documentation" begin
        text = install_page_text()
        @test occursin("juliacall", text)
        @test occursin("examples/python/fit_from_python.py", text)
        @test occursin("JUFITTER_RUN_PYTHON_INTEROP=1", text)
        @test occursin("experimental or deferred", text)
    end

    @testset "First-user path is executable and honest" begin
        home = home_page_text()
        quickstart = quickstart_page_text()
        install = install_page_text()
        readme = public_file_text(joinpath(ROOT, "README.md"))

        @test occursin("```@raw html\n<section class=\"jufitter-hero\">", home)
        @test occursin("data-jufitter-plot-group=\"home-first-fit\"", home)
        @test occursin("<code>FitResult</code> or <code>LikelihoodFitResult</code>", home)
        @test !occursin("collect(range", home)
        @test occursin("using CairoMakie", quickstart)
        @test !occursin("collect(range", quickstart)
        @test occursin("report=:both", quickstart)
        @test !occursin("println(report_text", quickstart)
        @test !occursin("Pkg.test()", install)
        @test occursin(r"CI workflow is\s+configured", install)
        @test occursin(r"must pass\s+before public release", install)
        @test !occursin("are exercised by the release test matrix", install)
        @test occursin(r"CI workflow is\s+configured", readme)
        @test occursin(r"must pass\s+before publication", readme)
        @test !occursin("canonical=", documenter_make_text())
    end

    @testset "Gallery overview is user-facing" begin
        gallery = gallery_page_text()
        @test occursin("Every gallery page follows the same pattern", gallery)
        @test occursin("## Run An Example", gallery)
        @test !occursin("## Editorial Standard", gallery)
        @test !occursin("docs_gallery_gate.jl", gallery)
        @test !occursin("docs_output_snapshots.jl", gallery)
        @test !occursin("09_docs_gallery_suite.jl", gallery)
    end

    @testset "Linear calibration is self-contained and dimensioned" begin
        calibration = public_file_text(joinpath(DOCS_SRC, "gallery", "linear_calibration.md"))
        @test occursin("controlled calibration record", calibration)
        @test occursin("filename=\"linear_calibration.pdf\"", calibration)
        @test occursin("\\mathrm{V\\,mm^{-1}}", calibration)
        @test !occursin("09_docs_gallery_suite.jl", calibration)
    end

    @testset "XY uncertainty tutorial isolates its statistical lesson" begin
        xy = public_file_text(joinpath(DOCS_SRC, "gallery", "xy_uncertainties.md"))
        @test occursin("sigma_x = fill(0.050", xy)
        @test occursin("sigma_U = fill(0.033", xy)
        @test occursin("filename=\"xy_uncertainties.pdf\"", xy)
        @test occursin("status = ok - no immediate issue", xy)
        @test occursin("\\mathrm{V\\,mm^{-1}}", xy)
        @test !occursin("asset generator", lowercase(xy))
    end

    @testset "Architecture page uses the package likelihood scale" begin
        architecture = architecture_page_text()
        @test occursin("appropriate ``-2\\log L`` objective", architecture)
        @test occursin("the ``-2\\log L`` minimum", architecture)
        @test !occursin("appropriate negative\nlog-likelihood", architecture)
        @test !occursin("plot_fit(result; report=", architecture)
    end

    @testset "Practitioner guidance follows statistical scale" begin
        guide = practitioner_page_text()

        @test occursin("there is no universal acceptable interval", guide)
        @test occursin("\\sqrt{2\\,\\mathrm{ndf}}", guide)
        @test occursin("r=(0.1,-0.1)", guide)
        @test occursin("WhiteningOperator", guide)
        @test occursin("cost=:auto", guide)
        @test occursin("profile_interval", guide)
        @test occursin("interval.profile_result", guide)
        @test occursin("cont = JuFitter.contour(", guide)
        @test !occursin(r"cont = contour\(", guide)
        @test !occursin("0.5 \\lesssim", guide)
    end

    @testset "Statistics guide makes profile geometry readable" begin
        guide = statistical_foundations_page_text()

        @test occursin("data-jufitter-plot-group=\"statistics-profile-matrix\"", guide)
        for style in ("screen", "article")
            @test occursin("data-jufitter-plot-style=\"$(style)\"", guide)
            @test occursin("saturation_profile_matrix_$(style)_light.png", guide)
            @test occursin("saturation_profile_matrix_$(style)_dark.png", guide)
        end
        @test occursin("filled one- and two-sigma profiled regions", guide)
        @test occursin("dashed local covariance ellipses", guide)
        @test occursin("correlation as a pointer, not a verdict", guide)
    end

    @testset "Plotting guide follows the public extension contract" begin
        guide = plotting_page_text()

        @test occursin("(result, figure)", guide)
        @test occursin("report=:plot", guide)
        @test occursin("show_stats=true", guide)
        @test occursin("appearance=:auto` currently resolves to the light", guide)
        @test occursin("plot_profile_matrix(matrix", guide)
        @test occursin("does not rerun the optimizer", guide)
        @test occursin("an x-y `FitResult`", guide)
        @test occursin("explicit `show_stats=true` or", guide)
        @test occursin("must provide `marginal_sigma`", guide)
        @test occursin(r"Pass a\s+`LaTeXString`", guide)
        @test occursin("pass a matching `xgrid`", guide)
        @test !occursin("being migrated", guide)
        @test !occursin("## Acceptance Tests", guide)
    end

    @testset "Statistical foundations states assumptions and limits" begin
        foundations = statistical_foundations_page_text()

        @test occursin("C(\\theta) = -2\\log L(\\theta)", foundations)
        @test occursin("\\chi^2_\\mathrm{common}=\\frac{2}{1+\\rho}=1.11", foundations)
        @test occursin("\\chi^2_\\mathrm{opposite}=\\frac{2}{1-\\rho}=10", foundations)
        @test occursin("normalized split-normal cost", foundations)
        @test occursin("cost=:gaussian_likelihood", foundations)
        @test occursin("stats.minus2loglik_min", foundations)
        @test occursin("`stats.cost_min` is the objective that was", foundations)
        @test occursin("not independent ``\\mathcal N(0,1)`` draws", foundations)
        @test occursin("plausible but unquantified bias", foundations)
        @test occursin("an auxiliary calibration ``g=1.00\\pm0.05``", foundations)
        @test occursin("linear Gaussian model with known, full-rank covariance", foundations)
        @test occursin("stored observation count `nobs`", foundations)
        @test occursin("complete, static observation covariance", foundations)
        @test occursin(r"Gaussian\s+auxiliary terms contribute both", foundations)
        @test occursin(
            "has likelihood meaning only when that objective follows the documented normalization",
            foundations,
        )
        @test occursin("fit_histogram_density", foundations)
        @test occursin("JuFitter therefore returns `NaN`", foundations)
        @test occursin("Wilks thresholds are asymptotic", foundations)
        @test occursin("same observations and event-selection domain", foundations)
        @test occursin("scale_covariance=:auto | :never | :always", foundations)
        @test !occursin("credible intervals", foundations)
    end

    @testset "Reference map preserves the public workflow boundaries" begin
        reference = reference_map_page_text()

        @test occursin("The data-generating process determines the objective", reference)
        @test occursin("named tuple `(result, figure)`", reference)
        @test occursin("actual additive covariance contributions", reference)
        @test occursin("A `WhiteningOperator` is the complete", reference)
        @test occursin("zero fitted", reference)
        @test occursin("complete parameter vector", reference)
        @test occursin("`plot_fit` does not accept `LikelihoodFitResult`", reference)
        @test occursin("profile_matrix_triage(matrix)", reference)
        @test occursin("does not rerun the", reference)
        @test occursin(r"do\s+not require Makie", reference)
        @test occursin("report=:plot | :console |", reference)
        @test !occursin("report=:plot,", reference)
    end
end
