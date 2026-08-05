using Test
using JuFitter

const ROOT = abspath(joinpath(@__DIR__, ".."))
const API_PAGES = [
    joinpath(ROOT, "docs", "src", "api.md"),
    joinpath(ROOT, "docs", "src", "api_fitting.md"),
    joinpath(ROOT, "docs", "src", "api_results.md"),
    joinpath(ROOT, "docs", "src", "api_plotting.md"),
]
const PUBLIC_API_DOC_EXEMPTIONS = Set([:JuFitter])
const FIT_ENTRYPOINTS = [
    :fit_model,
    :fit_custom,
    :fit_poisson_model,
    :fit_histogram_model,
    :fit_histogram_density,
    :fit_unbinned_model,
    :fit_extended_unbinned_model,
    :fit_indexed_model,
    :fit_multi_model,
]

function _public_exports()
    return sort!(setdiff(names(JuFitter; all=false), collect(PUBLIC_API_DOC_EXEMPTIONS)); by=string)
end

function _doc_text(name::Symbol)
    doc = @eval JuFitter (@doc $(name))
    return sprint(show, MIME("text/plain"), doc)
end

function _has_public_docstring(name::Symbol)
    text = strip(_doc_text(name))
    isempty(text) && return false
    text == "nothing" && return false
    occursin("No documentation found", text) && return false
    return true
end

function _api_page_text(path)
    isfile(path) || error("API reference page missing: $(relpath(path, ROOT))")
    return read(path, String)
end

function _api_reference_text()
    return join(_api_page_text.(API_PAGES), "\n")
end

@testset "Public API reference docstrings" begin
    exports = _public_exports()
    @test !isempty(exports)

    missing = Symbol[name for name in exports if !_has_public_docstring(name)]
    @test missing == Symbol[]

    overview_text = _api_page_text(first(API_PAGES))
    api_text = _api_reference_text()
    undocumented_on_page = Symbol[name for name in exports if !occursin(string(name), api_text)]
    @test undocumented_on_page == Symbol[]

    @testset "Reference page states the public contracts" begin
        required_sections = [
            "## Choose An Entry Point",
            "## Common Conventions",
            "## Gaussian Fits",
            "## Likelihood And Count Fits",
            "## Results",
            "## Profiles And Contours",
            "## Diagnostics And Reports",
            "# Plotting",
        ]
        @test all(section -> occursin(section, api_text), required_sections)

        required_contracts = [
            "one-based indices",
            "complete parameter vector",
            "inplace=true",
            "`whitening` is intentionally exclusive",
            "converged == false",
            "normalized ``-2\\log L``",
            "different uncertainty scales",
            "`gof(p)` is the data goodness-of-fit statistic",
            "on_failure=:throw",
            "A side that is not bracketed is returned as",
            "max_actions=5",
            "do not require Makie",
        ]
        @test all(contract -> occursin(contract, api_text), required_contracts)
        @test !occursin("ci_level", api_text)
        @test occursin("`:ok`, `:review`, or `:stop`", api_text)
        @test occursin("For `plot_fit`, use `show_stats`", api_text)
        @test occursin("zero fitted", api_text)
        @test occursin("[Fitting](api_fitting.md)", overview_text)
        @test occursin("[Results And Diagnostics](api_results.md)", overview_text)
        @test occursin("[Plotting](api_plotting.md)", overview_text)
    end

    @testset "Optional plotting boundary matches real methods" begin
        @test occursin("(result::FitResult, figure::Figure)", _doc_text(:fitplot))
        @test !occursin("plot_fit(model, x, y", _doc_text(:plot_fit))
        @test occursin("fit_axis(figure; index=1)", _doc_text(:fit_axis))
        @test occursin("plot_info_panel!(cell;", _doc_text(:plot_info_panel!))
    end

    @testset "Custom objectives state their inferential convention" begin
        custom_doc = _doc_text(:fit_custom)
        @test occursin("normalized `-2log(L)` cost", custom_doc)
        @test occursin("only arithmetic summaries", custom_doc)
        @test occursin("nobs", custom_doc)
    end

    @testset "Fit entry points state complete call contracts" begin
        for name in FIT_ENTRYPOINTS
            doc = _doc_text(name)
            @test occursin("->", doc)
            @test occursin("Result", doc)
            has_failure_contract = occursin("raise", lowercase(doc)) ||
                                   occursin("fail", lowercase(doc)) ||
                                   occursin("error", lowercase(doc))
            @test has_failure_contract
            @test occursin("$(name)", api_text)
        end

        @test occursin("# Example", _doc_text(:fit_model))
        @test occursin("# Example", _doc_text(:fit_custom))
        @test occursin("Minimal call shapes", overview_text)
        @test all(name -> occursin("`$(name)` | `$(name)(", overview_text), FIT_ENTRYPOINTS)
    end

    @testset "Constraint and statistic docstrings state their scope" begin
        @test occursin("complete parameter", _doc_text(:ConstraintSpec))
        likelihood_problem_doc = _doc_text(:LikelihoodFitProblem)
        @test occursin("complete", likelihood_problem_doc) &&
              occursin("fixed parameters", likelihood_problem_doc)
        @test occursin("custom loss", _doc_text(:FitStatistics))
        @test occursin(":model_relative", _doc_text(:ErrorComponent))
        @test occursin("zero fitted", _doc_text(:FixedParameter))
    end

    @testset "Fit options expose only effective controls" begin
        @test !hasfield(FitOptions, :ci_level)
        @test !occursin("ci_level", _doc_text(:FitOptions))
    end
end
