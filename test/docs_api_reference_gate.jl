using Test
using ScientificFitting
import REPL # Activate Julia's docstring parser in non-interactive test runs.

const ROOT = abspath(joinpath(@__DIR__, ".."))
const API_PAGES = [
    joinpath(ROOT, "docs", "src", "api.md"),
    joinpath(ROOT, "docs", "src", "api_fitting.md"),
    joinpath(ROOT, "docs", "src", "api_results.md"),
    joinpath(ROOT, "docs", "src", "api_plotting.md"),
    joinpath(ROOT, "docs", "src", "api_plotting_diagnostics.md"),
]
const PUBLIC_API_DOC_EXEMPTIONS = Set([:ScientificFitting])

function _public_exports()
    return sort!(setdiff(names(ScientificFitting; all=false), collect(PUBLIC_API_DOC_EXEMPTIONS)); by=string)
end

function _doc_text(name::Symbol)
    doc = @eval ScientificFitting (@doc $(name))
    # Check source markup, not terminal wrapping or stripped code/heading markers.
    return sprint(show, MIME("text/markdown"), doc)
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

@testset "Public API reference docstrings" begin
    exports = _public_exports()
    @test !isempty(exports)

    missing_docstrings = Symbol[name for name in exports if !_has_public_docstring(name)]
    @test missing_docstrings == Symbol[]

    api_text = join(_api_page_text.(API_PAGES), "\n")
    undocumented_on_page = Symbol[name for name in exports if !occursin(string(name), api_text)]
    @test undocumented_on_page == Symbol[]
end
