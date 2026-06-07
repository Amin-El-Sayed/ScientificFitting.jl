using Test
using JuFitter

const ROOT = abspath(joinpath(@__DIR__, ".."))
const API_PAGE = joinpath(ROOT, "docs", "src", "api.md")
const PUBLIC_API_DOC_EXEMPTIONS = Set([:JuFitter])

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

function _api_page_text()
    isfile(API_PAGE) || error("API reference page missing: $(relpath(API_PAGE, ROOT))")
    return read(API_PAGE, String)
end

@testset "Public API reference docstrings" begin
    exports = _public_exports()
    @test !isempty(exports)

    missing = Symbol[name for name in exports if !_has_public_docstring(name)]
    @test missing == Symbol[]

    api_text = _api_page_text()
    undocumented_on_page = Symbol[name for name in exports if !occursin(string(name), api_text)]
    @test undocumented_on_page == Symbol[]
end
