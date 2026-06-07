using Test
using JuFitter

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

@testset "Public API reference docstrings" begin
    exports = _public_exports()
    @test !isempty(exports)

    missing = Symbol[name for name in exports if !_has_public_docstring(name)]
    @test missing == Symbol[]
end
