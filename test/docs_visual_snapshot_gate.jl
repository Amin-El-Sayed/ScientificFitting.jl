using Test
using SHA

const ROOT = abspath(joinpath(@__DIR__, ".."))
const GALLERY_ASSETS = joinpath(ROOT, "docs", "src", "assets", "gallery")
const SNAPSHOT_MANIFEST = joinpath(@__DIR__, "docs_visual_snapshot_manifest.txt")

function _load_visual_manifest()
    entries = Dict{String,String}()
    isfile(SNAPSHOT_MANIFEST) || error("visual snapshot manifest missing: $SNAPSHOT_MANIFEST")

    for (lineno, raw_line) in enumerate(eachline(SNAPSHOT_MANIFEST))
        line = strip(split(raw_line, "#"; limit=2)[1])
        isempty(line) && continue
        parts = split(line, "="; limit=2)
        length(parts) == 2 || error("invalid manifest line $lineno: $raw_line")
        name = strip(parts[1])
        digest = lowercase(strip(parts[2]))
        occursin(r"^[0-9a-f]{64}$", digest) || error("invalid SHA-256 digest on line $lineno")
        haskey(entries, name) && error("duplicate visual snapshot entry: $name")
        entries[name] = digest
    end

    return entries
end

function _gallery_png_names()
    return sort!(String[name for name in readdir(GALLERY_ASSETS) if endswith(name, ".png")])
end

function _sha256_hex(path::AbstractString)
    return bytes2hex(sha256(read(path)))
end

@testset "Documentation visual snapshots" begin
    manifest = _load_visual_manifest()
    png_names = _gallery_png_names()

    @test !isempty(png_names)
    @test sort(collect(keys(manifest))) == png_names

    for name in png_names
        path = joinpath(GALLERY_ASSETS, name)
        @test _sha256_hex(path) == manifest[name]
    end
end
