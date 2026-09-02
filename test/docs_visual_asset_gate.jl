using Test

const ROOT = abspath(joinpath(@__DIR__, ".."))
const DOCS_SRC = joinpath(ROOT, "docs", "src")
const GALLERY_ASSETS = joinpath(DOCS_SRC, "assets", "gallery")

const REQUIRED_STYLE_PAIRS = Set(
    (style, appearance)
    for style in ("sans", "tex")
    for appearance in ("light", "dark")
)

const REQUIRED_PANEL_VARIANTS = Set(
    (style, panel, appearance)
    for style in ("sans", "tex")
    for panel in ("show", "hide")
    for appearance in ("light", "dark")
)

struct PngInfo
    width::Int
    height::Int
    bytes::Int
end

function docs_source_files()
    files = String[]
    for (dir, _, names) in walkdir(DOCS_SRC)
        for name in names
            if endswith(name, ".md") || endswith(name, ".css") || endswith(name, ".js")
                push!(files, joinpath(dir, name))
            end
        end
    end
    return sort(files)
end

function png_info(path::AbstractString)
    bytes = read(path)
    length(bytes) >= 24 || throw(ArgumentError("PNG is too small: $path"))
    bytes[1:8] == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a] ||
        throw(ArgumentError("invalid PNG signature: $path"))
    String(bytes[13:16]) == "IHDR" || throw(ArgumentError("missing IHDR chunk: $path"))

    width = parse_big_endian_uint32(bytes[17:20])
    height = parse_big_endian_uint32(bytes[21:24])
    return PngInfo(width, height, length(bytes))
end

function parse_big_endian_uint32(bytes)
    value = UInt32(0)
    for byte in bytes
        value = (value << 8) | UInt32(byte)
    end
    return Int(value)
end

function image_records()
    records = NamedTuple[]
    pattern = r"<img([^>]+)>"
    for page in docs_source_files()
        text = read(page, String)
        for img_match in eachmatch(pattern, text)
            tag = img_match.captures[1]
            src_match = match(r"src=\"([^\"]+\.png)\"", tag)
            isnothing(src_match) && continue

            group_match = match(r"data-scientificfitting-plot-group=\"([^\"]+)\"", tag)
            style_match = match(r"data-scientificfitting-plot-style=\"([^\"]+)\"", tag)
            panel_match = match(r"data-scientificfitting-plot-panel=\"([^\"]+)\"", tag)
            alt_match = match(r"alt=\"([^\"]+)\"", tag)
            src = src_match.captures[1]
            appearance = occursin("_dark", src) ? "dark" : occursin("_light", src) ? "light" : ""
            path = normpath(joinpath(dirname(page), src))
            push!(
                records,
                (
                    page=page,
                    src=src,
                    path=path,
                    group=isnothing(group_match) ? "" : group_match.captures[1],
                    style=isnothing(style_match) ? "" : style_match.captures[1],
                    panel=isnothing(panel_match) ? "" : panel_match.captures[1],
                    appearance=appearance,
                    alt=isnothing(alt_match) ? "" : alt_match.captures[1],
                ),
            )
        end
    end
    return records
end

function referenced_asset_names(records)
    return Set(basename(record.path) for record in records)
end

@testset "Documentation visual assets" begin
    records = image_records()
    @test !isempty(records)

    infos = Dict{String,PngInfo}()
    for record in records
        @test isfile(record.path)
        @test !isempty(strip(record.alt))
        info = png_info(record.path)
        infos[record.path] = info
        @test info.width >= 900
        @test info.height >= 500
        @test info.bytes >= 40_000
        # Stacked diagnostics legitimately use a slightly portrait canvas;
        # reject only shapes that become impractical in the documentation.
        @test 0.9 <= info.width / info.height <= 2.4
    end

    grouped = Dict{String,Vector{NamedTuple}}()
    for record in records
        isempty(record.group) && continue
        push!(get!(grouped, record.group, NamedTuple[]), record)
    end
    @test !isempty(grouped)

    for (group, group_records) in grouped
        observed = Set((record.style, record.appearance) for record in group_records)
        @test REQUIRED_STYLE_PAIRS ⊆ observed

        panels = Set(record.panel for record in group_records if !isempty(record.panel))
        if !isempty(panels)
            @test panels == Set(("show", "hide"))
            variants = Set(
                (record.style, record.panel, record.appearance) for record in group_records
            )
            @test REQUIRED_PANEL_VARIANTS ⊆ variants
        end

        # Light/dark assets must be interchangeable in-place. Panel state may
        # change the canvas because it changes content, while style may not.
        layout_keys = Set((record.style, record.panel) for record in group_records)
        for layout_key in layout_keys
            style_records = filter(
                record -> (record.style, record.panel) == layout_key,
                group_records,
            )
            dims = Set(
                (infos[record.path].width, infos[record.path].height)
                for record in style_records
            )
            @test length(dims) == 1
        end

        for record in group_records
            @test record.style in ("sans", "tex")
            @test record.panel in ("", "show", "hide")
            @test record.appearance in ("light", "dark")
        end
    end

    referenced = referenced_asset_names(records)
    for path in sort(readdir(GALLERY_ASSETS; join=true))
        endswith(path, ".png") || continue
        @test basename(path) in referenced
    end
end
