using Documenter
using JuFitter

makedocs(;
    modules=[JuFitter],
    sitename="JuFitter",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aminelsayed.github.io/JuFitter.jl",
        edit_link=nothing,
        repolink=nothing,
    ),
    pages=[
        "Home" => "index.md",
        "Roadmap" => "roadmap.md",
        "Research Landscape" => "research_landscape.md",
        "Plotting Design" => "plotting_design.md",
        "Overview" => "overview.md",
        "Statistical Foundations" => "statistical_foundations.md",
        "Backend Design" => "backend_design.md",
        "Mathematical Audit" => "mathematical_audit.md",
        "API Reference" => "api.md",
    ],
    checkdocs=:none,
    remotes=nothing,
)
