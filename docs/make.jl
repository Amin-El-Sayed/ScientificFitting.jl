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
        assets=["assets/jufitter.css", "assets/jufitter.js"],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => [
            "Install" => "install.md",
            "Quickstart" => "quickstart.md",
            "Gallery" => [
                "Overview" => "gallery.md",
                "Linear Calibration" => "gallery/linear_calibration.md",
                "Photoelectric Work Function" => "gallery/photoelectric_threshold.md",
                "Damped Oscillator" => "gallery/resonance_decay.md",
                "Full Covariance" => "gallery/full_covariance.md",
                "XY Uncertainties" => "gallery/xy_uncertainties.md",
                "Constraints and Profiles" => "gallery/constraints_profiles.md",
                "Poisson and Histograms" => "gallery/poisson_histogram.md",
                "Multi-Dataset Fit" => "gallery/multi_dataset.md",
            ],
        ],
        "Concepts" => [
            "Fitting for Practitioners" => "fitting_for_practitioners.md",
            "Statistical Foundations" => "statistical_foundations.md",
            "Plotting Design" => "plotting_design.md",
        ],
        "Reference" => [
            "API Reference" => "api.md",
            "Overview" => "overview.md",
        ],
        "Development Notes" => [
            "Backend Design" => "backend_design.md",
            "Performance" => "performance.md",
            "Roadmap" => "roadmap.md",
            "Documentation Plan" => "documentation_plan.md",
            "Maintenance Notes" => "maintenance.md",
            "Research Landscape" => "research_landscape.md",
        ],
    ],
    checkdocs=:none,
    remotes=nothing,
)
