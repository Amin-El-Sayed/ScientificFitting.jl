using Pkg

const DOCS_ROOT = @__DIR__
const PROJECT_ROOT = dirname(DOCS_ROOT)

if get(ENV, "JUFITTER_DOCS_SKIP_DEVELOP", "0") != "1"
    Pkg.develop(PackageSpec(path=PROJECT_ROOT))
    Pkg.instantiate()
end

using Documenter
using JuFitter
using CairoMakie

makedocs(;
    modules=[JuFitter],
    sitename="JuFitter",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link=nothing,
        repolink=nothing,
        assets=["assets/jufitter.css", "assets/jufitter.js"],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => [
            "Install" => "install.md",
            "Quickstart" => "quickstart.md",
            "How JuFitter Works" => "how_jufitter_works.md",
        ],
        "Gallery" => [
            "Overview" => "gallery.md",
            "Linear Calibration" => "gallery/linear_calibration.md",
            "XY Uncertainties" => "gallery/xy_uncertainties.md",
            "Full Covariance" => "gallery/full_covariance.md",
            "Damped Oscillator" => "gallery/resonance_decay.md",
            "Photoelectric Work Function" => "gallery/photoelectric_threshold.md",
            "Constraints and Profiles" => "gallery/constraints_profiles.md",
            "Poisson and Histograms" => "gallery/poisson_histogram.md",
            "Multi-Dataset Fit" => "gallery/multi_dataset.md",
        ],
        "Guides" => [
            "Fitting for Practitioners" => "fitting_for_practitioners.md",
            "Plotting Design" => "plotting_design.md",
        ],
        "Mathematics and Statistics" => [
            "Statistical Foundations" => "statistical_foundations.md",
        ],
        "Reference" => [
            "API Reference" => "api.md",
            "Reference Map" => "overview.md",
            "Technical Notes" => [
                "Backend Design" => "backend_design.md",
                "Performance" => "performance.md",
                "Maintenance Notes" => "maintenance.md",
            ],
        ],
    ],
    checkdocs=:none,
    remotes=nothing,
)
