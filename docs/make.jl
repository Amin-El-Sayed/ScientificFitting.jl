using Pkg

if get(ENV, "SCIENTIFICFITTING_DOCS_SKIP_DEVELOP", "0") != "1"
    Pkg.instantiate()
end

using Documenter
using ScientificFitting
using CairoMakie

makedocs(;
    modules=[ScientificFitting],
    sitename="ScientificFitting",
    format=Documenter.HTML(;
        # Keep raw gallery-card links identical for local and hosted static builds.
        prettyurls=false,
        edit_link=nothing,
        repolink=nothing,
        assets=[
            "assets/favicon.ico",
            "assets/scientificfitting.css",
            "assets/scientificfitting.js",
        ],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => [
            "Install" => "install.md",
            "Quickstart" => "quickstart.md",
            "How ScientificFitting Works" => "how_scientificfitting_works.md",
        ],
        "Gallery" => [
            "Overview" => "gallery.md",
            "Linear Calibration" => "gallery/linear_calibration.md",
            "XY Uncertainties" => "gallery/xy_uncertainties.md",
            "Full Covariance" => "gallery/full_covariance.md",
            "Damped Oscillator" => "gallery/resonance_decay.md",
            "Photoelectric Work Function" => "gallery/photoelectric_threshold.md",
            "Poisson and Histograms" => "gallery/poisson_histogram.md",
            "Constraints and Profiles" => "gallery/constraints_profiles.md",
            "Multi-Dataset Fit" => "gallery/multi_dataset.md",
        ],
        "Guides" => [
            "Fitting for Practitioners" => "fitting_for_practitioners.md",
            "Plotting and Customization" => "plotting_design.md",
        ],
        "Mathematics and Statistics" => [
            "Overview" => "statistical_foundations.md",
            "Gaussian Fits and Covariance" => "gaussian_models.md",
            "Parameters and Fit Quality" => "parameter_inference.md",
            "Profiles and Contours" => "profiles_contours.md",
            "Likelihoods and Model Comparison" => "likelihood_models.md",
        ],
        "Reference" => [
            "API Overview" => "api.md",
            "Fitting" => "api_fitting.md",
            "Results and Diagnostics" => "api_results.md",
            "Fit Plotting" => "api_plotting.md",
            "Diagnostic Plotting" => "api_plotting_diagnostics.md",
            "Citation" => "citation.md",
        ],
        "Developer" => [
            "Architecture" => "backend_design.md",
            "Performance" => "performance.md",
        ],
    ],
    checkdocs=:none,
    remotes=nothing,
)
