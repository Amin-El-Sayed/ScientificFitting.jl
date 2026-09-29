using Pkg

if get(ENV, "SCIENTIFICFITTING_DOCS_SKIP_DEVELOP", "0") != "1"
    Pkg.instantiate()
end

# Load documented extensions before Documenter inspects bindings and docstrings.
import BuildConstructors, DistributionsHEP, NativeMinuit
using Documenter
using ScientificFitting
using CairoMakie

makedocs(;
    modules=[ScientificFitting],
    sitename="ScientificFitting",
    format=Documenter.HTML(;
        # Keep raw gallery-card links identical for local and hosted static builds.
        prettyurls=false,
        collapselevel=1,
        edit_link=nothing,
        repolink=nothing,
        assets=[
            "assets/favicon.ico",
            "assets/scientificfitting.css",
        ],
    ),
    pages=[
        "ScientificFitting" => "index.md",
        "Getting Started" => [
            "Installation" => "install.md",
            "First Fit" => "quickstart.md",
            "Assessing a Fit" => "fitting_for_practitioners.md",
            "Migrating" => "migration.md",
        ],
        "Examples" => [
            "Overview" => "gallery.md",
            "X and Y Uncertainties" => "gallery/xy_uncertainties.md",
            "Counts and Histograms" => "gallery/poisson_histogram.md",
            "Profiles and Constraints" => "gallery/constraints_profiles.md",
            "Full Covariance" => "gallery/full_covariance.md",
            "Damped Oscillator" => "gallery/resonance_decay.md",
            "Photoelectric Work Function" => "gallery/photoelectric_threshold.md",
            "Multi-Dataset Fit" => "gallery/multi_dataset.md",
            "LHCb Mass Spectrum" => "gallery/lhcb_mass_spectrum.md",
        ],
        "Statistics" => [
            "Statistics Reference" => "statistics.md",
            "Validation" => "validation.md",
            "Glossary" => "glossary.md",
        ],
        "Plotting" => "plotting_design.md",
        "Python" => "python.md",
        "API Reference" => [
            "Overview" => "api.md",
            "Fitting" => "api_fitting.md",
            "Results and Diagnostics" => "api_results.md",
            "Fit Plotting" => "api_plotting.md",
            "Diagnostic Plotting" => "api_plotting_diagnostics.md",
        ],
        "Internals" => [
            "How ScientificFitting Works" => "how_scientificfitting_works.md",
            "Packages and Interfaces" => "interfaces.md",
            "Backend Design" => "backend_design.md",
        ],
        "About" => [
            "Scope and Alternatives" => "limits.md",
            "Citation and License" => "citation.md",
        ],
    ],
    checkdocs=:none,
    remotes=nothing,
)
