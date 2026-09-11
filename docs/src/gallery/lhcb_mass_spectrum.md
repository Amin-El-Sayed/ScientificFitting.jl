# A Peak In LHCb Open Data

**Question:** how many reconstructed three-kaon candidates belong to the
charged-B peak rather than the smooth background? Fit real collision counts
with Distributions.jl, BuildConstructors.jl and NativeMinuit.jl, then check how
the answer depends on the peak shape.

## Data And Selection

Source: **LHCb collaboration (2017)**, *Matter Antimatter Differences
(B meson decays to three hadrons) - Data Files*, CERN Open Data,
[DOI: 10.7483/OPENDATA.LHCB.AOF7.JH09](https://doi.org/10.7483/OPENDATA.LHCB.AOF7.JH09).
These are preselected 2011 pp collision candidates at 7 TeV, not simulation.
The data are CC0-1.0; CERN and LHCb do not endorse this analysis.

| Choice | Definition |
|:---|:---|
| File | `B2HHH_MagnetUp.root`, `DecayTree`; 3,420,295 candidates |
| Track selection | Every track: `ProbK > 0.5`, `ProbPi < 0.5`, `isMuon == 0`; 9,717 candidates remain |
| Charge | Both signs combined; only the MagnetUp sample |
| Mass | Sum three four-momenta with the kaon hypothesis, ``m_K=493.677\,\mathrm{MeV}/c^2`` |
| Fit window | ``5200\leq m_{KKK}<5600\,\mathrm{MeV}/c^2``; 7,368 candidates in 80 bins |

The cuts and momentum units follow the
[LHCb project notebook](https://github.com/lhcb/opendata-project/blob/master/LHCb_Open_Data_Project.ipynb).
The assigned kaon mass is from the
[Particle Data Group](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-charged-kaon-mass.pdf).
[Download the UnROOT preparation script](lhcb_prepare.jl) to rebuild the
histogram from the [original file](https://opendata.cern.ch/record/4900/files/B2HHH_MagnetUp.root).
It verifies the file checksum and accounts for candidates outside the window.
The full ROOT file is not downloaded during a documentation build.

```@setup lhcb
using ScientificFitting
cp(joinpath(dirname(pathof(ScientificFitting)), "..", "examples", "data", "lhcb_mass", "prepare.jl"), "lhcb_prepare.jl"; force=true)
cp(joinpath(dirname(pathof(ScientificFitting)), "..", "benchmarks", "lhcb_reference.py"), "lhcb_reference.py"; force=true)
```

```@example lhcb
using ScientificFitting, Distributions, DistributionsHEP, BuildConstructors, Printf
import NativeMinuit

edges = collect(5200.:5.:5600.)  # MeV/c^2; left-closed, right-open bins
counts = [
    12,16,19,19,24,26,40,64,67,92,152,228,333,419,571,683,
    738,747,711,570,434,291,219,147,88,57,35,36,25,19,16,15,
    13,12,12,3,16,17,13,10,11,11,9,7,10,13,9,12,
    11,6,10,6,8,8,12,13,10,8,8,12,9,7,14,13,
    11,9,7,9,12,7,9,4,15,3,10,9,7,4,3,3,
]
println("Candidates in fit window: ", sum(counts))
@assert length(edges) == length(counts)+1 && sum(counts) == 7368 # hide
```

## Model And Fit

Different track momenta and geometries need not have the same mass resolution.
Two Gaussians with a shared center give a compact empirical core-and-tail
model. They are **two resolution components of one peak**, not two particles.
An exponential represents the smooth background in this window.

```math
S(m)\propto f\,\mathcal N(m;\mu,\sigma)
 +(1-f)\,\mathcal N(m;\mu,r\sigma),\qquad B(m)\propto e^{-(m-5200)/\tau}.
```

Both shapes are normalized **inside the fit window**. The free yields ``N_s``
and ``N_b`` therefore count selected candidates there, not efficiency-corrected
decays. Independent Poisson bins have expectations
``\nu_i=N_s\int_i S(m)\,dm+N_b\int_i B(m)\,dm``; the adapter integrates the
distributions rather than sampling densities at bin centers.

```@example lhcb
@with_parameters(MassSpectrum;
    center::P, width::P, width_ratio::P, core_fraction::P,
    background_scale::P, signal_yield::P, background_yield::P, begin
    peak = MixtureModel(
        [Normal(center, width), Normal(center, width*width_ratio)],
        [core_fraction, 1-core_fraction])
    signal = truncated(peak, 5200., 5600.)
    # Exponential already has lower support zero; only truncate its upper end.
    background = 5200. + truncated(Exponential(background_scale); upper=400.)
    ExtendedMixtureModel([signal, background], [signal_yield, background_yield])
end)

constructor = ConstructorOfMassSpectrum(
    AdvancedParameter("center", 5284.; boundaries=(5250.,5310.)),
    AdvancedParameter("width", 14.; boundaries=(3.,30.)),
    AdvancedParameter("width_ratio", 2.; boundaries=(1.05,5.)),
    AdvancedParameter("core_fraction", 0.8; boundaries=(0.05,1.)),
    AdvancedParameter("background_scale", 400.; boundaries=(20.,5000.)),
    AdvancedParameter("signal_yield", 6500.; boundaries=(0.,15000.)),
    AdvancedParameter("background_yield", 1000.; boundaries=(0.,15000.)),
)
result = fit_distribution(constructor, edges, counts;
    solver=NativeMinuitSolver(), tol=1e-5)

for (name, value, error) in zip(keys(parameter_values(result)), result.params, result.param_stderr)
    @printf("%-18s = %.6g +/- %.3g\n", name, value, error)
end
@assert result.converged # hide
@assert isapprox(result.stats.chi2, 69.43671155; atol=1e-5) # hide
@assert maximum(abs.((result.params - [5284.7386633,14.7467591,1.67331788,0.63443312,531.8874376,6484.3962325,883.6045225]) ./ result.param_stderr)) < 0.002 # hide
@assert isapprox(result.param_stderr, [0.2432334,1.003249,0.1050925,0.139144,140.9726,93.58187,56.18607]; rtol=1e-3) # hide
```

The bounds separate the two widths and keep scales and yields physical; they
are not uncertainty estimates or priors. The errors above are local Hessian
errors, conditional on this model. `fitted_model(result)` returns the upstream
`ExtendedMixtureModel`; constructing a plot does not run the optimizer again.

## Inspect The Spectrum

```@example lhcb
include("lhcb_mass_plot.jl")  # downloadable Makie composition, no fitting code
figure = lhcb_mass_plot(result, constructor, edges, counts;
    theme=:sans, show_panel=true)
nothing # hide
```

```@setup lhcb
for style in (:sans, :tex), appearance in (:light, :dark), panel in (true, false)
    fig = lhcb_mass_plot(result, constructor, edges, counts; theme=style, appearance, show_panel=panel)
    save("lhcb_mass_$(style)_$(appearance)_$(panel).svg", fig)
end
```

```@raw html
<img class="scientificfitting-plot scientificfitting-plot-light" data-scientificfitting-plot-group="lhcb-mass" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="lhcb_mass_sans_light_true.svg" alt="LHCb mass spectrum, fitted signal and background, local mean band and deviance residuals, sans with panel">
<img class="scientificfitting-plot scientificfitting-plot-dark" data-scientificfitting-plot-group="lhcb-mass" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="lhcb_mass_sans_dark_true.svg" alt="LHCb mass spectrum and residuals, dark sans with panel">
<img class="scientificfitting-plot scientificfitting-plot-light" data-scientificfitting-plot-group="lhcb-mass" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="lhcb_mass_sans_light_false.svg" alt="LHCb mass spectrum and residuals, sans without panel">
<img class="scientificfitting-plot scientificfitting-plot-dark" data-scientificfitting-plot-group="lhcb-mass" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="lhcb_mass_sans_dark_false.svg" alt="LHCb mass spectrum and residuals, dark sans without panel">
<img class="scientificfitting-plot scientificfitting-plot-light" data-scientificfitting-plot-group="lhcb-mass" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="lhcb_mass_tex_light_true.svg" alt="LHCb mass spectrum and residuals, TeX with panel">
<img class="scientificfitting-plot scientificfitting-plot-dark" data-scientificfitting-plot-group="lhcb-mass" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="lhcb_mass_tex_dark_true.svg" alt="LHCb mass spectrum and residuals, dark TeX with panel">
<img class="scientificfitting-plot scientificfitting-plot-light" data-scientificfitting-plot-group="lhcb-mass" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="lhcb_mass_tex_light_false.svg" alt="LHCb mass spectrum and residuals, TeX without panel">
<img class="scientificfitting-plot scientificfitting-plot-dark" data-scientificfitting-plot-group="lhcb-mass" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="lhcb_mass_tex_dark_false.svg" alt="LHCb mass spectrum and residuals, dark TeX without panel">
```

The [plotting function](lhcb_mass_plot.jl) uses ordinary Makie axes. Bars are
68.3% Garwood intervals for individual Poisson means; their discreteness makes
coverage conservative. The local one-sigma band propagates covariance to bin means;
it is **not** a prediction interval for new counts. The lower axis shows
signed Poisson deviance residuals.

## Check Shape Dependence

Set the core fraction to one and fix the now-unused width ratio: the same
constructor becomes a single-Gaussian model. No likelihood rewrite is needed.

```@example lhcb
single = deepcopy(constructor)
update!(single, (core_fraction=1.,))
fix!(single, (:core_fraction, :width_ratio))
single_result = fit_distribution(single, edges, counts;
    solver=NativeMinuitSolver(), tol=1e-5)
for (label, fit) in (("one width", single_result), ("two widths", result))
    @printf("%-11s  Ns = %.1f +/- %.1f   D/ndf = %.2f/%d   p = %.4f\n",
        label, fit.params[6], fit.param_stderr[6], fit.stats.chi2, fit.stats.ndf, fit.stats.pvalue)
end
@printf("AIC(one) - AIC(two) = %.2f\n", single_result.stats.aic-result.stats.aic)
@assert single_result.converged # hide
@assert isapprox(single_result.stats.chi2, 101.9637663; atol=1e-5) # hide
```

The single width leaves excess structure and shifts the yield by about 148
candidates. This is a **model-sensitivity check**, not an independent systematic
error to add in quadrature. The quoted p-values use an asymptotic deviance
reference, approximate for sparse bins. A two-degree-of-freedom likelihood-ratio
significance would be unjustified: when one mixture weight vanishes, its width
is not identified. AIC can be compared here because data and likelihood
normalization are unchanged.

The width and mixture fraction are correlated. Profile the signal yield while
refitting the shape and background, rather than interpreting the local error alone:

```@example lhcb
scan = profile(result, 6; nsigma=2.5, npoints=31, on_failure=:throw)
interval = profile_interval(scan)
@printf("Ns profile interval at delta(-2 log L)=1: [%.1f, %.1f]\n",
    interval.lower, interval.upper)
@assert isfinite(interval.lower) && isfinite(interval.upper) # hide
@assert isapprox(interval.lower, 6391.4184; atol=1.) && isapprox(interval.upper, 6578.6286; atol=1.) # hide
```

This is an approximate one-parameter 68.3% interval under regular likelihood
conditions, conditional on the two-width model. It does not cover shape
misspecification or detector calibration. An independent SciPy-CDF/C++-Minuit2
implementation checks both minima and the signal-yield MINOS interval; the
documentation build checks those reference costs, parameters and interval
crossings, and regenerates every output above.
The [independent check](lhcb_reference.py) is kept in `benchmarks/lhcb_reference.py`
and runs from the repository checkout with NumPy, SciPy and iminuit.

**Scope:** this is a reproducible open-data fit, not a reproduction of an LHCb
publication. The fitted reconstructed centroid is not a precision measurement
of the B mass. No momentum-scale calibration, efficiency correction, production
asymmetry or detection asymmetry is estimated here; neither a branching fraction
nor a CP asymmetry follows from these yields.
