# Fit A Peak And Background: LHCb Data

**Question:** how many entries belong to a peak above a smooth background?
Use a Poisson likelihood to estimate its area, compare two peak shapes, and
check the signal uncertainty with a profile. DistributionsHEP supplies the
mixture, BuildConstructors names its parameters, and NativeMinuit minimizes
the same likelihood used for the subsequent profile.

## Data And Selection

Source: **LHCb collaboration (2017)**, *Matter Antimatter Differences
(B meson decays to three hadrons) - Data Files*, CERN Open Data,
[DOI: 10.7483/OPENDATA.LHCB.AOF7.JH09](https://doi.org/10.7483/OPENDATA.LHCB.AOF7.JH09).
These are real 2011 proton-proton collision data at 7 TeV, released under
CC0-1.0. Each candidate is a combination of three charged tracks. Assigning
each track the kaon mass gives one reconstructed parent mass ``m``. Decays
``B^\pm\to K^\pm K^+K^-`` produce a peak; unrelated track combinations form
a broad background. Detector resolution gives the peak its finite width.

LHCb supplied candidates that already pass trigger, momentum and vertex
selections ([preselection notebook](https://github.com/lhcb/opendata-project/blob/master/Background-Information-Notebooks/DataSelection.ipynb)).
We then apply the particle-identification cuts below to the **whole MagnetUp
file**, without random subsampling. MagnetDown is not used.

| Choice | Definition |
|:---|:---|
| File | `B2HHH_MagnetUp.root`, `DecayTree`; 3,420,295 candidates |
| Track selection | Every track: `ProbK > 0.5`, `ProbPi < 0.5`, `isMuon == 0`; 9,717 candidates remain |
| Charge | B-plus and B-minus candidates combined |
| Mass | Invariant mass of the three tracks, each assigned ``m_K=493.677\,\mathrm{MeV}/c^2`` |
| Fit window | ``5200\leq m_{KKK}<5600\,\mathrm{MeV}/c^2``; 7,368 candidates in 80 bins |

The thresholds follow the starting cuts in the
[LHCb project notebook](https://github.com/lhcb/opendata-project/blob/master/LHCb_Open_Data_Project.ipynb).
The kaon mass is from the
[Particle Data Group](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-charged-kaon-mass.pdf).
[Download the UnROOT preparation script](lhcb_prepare.jl) to rebuild the
histogram from the [original file](https://opendata.cern.ch/record/4900/files/B2HHH_MagnetUp.root).
It checks the file checksum, reconstructs masses and counts every selection
step. The histogram below is sufficient to run the fit without downloading ROOT data.

**Scope:** estimate the signal count in this selected sample and mass window.
The fit has **80 Poisson bins and seven free parameters**; reading 3.4 million
source candidates is a separate data-preparation step.

```@setup lhcb
using ScientificFitting
cp(joinpath(dirname(pathof(ScientificFitting)), "..", "examples", "data", "lhcb_mass", "prepare.jl"), "lhcb_prepare.jl"; force=true)
cp(joinpath(dirname(pathof(ScientificFitting)), "..", "benchmarks", "lhcb_reference.py"), "lhcb_reference.py"; force=true)
```

```@example lhcb
using ScientificFitting, Distributions, DistributionsHEP, BuildConstructors, Printf
import NativeMinuit  # keep ScientificFitting.profile unambiguous

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

Use two Gaussians with a common center for the peak, and an exponential for
the background. A narrow core plus a wider component approximates a mixture
of detector resolutions; it represents **one peak, not two particles**.
The exponential describes a smoothly falling background with one slope.

For ``W=[m_{\mathrm{lo}},m_{\mathrm{hi}})=[5200,5600)`` in ``\mathrm{MeV}/c^2``, define

```math
\begin{aligned}
g(m) &= f\,\mathcal N(m;\mu,\sigma)+(1-f)\,\mathcal N(m;\mu,r\sigma),\\
S(m) &= \frac{g(m)}{\int_W g(u)\,du},\qquad
B(m)=\frac{e^{-(m-m_{\mathrm{lo}})/\tau}}
 {\tau\,[1-e^{-(m_{\mathrm{hi}}-m_{\mathrm{lo}})/\tau}]}.
\end{aligned}
```

Both densities vanish outside ``W`` and integrate to one inside it.
``\mu,\sigma,r,f`` set the peak shape, ``\tau`` the background slope, and
``N_s,N_b`` the signal and background counts in the window. For each bin,

```math
\nu_i=N_s\int_{\mathrm{bin}\ i}S(m)\,dm+N_b\int_{\mathrm{bin}\ i}B(m)\,dm,
\qquad n_i\sim\operatorname{Poisson}(\nu_i).
```

`fit_distribution` minimizes ``-2\sum_i\log\operatorname{Poisson}(n_i;\nu_i)``.
It integrates the densities over each bin, rather than evaluating their heights
at bin centers. `AdvancedParameter` declares each name, start and bounds;
`::P` inserts its current value when BuildConstructors builds the distribution.

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
@printf("Poisson deviance / ndf = %.2f / %d; asymptotic p = %.4f\n",
    result.stats.chi2, result.stats.ndf, result.stats.pvalue)
@assert result.converged # hide
@assert isapprox(result.stats.chi2, 69.43671155; atol=1e-5) # hide
@assert maximum(abs.((result.params - [5284.7386633,14.7467591,1.67331788,0.63443312,531.8874376,6484.3962325,883.6045225]) ./ result.param_stderr)) < 0.002 # hide
@assert isapprox(result.param_stderr, [0.2432334,1.003249,0.1050925,0.139144,140.9726,93.58187,56.18607]; rtol=1e-3) # hide
```

The fitted signal is about **6,484 candidates with a local standard error of 94**.
The deviance is **69.44 for 73 degrees of freedom**, with an approximate
``p=0.60``: this check does not detect an overall lack of fit. Inspect the
residuals next, then test the peak-shape assumption below.

Bounds keep widths and yields physical; they are not priors. The reported
errors come from local curvature. `fitted_model(result)` returns an
`ExtendedMixtureModel`, so plotting or evaluating it does not refit.

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

The lower axis shows signed Poisson deviance residuals: look for runs of bins
that the fit consistently over- or underestimates. Bars are 68.3% Garwood
intervals for individual bin means, with conservative coverage from discrete
counts. The local one-sigma band propagates parameter covariance to the fitted
bin means, not to future counts. The [plotting function](lhcb_mass_plot.jl)
uses ordinary Makie axes.

## Check Shape Dependence

Set the core fraction to one and fix the now-unused width ratio: the same
constructor becomes a single-Gaussian model. No likelihood rewrite is needed.

```@example lhcb
single = deepcopy(constructor)
BuildConstructors.update!(single, (core_fraction=1.,))
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

The two-width model describes these data better: the deviance falls from
101.96 to 69.44, and the fitted signal increases by about 148 candidates.
This supports allowing a wider resolution component; it does not identify
a second physical signal. AIC is comparable because the data and likelihood
normalization are unchanged.

The yield shift measures sensitivity to the peak model, not an independent
error to add in quadrature. The p-values are approximate, especially in sparse
bins. Do not assign a standard two-parameter likelihood-ratio significance to
the improvement: at zero mixture weight the unused width is not identifiable.

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

The profile interval is close to the local ``N_s\pm94`` estimate here. It has
an approximate one-parameter 68.3% interpretation under regular likelihood
conditions and includes refitting the other parameters. Neither interval
includes uncertainty from choosing the wrong peak or background shape.

The [independent numerical check](lhcb_reference.py), run from the repository
with NumPy, SciPy and iminuit, compares both fits and the signal-yield MINOS
interval. The executed cells check those reference values; the agreement tests
the numerical implementation of this model.

## What The Yield Measures

``N_s`` counts the fitted peak in this selection and window, before efficiency
corrections. The [LHCb publication](https://arxiv.org/abs/1306.1246) instead fits
charges separately, uses a different selection and more detailed signal and
background shapes, then corrects detector and production effects to measure
CP asymmetry. Its ``22\,119\pm164`` yield is therefore not a target for this fit.

We do not apply its charm veto, so ``B\to DK``, ``D\to KK`` decays can also
contribute to the peak. This example estimates a selected yield, not a charmless
branching fraction or CP asymmetry. A precision mass measurement would also
require momentum-scale calibration.
