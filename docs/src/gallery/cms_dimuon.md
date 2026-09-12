# CMS Dimuon Mass Spectrum

**Question:** can one inclusive mass model describe 2.55 million muon pairs,
and how stable is its inferred signal yield?

The minimizer converges. The residuals still reject the model. Here we keep
the observations fixed, investigate two plausible modeling changes, and
check the kinematic dependence behind the mass peak.

## Data And Selection

Source: **Stefan Wunsch (2019)**, *DoubleMuParked dataset from 2012 in NanoAOD
format reduced on muons*, [CERN Open Data, DOI 10.7483/OPENDATA.CMS.LVG5.QT81](https://doi.org/10.7483/OPENDATA.CMS.LVG5.QT81).
The file contains CMS pp collision data from Runs 2012B/C at 8 TeV. It retains
muon kinematics and charges from the already-triggered primary datasets;
it is not an unbiased sample of all proton collisions. The release uses
validated runs but states that the reduced output received no further validation.
The data are CC0-1.0; CMS and CERN do not endorse this analysis.

| Stage | Selection | Events remaining |
|:---|:---|---:|
| Public ROOT file | `Run2012BC_DoubleMuParked_Muons.root` | 61,540,413 |
| Muon multiplicity | Exactly two reconstructed muons | 31,104,343 |
| Charge | Opposite signs | 24,067,843 |
| This fit | ``2.8\leq m_{\mu\mu}<3.4\ \mathrm{GeV}/c^2`` | 2,551,454 |

The first two cuts follow [ROOT's dimuon tutorial](https://root.cern.ch/doc/master/df102__NanoAODDimuonAnalysis_8C.html).
We add the mass window around ``J/\psi\to\mu^+\mu^-``. The parent mass is
reconstructed from the two measured four-momenta, using ``c=1``:

```math
E_i=\sqrt{m_i^2+p_{T,i}^2\cosh^2\eta_i},\qquad
m_{\mu\mu}^2=(E_1+E_2)^2-|\mathbf p_1+\mathbf p_2|^2.
```

The [preparation script](cms_prepare.jl) verifies the 2.1 GiB source file's
checksum and retains every selected mass and its muon kinematics. The
[downloadable CSV](cms_mass_histogram.csv) contains exact counts in 300 equal
2 MeV bins, plus three disjoint groups defined by
``\max(|\eta_1|,|\eta_2|)``. No random subsampling is used.

```@setup cms
using ScientificFitting
data_dir = joinpath(dirname(pathof(ScientificFitting)), "..", "examples", "data", "cms_dimuon")
cp(joinpath(data_dir, "prepare.jl"), "cms_prepare.jl"; force=true)
cp(joinpath(data_dir, "models.jl"), "cms_models.jl"; force=true)
cp(joinpath(data_dir, "mass_histogram.csv"), "cms_mass_histogram.csv"; force=true)
```

```@example cms
using ScientificFitting, Distributions, DistributionsHEP, Printf, LinearAlgebra
import NativeMinuit  # both packages export `profile`; keep that name unambiguous
include("cms_models.jl")  # model definitions, parameter bounds and MIGRAD steps

data = read_cms_histogram("cms_mass_histogram.csv")
edges, counts = data.edges, data.counts
println("Candidates: ", sum(counts), "; mass bins: ", length(counts))
println("Disjoint group counts: ", vec(sum(data.groups; dims=1)))
```

**Scope:** a model-validation study of selected candidates. The reduced file
does not contain trigger decisions, track-quality flags, decay vertices or
efficiency corrections. It cannot support a prompt-production cross section
or reproduce the selections of a CMS physics measurement.

## Model And Fit

The narrow ``J/\psi`` peak is broadened by momentum resolution. Energy loss
and reconstruction effects can produce asymmetric tails. Unrelated muons and
other processes contribute a continuum underneath it; this is background,
not an extra Gaussian measurement error.

We use `DoubleSidedBifurcatedCrystalBallDas` from DistributionsHEP: an
asymmetric Gaussian core, a left power-law tail and a right exponential tail.
Its parameters describe an empirical reconstructed-mass shape, not a decay
amplitude. All components are normalized inside the same mass window:

```math
\lambda(m)=N_s S(m)+N_b B(m),\qquad
\nu_i=\int_{e_i}^{e_{i+1}}\lambda(m)\,dm,\qquad
n_i\sim\operatorname{Poisson}(\nu_i).
```

`fit_distribution` integrates upstream distributions across each bin and
minimizes ``-2\log L``. Yields are parameterized in millions of candidates
for numerical scaling; bounds and optimizer step sizes are not priors.
The [complete model code](cms_models.jl) defines three comparisons:

| Model | Change | Reason to check it |
|:---|:---|:---|
| A | One core scale + exponential background | Baseline with smooth falling continuum |
| B | Keep the peak; use a positive quadratic background | The continuum need not have a constant logarithmic slope |
| C | Keep B's background; add a wider Gaussian at the same center | A mixture of resolutions need not have one core scale |

B uses ``z=(m-2.8)/0.6`` and a convex mixture of
``\operatorname{Beta}(1,3)``, ``\operatorname{Beta}(2,2)`` and
``\operatorname{Beta}(3,1)`` densities in ``z``: a positive, normalized
degree-two Bernstein polynomial. C uses
``S_C=f S_A+(1-f)\mathcal N_W(\mu,r\sigma)`` with ``r>1``; these are two
resolution components, not two particles. Related CMS
[scouting analyses](https://twiki.cern.ch/twiki/bin/view/CMSPublic/Run3DiMuonScouting)
use a Crystal Ball plus Gaussian and a Bernstein background. Their data and
polynomial degree differ; that precedent does not validate our choices.

### Check The Optimizer Before The Model

A single MIGRAD start can stop at a poor boundary solution. Its convergence
flag is not a global-minimum guarantee or a validation of the local covariance:

```@example cms
options = cms_baseline_options(sum(counts))
single = fit_distribution(cms_exponential, edges, counts;
    options..., tol=1e-4, maxiters=2000)
@printf("One start: converged=%s, D=%.2f, positive-definite covariance=%s\n",
    single.converged, single.stats.chi2, isposdef(Symmetric(single.param_covariance)))
```

For every model below, compare three explicit starts spanning 20–50 MeV core
scales, with the **same bounds, data and likelihood**. `initial_guesses` and
`multistart=3` select the converged candidate with the lowest cost, not the
model with the largest p-value. This checks starting-point sensitivity; it
does not prove global optimality.

```@example cms
results = fit_cms_models(edges, counts)
for (label, fit) in zip(("A", "B", "C"), results)
    i = findfirst(==("signal_million"), fit.problem.parameter_names)
    @printf("%s  converged=%s  D/ndf=%.2f/%d  asymptotic p=%.3g\n",
        label, fit.converged, fit.stats.chi2, fit.stats.ndf, fit.stats.pvalue)
    @printf("   signal = %.0f +/- %.0f candidates (local error)\n",
        1e6fit.params[i], 1e6fit.param_stderr[i])
end
```

```@setup cms
for (fit, reference) in zip(results, (791.29897, 656.202, 417.422))
    @assert fit.converged && all(isfinite, fit.params) && all(>(0), fit.param_stderr)
    @assert all(isfinite, fit.param_covariance)
    @assert isposdef(Symmetric(fit.param_covariance))
    @assert isapprox(fit.stats.chi2, reference; atol=.05)
    means = cms_binmeans(fit, edges)
    @assert all(>(0), means)
    @assert isapprox(sum(means), sum(DistributionsHEP.yields(fitted_model(fit))); rtol=1e-12)
    @assert isapprox(sum(abs2, cms_residuals(counts, means)), fit.stats.chi2; rtol=1e-8)
end
```

## Inspect The Residuals

The peak curves nearly coincide, but the signed deviance residuals reveal
structure. For each model,

```math
r_i=\operatorname{sign}(n_i-\nu_i)\sqrt{2[\nu_i-n_i+n_i\log(n_i/\nu_i)]},
\qquad D=\sum_i r_i^2.
```

At these large counts they are approximately standard normal under a correct
model, apart from the dependence introduced by fitting.

```@example cms
include("cms_mass_plot.jl")  # ordinary Makie composition; no refitting
fig = cms_mass_plot(results, edges, counts; theme=:sans, show_panel=true)
nothing # hide
```

```@setup cms
for style in (:sans, :tex), appearance in (:light, :dark), panel in (true, false)
    save("cms_mass_$(style)_$(appearance)_$(panel).svg",
        cms_mass_plot(results, edges, counts; theme=style, appearance, show_panel=panel))
end
```

```@raw html
<img class="scientificfitting-plot scientificfitting-plot-light" data-scientificfitting-plot-group="cms-mass" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="cms_mass_sans_light_true.svg" alt="CMS mass spectrum and three residual panels; all candidate counts unchanged; sans with panel">
<img class="scientificfitting-plot scientificfitting-plot-dark" data-scientificfitting-plot-group="cms-mass" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="cms_mass_sans_dark_true.svg" alt="CMS mass comparison, dark sans with panel">
<img class="scientificfitting-plot scientificfitting-plot-light" data-scientificfitting-plot-group="cms-mass" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="cms_mass_sans_light_false.svg" alt="CMS mass comparison, sans without panel">
<img class="scientificfitting-plot scientificfitting-plot-dark" data-scientificfitting-plot-group="cms-mass" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="cms_mass_sans_dark_false.svg" alt="CMS mass comparison, dark sans without panel">
<img class="scientificfitting-plot scientificfitting-plot-light" data-scientificfitting-plot-group="cms-mass" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="cms_mass_tex_light_true.svg" alt="CMS mass comparison, TeX with panel">
<img class="scientificfitting-plot scientificfitting-plot-dark" data-scientificfitting-plot-group="cms-mass" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="cms_mass_tex_dark_true.svg" alt="CMS mass comparison, dark TeX with panel">
<img class="scientificfitting-plot scientificfitting-plot-light" data-scientificfitting-plot-group="cms-mass" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="cms_mass_tex_light_false.svg" alt="CMS mass comparison, TeX without panel">
<img class="scientificfitting-plot scientificfitting-plot-dark" data-scientificfitting-plot-group="cms-mass" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="cms_mass_tex_dark_false.svg" alt="CMS mass comparison, dark TeX without panel">
```

Bars on the counts use the large-count ``\sqrt n`` approximation. The shaded
``\pm2`` residual range is a visual reference, not a simultaneous confidence
band. Even C gives ``D\simeq417`` for 288 degrees of freedom: **the remaining
discrepancy is not explained by Poisson fluctuations under that model**.

## Check The Kinematic Mixture

Pseudorapidity ``\eta=-\log\tan(\theta/2)`` describes the muon's angle to
the beam. Different angles traverse different detector regions; different
momenta also change the curvature measurement. Mixing such candidates can
mix their mass resolutions. In the
[CMS 2011 quarkonium measurement](https://arxiv.org/abs/1502.04155), the fits
use narrow rapidity/transverse-momentum bins where resolution varies little.
That is a different dataset at 7 TeV, not a result we reproduce here.

Refit A in three disjoint muon-angle groups. They preserve every selected
candidate and use the same 300 bins. This is an exploratory check, not a
detector calibration or a selection optimized to improve fit quality.

```@example cms
group_results = fit_cms_groups(edges, data.groups, first(results))
for (label, fit) in zip(("[0,0.9)", "[0.9,1.4)", "[1.4,Inf)"), group_results)
    @printf("max |eta| %-10s  core scale = %.2f +/- %.2f MeV/c^2  D/ndf=%.2f/%d\n",
        label, 1000fit.params[2], 1000fit.param_stderr[2], fit.stats.chi2, fit.stats.ndf)
end
```

```@setup cms
for (fit, scale, deviance) in zip(group_results, (.0246936, .0365539, .0464352), (740.,436.85,380.05))
    @assert fit.converged && isapprox(fit.params[2], scale; atol=1e-5)
    @assert isapprox(fit.stats.chi2, deviance; atol=.1)
end
for style in (:sans, :tex), appearance in (:light, :dark), panel in (true, false)
    save("cms_kinematics_$(style)_$(appearance)_$(panel).svg",
        cms_kinematic_plot(group_results, edges, data.groups; theme=style, appearance, show_panel=panel))
end
```

```@raw html
<img class="scientificfitting-plot scientificfitting-plot-light" data-scientificfitting-plot-group="cms-kinematics" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="cms_kinematics_sans_light_true.svg" alt="Normalized mass histograms by muon angle and fitted core scales; sans with panel">
<img class="scientificfitting-plot scientificfitting-plot-dark" data-scientificfitting-plot-group="cms-kinematics" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="cms_kinematics_sans_dark_true.svg" alt="Muon angle comparison, dark sans with panel">
<img class="scientificfitting-plot scientificfitting-plot-light" data-scientificfitting-plot-group="cms-kinematics" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="cms_kinematics_sans_light_false.svg" alt="Muon angle comparison, sans without panel">
<img class="scientificfitting-plot scientificfitting-plot-dark" data-scientificfitting-plot-group="cms-kinematics" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="cms_kinematics_sans_dark_false.svg" alt="Muon angle comparison, dark sans without panel">
<img class="scientificfitting-plot scientificfitting-plot-light" data-scientificfitting-plot-group="cms-kinematics" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="cms_kinematics_tex_light_true.svg" alt="Muon angle comparison, TeX with panel">
<img class="scientificfitting-plot scientificfitting-plot-dark" data-scientificfitting-plot-group="cms-kinematics" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="cms_kinematics_tex_dark_true.svg" alt="Muon angle comparison, dark TeX with panel">
<img class="scientificfitting-plot scientificfitting-plot-light" data-scientificfitting-plot-group="cms-kinematics" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="cms_kinematics_tex_light_false.svg" alt="Muon angle comparison, TeX without panel">
<img class="scientificfitting-plot scientificfitting-plot-dark" data-scientificfitting-plot-group="cms-kinematics" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="cms_kinematics_tex_dark_false.svg" alt="Muon angle comparison, dark TeX without panel">
```

The effective core scale changes from roughly 25 to 46 MeV. Its plotted
one-sigma error is a **local model-conditional error**, not an uncertainty on
a validated detector resolution. The groups themselves still have poor fit
quality. C's components therefore cannot be identified with these groups.

## Interpretation

- Changing only the background moves the signal yield by about **20,000
  candidates**, versus local errors of about 1,800–2,400. This is sensitivity
  to the assumed model, not a calibrated systematic uncertainty from three fits.
- C improves the residuals but does not resolve them. A different minimizer
  or a narrower local error bar would not fix that model discrepancy.
- The comparisons use the same data, window and likelihood. A and B are not
  nested; introducing C's extra component has boundary/identifiability issues.
  Do not interpret their deviance differences with a naive one-parameter
  likelihood-ratio test.

These models were investigated on this same spectrum, not chosen in an
independent validation sample. The comparisons are exploratory model checks,
not a calibrated model-selection or discovery test.

## What To Do Before Reporting A Yield

Check momentum and angular dependence, source selection and reconstruction
quality before expanding the mass model. The reduced file contains 6,893
selected pairs with at least one muon beyond ``|\eta|=2.4``; they remain in
this example. Missing quality and trigger fields prevent us from deciding
which reconstruction or selection effects explain the residual structure.
A measurement would need richer event information, detector validation and
an uncertainty study, not successively higher polynomial orders until a
p-value crosses a threshold.

### Run The Individual-Event Fit

The histogram check and the large-event throughput check are different tasks.
To use all **2,551,454 individual masses**, download the original file from
the source record, install UnROOT, and run:

```julia
using UnROOT
include("cms_prepare.jl")
events = prepare_cms("Run2012BC_DoubleMuParked_Muons.root")
write_cms_histogram("rebuilt_histogram.csv", events)  # reproduces the CSV above

# Keep the fitted starting values, model, bounds and solver; only change the data form.
options = cms_baseline_options(length(events.mass))
unbinned = fit_distribution(cms_exponential, events.mass;
    options..., p0=first(results).params, tol=1e-4, maxiters=1000)
println(report_text(unbinned))
```

The full unbinned baseline fit was executed: about **55 s**, including remaining
compilation, with approximately **3.1 GiB process peak resident memory** on the
development Mac. ROOT preparation is separate. These are environment-dependent
observations, not a cross-package benchmark. Event likelihoods retain every
mass; the binned residual checks above remain necessary for assessing the shape.
Convergence alone does not validate either analysis.
The [integration and event-storage repairs](../backend_design.md#v03-ecosystem)
are separate numerical checks; they do not change the selected events or
remove the residual model discrepancies.
