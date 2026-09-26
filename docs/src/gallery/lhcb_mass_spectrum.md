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
CC0-1.0. Each candidate combines three charged tracks; assigning each track
the kaon mass gives one reconstructed parent mass ``m``. Decays
``B^\pm\to K^\pm K^+K^-`` produce a peak, unrelated track combinations a broad
background, and detector resolution gives the peak its finite width.

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
    solver=NativeMinuitSolver(), tol=1e-3)  # also used for subsequent profile refits

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

Bounds keep widths and yields physical; they are
[not priors](../statistics.md#Fixed-Parameters-And-Bounds). The reported errors
come from [local curvature](../statistics.md#Local-Parameter-Covariance).
`fitted_model(result)` returns an `ExtendedMixtureModel`, so plotting or
evaluating it does not refit.

**What SF adds here:** the model itself comes from DistributionsHEP and
BuildConstructors; MIGRAD comes from NativeMinuit.

| Step | Handled by SF |
|:---|:---|
| Model to likelihood | Integrate component distributions over bins, apply their yields, validate counts and evaluate the Poisson likelihood. |
| Constructor to solver | Transfer names, starts, bounds and fixed values; prepare derivatives and the correct objective scale. |
| Fit to inference | Return covariance, deviance, AIC and named values; retain the model and solver for profile refits and plots. |

[NativeMinuit](https://github.com/fkguo/NativeMinuit.jl) also provides binned
likelihoods, HESSE, MINOS and contours; a direct implementation can reach the
same result. SF supplies the adapters above and a common result and diagnostics
API across solvers. The custom spectrum drawing below remains ordinary Makie
code.

## Inspect The Spectrum

Evaluate each fitted component over the bin edges. The mean band uses
``\sigma_{\nu_i}^2 = J_i\operatorname{Cov}(\hat p)J_i^\mathsf{T}``, where
``J_{ij}=\partial\nu_i/\partial p_j``; these are post-fit calculations.

For the data points, ``n\pm\sqrt n`` is a poor guide at small counts: an empty
bin would even get a zero-width bar. Instead use **Garwood intervals** for
each Poisson mean ``\nu``. For confidence level ``1-\alpha=0.6827``, they are

```math
\nu_{\mathrm{lo}} =
\begin{cases}
0, & n=0,\\
\tfrac12\chi^2_{2n,\,\alpha/2}, & n>0,
\end{cases}
\qquad
\nu_{\mathrm{hi}}=\tfrac12\chi^2_{2(n+1),\,1-\alpha/2}.
```

Here ``\chi^2_{k,q}`` is the ``q``-quantile with ``k`` degrees of freedom.
The formula inverts the two Poisson tail tests; for ``n=0`` it gives
``[0,1.84]``, rather than ``[0,0]``. Coverage is at least the nominal level
because counts are discrete. These bars describe individual bins, not the
fitted signal-yield error, and **are not weights in the fit**.
See the [central Poisson interval formula](https://docs.astropy.org/en/stable/api/astropy.stats.poisson_conf_interval.html).

```@example lhcb
using ForwardDiff

names = keys(parameter_values(result))
build(p) = build_model(constructor, NamedTuple{names}(Tuple(p)))
component_counts(model) = [n .* diff(cdf.(d, edges))
    for (d, n) in zip(components(model), DistributionsHEP.yields(model))]
signal, background = component_counts(fitted_model(result))
means = signal + background
J = ForwardDiff.jacobian(p -> sum(component_counts(build(p))), result.params)
sigma = sqrt.(vec(sum((J * result.param_covariance) .* J; dims=2)))

# Signed Poisson deviance residuals; a zero count contributes no n*log(n/mu) term.
residual = [sign(n-m)*sqrt(max(2*(m-n+(n==0 ? 0 : n*log(n/m))), 0))
    for (n, m) in zip(counts, means)]

# Central 68.3% Garwood intervals for each bin mean, including empty bins.
tail = (1 - 0.682689492137)/2
lower = [n==0 ? 0. : quantile(Chisq(2n), tail)/2 for n in counts]
upper = [quantile(Chisq(2(n+1)), 1-tail)/2 for n in counts]
@assert isapprox(quantile(Chisq(2), 1-tail)/2, 1.841021645; atol=1e-8) # hide
centers = (edges[1:end-1] + edges[2:end])/2
step_x = repeat(edges; inner=2)[2:end-1]  # vertical steps at the actual bin edges
@printf("Expected count: %.2f; sum of squared deviance residuals: %.2f\n",
    sum(means), sum(abs2, residual))
@assert all(isfinite, sigma) && all(sigma .> 0) # hide
@assert isapprox(sum(means), sum(DistributionsHEP.yields(fitted_model(result))); rtol=1e-12) # hide
@assert isapprox(sum(abs2, residual), result.stats.chi2; atol=1e-8) # hide
```

Compose two ordinary Makie axes from these arrays; `theme`, `appearance` and
`show_panel` are independent options. Change the Makie calls directly to add
or restyle elements.

```@example lhcb
using CairoMakie, LaTeXStrings

"""Draw the fitted bin means and residuals computed above; do not refit."""
function spectrum_figure(; theme=:sans, appearance=:light, show_panel=true)
    pal = plot_palette(theme; appearance)
    return with_theme(plot_theme(theme; appearance)) do
        fig = Figure(size=(show_panel ? 1180 : 900, 740))
        ax = Axis(fig[1,1]; title="LHCb open data: three-kaon mass",
            ylabel="candidates / bin")
        band!(ax, step_x, repeat(means-sigma; inner=2), repeat(means+sigma; inner=2);
            color=(pal.band_color, 0.28), label="local 1-sigma mean band")
        lines!(ax, step_x, repeat(means; inner=2);
            color=pal.fit_color, label="two-width peak + background")
        lines!(ax, step_x, repeat(background; inner=2);
            color=pal.reference_color, linestyle=:dash, label="background")
        errorbars!(ax, centers, counts, counts-lower, upper-counts;
            color=pal.yerr_color, whiskerwidth=pal.error_whiskerwidth)
        scatter!(ax, centers, counts; color=pal.data_color,
            markersize=pal.data_markersize, label="data (68% Poisson intervals)")
        hidexdecorations!(ax; grid=false)

        rx = Axis(fig[2,1]; xlabel=theme==:tex ? L"m_{KKK}\,/(\mathrm{MeV}\,c^{-2})" :
            "three-kaon mass / (MeV/c^2)", ylabel="deviance residual")
        barplot!(rx, centers, residual; width=0.8 .* diff(edges), color=pal.fit_color)
        hlines!(rx, [0.]; color=pal.stats_color)
        hlines!(rx, [-2., 2.]; color=pal.reference_color, linestyle=:dash)
        linkxaxes!(ax, rx)
        xlims!(ax, first(edges), last(edges))
        ylims!(ax, 0, 1.08maximum(upper))
        rowsize!(fig.layout, 1, Auto(3))
        rowsize!(fig.layout, 2, Auto(1))

        if show_panel
            p, e = result.params, result.param_stderr
            labels = theme==:tex ? Any[
                LaTeXString(@sprintf("N_s = %.0f \\pm %.0f", p[6], e[6])),
                LaTeXString(@sprintf("N_b = %.0f \\pm %.0f", p[7], e[7])),
                LaTeXString(@sprintf("\\mu = %.2f \\pm %.2f\\;\\mathrm{MeV}/c^2", p[1], e[1])),
                LaTeXString(@sprintf("\\sigma = %.2f \\pm %.2f\\;\\mathrm{MeV}/c^2", p[2], e[2])),
            ] : Any[@sprintf("signal = %.0f +/- %.0f candidates", p[6], e[6]),
                @sprintf("background = %.0f +/- %.0f candidates", p[7], e[7]),
                @sprintf("centroid = %.2f +/- %.2f MeV/c^2", p[1], e[1]),
                @sprintf("core width = %.2f +/- %.2f MeV/c^2", p[2], e[2])]
            plot_info_panel!(fig[1:2,2]; theme, appearance, legend_source=ax,
                title="Window-conditional fit", parameter_lines=labels,
                statistic_lines=[@sprintf("D / ndf = %.2f / %d", result.stats.chi2, result.stats.ndf),
                    @sprintf("asymptotic p = %.3f", result.stats.pvalue)])
        else
            Legend(fig[3,1], ax; nbanks=2, tellwidth=false)
        end
        resize_plot_to_layout!(fig; minimum_axis_size=(600, nothing))
        fig
    end
end

figure = spectrum_figure(; theme=:sans, show_panel=true)
nothing # hide
```

```@setup lhcb
fig = spectrum_figure(; theme=:sans, appearance=:light, show_panel=true)
save("lhcb_mass_sans_light_true.svg", fig)
```

```@raw html
<img class="scientificfitting-plot" src="lhcb_mass_sans_light_true.svg" alt="LHCb mass spectrum, fitted signal and background, local mean band and deviance residuals, sans with panel">
```

The lower axis shows signed Poisson deviance residuals: look for runs of bins
that the fit consistently over- or underestimates. The shaded mean band
contains parameter uncertainty, not the additional fluctuation of future counts.

## Check Shape Dependence

Set the core fraction to one and fix the now-unused width ratio: the same
constructor becomes a single-Gaussian model.

```@example lhcb
single = deepcopy(constructor)
BuildConstructors.update!(single, (core_fraction=1.,))
fix!(single, (:core_fraction, :width_ratio))
single_result = fit_distribution(single, edges, counts;
    solver=NativeMinuitSolver(), tol=1e-3)
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

## See The Yield Uncertainty

Could a different peak width or background level change the signal count?
Fix ``N_s`` at a series of trial values and refit the other six parameters.
Plot the resulting cost increase against the local covariance parabola:

```@example lhcb
scan = profile(result, 6; nsigma=2.5, npoints=31, on_failure=:throw)
interval = profile_interval(scan)
@printf("Ns profile interval at delta(-2 log L)=1: [%.1f, %.1f]\n",
    interval.lower, interval.upper)
profile_options = (local_sigma=result.param_stderr[6],
    title="Uncertainty of the signal count", xlabel="signal count Ns",
    ylabel="Delta (-2 log L)")
profile_figure = plot_profile(scan; profile_options...)
@assert isfinite(interval.lower) && isfinite(interval.upper) # hide
@assert isapprox(interval.lower, 6391.4184; atol=1.) && isapprox(interval.upper, 6578.6286; atol=1.) # hide
# Regression: these two outer refits once reported success above the minimum. # hide
@assert isapprox(scan.delta_cost[[5, 27]], [3.43902, 3.27872]; atol=0.005) # hide
nothing # hide
```

```@setup lhcb
fig = plot_profile(scan; theme=:sans, appearance=:light, profile_options...)
save("lhcb_yield_profile_sans_light.svg", fig)
```

```@raw html
<img class="scientificfitting-plot" src="lhcb_yield_profile_sans_light.svg" alt="Signal-yield profile and local covariance parabola, sans style">
```

The curve is close to a parabola in the relevant range. Its crossings at
``\Delta(-2\log L)=1`` give about ``[6391,6579]`` candidates, consistent with
``N_s\pm94``. This is an approximate **one-parameter 68.3% interval**, with the
shape and background refitted rather than frozen
([why a profile is not a slice](../statistics.md#Why-A-Profile-Is-Not-A-Slice)).

To see the trade-off between signal and background, vary both yields together.
At each grid point, refit the five shape parameters:

```@example lhcb
joint = ScientificFitting.contour(result, 6, 7;
    npoints=25, nsigma=3, on_failure=:throw)
contour_options = (local_covariance=result.param_covariance,
    local_center=result.params[[6, 7]], title="Signal and background counts",
    xlabel="signal count Ns", ylabel="background count Nb")
contour_figure = plot_contour(joint; contour_options...)
@printf("Local correlation of signal and background yields: %.3f\n",
    result.param_correlation[6,7])
@assert all(isfinite, joint.delta_cost) # hide
@assert result.param_correlation[6,7] < 0 # hide
@assert minimum(vcat(joint.delta_cost[1,:], joint.delta_cost[end,:], joint.delta_cost[:,1], joint.delta_cost[:,end])) > maximum(joint.levels) # hide
nothing # hide
```

```@setup lhcb
fig = plot_contour(joint; theme=:sans, appearance=:light, contour_options...)
save("lhcb_yield_contour_sans_light.svg", fig)
```

```@raw html
<img class="scientificfitting-plot" src="lhcb_yield_contour_sans_light.svg" alt="Joint profiled signal and background confidence regions with local covariance ellipses, sans style">
```

The filled regions use ``\Delta(-2\log L)=2.30,6.18`` for approximate **joint
68.3% and 95.45% coverage**; dashed ellipses show the local covariance
approximation. Their tilt shows how increasing one yield can be compensated
by decreasing the other; one- and two-parameter thresholds differ
([Profiles And Contours](../statistics.md#Profiles-And-Contours)). Neither
calculation includes uncertainty from choosing the wrong peak or background
shape.

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
contribute to the peak. A precision mass measurement would also require
momentum-scale calibration.
