using ScientificFitting, Distributions, DistributionsHEP
import NativeMinuit  # avoid its separate `profile` binding

"""Read the checked-in mass histogram; integer counts retain every selected pair."""
function read_cms_histogram(path)
    rows = split.(readlines(path)[2:end], ',')
    edges = [parse.(Float64, first.(rows)); parse(Float64, last(rows)[2])]
    groups = reduce(vcat, [permutedims(parse.(Int, row[4:6])) for row in rows])
    counts = parse.(Int, getindex.(rows, 3))
    @assert all(>=(0), groups) && vec(sum(groups; dims=2)) == counts
    @assert sum(counts) == 2_551_454 && length(edges) == 301
    @assert all(isapprox.(diff(edges), 0.002; atol=1e-14))
    @assert edges[2:end] == parse.(Float64, getindex.(rows, 2))
    return (; edges, counts, groups)
end

"""
Window-normalized asymmetric core, left power tail and right exponential tail.
First six parameters: center, core scale, alpha_left, n_left, alpha_right, asymmetry.
Masses and scales use GeV/c^2. This is an empirical reconstructed-mass shape.
"""
cms_peak(p) = truncated(DoubleSidedBifurcatedCrystalBallDas(
    p[1], p[2], p[6], p[3], p[4], p[5]), 2.8, 3.4)

"""Model A: one core scale plus exponential background; p[8:9] are yields in millions."""
function cms_exponential(p)
    background = 2.8 + truncated(Exponential(p[7]); upper=0.6)
    return ExtendedMixtureModel([cms_peak(p), background], 1e6 .* p[8:9])
end

# Each scaled Beta is a normalized degree-two Bernstein basis density.
const CMS_BERNSTEIN = [2.8 + 0.6Beta(k+1, 3-k) for k in 0:2]

"""Positive quadratic background; two log weight ratios remove the redundant scale."""
function cms_background(p)
    weights = [exp(p[7]), exp(p[8]), one(p[7])]
    return MixtureModel(CMS_BERNSTEIN, weights ./ sum(weights))
end

"""Model B: same peak, quadratic background; p[9:10] are yields in millions."""
cms_quadratic(p) = ExtendedMixtureModel(
    [cms_peak(p), cms_background(p)], 1e6 .* p[9:10])

"""
Model C: add a wider Gaussian with the same center. p[11] is the core fraction;
p[12] is the ratio of Gaussian width to core scale. These are resolution
components, not distinct particles or identified kinematic populations.
"""
function cms_two_width(p)
    broad = truncated(Normal(p[1], p[2]*p[12]), 2.8, 3.4)
    signal = MixtureModel([cms_peak(p), broad], [p[11], 1-p[11]])
    return ExtendedMixtureModel([signal, cms_background(p)], 1e6 .* p[9:10])
end

"""MIGRAD bounds and numerical steps for model A; steps are not prior uncertainties."""
function cms_baseline_options(n)
    p0 = [3.097, .025, 1.5, 3., 1.5, 0., 1., .85n/1e6, .15n/1e6]
    bounds = ([3.07,.005,.5,1.05,.5,-.5,.05,0.,0.],
              [3.12,.12,4.,40.,4.,.5,30.,6.,6.])
    steps = [.001,.001,.1,.2,.1,.02,.1,.001,.001]
    parameter_names = ["center", "core_scale", "alpha_left", "n_left", "alpha_right",
        "asymmetry", "background_scale", "signal_million", "background_million"]
    return (; p0, bounds, parameter_names, solver=NativeMinuitSolver(; steps))
end

"""Two additional starts spanning narrow/broad cores, without changing the model or cuts."""
function cms_start_guesses(p0)
    return map(((.020, 1.2, 3.), (.050, 2., 5.))) do (width, alpha, n)
        p = copy(p0)
        p[2], p[3], p[4], p[5] = width, alpha, n, alpha
        p
    end
end

"""Fit three hypotheses to the same bins; compare three starts per model by likelihood."""
function fit_cms_models(edges, counts)
    options = cms_baseline_options(sum(counts))
    a = fit_distribution(cms_exponential, edges, counts;
        options..., initial_guesses=cms_start_guesses(options.p0), multistart=3,
        tol=1e-4, maxiters=2000)
    p = a.params
    # Initialize the polynomial from the exponential's slope, not fitted answers.
    p0 = [p[1:6]; .6/p[7]; .3/p[7]; p[8:9]]
    bounds = ([options.bounds[1][1:6]; -5.; -5.; 0.; 0.],
              [options.bounds[2][1:6]; 5.; 5.; 6.; 6.])
    steps = [.001,.001,.1,.2,.1,.02,.05,.05,.001,.001]
    names = ["center", "core_scale", "alpha_left", "n_left", "alpha_right", "asymmetry",
        "background_logweight_left", "background_logweight_mid", "signal_million", "background_million"]
    b = fit_distribution(cms_quadratic, edges, counts; p0, bounds,
        parameter_names=names, solver=NativeMinuitSolver(; steps),
        initial_guesses=cms_start_guesses(p0), multistart=3, tol=1e-4, maxiters=2000)
    p0 = [b.params; .7; 1.8]
    p0[2] = .023  # start with distinct narrow and broad scales
    c = fit_distribution(cms_two_width, edges, counts; p0,
        bounds=([bounds[1]; .01; 1.01], [bounds[2]; .99; 5.]),
        parameter_names=[names; "core_fraction"; "width_ratio"],
        solver=NativeMinuitSolver(steps=[steps; .02; .05]),
        initial_guesses=cms_start_guesses(p0), multistart=3, tol=1e-4, maxiters=4000)
    return (a, b, c)
end

"""Refit model A in disjoint eta groups; group counts must sum to the inclusive histogram."""
function fit_cms_groups(edges, groups, baseline)
    return map(eachcol(groups)) do counts
        options = cms_baseline_options(sum(counts))
        p0 = copy(baseline.params)
        p0[8:9] .*= sum(counts)/sum(groups)
        fit_distribution(cms_exponential, edges, counts; options..., p0,
            initial_guesses=cms_start_guesses(p0), multistart=3, tol=1e-4, maxiters=2000)
    end
end

"""Bin expectations from upstream CDFs, for plotting an already-fitted extended model."""
cms_binmeans(result, edges) = sum(n .* diff(cdf.(d, edges))
    for (d, n) in zip(components(fitted_model(result)),
        DistributionsHEP.yields(fitted_model(result))))

"""Signed Poisson deviance residuals; the zero-count term is defined by continuity."""
cms_residuals(counts, means) = [sign(n-m)*sqrt(max(2*(m-n+(n==0 ? 0 : n*log(n/m))), 0))
    for (n, m) in zip(counts, means)]
