# Scope And Alternatives

ScientificFitting is a frequentist fitting package for measured data with
explicit uncertainty models.

## Not Implemented

| Task | Status | Use instead |
| --- | --- | --- |
| Posterior sampling, credible intervals | out of scope | [Turing.jl](https://turinglang.org); a ScientificFitting likelihood can be reused there via the [external-likelihood interface](interfaces.md) |
| ODE/PDE parameter estimation | out of scope | [SciML](https://docs.sciml.ai) (DiffEqParamEstim, SciMLSensitivity) |
| Robust M-estimators (Huber, Tukey) | not implemented | model heavy tails explicitly with `fit_likelihood_model` (for example Laplace or Student-t errors), which keeps the likelihood interpretable |
| Multidimensional x (surfaces, fields) | not implemented | `fit_indexed_model` fits observations indexed 1..n that have no meaningful x coordinate; general regression on multidimensional x needs another tool |
| Errors-in-variables beyond linearization | not implemented | the effective-variance propagation assumes small x errors and smooth models; latent-variable EIV needs an explicit model, for example in Turing.jl |
| Global optimization guarantees | out of scope | all solvers are local; `multistart` and explicit `initial_guesses` mitigate, not guarantee |
| Effect-size or hypothesis-testing frameworks | out of scope | [HypothesisTests.jl](https://github.com/JuliaStats/HypothesisTests.jl) |

## When A Simpler Tool Suffices

For unweighted or diagonally weighted least squares without diagnostics,
[LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl) (which powers this
package's fast path) or [CurveFit.jl](https://docs.sciml.ai/CurveFit/stable/)
are smaller dependencies. Use ScientificFitting when uncertainties are part
of the question: correlated or x errors, count likelihoods, priors, profile
intervals, and reviewable diagnostics.

## Known Approximations

- **Effective variance for x errors** is a first-order linearization; its
  determinant term makes results differ slightly but reproducibly from
  ODR-convention tools ([details](statistics.md#Relation-To-ODR-And-The-York-Method)).
- **Local covariance** is curvature at one point; profiles exist because it
  can fail ([details](statistics.md#Local-Parameter-Covariance)).
- **Wilks thresholds** for profiles and contours are asymptotic; non-regular
  problems need simulation ([details](statistics.md#Profiles-And-Contours)).
- **Prediction and confidence bands** are pointwise intervals (68.27% at the
  default `nsigma=1`) under approximate normality of the estimator, not
  simultaneous bands.
