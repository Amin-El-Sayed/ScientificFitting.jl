# Related Packages And Citation

## Related Packages

ScientificFitting brings measurement-error models, parameter constraints,
diagnostics, profiles, and editable plots into one fitting workflow. Its role
is to connect statistical analysis to numerical tools, not to replace their
specialized algorithms. These Julia packages cover complementary parts of that
workflow or offer a different modeling approach:

| Package | Focus and relationship to ScientificFitting |
|---|---|
| [LsqFit.jl](https://julianlsolvers.github.io/LsqFit.jl/latest/) | Nonlinear least squares with Levenberg-Marquardt; already powers our compatible Gaussian least-squares path. |
| [Optimization.jl](https://docs.sciml.ai/Optimization/stable/) and [Optim.jl](https://docs.sciml.ai/Optimization/stable/optimization_packages/optim/) | A solver interface and Julia optimization algorithms, respectively. ScientificFitting already uses them for general objectives and constrained fits. |
| [NonlinearSolve.jl](https://docs.sciml.ai/NonlinearSolve/stable/solvers/nonlinear_least_squares_solvers/) | Residual-based solvers, including Gauss-Newton, Levenberg-Marquardt, and trust-region methods. A candidate for an additional least-squares path, not currently integrated. |
| [NativeMinuit.jl](https://github.com/fkguo/NativeMinuit.jl) / [Minuit2.jl](https://github.com/JuliaHEP/Minuit2.jl) | A Julia port / Julia bindings to C++ Minuit2. Both expose minimization and likelihood error analysis; neither is currently a selectable ScientificFitting backend. |
| [RooFitLite.jl](https://github.com/JuliaHEP/RooFitLite.jl) | RooFit-style model construction for Julia HEP analyses, with optional Minuit2 fitting. An alternative model-building workflow; there is no direct adapter yet. |
| [Distributions.jl](https://juliastats.org/Distributions.jl/stable/) | Probability distributions and density/mass evaluation. Its `logpdf` methods can already be used in our likelihood callbacks. |
| [NumericalDistributions.jl](https://github.com/mmikhasenko/NumericalDistributions.jl) / [DistributionsHEP.jl](https://github.com/JuliaHEP/DistributionsHEP.jl) | Numerically normalized distributions / HEP-specific shapes and extended mixtures. Targets for tested distribution-object integration in v0.3. |
| [BuildConstructors.jl](https://github.com/RUB-EP1/BuildConstructors.jl) | Constructs domain models from named parameters and their metadata. A solver-independent integration target for v0.3, not another minimizer. |
| [GLM.jl](https://juliastats.org/GLM.jl/stable/) | Linear and generalized linear regression with formulas, tables, and link functions. A focused alternative when that model structure fits the analysis. |
| [Turing.jl](https://turinglang.org/docs/core-functionality/) | Probabilistic modeling and posterior inference. Appropriate when the target is a posterior distribution rather than the likelihood-based estimates and intervals provided here. |

The [v0.3 integration scope](backend_design.md#v03-ecosystem) separates planned
adapters from currently supported APIs. A related package is not necessarily
a drop-in backend: probability models, minimizers, and posterior samplers have
different contracts. Existing callbacks remain available for custom models.

## Inspiration And Attribution

Part of the motivation came from kafe2's treatment of measurement uncertainties,
fit reports, and profile/contour diagnostics:

- J. Gäßler, G. Quast, D. Savoiu, and C. Verstege,
  [*kafe2 -- a Modern Tool for Model Fitting in Physics Lab Courses*](https://arxiv.org/abs/2210.12768),
  arXiv:2210.12768 (2022).

kafe2 offers a Python/YAML fitting workflow with Matplotlib. ScientificFitting
has a Julia numerical core, optional Makie plotting, and a Python interface
with native Matplotlib plots. Neither requires adopting the other's workflow.

ScientificFitting also uses Julia packages for differentiation, quadrature, and
plotting. When an upstream method contributes materially to a publication,
follow that project's citation guidance as well.

## Citation

ScientificFitting includes a root-level `CITATION.cff` with the author, package
version, and license. If the software contributes to a scientific result,
citation is requested. Record the exact release or commit together with the
Julia version and the uncertainty or likelihood model used in the analysis.

Citation is a request for scientific attribution, not a condition of use.

## License

ScientificFitting is distributed under the
[MIT License](https://opensource.org/license/mit). The copyright and permission
notice must remain with copies or substantial portions of the software.

ScientificFitting was made with care and with assistance from AI tools; every
released change remains subject to human review and approval.
