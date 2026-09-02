# Citation And Related Work

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

## Related Work

Part of the motivation for ScientificFitting came from kafe2, especially its
practical treatment of measurement uncertainties, readable fit reports, and
profile/contour diagnostics:

- J. Gäßler, G. Quast, D. Savoiu, and C. Verstege,
  [*kafe2 -- a Modern Tool for Model Fitting in Physics Lab Courses*](https://arxiv.org/abs/2210.12768),
  arXiv:2210.12768 (2022).

kafe2 provides Python and YAML interfaces and uses Matplotlib for plotting.
ScientificFitting provides native Julia problem and result types, with plotting
supplied separately through CairoMakie so that fitting, diagnostics, profiles,
and text reports remain usable without loading a plotting backend.

ScientificFitting also uses established Julia packages for optimization,
automatic differentiation, probability distributions, quadrature, and
plotting. When a publication depends materially on a particular numerical
method, follow that upstream project's citation guidance as well.

ScientificFitting was made with care and with assistance from AI tools; every
released change remains subject to human review and approval.
