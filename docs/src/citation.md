# Citation And License

See [Packages and Interfaces](interfaces.md) for the technical package overview
and executable integration examples.

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
