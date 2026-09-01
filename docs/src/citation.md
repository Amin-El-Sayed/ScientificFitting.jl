# Citation, License, And Related Work

## Cite The Exact Software Release

Scientific results should identify the software version that produced them.
JuFitter ships a root-level `CITATION.cff`; GitHub can render it as a citation
and export BibTeX. The first archived public release should also receive a
Zenodo DOI, which will become the preferred stable citation target.

At minimum, record:

- the JuFitter version or exact Git commit;
- the Julia version;
- the fit entry point and cost or likelihood;
- the uncertainty, covariance, and parameter-constraint model;
- whether reported intervals came from local covariance or profiles.

This information is more useful for reproducibility than citing an unversioned
repository homepage.

## MIT License And Citation

JuFitter is released under the [MIT License](https://opensource.org/license/mit).
The license permits use, modification, and redistribution subject to preserving
its copyright and permission notice. Academic citation is requested as
scientific attribution, but it is not an additional license condition.

## Related Work

JuFitter is an independent Julia implementation. Its emphasis on readable fit
reports, practical diagnostics, and profile-versus-local-covariance comparisons
was informed in part by kafe2 and its documentation:

- J. Gäßler, G. Quast, D. Savoiu, and C. Verstege,
  [*kafe2 -- a Modern Tool for Model Fitting in Physics Lab Courses*](https://arxiv.org/abs/2210.12768),
  arXiv:2210.12768 (2022).

JuFitter does not contain kafe2 code. Its numerical implementation builds on
Julia packages including LsqFit, Optimization.jl, DifferentiationInterface,
ForwardDiff, Distributions, QuadGK, and the optional Makie plotting ecosystem.
A publication that depends materially on a particular upstream algorithm should
also follow that project's citation guidance; citing JuFitter does not replace
method-specific attribution.

## Software Paper

A peer-reviewed software paper is useful once the public package has external
users, a stable interface, and an open development record. Until then, a
versioned software DOI provides an honest citation target without implying that
the scientific and software review has already happened.
