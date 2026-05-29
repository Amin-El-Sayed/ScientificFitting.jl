# Research Landscape

The fitting ecosystem already contains strong tools. JuFitter should not copy a
single package; it should combine the best ideas in a Julia-native design.

## Relevant Tools

- **kafe2** emphasizes statistically explicit fits, uncertainties in x and y,
  covariance matrices, constraints, profiles, contours, and report generation.
- **iminuit/Minuit** is a reference point for robust numerical minimization,
  covariance estimates, profile scans, and contour-based uncertainty analysis.
- **SciPy** provides widely used low-level optimization and curve-fitting
  primitives with simple entry points.
- **lmfit** adds user-friendly parameter objects, bounds, constraints, reports,
  and confidence intervals on top of SciPy-style workflows.
- **ROOT fitting tools** are common in particle and nuclear physics, especially
  where Minuit-style likelihood fits and contour diagnostics are expected.
- **Julia packages such as LsqFit, Optim, Optimization, Turing, Measurements,
  and Makie** cover important parts of the stack, but not the full
  out-of-the-box scientific fitting experience JuFitter targets.

## What Users Need

Scientific users typically need:

- minimal syntax for standard fits,
- correct covariance and uncertainty propagation,
- clear treatment of x-errors, y-errors, correlations, priors, and constraints,
- profile-likelihood intervals when local covariance is not enough,
- plots that are immediately suitable for reports or papers,
- diagnostics that explain when a result should not be trusted,
- examples that map directly to laboratory and analysis workflows,
- reproducibility through documented model assumptions and exported reports.

## JuFitter Differentiation

JuFitter should become the standard tool by making the correct workflow also
the easiest workflow:

- Defaults should choose statistically defensible costs and plot layouts.
- Advanced controls should remain visible and composable, not hidden behind a
  black-box interface.
- Every result should expose enough structure for reports, notebooks, and
  automated analysis pipelines.
- Documentation should teach statistical fitting, not only describe functions.
