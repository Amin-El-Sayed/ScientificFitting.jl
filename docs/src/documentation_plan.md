# Documentation Strategy

This page is the working contract for the documentation rewrite. It exists to
prevent the site from becoming a loose collection of examples, API fragments,
and generated-looking explanations.

## Reader Paths

The documentation must serve four different readers without forcing them
through the same path.

```@raw html
<div class="jufitter-flow">
  <div class="jufitter-flow-step"><strong>New User</strong><span>Install, quickstart, first plot, first diagnosis.</span></div>
  <div class="jufitter-flow-step"><strong>Example-Driven Scientist</strong><span>Start from a realistic workflow and adapt it.</span></div>
  <div class="jufitter-flow-step"><strong>Practitioner</strong><span>Decide what to do when a fit looks suspicious.</span></div>
  <div class="jufitter-flow-step"><strong>Statistics Reader</strong><span>Understand and justify the mathematical method.</span></div>
  <div class="jufitter-flow-step"><strong>API User</strong><span>Find every argument, default, return value, and failure mode.</span></div>
</div>
```

The site order follows that priority:

1. Getting Started: install, first fit, and how the fit machinery works.
2. Gallery: a gradual progression from basic least squares to advanced
   diagnostics and multi-dataset workflows.
3. Guides.
4. Mathematics and Statistics.
5. Reference.
6. Engineering Notes for maintainers and contributors.

The Gallery comes before Guides because many scientists learn fitting by
recognizing a workflow close to their own experiment. Guides explain general
rules after the reader has seen concrete examples.

Current first-user status: the landing page, installation page, quickstart,
`How JuFitter Works`, and `Fitting for Practitioners` have completed this
page-level pass. They present the
pre-release installation state honestly, show a real plot immediately, use
visible explicit arrays as the only numerical source, and keep every displayed
report synchronized with the executable snapshot gate. The architecture page
now follows the actual typed-problem, validation, objective, solver, result, and
post-fit flow without implying that bounds are additive penalties or that
profile scans are the primary solver. The practitioner guide now starts from
the measurement process, uses ndf-dependent goodness-of-fit reasoning, and
turns profile/contour geometry into concrete reporting decisions. The next
plotting guide now documents the actual optional-Makie boundary, layout and
report controls, style precedence, post-fit composition helpers, uncertainty
bands, and the compute-once profile-matrix workflow. The next editorial unit is
`Statistical Foundations`, followed by the API reference.

## Core Narrative

The documentation should teach one repeated workflow:

1. State the scientific question.
2. Inspect the measured quantities and uncertainties.
3. Choose a model and explain its parameters.
4. Choose the statistical cost and justify it.
5. Fit.
6. Diagnose the result.
7. Report the parameter values, uncertainty, and limitations.
8. Decide whether the model is good enough or needs revision.

Every tutorial and gallery page should follow this pattern unless there is a
specific reason not to.

## Gallery Standard

Gallery pages are not decorative examples. Each one must be a complete
scientific workflow.

Required sections:

- **Question:** what physical, engineering, or statistical quantity is being
  determined?
- **Data:** column meanings, units, uncertainty sources, and why the dataset is
  realistic.
- **Model:** formula, parameter meanings, assumptions, and expected failure
  modes.
- **Fit:** full executable code, not just a fragment.
- **Diagnostics:** residuals, chi-square, p-value, dashboard, profiles, or
  contours where relevant.
- **Interpretation:** final result with uncertainty and what it means.
- **What can go wrong:** at least one realistic warning sign or limitation.

Synthetic data are acceptable only for controlled demonstrations of a method.
The main gallery should use real or realistically messy datasets with visible
uncertainties, imperfect residuals, and meaningful interpretation.

Across the gallery, uncertainty examples must deliberately cover diagonal and
correlated covariance, absolute and relative components, x and y uncertainty,
heteroskedastic point-wise uncertainty, confidence bands, prediction bands, and
derived-quantity propagation. Individual pages should use the uncertainty model
that belongs to their scientific question rather than adding complexity only to
tick a coverage box.

Current coverage and remaining gaps:

| Workflow | Uncertainty semantics | Status |
| --- | --- | --- |
| Linear calibration | heteroskedastic absolute y errors; 1-sigma prediction band | covered |
| Photoelectric threshold | heteroskedastic x/y errors; two 1-sigma fit bands; propagated line intersection | covered |
| Damped oscillator | x/y errors in a nonlinear model | covered |
| Full-covariance decay | dense correlated y covariance | covered |
| XY calibration | effective-variance x/y propagation | covered |
| Constraints and profiles | data uncertainty, prior information, local covariance versus profile regions | covered |
| Poisson and histogram fits | likelihood-defined count uncertainty rather than Gaussian error bars | covered |
| Multi-dataset fit | dataset-specific absolute y errors with shared parameters | covered |
| Relative uncertainty components | instrument scale or gain error | missing gallery workflow |
| Correlated x uncertainty | shared calibration uncertainty on the independent variable | missing gallery workflow |

This table is a coverage audit, not a reason to make every figure visually busy.

Visible tutorial code should look like a real notebook. Prefer explicit curated
arrays or small CSV reads. Synthetic or controlled data generation belongs in
asset generator scripts, not in beginner-facing code blocks, unless the page is
explicitly teaching simulation.

Documentation plot style switching must use real Makie-rendered assets for
each supported style and appearance. CSS recoloring or image inversion is not
acceptable because it changes typography, contrast, and scientific hierarchy
outside the plotting backend.

## Mathematics and Statistics Standard

The mathematics section is not a library manual. It must justify the methods in
a scientific context.

Each methods page should answer:

- What assumptions does this method make?
- What quantity is minimized or estimated?
- Which terms are constants, and which affect the fitted parameters?
- When are local covariance errors valid?
- When do profiles or contours replace local symmetric errors?
- What does the p-value or deviance mean, and what does it not mean?
- Which approximation is being used, and when does it fail?
- What should a scientist cite, report, or check in a critical analysis?

API calls may appear at the end of a methods page, but the main body should
explain the method itself.

The minimum statistics curriculum is:

- weighted chi-square and whitening,
- full Gaussian likelihood and log determinants,
- diagonal versus dense covariance,
- x-uncertainty and effective variance approximations,
- Poisson counts and deviance,
- histogram, unbinned, and extended likelihoods,
- priors, bounds, fixed parameters, and constraints,
- degrees of freedom and goodness-of-fit,
- local Hessian/Jacobian covariance,
- profile likelihoods, Wilks thresholds, and contour interpretation,
- AIC, BIC, model comparison, and their limitations.

## Guides Standard

Guides are for decisions during use. They should be shorter than the math
reference and more operational than the gallery.

Examples:

- Which uncertainty input should I use?
- My chi-square is terrible. What do I inspect first?
- When should I use a likelihood fit instead of least squares?
- How do I decide whether local errors are trustworthy?
- How do I make a plot/report for a lab notebook or paper?

Guides should link to gallery pages for full workflows and to mathematics pages
for justification.

## API Reference Standard

Autodocs alone are not sufficient. Every public export must eventually document:

- required arguments,
- keyword arguments and defaults,
- units and statistical semantics where relevant,
- return type,
- failure modes,
- diagnostics emitted,
- one minimal example,
- one realistic example or link to a gallery page.

## Page Completion Gate

A page is not ready unless:

- It has a clear reader and purpose.
- It uses concrete data, formulas, or code rather than generic prose.
- Every plot explains bands, errors, markers, and reports.
- Every code block is executable or explicitly marked as a fragment.
- No private paths, course-internal language, or placeholder text remain.
- The page links to the next useful page in the reader path.
- `julia --project=docs docs/make.jl` builds after the change.

## Automation Rules

Automated documentation work should proceed page by page. A heartbeat or agent
run must:

1. Read this strategy, `ROADMAP.md`, `RELEASE_AUDIT.md`, and `git status`.
2. Pick one page or one tightly scoped navigation change.
3. Improve it to the completion gate above.
4. Run the docs build.
5. Leave the worktree clean or stop with a precise blocker.

Do not spread shallow edits across many pages. One excellent page is more useful
than ten plausible pages.
