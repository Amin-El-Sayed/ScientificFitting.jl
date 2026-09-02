# Contributing

ScientificFitting welcomes bug reports, scientific examples, documentation corrections,
performance evidence, and focused code contributions. The package treats
statistical meaning and failure behavior as part of the public API.

## Report A Problem

Please include:

- a minimal executable Julia example;
- Julia and ScientificFitting versions;
- the fit entry point, uncertainty model, and solver controls used;
- the complete error or diagnostic output;
- small redistributable data, or a synthetic reproduction when the original
  data cannot be shared;
- the expected behavior and why it is scientifically justified.

Do not post confidential measurements, credentials, or private file paths.

## Propose A Change

Open an issue before a large architectural change. A numerical or statistical
change should state its mathematical convention, assumptions, failure cases,
and an analytic or independent numerical reference. A performance claim needs
a reproducible benchmark on named hardware; a screenshot alone is not evidence.

Keep pull requests narrow. New public behavior needs a docstring, focused test,
and documentation update. Do not commit generated `docs/build/`, example output,
local benchmark output, or unrelated formatting changes.

## Development Workflow

Read [DEVELOPMENT.md](DEVELOPMENT.md) for branches and test commands and
[MAINTAINERS.md](MAINTAINERS.md) for numerical and plotting invariants. Run the
focused gate for the changed subsystem and `git diff --check` before opening a
pull request. The maintainer may request the broader package gate for changes
that affect statistical semantics or shared architecture.

By contributing, you agree that your contribution is licensed under the
repository's MIT License.
