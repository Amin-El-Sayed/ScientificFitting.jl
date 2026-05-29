# JuFitter Test Layout

The test suite is split by purpose, not by implementation file.

- `regression/`: current API and behavior coverage kept stable while the package is refactored.
- `statistics/`: analytic reference cases, likelihood conventions, coverage checks, and covariance semantics.
- `numerics/`: conditioning, fallback paths, optimizer behavior, and derivative accuracy.
- `api/`: user-facing constructor, keyword, validation, and error-message behavior.
- `plots/`: plot generation, layout robustness, export formats, and visual-regression checks.

New tests should be narrow and deterministic. Expensive Monte Carlo or benchmark checks belong in
`benchmarks/` or a dedicated validation job, not in the default unit-test path.
