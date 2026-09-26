# Glossary

**Pull** — a residual divided by its standard uncertainty,
``z_i=(d_i-m_i)/\sigma_i``. Under a correct model, pulls scatter around zero
with scale one; with non-diagonal covariance, whitened components lose their
per-point meaning.

**Whitening** — re-expressing correlated residuals as
``z=L^{-1}r`` with ``V=LL^T``, so that ``\chi^2=z^Tz``. It changes
coordinates, not the data or the model.

**Effective variance** — the propagated per-point variance
``\sigma_y^2+(\partial f/\partial x)^2\sigma_x^2`` that turns x uncertainty
into a parameter-dependent y covariance.

**ndf** — degrees of freedom: observations (including Gaussian auxiliary
terms) minus free parameters. Reduced statistics and p-values require
``\mathrm{ndf}>0``.

**Deviance** — the Poisson likelihood-ratio statistic against the saturated
model; it fills the chi-square fields for count fits and shares their
asymptotic calibration.

**Profile** — the cost minimized over all nuisance parameters at each fixed
value of one parameter of interest. A slice (nuisances frozen) overstates the
information; the refit is the definition.

**Profile interval** — the parameter range where the profiled
``\Delta C`` stays below a threshold; `threshold=1` corresponds to 68.27%
asymptotically and is reported on the same scale as `param_stderr`.

**Confidence band** — the pointwise 68.27% (at `nsigma=1`) interval for the
fitted mean curve from the parameter covariance.

**Prediction band** — the confidence band widened by the observation
uncertainty: where a new measurement would plausibly land.

**Covariance scaling** — multiplying the parameter covariance by
``\chi^2/\mathrm{ndf}`` when the residual scale was estimated from the data
(`scale_covariance`); profiles and contours follow the same scale.

**Finding** — one structured diagnostic (severity, stable code, evidence,
recommendation) produced by `diagnose`; the dashboard aggregates findings
into `ok`/`review`/`stop`.
