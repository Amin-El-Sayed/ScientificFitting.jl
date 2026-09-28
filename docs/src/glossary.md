# Glossary

**Pull** — a residual divided by its standard uncertainty,
``z_i=(d_i-m_i)/\sigma_i``, with data value ``d_i``, model prediction
``m_i``, and standard uncertainty ``\sigma_i``. Under a correct model, pulls
scatter around zero with scale one. With a non-diagonal covariance the
reported pulls are whitened residuals (see Whitening); they keep the
scale-one property but no longer refer to individual data points.

**Whitening** — re-expressing the residual vector ``r``, whose covariance is
``V``, as ``z=L^{-1}r`` with ``V=LL^T``, so that ``\chi^2=z^Tz``. It changes
coordinates, not the data or the model.

**Effective variance** — the propagated per-point variance
``\sigma_{y,i}^2+(\partial f/\partial x)^2|_{x_i}\,\sigma_{x,i}^2``, where
``f`` is the model curve evaluated at the current parameters; it turns
x uncertainty into a parameter-dependent y covariance.

**ndf** — degrees of freedom: observations (including Gaussian priors and
constraints, which count as extra observations) minus free parameters.
Reduced statistics and p-values require ``\mathrm{ndf}>0``.

**Deviance** — the Poisson likelihood-ratio statistic against the saturated
model, the reference model with one free expectation per count that
reproduces every observation exactly; it fills the chi-square fields for
count fits and shares their asymptotic calibration.

**Profile** — the cost minimized over all remaining parameters (the nuisance
parameters) at each fixed value of one parameter of interest. Freezing the
others at their best-fit values (a slice) yields intervals that are too
narrow; the re-minimization at each fixed value is what defines the profile.

**Profile interval** — the parameter range where the profiled cost exceeds
the cost at the global minimum by at most the threshold,
``\Delta C = C_\mathrm{profile} - C_\mathrm{min} \le`` `threshold`.
`threshold=1` corresponds to 68.27% asymptotically for one parameter;
``\Delta C`` is divided by the same covariance scale as `param_stderr`, so
both report matching uncertainties.

**Confidence band** — the pointwise 68.27% (at `nsigma=1`) interval for the
fitted mean curve from the parameter covariance.

**Prediction band** — the confidence band widened by the observation
uncertainty: where a new measurement would plausibly land.

**Covariance scaling** — multiplying the parameter covariance by
``\chi^2/\mathrm{ndf}`` when no y or x uncertainties were supplied, so the
residual scale is estimated from the data (`scale_covariance=:auto`);
profiles and contours follow the same scale.

**Finding** — one structured diagnostic (severity, stable code, evidence,
recommendation) produced by `diagnose`; the dashboard aggregates findings
into `ok`/`review`/`stop`.
