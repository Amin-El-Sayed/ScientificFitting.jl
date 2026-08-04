# Damped Oscillator Data

`pohl_wheel_free_decay.csv` is a curated laboratory dataset for a freely
decaying torsion oscillator. It contains every tenth sample from a 50 Hz
acquisition between 20.18 s and 79.98 s, giving 300 unfiltered angle samples at
approximately 5 Hz. No smoothing or interpolation was applied.

Columns:

- `time_s`: time in seconds.
- `phi_rad`: angular displacement in radians.
- `sigma_phi_rad`: assigned standard uncertainty of the angular displacement.

The angle was obtained from the measured path displacement at a radius of
91.9 mm. The acquisition analysis assigned a 1 mm path-length uncertainty, so
the corresponding angular uncertainty is

```text
0.001 m / 0.0919 m = 0.0108814 rad.
```

This value is an instrument-based uncertainty assignment, not an empirical
standard deviation estimated from repeated decay records. The example keeps it
unchanged and lets the goodness-of-fit diagnostics reveal whether the resulting
independent-Gaussian model is consistent with the observed residuals.

The file is intended for documentation, plotting, and regression examples; raw
acquisition notes and unrelated laboratory files are intentionally not
distributed.
