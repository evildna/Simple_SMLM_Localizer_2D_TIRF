# Changelog

## 1.0.0 — 2026-11 (first public release)

Rewritten from the unpublished 0.1.0 draft, which had never been compiled or run.

- All spatial parameters in nm (DoG σ, candidate separation, fit window, σ gate), converted with the pixel size using
  the same rules as the manuscript's analysis; defaults changed to the manuscript values (σ gate 86–235 nm instead of
  0.7–2.5 px; fit window and candidate neighbourhood now scale with the pixel size).
- Fit: bounded Levenberg–Marquardt with the reference start values and bounds (was unbounded with a post-hoc
  position filter).
- Detection threshold k × robust SD; MAD over the whole frame (was a subsample); Gaussian filter with scipy-identical
  kernel and mirror boundaries.
- Intensities labelled ADU unless a photon conversion is given; uncertainty only reported when calibrated (was
  labelled "photon" in all cases).
- "RCC" drift correction renamed to localisation-based drift correction (it correlates with the first segment) and
  off by default.
- Render: histogram + Gaussian blur (was per-localisation Gaussians scaled by an uncalibrated uncertainty); optional
  time-binned stack.
- New: plugins.config (menu entry for the jar; the 0.1.0 jar would not have appeared in the menu), Save_CSV option for
  batch/headless use, headless-safe output, locale-independent logging, `t [s]` column, static API for scripting.
- Groovy script generated from the Java source; now compatible with Groovy 2.4 (older Fiji) and statically compiled
  (≈40× faster under Groovy 4).
- Tests, sample data, validation against the reference implementation, CI workflow.
