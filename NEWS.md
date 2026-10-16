# lwqs 0.5.1

Maintenance release to restore the package to CRAN (the 0.5.0 version was
archived on 2026-02-10 because the maintainer email was undeliverable) and to
update for current dependency versions.

## Maintainer

* Updated maintainer email address.

## Bug fixes

* `lwqs()` now extracts the time-varying mixture effect and its standard error
  through the documented `predict()` interface for `mgcv`/`gamm4` fits, instead
  of reading positionally from the list invisibly returned by `plot.gam()`. The
  previous approach depended on undocumented element ordering and also drew a
  plot as a side effect.

## Improvements

* Added a `ci_level` argument to `lwqs()` (default 0.95); the confidence band is
  now built from the corresponding normal critical value rather than a
  hard-coded 1.96.
* Removed the hard-coded `ylim(-1, 1)` on the effect plot, which could silently
  clip effects and confidence bands that fell outside that range.
* Added input validation for `data`, `timevar`, `outcome`, and `ID` with
  informative error messages.
* The per-time-point WQS results are now assembled by named column rather than
  by position, removing a potential silent column-misalignment.

## Compatibility

* Updated examples and vignette to use the `gWQS` argument `b_constr`
  (`b1_constr` is deprecated in current `gWQS`).
* Replaced deprecated `ggplot2` `size` aesthetic with `linewidth` for line
  geoms.
* `geom_smooth()` is now called with an explicit `formula = y ~ x`.
