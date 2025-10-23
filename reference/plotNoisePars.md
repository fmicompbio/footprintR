# Plot background noise model fit diagnostics.

Plot background noise model fit diagnostics.

## Usage

``` r
plotNoisePars(fit)
```

## Arguments

- fit:

  A list returned by
  [`estimateNoiseParsWindows()`](https://fmicompbio.github.io/footprintR/reference/estimateNoiseParsWindows.md)
  with `return_data = TRUE`. Must contain:

  - `coefficients`: named numeric vector
    `c(intercept, slopeMean, slopeInvDepth)` for the fitted model.

  - `data`: `data.frame` with columns `mean`, `noise`, `depth`,
    `sample`, holding the filtered per-window points used for the fit.

## Value

A ggplot object.
