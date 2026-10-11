# csmGmmCorr

`csmGmmCorr` fits a correlated conditionally symmetric multidimensional
Gaussian mixture for testing composite null hypotheses. It extends the csmGmm
model by replacing the identity covariance with a common correlation matrix
whose smallest eigenvalue is bounded away from zero. The remaining structure
is retained: all signed configurations, support-specific weights and
magnitude components, and the magnitude-ordering restriction.

One fitted model supplies local false discovery rates for three targets:

- pleiotropy: more than `t` nonzero effects;
- replication: more than `t` nonzero effects with one common sign;
- mediation: more than `t` nonzero coordinates and at least one complete
  coordinate pair.

## Installation

```r
install.packages("remotes")
remotes::install_github("ZLiu20/csmGmmCorr")
```

The package imports only `stats`. A C compiler is required for installation.

## Basic analysis

```r
library(csmGmmCorr)

spec <- csm_spec(K = 4, b_mu = 15)
fit <- csm_fit(Z, spec, starts = c("null", "shared"))

lfdr <- csm_lfdr(fit, Z, target = "pleiotropy", t = 1)
discoveries <- csm_stepup(lfdr, q = 0.1)$rejected
ordered <- csm_adjust(Z, lfdr, target = "pleiotropy", t = 1, q = 0.1)

csm_test(fit, Z, target = "replication", t = 2, q = 0.1)
```

Rows of `Z` are features and columns are correlated z-statistics.
`csm_control(fix_sigma = TRUE)` recovers the independence working model;
passing a correlation matrix to `fix_sigma` produces a plug-in correlated fit.

The estimator uses constrained ECM updates from several starts. The E-step
evaluates sign pairs, the magnitude step solves constrained quadratic
subproblems, and the correlation step remains inside the set of correlation
matrices with eigenvalues at least `kappa`. SQUAREM acceleration is enabled by
default. Set both penalties to zero for maximum-likelihood fitting.

`tests/run-tests.R` contains numerical reference checks, constrained-update
checks, fitting checks, and round-trip conversion tests. The companion
`csmGmmCorr_reproduce` repository contains the manuscript simulations,
analyses, tables, and figures.

## Citation

Liu Z. and Sun R. (2026). Correlated csmGmm: composite null testing with
correlated statistics. Manuscript.

Sun R., McCaw Z. R., and Lin X. (2025). Testing a large number of composite
null hypotheses using conditionally symmetric multidimensional Gaussian
mixtures in genome-wide studies. *Journal of the American Statistical
Association*, 120, 605-617.

## License

GPL (>= 3).
