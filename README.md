# csmGmmCorr

R package for the correlated conditionally symmetric multidimensional Gaussian
mixture (correlated csmGmm) of Liu and Sun. The package tests composite null
hypotheses of the form "at most *t* nonzero effects" when the *K* statistics of a
feature are correlated, for example statistics of several traits measured in the same
individuals.

The model is csmGmm (Sun, McCaw and Lin, 2025, *JASA* 120:605-617) with the identity
covariance replaced by a common correlation matrix Σ with smallest eigenvalue at least κ.
Everything else in csmGmm is kept:

- all 3^K signed configurations;
- support-specific weights;
- support-specific magnitude components;
- the magnitude-ordering constraint.

One fit serves three targets and every threshold *t*:

- **pleiotropy:** more than *t* nonzero effects;
- **replication:** more than *t* nonzero effects, all with the same sign;
- **mediation** (*K* = 2G, coordinates ordered (α₁, β₁, α₂, β₂, ...)): more than *t*
  nonzero effects, including at least one complete pair (α_g, β_g).

## Installation

```r
# install.packages("remotes")
remotes::install_github("<github-user>/csmGmmCorr")
```

The package depends only on base R (`stats`). It contains a small C kernel, so a C
compiler is needed (Rtools on Windows, the Xcode command line tools on macOS).

## Example

```r
library(csmGmmCorr)
spec <- csm_spec(K = 4, b_mu = 15)               # 2 magnitude components for supports of
                                                 # size one, 1 for larger supports
fit  <- csm_fit(Z, spec, starts = c("null", "shared"))   # penalized fit (default)
lf   <- csm_lfdr(fit, Z, target = "pleiotropy", t = 1)   # local false discovery rates
rej  <- csm_stepup(lf, q = 0.1)$rejected                 # step-up rule
adj  <- csm_adjust(Z, lf, "pleiotropy", 1, q = 0.1)      # monotone adjustment
csm_test(fit, Z, "replication", t = 2, q = 0.1)          # all steps in one call
```

`Z` is an N x K matrix of z-statistics, one row per feature.

## Functions and the manuscript

| Manuscript | Function |
|---|---|
| Model, magnitude components and ordering constraint (Section 3.1) | `csm_spec`, `csm_structure`, `csm_params` |
| Penalized estimator, starting values, selection of the start (Section 3.2, Algorithm 1) | `csm_fit`, `csm_control` |
| Target-specific local false discovery rates (Section 3.3) | `csm_lfdr`, `predict` |
| Step-up rule (Section 3.3, Equation 11) | `csm_stepup` |
| Extrapolation check (Remark 1) | `csm_extrapolation` |
| Monotone adjustment (Section 4.3, Theorem 3) | `csm_adjust` |
| Incongruity counts (Supplementary Section 7) | `csm_incongruity` |
| Null configurations of each target | `csm_null_components` |
| Simulation from a fitted or specified model | `csm_simulate` |
| Comparison of magnitude structures by BIC | `csm_select` |
| csmGmm of Sun et al. (Σ fixed at I) and plug-in Σ | `csm_control(fix_sigma = TRUE)`, `csm_control(fix_sigma = R)` |
| Conversion to and from the csmGmm list format | `csm_params`, `as_sun_lists` |

Defaults follow the manuscript: penalty constants λ_μ = 30 and λ_w = 10
(`csm_control(penalty = 30, weight_penalty = 10)`; `0, 0` gives maximum likelihood),
lower bound of the nonzero magnitudes a_μ = 0.5, κ = 0.05, and two magnitude components
for supports of size one and one otherwise. The manuscript uses an upper bound
b_μ = 15 (`csm_spec(b_mu = 15)`).

## Computation

The E-step works on sign pairs (h, −h), which share a weight and a quadratic term, and
visits the sign patterns of a support in Gray-code order. It runs in several threads
(`csm_control(cores = )`). The magnitude update solves a box-constrained quadratic
program per support, using ADMM when the ordering constraint binds. The correlation
update is a projected-gradient step onto {diag = 1, eigenvalues ≥ κ}. Each block is
accepted only if it does not decrease its conditional objective, and the iterations are
accelerated by SQUAREM. The cost per iteration is of order N times the number of
signed configurations; the full model is practical for K up to about 6.

## Tests

`tests/run-tests.R` checks the exact lfdr against brute force, the numerical
examples of the manuscript (Example 1), the M-step blocks and the monotone adjustment.

## Reproducing the manuscript

All simulations, analyses, figures and tables of the manuscript are reproduced by the
companion repository `csmGmmCorr_reproduce`.

## Citation

Liu Z. and Sun R. Correlated csmGmm: composite null testing with correlated statistics.
Manuscript, 2026.

Sun R., McCaw Z. R. and Lin X. (2025). Testing a large number of composite null
hypotheses using conditionally symmetric multidimensional Gaussian mixtures in
genome-wide studies. *Journal of the American Statistical Association*, 120, 605-617.

## License

GPL (>= 3).
