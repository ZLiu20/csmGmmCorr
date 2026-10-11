# csmGmmCorr 1.3.0

- Added `fix_sigma` to fit the independence working model or a supplied
  plug-in correlation matrix with the same mixture constraints.
- Retained maximum-likelihood operation when both penalties are zero.

# csmGmmCorr 1.2.0

- Changed the default component structure to two components for supports of
  size one and one component for larger supports.

# csmGmmCorr 1.1.0

- Added `csm_select()` for BIC comparison of magnitude-component structures.
- Added `csm_structure()` to construct the corresponding component counts.

# csmGmmCorr 1.0.1

- Stabilized zero-mass support sizes under the weight penalty.
- Added the sign-paired native E-step and skipped numerically negligible
  non-null groups while retaining null groups for local-fdr evaluation.

# csmGmmCorr 1.0.0

- Introduced the correlated support-specific csmGmm and constrained ECM
  estimator.
- Added pleiotropy, same-sign replication, and paired-mediation targets.
- Added the local-fdr step-up rule, monotone adjustment, extrapolation and
  incongruity diagnostics, and csmGmm list-format conversion.
- Added SQUAREM acceleration, optional warm starts, and penalized estimation.
