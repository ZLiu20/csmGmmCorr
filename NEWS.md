# csmGmmCorr 1.3.0

* New `csm_control(fix_sigma = )`: `TRUE` fixes the correlation matrix at the identity, which gives the csmGmm of Sun, McCaw and Lin (2025) (independence working model, `symm_fit_ind_EM`) with the same weights, magnitudes and ordering constraint; a correlation matrix fixes it at that value (plug-in correlation, as in their `symm_fit_cor_EM`, for any K). With `penalty = 0, weight_penalty = 0` these are maximum likelihood estimators. Used as comparison methods.

# csmGmmCorr 1.2.0

* The default number of magnitude components in `csm_spec()` is now `csm_structure(K, 2, 1)`: two components for each support of size one and one for each larger support, the structure of the example in Sun, McCaw and Lin (2025). With one component per support, the spread of effect sizes among single-signal features can be absorbed by components of size-two supports with a small second magnitude, which inflates the false discovery rate for mediation. `M = 1` gives the previous default.

# csmGmmCorr 1.1.0

* New `csm_select()`: fits several numbers of magnitude components per support (by default 1/1, 2/1 and 2/2 for supports of size one / two or more) and keeps the fit with the smallest BIC. With effect sizes that vary, one component per support can let a support absorb the spread of another and bias the local false discovery rates. `csm_structure()` builds the corresponding `M` for `csm_spec()`.

# csmGmmCorr 1.0.1

* Fixed: with the weight penalty, a support size whose total posterior mass underflows to zero (possible for large K) gave a non-finite objective and stopped the fit with an error. Its within-size shares are now defined by their limit (1/n_s), and a non-finite objective ends a run with status "nonfinite" instead of an error.
* Faster E-step kernel, same results to rounding (relative differences below 1e-13): exp(+-g) is updated along the Gray code by one multiplication per pair; groups whose every term is below exp(-50) relative to the largest term of a feature are skipped (null groups are always kept when computing local false discovery rates); and only the group sums of the D statistics used by the M-step are accumulated. About 3x faster at fitted parameters for K = 10.

# csmGmmCorr 1.0.0

* Complete rewrite to match the revised manuscript: the correlated csmGmm with support-specific weights, support-specific magnitude components and the ordering constraint.
* Constrained ECM on the closed simplex (no weight floor), from several starting values, including an embedded shared-magnitude fit.
* Targets: pleiotropy, replication (same sign) and mediation, for any t.
* Added the step-up rule, the monotone adjustment with certified decision mode, the extrapolation diagnostic and the incongruity diagnostic.
* Compatibility with the csmGmm list format.
* The pooled-size and shared-magnitude restrictions are not offered as models; the only restriction is the largest support size S.
* SQUAREM acceleration of the ECM algorithm (on by default; `csm_control(accelerate = FALSE)` for plain ECM).
* E-step in C on sign pairs, enumerated in Gray-code order, multithreaded (`csm_control(cores = )`); optional warm start on a subset of features (`warm_n`, off by default); Aitken stopping rule with tolerance 0.01/N.
* Removed the Python and interval-arithmetic backend and all compiled code.
* Penalized estimator (on by default): `csm_control(penalty = 30, weight_penalty = 10)` shrinks each support's magnitudes toward the mean of supports of the same size and adds pseudo-features to its weight share within that size. Both terms are of order 1/N, so the model and the asymptotic theory are unchanged; `penalty = 0, weight_penalty = 0` gives maximum likelihood.
