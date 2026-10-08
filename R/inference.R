# Decision rule, reporting and simulation (main-text Section 3.3).

# ---- Step-up rule on the average local false discovery rate ----
csm_stepup <- function(lfdr, q = 0.1) {
  if (!is.numeric(lfdr) || anyNA(lfdr) || any(lfdr < 0 | lfdr > 1)) {
    stop("'lfdr' must be numeric values in [0, 1] without missing values.", call. = FALSE)
  }
  .check_scalar(q, "q", 0, 1)
  J <- length(lfdr)
  o <- order(lfdr, seq_len(J))
  avg <- cumsum(lfdr[o]) / seq_len(J)
  ok <- which(avg <= q)
  r <- if (length(ok)) max(ok) else 0L
  rejected <- logical(J)
  if (r > 0L) rejected[o[seq_len(r)]] <- TRUE
  list(rejected = rejected, n_rejected = r,
       cutoff = if (r > 0L) lfdr[o[r]] else NA_real_,
       average = if (r > 0L) avg[r] else NA_real_)
}

# ---- Extrapolation diagnostic ----
csm_extrapolation <- function(fit, Z, target = c("pleiotropy", "replication", "mediation"),
                              t = 1L, rejected, cutoff, pairs = NULL, margin = 3) {
  target <- match.arg(target)
  Z <- .check_z(Z, fit$spec$K)
  if (!is.logical(rejected) || length(rejected) != nrow(Z)) {
    stop("'rejected' must be a logical vector with one entry per row of 'Z'.", call. = FALSE)
  }
  p <- fit$params
  used <- p$varpi > 0
  cap <- vapply(seq_len(fit$spec$K), function(k) {
    v <- p$mu[used & fit$spec$grp_A[, k], k]
    if (length(v)) max(v) else fit$spec$a_mu
  }, 0) + margin
  flag <- rep(NA, nrow(Z))
  lw <- rep(NA_real_, nrow(Z))
  if (any(rejected)) {
    Zr <- Z[rejected, , drop = FALSE]
    Zw <- sign(Zr) * pmin(abs(Zr), rep(cap, each = nrow(Zr)))
    lw[rejected] <- csm_lfdr(fit, Zw, target, t, pairs)
    flag[rejected] <- lw[rejected] > cutoff
  }
  list(flag = flag, lfdr_winsorized = lw, caps = cap)
}

# ---- Fit-based inference for one target ----
csm_test <- function(fit, Z, target = c("pleiotropy", "replication", "mediation"), t = 1L,
                     q = 0.1, pairs = NULL, adjust = FALSE, extrapolation = TRUE,
                     allow_unconverged = FALSE) {
  target <- match.arg(target)
  if (!inherits(fit, "csm_fit")) stop("'fit' must come from csm_fit().", call. = FALSE)
  Z <- .check_z(Z, fit$spec$K)
  J <- nrow(Z)
  out <- list(target = target, t = t, q = q, J = J)
  if ((!fit$converged && !allow_unconverged) || is.null(fit$params)) {
    out$status <- "failed"
    out$lfdr <- rep(NA_real_, J)
    out$rejected <- logical(J)
    out$n_rejected <- 0L
    return(structure(out, class = "csm_test"))
  }
  out$status <- if (fit$converged) "ok" else "unconverged"
  lfdr <- csm_lfdr(fit, Z, target, t, pairs, cores = fit$control$cores)
  su <- csm_stepup(lfdr, q)
  out$lfdr <- lfdr
  out$rejected <- su$rejected
  out$n_rejected <- su$n_rejected
  out$cutoff <- su$cutoff
  if (extrapolation) {
    ex <- csm_extrapolation(fit, Z, target, t, su$rejected, su$cutoff, pairs)
    out$extrapolation <- ex$flag
    out$lfdr_winsorized <- ex$lfdr_winsorized
  }
  if (adjust) out$adjusted <- csm_adjust(Z, lfdr, target, t, q)
  structure(out, class = "csm_test")
}

predict.csm_fit <- function(object, newdata, target = c("pleiotropy", "replication", "mediation"),
                            t = 1L, pairs = NULL, ...) {
  csm_lfdr(object, newdata, match.arg(target), t, pairs)
}

print.csm_test <- function(x, ...) {
  cat(sprintf("Correlated csmGmm test: %s, t = %d, q = %g, %d features\n", x$target,
              as.integer(x$t), x$q, x$J))
  if (x$status == "failed") {
    cat("  Fit failed: no discoveries.\n")
    return(invisible(x))
  }
  if (x$status == "unconverged") cat("  WARNING: based on an unconverged fit.\n")
  cat(sprintf("  discoveries (unadjusted lfdr): %d\n", x$n_rejected))
  if (!is.null(x$extrapolation) && x$n_rejected > 0L) {
    cat(sprintf("  flagged as extrapolation-driven: %d\n", sum(x$extrapolation, na.rm = TRUE)))
  }
  if (!is.null(x$adjusted)) cat(sprintf("  discoveries (monotone adjustment): %d\n", x$adjusted$n_rejected))
  invisible(x)
}

# ---- Simulate from a correlated csmGmm ----
csm_simulate <- function(n, params, spec) {
  if (!inherits(spec, "csm_spec")) stop("'spec' must come from csm_spec().", call. = FALSE)
  .check_scalar(n, "n", 1)
  cp <- spec$comp
  w <- params$varpi[cp$grp] / 2^cp$size
  r <- sample.int(length(w), n, replace = TRUE, prob = w)
  Mean <- cp$H[r, , drop = FALSE] * params$mu[cp$grp[r], , drop = FALSE]
  Z <- Mean + matrix(stats::rnorm(n * spec$K), n) %*% chol(params$Sigma)
  list(Z = Z, H = cp$H[r, , drop = FALSE], group = cp$grp[r])
}
