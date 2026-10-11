# Decisions, diagnostics, reporting, and simulation.

#' Step-up rule based on average local false discovery rates
#' @export
csm_stepup <- function(lfdr, q = 0.1) {
  if (!is.numeric(lfdr) || anyNA(lfdr) || any(lfdr < 0 | lfdr > 1)) {
    stop("'lfdr' must contain values in [0, 1] without missing values.", call. = FALSE)
  }
  q <- .check_numeric_scalar(q, "q", 0, 1)
  n <- length(lfdr)
  ordering <- order(lfdr, seq_len(n))
  average <- cumsum(lfdr[ordering]) / seq_len(n)
  eligible <- which(average <= q)
  r <- if (length(eligible)) max(eligible) else 0L
  rejected <- logical(n)
  if (r) rejected[ordering[seq_len(r)]] <- TRUE
  list(
    rejected = rejected,
    n_rejected = r,
    cutoff = if (r) lfdr[ordering[r]] else NA_real_,
    average = if (r) average[r] else NA_real_
  )
}

#' Diagnose discoveries driven by extrapolation beyond fitted magnitudes
#' @export
csm_extrapolation <- function(fit, z,
                              target = c("pleiotropy", "replication", "mediation"),
                              t = 1L, rejected, cutoff, pairs = NULL, margin = 3) {
  target <- match.arg(target)
  z <- .check_z(z, fit$spec$K)
  if (!is.logical(rejected) || length(rejected) != nrow(z)) {
    stop("'rejected' must have one logical value per row of 'z'.", call. = FALSE)
  }
  par <- fit$params
  used <- par$varpi > 0
  caps <- vapply(seq_len(fit$spec$K), function(k) {
    values <- par$mu[used & fit$spec$grp_A[, k], k]
    (if (length(values)) max(values) else fit$spec$a_mu) + margin
  }, numeric(1))
  flag <- rep(NA, nrow(z))
  winsorized <- rep(NA_real_, nrow(z))
  if (any(rejected)) {
    zr <- z[rejected, , drop = FALSE]
    zw <- sign(zr) * pmin(abs(zr), rep(caps, each = nrow(zr)))
    winsorized[rejected] <- csm_lfdr(fit, zw, target, t, pairs)
    flag[rejected] <- winsorized[rejected] > cutoff
  }
  list(flag = flag, lfdr_winsorized = winsorized, caps = caps)
}

#' Complete target-specific analysis from a fitted model
#' @export
csm_test <- function(fit, z,
                     target = c("pleiotropy", "replication", "mediation"),
                     t = 1L, q = 0.1, pairs = NULL, adjust = FALSE,
                     extrapolation = TRUE, allow_unconverged = FALSE) {
  target <- match.arg(target)
  if (!inherits(fit, "csm_fit")) stop("'fit' must come from csm_fit().", call. = FALSE)
  z <- .check_z(z, fit$spec$K)
  out <- list(target = target, t = t, q = q, J = nrow(z))
  if ((!isTRUE(fit$converged) && !allow_unconverged) || is.null(fit$params)) {
    out$status <- "failed"
    out$lfdr <- rep(NA_real_, nrow(z))
    out$rejected <- logical(nrow(z))
    out$n_rejected <- 0L
    return(structure(out, class = "csm_test"))
  }
  out$status <- if (isTRUE(fit$converged)) "ok" else "unconverged"
  out$lfdr <- csm_lfdr(fit, z, target, t, pairs, cores = fit$control$cores)
  decision <- csm_stepup(out$lfdr, q)
  out$rejected <- decision$rejected
  out$n_rejected <- decision$n_rejected
  out$cutoff <- decision$cutoff
  if (extrapolation) {
    diagnostic <- csm_extrapolation(fit, z, target, t, out$rejected,
                                    decision$cutoff, pairs)
    out$extrapolation <- diagnostic$flag
    out$lfdr_winsorized <- diagnostic$lfdr_winsorized
  }
  if (adjust) out$adjusted <- csm_adjust(z, out$lfdr, target, t, q)
  structure(out, class = "csm_test")
}

#' @export
predict.csm_fit <- function(object, newdata,
                            target = c("pleiotropy", "replication", "mediation"),
                            t = 1L, pairs = NULL, ...) {
  csm_lfdr(object, newdata, match.arg(target), t, pairs)
}

#' @export
print.csm_test <- function(x, ...) {
  cat(sprintf("Correlated csmGmm test: %s, t = %d, q = %g, %d features\n",
              x$target, as.integer(x$t), x$q, x$J))
  if (x$status == "failed") {
    cat("  Fit failed: no discoveries.\n")
    return(invisible(x))
  }
  if (x$status == "unconverged") cat("  WARNING: based on an unconverged fit.\n")
  cat(sprintf("  discoveries (unadjusted lfdr): %d\n", x$n_rejected))
  if (!is.null(x$extrapolation) && x$n_rejected) {
    cat(sprintf("  flagged as extrapolation-driven: %d\n",
                sum(x$extrapolation, na.rm = TRUE)))
  }
  if (!is.null(x$adjusted)) {
    cat(sprintf("  discoveries (monotone adjustment): %d\n", x$adjusted$n_rejected))
  }
  invisible(x)
}

#' Simulate from a correlated csmGmm
#' @export
csm_simulate <- function(n, params, spec) {
  if (!inherits(spec, "csm_spec")) stop("'spec' must come from csm_spec().", call. = FALSE)
  n <- .check_integer_scalar(n, "n", 1L, .Machine$integer.max)
  pair <- spec$comp
  weights <- params$varpi[pair$group] / 2^pair$size
  state <- sample.int(nrow(pair$H), n, replace = TRUE, prob = weights)
  means <- pair$H[state, , drop = FALSE] * params$mu[pair$group[state], , drop = FALSE]
  noise <- matrix(stats::rnorm(n * spec$K), nrow = n) %*% chol(params$Sigma)
  list(Z = means + noise, H = pair$H[state, , drop = FALSE], group = pair$group[state])
}
