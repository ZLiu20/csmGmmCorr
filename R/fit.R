# Constrained ECM fitting from several starting values (main-text Section 3.2
# and Supplementary Section S2).

# ---- Control parameters for csm_fit ----
csm_control <- function(max_iter = 5000L, tol_loglik = NULL, tol_param = Inf,
                        accelerate = TRUE, warm_n = 0L, penalty = 30, weight_penalty = 10,
                        block_cells = 5e6,
                        sigma_steps = 50L, fix_sigma = FALSE, cores = 1L, verbose = FALSE) {
  .check_scalar(penalty, "penalty", 0)
  .check_scalar(weight_penalty, "weight_penalty", 0)
  .check_scalar(warm_n, "warm_n", 0)
  .check_scalar(cores, "cores", 1)
  .check_scalar(max_iter, "max_iter", 1)
  if (!is.null(tol_loglik)) .check_scalar(tol_loglik, "tol_loglik", 0)
  if (!is.numeric(tol_param) || length(tol_param) != 1L || is.na(tol_param) || tol_param < 0) {
    stop("'tol_param' must be a nonnegative number or Inf.", call. = FALSE)
  }
  .check_scalar(block_cells, "block_cells", 1e3)
  .check_scalar(sigma_steps, "sigma_steps", 0)
  list(max_iter = as.integer(max_iter), tol_loglik = tol_loglik, tol_param = tol_param,
       warm_n = as.integer(warm_n), penalty = as.numeric(penalty),
       weight_penalty = as.numeric(weight_penalty),
       accelerate = isTRUE(accelerate), block_cells = block_cells,
       sigma_steps = as.integer(sigma_steps), fix_sigma = .check_fix_sigma(fix_sigma), cores = as.integer(cores),
       verbose = isTRUE(verbose))
}

# fix_sigma: FALSE (estimate Sigma), TRUE (Sigma fixed at the identity) or a correlation
# matrix at which Sigma is fixed (e.g. a plug-in estimate).
.check_fix_sigma <- function(x) {
  if (is.matrix(x)) {
    if (!is.numeric(x) || nrow(x) != ncol(x) || any(!is.finite(x)) || max(abs(x - t(x))) > 1e-10 ||
        max(abs(diag(x) - 1)) > 1e-10 || min(eigen(x, TRUE, TRUE)$values) <= 0) {
      stop("'fix_sigma' must be TRUE, FALSE or a positive definite correlation matrix.", call. = FALSE)
    }
    return((x + t(x)) / 2)
  }
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("'fix_sigma' must be TRUE, FALSE or a positive definite correlation matrix.", call. = FALSE)
  }
  x
}
.fixed_sigma <- function(control, K) {
  fs <- control$fix_sigma
  if (is.matrix(fs)) {
    if (nrow(fs) != K) stop("'fix_sigma' must be a K x K matrix.", call. = FALSE)
    return(fs)
  }
  if (isTRUE(fs)) diag(K) else NULL
}

# One ECM map: with the responsibilities of `es` (computed at p), update the
# weights, the magnitudes and the correlation matrix in turn.
.ecm_map <- function(es, p, spec, XtX, control, N) {
  pn <- p
  pn$varpi <- .update_weights(es$T, spec, control$weight_penalty / N)
  st <- .mu_stats(es, spec)
  if (control$penalty > 0) st <- .penalize_mu_stats(st, p$mu, spec, es$Omega, control$penalty / N)
  pn$mu <- .update_mu(st, spec, p$mu, es$Omega)
  if (is.null(.fixed_sigma(control, spec$K))) {   # fixed Sigma: stays at its starting value
    S <- .scatter(es, pn$mu, spec, XtX)
    pn$Sigma <- .update_sigma(S, p$Sigma, spec$kappa, control$sigma_steps)
  }
  pn
}

# Combination a * p0 + b * p1 + c * p2 of parameters, mapped back to the
# feasible region. Weights are combined on the log scale and renormalized, so
# extrapolation never sets a positive weight to zero (an ECM step cannot revive
# a zero weight); magnitudes are projected onto the box and the ordering
# constraint, and the correlation matrix onto {diag = 1, eigenvalues >= kappa}.
.combine <- function(p0, p1, p2, a, b, c, spec) {
  q <- p0
  pos <- p0$varpi > 0 & p1$varpi > 0 & p2$varpi > 0
  lw <- rep(-Inf, length(pos))
  lw[pos] <- a * log(p0$varpi[pos]) + b * log(p1$varpi[pos]) + c * log(p2$varpi[pos])
  # A weight may shrink by at most a factor 1000 relative to the last ECM step,
  # so extrapolation cannot push small weights to numerically dead values.
  lw[pos] <- pmax(lw[pos], log(p2$varpi[pos]) - log(1000))
  lw[pos] <- lw[pos] - max(lw[pos])
  q$varpi <- exp(lw) / sum(exp(lw))
  q$mu <- .project_mu(a * p0$mu + b * p1$mu + c * p2$mu, spec)
  q$Sigma <- .project_cor(a * p0$Sigma + b * p1$Sigma + c * p2$Sigma, spec$kappa)
  q
}

.pdiff <- function(a, b) c(a$varpi - b$varpi, a$mu - b$mu, a$Sigma - b$Sigma)

# Coordinates used for extrapolation: log weights (positive ones), magnitudes
# and the correlation matrix.
.ptrans <- function(p, pos) c(log(p$varpi[pos]), p$mu, p$Sigma)

# Constrained ECM with SQUAREM acceleration (Varadhan and Roland, 2008). Each
# cycle takes two ECM steps p -> p1 -> p2, extrapolates along them, projects
# the extrapolated point onto the feasible region and applies one more ECM
# step. The result is kept only if its log-likelihood is at least that of p2,
# so the log-likelihood never decreases and every accepted point is feasible.
# Acceleration changes the path to the maximum, not the model or the estimator.
.fit_one <- function(X, spec, p0, control) {
  # Default stopping tolerance: a remaining average log-likelihood increase of
  # 0.01 / N, well within the O(1/N) deficit allowed by the FDR theory.
  if (is.null(control$tol_loglik)) control$tol_loglik <- 0.01 / nrow(X)
  XtX <- crossprod(X) / nrow(X)
  est <- function(p) .estep(X, spec, p, block_cells = control$block_cells, cores = control$cores)
  N <- nrow(X)
  lam <- control$penalty
  alw <- control$weight_penalty
  map <- function(es, p) .ecm_map(es, p, spec, XtX, control, N)
  # Objective: average log-likelihood minus (penalty / N) times the
  # within-size dispersion of the magnitudes (zero when penalty = 0).
  est0 <- est
  est <- function(p) {
    e <- est0(p)
    e$J <- e$F - (if (lam > 0) lam * .mag_penalty(p$mu, spec) / N else 0) -
      (if (alw > 0) alw * .weight_penalty(p$varpi, spec) / N else 0)
    e
  }
  p <- p0
  es0 <- est(p)
  trace <- es0$J
  status <- "max_iter"
  maps <- 0L
  eps <- function(F) 1e-10 * (1 + abs(F))
  next_report <- 50L
  cycle_F <- es0$J
  step_max <- 1
  while (maps < control$max_iter) {
    p1 <- map(es0, p)
    es1 <- est(p1)
    p2 <- map(es1, p1)
    es2 <- est(p2)
    maps <- maps + 2L
    if (!is.finite(es1$J) || !is.finite(es2$J)) {
      status <- "nonfinite"
      break
    }
    if (es1$J < es0$J - eps(es0$J) || es2$J < es1$J - eps(es1$J)) {
      status <- "likelihood_decrease"
      break
    }
    trace <- c(trace, es1$J, es2$J)
    # Converged when the Aitken estimate from the two plain ECM steps and, with
    # acceleration, the total gain over the last three cycles are both below
    # the tolerance (a single quiet pair of steps can be misleading).
    gap <- .aitken_gap(c(es0$J, es1$J, es2$J))
    if (control$accelerate) {
      # also require a small total gain over the last three cycles
      hist <- c(cycle_F, es2$J)
      gap <- max(gap, if (length(hist) > 3L) hist[length(hist)] - hist[length(hist) - 3L] else Inf)
    }
    dpar <- max(abs(.pdiff(p2, p1)))
    if (gap <= control$tol_loglik && dpar <= control$tol_param) {
      p <- p2
      es0 <- es2
      status <- "converged"
      break
    }
    acc_p <- p2
    acc_es <- es2
    if (control$accelerate) {
      pos <- p$varpi > 0 & p1$varpi > 0 & p2$varpi > 0
      r <- .ptrans(p1, pos) - .ptrans(p, pos)
      v <- .ptrans(p2, pos) - .ptrans(p1, pos) - r
      nv <- sqrt(sum(v^2))
      if (nv > 0) {
        alpha <- max(-step_max, min(-1, -sqrt(sum(r^2)) / nv))
        pq <- .combine(p, p1, p2, (1 + alpha)^2, -2 * alpha * (1 + alpha), alpha^2, spec)
        cand <- tryCatch({
          esq <- est(pq)
          p3 <- map(esq, pq)
          list(p = p3, es = est(p3))
        }, error = function(e) NULL)
        maps <- maps + 1L
        if (!is.null(cand) && is.finite(cand$es$J) && cand$es$J >= es2$J) {
          acc_p <- cand$p
          acc_es <- cand$es
          trace <- c(trace, acc_es$J)
          # Adaptive maximum step length, as in the SQUAREM package.
          if (alpha <= -step_max) step_max <- 4 * step_max
        } else {
          step_max <- max(1, step_max / 4)
        }
      }
    }
    cycle_F <- c(cycle_F, acc_es$J)
    p <- acc_p
    es0 <- acc_es
    if (control$verbose && maps >= next_report) {
      message(sprintf("  ECM steps %d: average log-likelihood %.10f", maps, es0$J))
      next_report <- next_report + 50L
    }
  }
  list(params = structure(p, class = "csm_params"), loglik = es0$F, objective = es0$J,
       iterations = maps, status = status, trace = trace)
}

.aitken_gap <- function(trace) {
  n <- length(trace)
  if (n < 3L) return(Inf)
  d1 <- trace[n] - trace[n - 1L]
  d0 <- trace[n - 1L] - trace[n - 2L]
  if (d1 == 0 && d0 == 0) return(0)
  if (d1 <= 0) return(Inf)
  if (d0 <= 0) return(Inf)
  a <- d1 / d0
  if (a >= 1) return(Inf)
  d1 / (1 - a)
}

# Warm start: when N exceeds 3 * warm_n, the ECM algorithm is first run on a
# fixed subset of warm_n features and then continued on all N features from
# that point, with the full-data stopping rule. Only the starting point of the
# full-data run changes; the estimator is still the full-data maximum.
.warm_rows <- function(N, m) {
  # deterministic, well-spread subset (golden-ratio scrambling, no RNG use)
  sort(order((as.numeric(seq_len(N)) * 0.6180339887498949) %% 1)[seq_len(m)])
}

.fit_warm <- function(X, spec, p0, control) {
  N <- nrow(X)
  m <- control$warm_n
  if (m < 1L || N <= 3L * m) return(.fit_one(X, spec, p0, control))
  sub <- .fit_one(X[.warm_rows(N, m), , drop = FALSE], spec, p0, control)
  full <- .fit_one(X, spec, unclass(sub$params), control)
  full$warm_iterations <- sub$iterations
  full
}

# ---- Starting values --------------------------------------------------------
.base_magnitudes <- function(X, spec) {
  b <- apply(abs(X), 2, function(v) {
    w <- v[v > 2]
    if (length(w) >= 5L) mean(w) else 2.5
  })
  pmin(spec$b_mu, pmax(spec$a_mu, b))
}

.magnitudes_from_base <- function(base, spec) {
  G <- length(spec$grp_b)
  K <- spec$K
  nb <- 2L^K
  MU <- matrix(0, G, K)
  for (g in which(spec$grp_b != 0L)) {
    b <- spec$grp_b[g]
    Mb <- spec$M[b + 1L]
    m <- spec$grp_m[g]
    f <- if (b == nb - 1L) 1 + 0.5 * (m - 1) / max(1, Mb - 1) else 0.5 + 0.5 * m / Mb
    a <- spec$grp_A[g, ]
    MU[g, a] <- base[a] * f
  }
  .project_mu(MU, spec)
}

.start_params <- function(X, spec, type) {
  G <- length(spec$grp_b)
  base <- .base_magnitudes(X, spec)
  varpi <- switch(type,
    equal = rep(1 / G, G),
    null = ifelse(spec$grp_b == 0L, 0.5, 0.5 / (G - 1L))
  )
  list(varpi = varpi, mu = .magnitudes_from_base(base, spec), Sigma = diag(spec$K))
}

.random_start <- function(X, spec, base_p) {
  G <- length(spec$grp_b)
  K <- spec$K
  v <- base_p$varpi * exp(stats::rnorm(G, 0, 0.5))
  v <- v / sum(v)
  if (spec$magnitudes == "shared") {
    mu <- .base_magnitudes(X, spec) * exp(stats::rnorm(K, 0, 0.25))
    mu <- pmin(spec$b_mu, pmax(spec$a_mu, mu))
    MU <- spec$grp_A * rep(mu, each = G)
  } else {
    MU <- base_p$mu * exp(stats::rnorm(G * K, 0, 0.25))
    MU <- .project_mu(MU, spec)
  }
  rho <- stats::runif(1, 0, 0.3)
  Sigma <- .project_cor((1 - rho) * diag(K) + rho, spec$kappa)
  list(varpi = v, mu = MU, Sigma = Sigma)
}

# Fit of the shared-magnitude restriction, embedded in the full model.
.shared_start <- function(X, spec, control) {
  sh <- .make_spec(spec$K, 1L, spec$S, "shared", spec$a_mu, spec$b_mu, spec$kappa)
  ctl <- control
  ctl$penalty <- 0
  ctl$weight_penalty <- 0
  st <- .start_params(X, sh, "null")
  fs <- .fixed_sigma(ctl, spec$K)
  if (!is.null(fs)) st$Sigma <- fs
  f <- .fit_warm(X, sh, st, ctl)
  mu_shared <- vapply(seq_len(spec$K), function(k) {
    max(f$params$mu[sh$grp_A[, k], k])
  }, 0)
  vb <- f$params$varpi
  names(vb) <- sh$grp_b
  varpi <- vb[as.character(spec$grp_b)] / spec$M[spec$grp_b + 1L]
  varpi <- as.vector(varpi) / sum(varpi)
  list(varpi = varpi, mu = .magnitudes_from_base(mu_shared, spec), Sigma = f$params$Sigma)
}

# ---- Fit a correlated csmGmm ----
csm_fit <- function(X, spec, starts = c("null", "equal", "shared"), n_random = 0L,
                    init = NULL, control = csm_control(), seed = NULL) {
  if (!inherits(spec, "csm_spec")) stop("'spec' must come from csm_spec().", call. = FALSE)
  X <- .check_z(X, spec$K, "X")
  if (nrow(X) < 2L) stop("'X' must have at least two rows.", call. = FALSE)
  if (is.null(spec$b_mu)) spec$b_mu <- max(spec$a_mu + 1, ceiling(max(abs(X))) + 1)
  starts <- unique(starts)
  if (!all(starts %in% c("null", "equal", "shared"))) {
    stop("'starts' may contain \"null\", \"equal\" and \"shared\".", call. = FALSE)
  }
  if (spec$magnitudes == "shared") starts <- setdiff(starts, "shared")
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(init) && inherits(init, "csm_params")) init <- list(init)
  sl <- list()
  for (s in starts) {
    sl[[s]] <- if (s == "shared") {
      tryCatch(.shared_start(X, spec, control), error = function(e) e)
    } else .start_params(X, spec, s)
  }
  for (i in seq_len(n_random)) {
    sl[[paste0("random", i)]] <- .random_start(X, spec, .start_params(X, spec, "null"))
  }
  # Fixed Sigma: the identity gives the csmGmm of Sun, McCaw and Lin (independence working
  # model); a supplied correlation matrix gives a plug-in correlated fit.
  fs <- .fixed_sigma(control, spec$K)
  if (!is.null(fs)) {
    if (min(eigen(fs, TRUE, TRUE)$values) < spec$kappa) {
      stop("'fix_sigma' has an eigenvalue below kappa; project it first.", call. = FALSE)
    }
    for (s in names(sl)) if (!inherits(sl[[s]], "error")) sl[[s]]$Sigma <- fs
  }
  for (i in seq_along(init)) {
    p <- init[[i]]
    bad <- .feasibility(p, spec)
    if (length(bad)) stop(sprintf("init[[%d]] is infeasible: %s.", i, paste(bad, collapse = "; ")), call. = FALSE)
    sl[[paste0("init", i)]] <- unclass(p)
  }
  if (!length(sl)) stop("No starting values requested.", call. = FALSE)
  runs <- vector("list", length(sl))
  names(runs) <- names(sl)
  for (s in names(sl)) {
    if (control$verbose) message("Start: ", s)
    runs[[s]] <- if (inherits(sl[[s]], "error")) {
      list(status = "error", message = conditionMessage(sl[[s]]), loglik = NA_real_, objective = NA_real_,
           iterations = 0L)
    } else {
      tryCatch(.fit_warm(X, spec, sl[[s]], control), error = function(e) {
        list(status = "error", message = conditionMessage(e), loglik = NA_real_, objective = NA_real_,
             iterations = 0L)
      })
    }
  }
  summary <- data.frame(
    start = names(runs),
    status = vapply(runs, function(r) r$status, ""),
    iterations = vapply(runs, function(r) as.integer(r$iterations), 0L),
    warm_iterations = vapply(runs, function(r) if (is.null(r$warm_iterations)) NA_integer_
                             else as.integer(r$warm_iterations), 0L),
    loglik = vapply(runs, function(r) r$loglik, 0),
    objective = vapply(runs, function(r) if (is.null(r$objective)) NA_real_ else r$objective, 0),
    stringsAsFactors = FALSE, row.names = NULL
  )
  # Fits are compared by the objective they maximize: the penalized
  # log-likelihood (equal to the log-likelihood when both penalties are 0).
  ok <- which(summary$status == "converged")
  converged <- length(ok) > 0L
  pick <- if (converged) ok[which.max(summary$objective[ok])] else {
    fin <- which(is.finite(summary$objective))
    if (length(fin)) fin[which.max(summary$objective[fin])] else NA_integer_
  }
  structure(list(
    params = if (is.na(pick)) NULL else runs[[pick]]$params,
    loglik = if (is.na(pick)) NA_real_ else runs[[pick]]$loglik,
    objective = if (is.na(pick)) NA_real_ else runs[[pick]]$objective,
    selected = if (is.na(pick)) NA_character_ else summary$start[pick],
    converged = converged,
    runs = summary,
    traces = lapply(runs, function(r) r$trace),
    spec = spec, n = nrow(X), control = control
  ), class = "csm_fit")
}

print.csm_fit <- function(x, ...) {
  cat("Correlated csmGmm fit\n")
  cat(sprintf("  K = %d, N = %d, largest support size %d, %d support-component pairs\n",
              x$spec$K, x$n, x$spec$S, length(x$spec$grp_b)))
  if (x$converged) {
    cat(sprintf("  selected start: %s, average log-likelihood %.6f\n", x$selected, x$loglik))
  } else {
    cat("  NO START CONVERGED: the analysis is recorded as a failure and makes no discoveries.\n")
  }
  print(x$runs, row.names = FALSE)
  invisible(x)
}

logLik.csm_fit <- function(object, ...) {
  sp <- object$spec
  G <- length(sp$grp_b)
  n_w <- G - 1L
  n_mu <- if (sp$magnitudes == "shared") sp$K else sum(sp$grp_A)
  n_sig <- if (is.null(.fixed_sigma(object$control, sp$K))) sp$K * (sp$K - 1L) / 2 else 0
  structure(object$loglik * object$n, df = n_w + n_mu + n_sig,
            nobs = object$n, class = "logLik")
}
