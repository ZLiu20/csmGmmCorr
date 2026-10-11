# Model specification, parameter objects and component enumeration.
#
# Notation follows the manuscript. A support is indexed by its binary label
# b = sum_k 2^(K-k) |h_k| (coordinate 1 is the most significant bit), as in
# Sun, McCaw and Lin (2025). A "group" is a pair (b, m) of a support and one of
# its M_b magnitude components; a "component" is a signed configuration h on
# that support together with m. Parameters are stored at the group level:
#   varpi[g] = total probability of group g (= 2^|A_b| * pi_{bm}),
#   mu[g, ]  = magnitudes of group g (positive on A_b, zero elsewhere),
#   Sigma    = common correlation matrix.

.support_table <- function(K) {
  b <- seq_len(2L^K) - 1L
  A <- matrix(FALSE, length(b), K)
  for (k in seq_len(K)) A[, k] <- bitwAnd(b, 2L^(K - k)) > 0L
  list(b = b, A = A, size = rowSums(A))
}

.sign_patterns <- function(s) {
  if (s == 0L) return(matrix(0L, 1L, 0L))
  unname(as.matrix(expand.grid(rep(list(c(1L, -1L)), s))))
}

.components <- function(K, grp_A) {
  G <- nrow(grp_A)
  sizes <- rowSums(grp_A)
  n_r <- 2L^sizes
  H <- matrix(0L, sum(n_r), K)
  pats <- lapply(0:K, .sign_patterns)
  start <- c(0L, cumsum(n_r))[seq_len(G)]
  for (g in seq_len(G)) {
    s <- sizes[g]
    if (s > 0L) H[start[g] + seq_len(n_r[g]), which(grp_A[g, ])] <- pats[[s + 1L]]
  }
  grp <- rep(seq_len(G), n_r)
  list(H = H, grp = grp, size = sizes[grp])
}

# Sign pairs: the components h and -h of a group share their weight and the
# magnitude of their mean, so the E-step works with one representative per
# pair (the pattern whose last active coordinate is positive). Within each
# group, the sign patterns are ordered so that row j and row 2^s + 1 - j are
# negatives of each other. The null component is its own pair.
.pairs <- function(comp, G) {
  n_r <- tabulate(comp$grp, G)
  start <- match(seq_len(G), comp$grp)
  keep <- logical(nrow(comp$H))
  for (g in seq_len(G)) keep[start[g] + seq_len(max(1L, n_r[g] %/% 2L)) - 1L] <- TRUE
  list(H = comp$H[keep, , drop = FALSE], grp = comp$grp[keep], size = comp$size[keep],
       self = comp$size[keep] == 0L)
}

.check_scalar <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower || x > upper) {
    stop(sprintf("'%s' must be a finite number in [%g, %g].", name, lower, upper), call. = FALSE)
  }
  invisible(as.numeric(x))
}

# ---- Specify a correlated csmGmm ----
csm_spec <- function(K, M = csm_structure(K, 2L, 1L), S = K, a_mu = 0.5, b_mu = NULL, kappa = 0.05) {
  .make_spec(K, M, S, "support", a_mu, b_mu, kappa)
}

# Internal constructor. magnitudes = "shared" (one magnitude per coordinate,
# used by every support) is used only to compute a starting value for the
# correlated csmGmm; it is not offered as a model.
.make_spec <- function(K, M, S, magnitudes, a_mu, b_mu, kappa) {
  if (!is.numeric(K) || length(K) != 1L || K != round(K) || K < 2 || K > 12) {
    stop("'K' must be an integer between 2 and 12.", call. = FALSE)
  }
  K <- as.integer(K)
  if (!is.numeric(S) || length(S) != 1L || S != round(S) || S < 1 || S > K) {
    stop("'S' must be an integer between 1 and K.", call. = FALSE)
  }
  S <- as.integer(S)
  nb <- 2L^K
  if (!is.numeric(M) || anyNA(M) || any(M != round(M)) || any(M < 1)) {
    stop("'M' must contain positive integers.", call. = FALSE)
  }
  M <- as.integer(M)
  if (length(M) == 1L) M <- rep(M, nb)
  if (length(M) != nb) stop("'M' must have length 1 or 2^K.", call. = FALSE)
  M[1L] <- 1L
  if (magnitudes == "shared" && any(M != 1L)) stop("Shared magnitudes require M = 1.", call. = FALSE)
  .check_scalar(a_mu, "a_mu", 1e-8)
  if (!is.null(b_mu)) .check_scalar(b_mu, "b_mu", a_mu * (1 + 1e-8))
  .check_scalar(kappa, "kappa", 1e-8, 1 - 1e-8)

  sup <- .support_table(K)
  keep_b <- sup$b[sup$size <= S]
  grp_b <- rep(keep_b, M[keep_b + 1L])
  grp_m <- sequence(M[keep_b + 1L])
  grp_A <- sup$A[grp_b + 1L, , drop = FALSE]
  grp_size <- rowSums(grp_A)
  full <- which(grp_b == nb - 1L)
  other <- lapply(seq_len(K), function(k) which(grp_b != 0L & grp_b != nb - 1L & grp_A[, k]))
  comp <- .components(K, grp_A)
  pair <- .pairs(comp, length(grp_b))
  structure(list(
    K = K, S = S, M = M, magnitudes = magnitudes,
    a_mu = as.numeric(a_mu), b_mu = if (is.null(b_mu)) NULL else as.numeric(b_mu),
    kappa = as.numeric(kappa),
    grp_b = grp_b, grp_m = grp_m, grp_A = grp_A, grp_size = grp_size,
    full = full, other = other, comp = comp, pair = pair
  ), class = "csm_spec")
}

print.csm_spec <- function(x, ...) {
  cat("Correlated csmGmm specification\n")
  cat(sprintf("  K = %d, largest support size S = %d\n", x$K, x$S))
  cat(sprintf("  groups (support x magnitude component): %d, signed components: %d\n",
              length(x$grp_b), nrow(x$comp$H)))
  cat(sprintf("  region: a_mu = %g, b_mu = %s, kappa = %g\n", x$a_mu,
              if (is.null(x$b_mu)) "set from data at fit" else format(x$b_mu), x$kappa))
  invisible(x)
}

.ordering_ok <- function(MU, spec, tol = 1e-10) {
  if (length(spec$full) == 0L) return(TRUE)
  for (k in seq_len(spec$K)) {
    o <- spec$other[[k]]
    if (length(o) && min(MU[spec$full, k]) < max(MU[o, k]) - tol) return(FALSE)
  }
  TRUE
}

.feasibility <- function(p, spec, tol = 1e-8) {
  G <- length(spec$grp_b)
  K <- spec$K
  msg <- character(0)
  if (!is.numeric(p$varpi) || length(p$varpi) != G || any(!is.finite(p$varpi)) ||
      any(p$varpi < -tol) || abs(sum(p$varpi) - 1) > 1e-6) msg <- c(msg, "weights must lie in the simplex")
  if (!is.matrix(p$mu) || !identical(dim(p$mu), c(G, K)) || any(!is.finite(p$mu))) {
    msg <- c(msg, "magnitude matrix has the wrong dimension")
  } else {
    act <- spec$grp_A
    if (any(p$mu[!act] != 0)) msg <- c(msg, "magnitudes must be zero outside the support")
    if (any(p$mu[act] < spec$a_mu - tol) ||
        (!is.null(spec$b_mu) && any(p$mu[act] > spec$b_mu + tol))) msg <- c(msg, "magnitudes outside [a_mu, b_mu]")
    if (!.ordering_ok(p$mu, spec, tol)) msg <- c(msg, "the ordering constraint is violated")
    if (spec$magnitudes == "shared") {
      for (k in seq_len(K)) if (diff(range(p$mu[act[, k], k])) > tol) {
        msg <- c(msg, "magnitudes are not shared across supports"); break
      }
    }
  }
  S <- p$Sigma
  if (!is.matrix(S) || !identical(dim(S), c(K, K)) || any(!is.finite(S)) ||
      max(abs(S - t(S))) > 1e-10 || max(abs(diag(S) - 1)) > 1e-10) {
    msg <- c(msg, "Sigma must be a K x K correlation matrix")
  } else if (min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) < spec$kappa - 1e-8) {
    msg <- c(msg, "Sigma has an eigenvalue below kappa")
  }
  msg
}

# ---- Construct a parameter object ----
csm_params <- function(spec, varpi = NULL, mu, Sigma = diag(spec$K), pi = NULL) {
  if (!inherits(spec, "csm_spec")) stop("'spec' must come from csm_spec().", call. = FALSE)
  G <- length(spec$grp_b)
  K <- spec$K
  if (is.null(varpi) == is.null(pi)) stop("Give exactly one of 'varpi' and 'pi'.", call. = FALSE)
  from_list <- function(x, what) {
    # Sun et al. list format: element b + 1 holds the values for support b.
    if (length(x) != 2L^K) stop(sprintf("'%s' as a list must have 2^K elements.", what), call. = FALSE)
    x
  }
  if (!is.null(pi)) {
    if (is.list(pi)) {
      pi <- from_list(pi, "pi")
      varpi <- vapply(seq_len(G), function(g) pi[[spec$grp_b[g] + 1L]][spec$grp_m[g]], 0)
    } else varpi <- pi
    if (length(varpi) != G) stop("'pi' has the wrong length.", call. = FALSE)
    varpi <- varpi * 2^spec$grp_size
  } else if (is.list(varpi)) {
    varpi <- from_list(varpi, "varpi")
    varpi <- vapply(seq_len(G), function(g) varpi[[spec$grp_b[g] + 1L]][spec$grp_m[g]], 0)
  }
  if (is.list(mu)) {
    mu <- from_list(mu, "mu")
    MU <- matrix(0, G, K)
    for (g in seq_len(G)) {
      b <- spec$grp_b[g]
      if (b > 0L) {
        v <- as.matrix(mu[[b + 1L]])
        if (nrow(v) != K) stop("Each element of 'mu' must be a K x M_b matrix.", call. = FALSE)
        MU[g, ] <- v[, spec$grp_m[g]] * spec$grp_A[g, ]
      }
    }
    mu <- MU
  } else if (is.numeric(mu) && length(mu) == K && !is.matrix(mu)) {
    mu <- spec$grp_A * rep(mu, each = G)
  }
  p <- list(varpi = as.numeric(varpi), mu = mu, Sigma = Sigma)
  if (is.null(spec$b_mu)) spec$b_mu <- Inf
  bad <- .feasibility(p, spec)
  if (length(bad)) stop("Infeasible parameter: ", paste(bad, collapse = "; "), ".", call. = FALSE)
  p$varpi <- pmax(p$varpi, 0) / sum(pmax(p$varpi, 0))
  structure(p, class = "csm_params")
}

# ---- Convert parameters to the list format of the csmGmm package ----
as_sun_lists <- function(object, spec = NULL) {
  if (inherits(object, "csm_fit")) {
    spec <- object$spec
    object <- object$params
  }
  if (is.null(spec)) stop("'spec' is required.", call. = FALSE)
  K <- spec$K
  nb <- 2L^K
  sup <- .support_table(K)
  muList <- vector("list", nb)
  piList <- vector("list", nb)
  for (b in seq_len(nb) - 1L) {
    g <- which(spec$grp_b == b)
    Mb <- spec$M[b + 1L]
    if (length(g)) {
      muList[[b + 1L]] <- t(object$mu[g, , drop = FALSE])
      piList[[b + 1L]] <- object$varpi[g] / 2^sup$size[b + 1L]
    } else {
      muList[[b + 1L]] <- matrix(0, K, Mb)
      piList[[b + 1L]] <- rep(0, Mb)
    }
  }
  list(muList = muList, piList = piList, Sigma = object$Sigma)
}
