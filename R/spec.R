# Model specification, parameter validation, and signed-state enumeration.

.check_integer_scalar <- function(x, name, lower, upper) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x) &&
    x >= lower && x <= upper
  if (!ok) {
    stop(sprintf("'%s' must be an integer between %d and %d.", name, lower, upper),
         call. = FALSE)
  }
  as.integer(x)
}

.check_numeric_scalar <- function(x, name, lower = -Inf, upper = Inf) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= lower && x <= upper
  if (!ok) {
    stop(sprintf("'%s' must be a finite number in [%g, %g].", name, lower, upper),
         call. = FALSE)
  }
  as.numeric(x)
}

.support_table <- function(K) {
  labels <- seq_len(2^K) - 1L
  membership <- matrix(FALSE, nrow = length(labels), ncol = K)
  for (k in seq_len(K)) {
    membership[, k] <- bitwAnd(labels, bitwShiftL(1L, K - k)) != 0L
  }
  list(b = labels, A = membership, size = rowSums(membership))
}

.sign_patterns <- function(s) {
  if (s == 0L) return(matrix(0L, nrow = 1L, ncol = 0L))
  grid <- expand.grid(rep(list(c(1L, -1L)), s), KEEP.OUT.ATTRS = FALSE)
  unname(as.matrix(grid))
}

.block_starts <- function(counts, first = 1L) {
  if (length(counts) <= 1L) return(as.integer(first))
  as.integer(cumsum(c(first, counts[-length(counts)])))
}

.enumerate_components <- function(K, group_support) {
  sizes <- rowSums(group_support)
  per_group <- 2L^sizes
  first <- .block_starts(per_group)
  H <- matrix(0L, nrow = sum(per_group), ncol = K)
  for (g in seq_len(nrow(group_support))) {
    active <- which(group_support[g, ])
    if (length(active)) {
      rows <- first[g] + seq_len(per_group[g]) - 1L
      H[rows, active] <- .sign_patterns(length(active))
    }
  }
  group <- rep(seq_len(nrow(group_support)), per_group)
  list(H = H, group = group, size = sizes[group])
}

.pair_components <- function(components, n_groups) {
  counts <- tabulate(components$group, nbins = n_groups)
  first <- .block_starts(counts)
  keep <- logical(nrow(components$H))
  for (g in seq_len(n_groups)) {
    n_keep <- max(1L, counts[g] %/% 2L)
    keep[first[g] + seq_len(n_keep) - 1L] <- TRUE
  }
  list(
    H = components$H[keep, , drop = FALSE],
    group = components$group[keep],
    grp = components$group[keep],
    size = components$size[keep],
    self = components$size[keep] == 0L
  )
}

#' Number of magnitude components by binary support
#'
#' @param K Number of coordinates.
#' @param single Number of components for supports of size one.
#' @param multiple Number of components for supports of size two or larger.
#' @return Integer vector of length `2^K`; the first element is always one.
#' @export
csm_structure <- function(K, single = 2L, multiple = 1L) {
  K <- .check_integer_scalar(K, "K", 2L, 12L)
  single <- .check_integer_scalar(single, "single", 1L, .Machine$integer.max)
  multiple <- .check_integer_scalar(multiple, "multiple", 1L, .Machine$integer.max)
  supports <- .support_table(K)
  out <- rep.int(as.integer(multiple), 2^K)
  out[supports$size == 1L] <- as.integer(single)
  out[1L] <- 1L
  out
}

#' Specify a correlated csmGmm
#' @export
csm_spec <- function(K, M = csm_structure(K, 2L, 1L), S = K,
                     a_mu = 0.5, b_mu = NULL, kappa = 0.05) {
  K <- .check_integer_scalar(K, "K", 2L, 12L)
  S <- .check_integer_scalar(S, "S", 1L, K)
  n_supports <- 2^K
  if (!is.numeric(M) || anyNA(M) || any(!is.finite(M)) ||
      any(M != round(M)) || any(M < 1)) {
    stop("'M' must contain positive integers.", call. = FALSE)
  }
  if (length(M) == 1L) M <- rep.int(as.integer(M), n_supports)
  if (length(M) != n_supports) {
    stop("'M' must have length one or 2^K.", call. = FALSE)
  }
  M <- as.integer(M)
  M[1L] <- 1L
  a_mu <- .check_numeric_scalar(a_mu, "a_mu", 1e-8, Inf)
  if (!is.null(b_mu)) b_mu <- .check_numeric_scalar(b_mu, "b_mu", a_mu * (1 + 1e-8), Inf)
  kappa <- .check_numeric_scalar(kappa, "kappa", 1e-8, 1 - 1e-8)

  supports <- .support_table(K)
  kept_supports <- supports$b[supports$size <= S]
  group_b <- rep.int(kept_supports, M[kept_supports + 1L])
  group_m <- sequence(M[kept_supports + 1L])
  group_A <- supports$A[group_b + 1L, , drop = FALSE]
  group_size <- rowSums(group_A)
  components <- .enumerate_components(K, group_A)
  components$grp <- components$group
  pairs <- .pair_components(components, length(group_b))
  full <- which(group_b == n_supports - 1L)
  other <- lapply(seq_len(K), function(k) {
    which(group_b != 0L & group_b != n_supports - 1L & group_A[, k])
  })

  structure(list(
    K = K,
    S = S,
    M = M,
    magnitudes = "support",
    a_mu = a_mu,
    b_mu = if (is.null(b_mu)) NULL else as.numeric(b_mu),
    kappa = kappa,
    grp_b = group_b,
    grp_m = group_m,
    grp_A = group_A,
    grp_size = group_size,
    full = full,
    other = other,
    comp = components,
    pair = pairs
  ), class = "csm_spec")
}

.make_spec <- function(K, M, S, magnitudes = c("support", "shared"),
                       a_mu = 0.5, b_mu = NULL, kappa = 0.05) {
  magnitudes <- match.arg(magnitudes)
  spec <- csm_spec(K, M, S, a_mu, b_mu, kappa)
  spec$magnitudes <- magnitudes
  spec
}

#' @export
print.csm_spec <- function(x, ...) {
  cat("Correlated csmGmm specification\n")
  cat(sprintf("  K = %d, largest support size S = %d\n", x$K, x$S))
  cat(sprintf("  groups: %d; signed components: %d\n",
              length(x$grp_b), nrow(x$comp$H)))
  upper <- if (is.null(x$b_mu)) "infinity" else format(x$b_mu)
  cat(sprintf("  magnitude interval: [%g, %s]; kappa = %g\n",
              x$a_mu, upper, x$kappa))
  invisible(x)
}

.ordering_ok <- function(mu, spec, tolerance = 1e-10) {
  if (!length(spec$full)) return(TRUE)
  for (k in seq_len(spec$K)) {
    competitors <- spec$other[[k]]
    if (length(competitors) &&
        min(mu[spec$full, k]) < max(mu[competitors, k]) - tolerance) {
      return(FALSE)
    }
  }
  TRUE
}

.parameter_violations <- function(par, spec, tolerance = 1e-8) {
  G <- length(spec$grp_b)
  K <- spec$K
  problems <- character()
  if (!is.numeric(par$varpi) || length(par$varpi) != G ||
      any(!is.finite(par$varpi)) || any(par$varpi < -tolerance) ||
      abs(sum(par$varpi) - 1) > tolerance) {
    problems <- c(problems, "mixture weights must lie in the simplex")
  }
  if (!is.matrix(par$mu) || !identical(dim(par$mu), c(G, K)) ||
      any(!is.finite(par$mu))) {
    problems <- c(problems, "magnitude matrix has the wrong dimension")
  } else {
    inactive <- !spec$grp_A
    if (any(abs(par$mu[inactive]) > tolerance)) {
      problems <- c(problems, "magnitudes must be zero outside the support")
    }
    active <- spec$grp_A
    if (any(par$mu[active] < spec$a_mu - tolerance) ||
        (!is.null(spec$b_mu) && any(par$mu[active] > spec$b_mu + tolerance))) {
      problems <- c(problems, "active magnitudes fall outside the specified interval")
    }
    if (!.ordering_ok(par$mu, spec, tolerance)) {
      problems <- c(problems, "magnitude ordering is violated")
    }
  }
  Sigma <- par$Sigma
  if (!is.matrix(Sigma) || !identical(dim(Sigma), c(K, K)) ||
      any(!is.finite(Sigma)) || max(abs(Sigma - t(Sigma))) > 1e-10 ||
      max(abs(diag(Sigma) - 1)) > 1e-10) {
    problems <- c(problems, "Sigma must be a K by K correlation matrix")
  } else if (min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values) <
             spec$kappa - tolerance) {
    problems <- c(problems, "Sigma violates the eigenvalue lower bound")
  }
  problems
}

.feasibility <- function(par, spec, tolerance = 1e-8) {
  .parameter_violations(par, spec, tolerance)
}

.list_argument <- function(x, label, expected) {
  if (!is.list(x) || length(x) != expected) {
    stop(sprintf("'%s' must be a list of length %d.", label, expected), call. = FALSE)
  }
  x
}

#' Construct parameters of a correlated csmGmm
#' @export
csm_params <- function(spec, varpi = NULL, mu, Sigma = diag(spec$K), pi = NULL) {
  if (!inherits(spec, "csm_spec")) stop("'spec' must come from csm_spec().", call. = FALSE)
  if (is.null(varpi) == is.null(pi)) {
    stop("Give exactly one of 'varpi' and 'pi'.", call. = FALSE)
  }
  G <- length(spec$grp_b)
  K <- spec$K
  n_supports <- 2^K

  if (!is.null(pi)) {
    if (is.list(pi)) {
      pi <- .list_argument(pi, "pi", n_supports)
      varpi <- vapply(seq_len(G), function(g) {
        pi[[spec$grp_b[g] + 1L]][spec$grp_m[g]]
      }, numeric(1))
    } else {
      varpi <- pi
    }
    if (length(varpi) != G) stop("'pi' has the wrong length.", call. = FALSE)
    # csmGmm stores the total mass of a support-component block.  The
    # conditional sign symmetry is applied later when the block mass is
    # distributed over its 2^s signed configurations.
    varpi <- as.numeric(varpi)
  } else {
    if (is.list(varpi)) {
      varpi <- .list_argument(varpi, "varpi", n_supports)
      varpi <- vapply(seq_len(G), function(g) {
        varpi[[spec$grp_b[g] + 1L]][spec$grp_m[g]]
      }, numeric(1))
    }
    if (length(varpi) != G) stop("'varpi' has the wrong length.", call. = FALSE)
  }

  if (is.list(mu)) {
    mu <- .list_argument(mu, "mu", n_supports)
    magnitude <- matrix(0, nrow = G, ncol = K)
    for (g in seq_len(G)) {
      b <- spec$grp_b[g]
      if (b > 0L) {
        block <- as.matrix(mu[[b + 1L]])
        if (nrow(block) != K) {
          stop("Every non-null element of 'mu' must be a K by M_b matrix.", call. = FALSE)
        }
        magnitude[g, ] <- block[, spec$grp_m[g]] * spec$grp_A[g, ]
      }
    }
    mu <- magnitude
  } else if (is.numeric(mu) && !is.matrix(mu) && length(mu) == K) {
    mu <- spec$grp_A * rep(as.numeric(mu), each = G)
  }

  out <- list(varpi = as.numeric(varpi), mu = as.matrix(mu), Sigma = as.matrix(Sigma))
  violations <- .parameter_violations(out, spec)
  if (length(violations)) {
    stop(sprintf("Infeasible parameter: %s.", paste(violations, collapse = "; ")), call. = FALSE)
  }
  out$varpi <- pmax(out$varpi, 0)
  out$varpi <- out$varpi / sum(out$varpi)
  structure(out, class = "csm_params")
}

#' Convert parameters to the list representation used by csmGmm
#' @export
as_sun_lists <- function(object, spec = NULL) {
  if (inherits(object, "csm_fit")) {
    spec <- object$spec
    object <- object$params
  }
  if (is.null(spec) || !inherits(spec, "csm_spec")) {
    stop("'spec' is required.", call. = FALSE)
  }
  K <- spec$K
  n_supports <- 2^K
  mu_list <- vector("list", n_supports)
  pi_list <- vector("list", n_supports)
  for (b in 0:(n_supports - 1L)) {
    groups <- which(spec$grp_b == b)
    Mb <- spec$M[b + 1L]
    if (length(groups)) {
      mu_list[[b + 1L]] <- t(object$mu[groups, , drop = FALSE])
      pi_list[[b + 1L]] <- object$varpi[groups]
    } else {
      mu_list[[b + 1L]] <- matrix(0, nrow = K, ncol = Mb)
      pi_list[[b + 1L]] <- rep(0, Mb)
    }
  }
  list(muList = mu_list, piList = pi_list, Sigma = object$Sigma)
}
