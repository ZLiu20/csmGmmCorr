# Gaussian mixture evaluation and target-specific local false discovery rates.

.check_z <- function(z, K, name = "z") {
  if (is.data.frame(z)) z <- as.matrix(z)
  if (is.numeric(z) && is.null(dim(z)) && length(z) == K) {
    z <- matrix(z, nrow = 1L)
  }
  if (!is.matrix(z) || !is.numeric(z) || ncol(z) != K) {
    stop(sprintf("'%s' must be a numeric matrix with %d columns.", name, K), call. = FALSE)
  }
  if (any(!is.finite(z))) {
    stop(sprintf("'%s' contains missing or infinite values.", name), call. = FALSE)
  }
  storage.mode(z) <- "double"
  z
}

.mixture_pieces <- function(spec, par) {
  pair <- spec$pair
  active <- which(par$varpi[pair$group] > 0)
  U <- chol(par$Sigma)
  precision <- chol2inv(U)
  means <- pair$H[active, , drop = FALSE] * par$mu[pair$group[active], , drop = FALSE]
  mean_precision <- means %*% precision
  log_weight <- log(par$varpi[pair$group[active]]) - pair$size[active] * log(2)
  quadratic <- rowSums(mean_precision * means)
  active_groups <- unique(pair$group[active])
  active_coordinates <- lapply(active_groups, function(g) which(spec$grp_A[g, ]))
  row_group <- rep(active_groups, lengths(active_coordinates))
  row_coordinate <- as.integer(unlist(active_coordinates, use.names = FALSE))
  group_sizes <- spec$grp_size[active_groups]
  group_counts <- ifelse(group_sizes == 0L, 1L, 2L^(group_sizes - 1L))
  group_starts <- .block_starts(group_counts, 0L)
  group_offsets <- .block_starts(lengths(active_coordinates), 0L)
  list(
    pair = pair,
    active = active,
    means = means,
    mean_precision = mean_precision,
    precision = precision,
    self = pair$self[active],
    row_group = row_group,
    row_coordinate = row_coordinate,
    group_start = as.integer(group_starts),
    group_size = as.integer(group_sizes),
    group_offset = as.integer(group_offsets),
    group_coordinate = as.integer(row_coordinate - 1L),
    group_mu = as.numeric(unlist(lapply(seq_along(active_groups), function(j) {
      par$mu[active_groups[j], active_coordinates[[j]]]
    }), use.names = FALSE)),
    log_constant = log_weight - 0.5 * quadratic,
    gaussian_constant = -0.5 * spec$K * log(2 * pi) - sum(log(diag(U)))
  )
}

.block_rows <- function(n_components, block_cells) {
  if (is.null(block_cells)) block_cells <- 5e6
  max(1L, as.integer(floor(block_cells / max(1L, n_components))))
}

.paired_block <- function(x, pieces) {
  linear <- tcrossprod(x, pieces$mean_precision)
  plus <- sweep(linear, 2L, pieces$log_constant, "+")
  minus <- sweep(-linear, 2L, pieces$log_constant, "+")
  if (any(pieces$self)) minus[, pieces$self] <- -Inf
  row_max <- pmax(apply(plus, 1L, max), apply(minus, 1L, max))
  exp_plus <- exp(plus - row_max)
  exp_minus <- exp(minus - row_max)
  denominator <- rowSums(exp_plus) + rowSums(exp_minus)
  row_constant <- -0.5 * rowSums((x %*% pieces$precision) * x) +
    pieces$gaussian_constant
  list(
    plus = exp_plus,
    minus = exp_minus,
    denominator = denominator,
    log_density = row_max + log(denominator) + row_constant
  )
}

.group_sums <- function(D, spec, pieces) {
  out <- matrix(0, nrow = length(pieces$row_group), ncol = spec$K)
  pair_h <- pieces$pair$H[pieces$active, , drop = FALSE]
  groups <- pieces$pair$group[pieces$active]
  for (r in seq_along(pieces$row_group)) {
    g <- pieces$row_group[r]
    k <- pieces$row_coordinate[r]
    idx <- which(groups == g & pair_h[, k] != 0)
    if (length(idx)) {
      out[r, ] <- colSums(D[idx, , drop = FALSE] * pair_h[idx, k])
    }
  }
  out
}

.estep_r <- function(x, spec, par, stats = TRUE, block_cells = 5e6) {
  x <- .check_z(x, spec$K, "x")
  pieces <- .mixture_pieces(spec, par)
  n <- nrow(x)
  n_pairs <- nrow(spec$pair$H)
  active_n <- length(pieces$active)
  block_n <- .block_rows(active_n, block_cells)
  pair_mass <- numeric(active_n)
  signed_first_moment <- matrix(0, nrow = active_n, ncol = spec$K)
  loglik <- 0

  for (from in seq.int(1L, n, by = block_n)) {
    idx <- from:min(n, from + block_n - 1L)
    xb <- x[idx, , drop = FALSE]
    b <- .paired_block(xb, pieces)
    loglik <- loglik + sum(b$log_density)
    if (stats) {
      normalized_plus <- b$plus / b$denominator
      normalized_minus <- b$minus / b$denominator
      pair_mass <- pair_mass + colSums(normalized_plus + normalized_minus)
      signed_first_moment <- signed_first_moment +
        crossprod(normalized_plus - normalized_minus, xb)
    }
  }

  out <- list(F = loglik / n, Omega = pieces$precision)
  if (stats) {
    out$T <- numeric(n_pairs)
    out$T[pieces$active] <- pair_mass / n
    out$D <- matrix(0, nrow = n_pairs, ncol = spec$K)
    out$D[pieces$active, ] <- signed_first_moment / n
    out$S <- .group_sums(out$D[pieces$active, , drop = FALSE], spec, pieces)
    out$Sg <- pieces$row_group
    out$Sk <- pieces$row_coordinate
  }
  out
}

# The compiled implementation is an optimization only. Keeping a single
# reference path here guarantees identical behavior when native code is absent.
.estep <- function(x, spec, par, stats = TRUE, block_cells = NULL, cores = 1L) {
  x <- .check_z(x, spec$K, "x")
  if (!is.loaded("csm_kernel", PACKAGE = "csmGmmCorr")) {
    return(.estep_r(x, spec, par, stats = stats,
                    block_cells = if (is.null(block_cells)) 5e6 else block_cells))
  }
  pieces <- .mixture_pieces(spec, par)
  row_constant <- -0.5 * rowSums((x %*% pieces$precision) * x) +
    pieces$gaussian_constant
  native <- .Call("csm_kernel", x, pieces$precision, pieces$log_constant,
                  pieces$self, row_constant, as.integer(stats), logical(0),
                  pieces$group_start, pieces$group_size, pieces$group_offset,
                  pieces$group_coordinate, pieces$group_mu,
                  as.integer(max(1L, cores)), PACKAGE = "csmGmmCorr")
  out <- list(F = sum(native[[1L]]) / nrow(x), Omega = pieces$precision)
  if (stats) {
    out$T <- numeric(nrow(spec$pair$H))
    out$T[pieces$active] <- native[[2L]] / nrow(x)
    out$D <- matrix(0, nrow = nrow(spec$pair$H), ncol = spec$K)
    out$S <- native[[3L]] / nrow(x)
    out$Sg <- pieces$row_group
    out$Sk <- pieces$row_coordinate
  }
  out
}

.check_pairs <- function(pairs, K) {
  if (is.null(pairs)) {
    if (K %% 2L) stop("Mediation with odd K requires 'pairs'.", call. = FALSE)
    return(matrix(seq_len(K), ncol = 2L, byrow = TRUE))
  }
  pairs <- as.matrix(pairs)
  ok <- ncol(pairs) == 2L && all(pairs %in% seq_len(K)) &&
    !anyDuplicated(as.vector(pairs))
  if (!ok) stop("'pairs' must be a two-column matrix of distinct coordinates.", call. = FALSE)
  storage.mode(pairs) <- "integer"
  pairs
}

#' Identify null signed configurations
#' @export
csm_null_components <- function(H,
                                target = c("pleiotropy", "replication", "mediation"),
                                t = 1L, pairs = NULL) {
  target <- match.arg(target)
  H <- as.matrix(H)
  K <- ncol(H)
  if (!is.numeric(t) || length(t) != 1L || !is.finite(t) ||
      t != round(t) || t < 1L || t >= K) {
    stop("'t' must be an integer satisfying 1 <= t < K.", call. = FALSE)
  }
  nz <- rowSums(H != 0)
  if (target == "pleiotropy") return(nz <= t)
  if (target == "replication") {
    same_sign_alt <- nz > t & abs(rowSums(H)) == nz
    return(!same_sign_alt)
  }
  pairs <- .check_pairs(pairs, K)
  complete <- rowSums(H[, pairs[, 1L], drop = FALSE] != 0 &
                      H[, pairs[, 2L], drop = FALSE] != 0)
  complete == 0L | nz <= t
}

.lfdr_core_r <- function(z, spec, par, null_pair, block_cells = 5e6) {
  pieces <- .mixture_pieces(spec, par)
  null_active <- as.logical(null_pair[pieces$active])
  if (!any(null_active)) return(numeric(nrow(z)))
  if (all(null_active)) return(rep(1, nrow(z)))
  out <- numeric(nrow(z))
  block_n <- .block_rows(length(pieces$active), block_cells)
  for (from in seq.int(1L, nrow(z), by = block_n)) {
    idx <- from:min(nrow(z), from + block_n - 1L)
    b <- .paired_block(z[idx, , drop = FALSE], pieces)
    numer <- rowSums((b$plus + b$minus)[, null_active, drop = FALSE])
    out[idx] <- numer / b$denominator
  }
  pmin(1, pmax(0, out))
}

.lfdr_core <- function(z, spec, par, null_pair, block_cells = 5e6, cores = 1L) {
  if (!is.loaded("csm_kernel", PACKAGE = "csmGmmCorr")) {
    return(.lfdr_core_r(z, spec, par, null_pair, block_cells))
  }
  pieces <- .mixture_pieces(spec, par)
  active_null <- as.logical(null_pair[pieces$active])
  if (!any(active_null)) return(numeric(nrow(z)))
  if (all(active_null)) return(rep(1, nrow(z)))
  row_constant <- -0.5 * rowSums((z %*% pieces$precision) * z) +
    pieces$gaussian_constant
  native <- .Call("csm_kernel", z, pieces$precision, pieces$log_constant,
                  pieces$self, row_constant, 2L, active_null,
                  pieces$group_start, pieces$group_size, pieces$group_offset,
                  pieces$group_coordinate, pieces$group_mu,
                  as.integer(max(1L, cores)), PACKAGE = "csmGmmCorr")
  pmin(1, pmax(0, native[[4L]]))
}

#' Local false discovery rates for a composite null
#' @export
csm_lfdr <- function(object, z,
                     target = c("pleiotropy", "replication", "mediation"),
                     t = 1L, pairs = NULL, spec = NULL,
                     block_cells = 5e6, cores = 1L) {
  target <- match.arg(target)
  if (inherits(object, "csm_fit")) {
    spec <- object$spec
    par <- object$params
  } else {
    par <- object
    if (!inherits(spec, "csm_spec")) {
      stop("'spec' is required when 'object' is a parameter object.", call. = FALSE)
    }
  }
  z <- .check_z(z, spec$K)
  null_pair <- csm_null_components(spec$pair$H, target, t, pairs)
  .lfdr_core(z, spec, par, null_pair, block_cells, cores)
}
