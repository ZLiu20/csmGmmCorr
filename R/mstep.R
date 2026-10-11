# Constrained conditional maximization for weights, magnitudes, and correlation.

.update_weights <- function(T, spec, pseudo = 0) {
  G <- length(spec$grp_b)
  group_mass <- numeric(G)
  totals <- rowsum(T, spec$pair$group, reorder = FALSE)
  group_mass[as.integer(rownames(totals))] <- totals[, 1L]
  group_mass <- pmax(group_mass, 0)

  if (pseudo > 0) {
    for (s in unique(spec$grp_size)) {
      idx <- which(spec$grp_size == s)
      total <- sum(group_mass[idx])
      if (total > 0) {
        group_mass[idx] <- total * (group_mass[idx] + pseudo) /
          (total + length(idx) * pseudo)
      }
    }
  }
  if (!is.finite(sum(group_mass)) || sum(group_mass) <= 0) {
    stop("The posterior component masses are not finite and positive.", call. = FALSE)
  }
  group_mass / sum(group_mass)
}

.weight_penalty <- function(varpi, spec) {
  value <- 0
  for (s in unique(spec$grp_size)) {
    idx <- which(spec$grp_size == s)
    total <- sum(varpi[idx])
    shares <- if (total > 0) varpi[idx] / total else rep(1 / length(idx), length(idx))
    value <- value - sum(log(pmax(shares, .Machine$double.xmin)))
  }
  value
}

.mu_stats <- function(es, spec) {
  G <- length(spec$grp_b)
  K <- spec$K
  C <- matrix(0, nrow = G, ncol = K)
  if (length(es$Sg)) {
    C[cbind(es$Sg, es$Sk)] <- rowSums(es$S * es$Omega[es$Sk, , drop = FALSE])
  }
  W <- array(0, dim = c(G, K, K))
  H <- spec$pair$H
  groups <- spec$pair$group
  for (k in seq_len(K)) {
    for (l in seq_len(K)) {
      values <- H[, k] * H[, l] * es$T
      sums <- rowsum(values, groups, reorder = FALSE)
      W[as.integer(rownames(sums)), k, l] <- sums[, 1L]
    }
  }
  Tg <- numeric(G)
  totals <- rowsum(es$T, groups, reorder = FALSE)
  Tg[as.integer(rownames(totals))] <- totals[, 1L]
  list(c = unname(C), C = unname(C), W = W, Tg = Tg)
}

.q_mu <- function(mu, stats, precision) {
  value <- sum(stats$C * mu)
  for (g in seq_len(nrow(mu))) {
    value <- value - 0.5 * sum(precision * stats$W[g, , ] * tcrossprod(mu[g, ]))
  }
  as.numeric(value)
}

.box_qp <- function(R, cvec, lower, upper, start, tolerance = 1e-12,
                    max_sweeps = 10000L) {
  unconstrained <- tryCatch(solve(R, cvec), error = function(e) NULL)
  if (!is.null(unconstrained) && all(is.finite(unconstrained)) &&
      all(unconstrained >= lower) && all(unconstrained <= upper)) {
    return(as.numeric(unconstrained))
  }
  x <- pmin(upper, pmax(lower, start))
  diagonal <- diag(R)
  diagonal[diagonal <= 0 | !is.finite(diagonal)] <- 1
  for (iteration in seq_len(max_sweeps)) {
    largest <- 0
    for (k in seq_along(x)) {
      candidate <- x[k] + (cvec[k] - sum(R[k, ] * x)) / diagonal[k]
      candidate <- min(upper[k], max(lower[k], candidate))
      largest <- max(largest, abs(candidate - x[k]))
      x[k] <- candidate
    }
    if (largest <= tolerance * (1 + max(abs(x)))) break
  }
  x
}

.ordering_threshold <- function(full_values, other_values) {
  lower <- min(full_values)
  upper <- max(other_values)
  if (lower >= upper) return(lower)
  balance <- function(x) {
    sum(pmax(x - full_values, 0)) - sum(pmax(other_values - x, 0))
  }
  for (iteration in seq_len(200L)) {
    midpoint <- (lower + upper) / 2
    if (balance(midpoint) > 0) upper <- midpoint else lower <- midpoint
    if (upper - lower <= 1e-14 * (1 + abs(midpoint))) break
  }
  (lower + upper) / 2
}

.project_mu <- function(mu, spec) {
  mu <- as.matrix(mu)
  mu[!spec$grp_A] <- 0
  upper <- if (is.null(spec$b_mu)) Inf else spec$b_mu
  mu[spec$grp_A] <- pmin(upper, pmax(spec$a_mu, mu[spec$grp_A]))
  if (!length(spec$full)) return(mu)

  for (k in seq_len(spec$K)) {
    other <- spec$other[[k]]
    if (!length(other)) next
    full_values <- mu[spec$full, k]
    other_values <- mu[other, k]
    if (min(full_values) >= max(other_values)) next
    threshold <- .ordering_threshold(full_values, other_values)
    threshold <- min(upper, max(spec$a_mu, threshold))
    mu[spec$full, k] <- pmax(full_values, threshold)
    mu[other, k] <- pmin(other_values, threshold)
  }
  mu
}

.mu_admm <- function(stats, spec, start, precision,
                     max_iterations = 3000L, tolerance = 1e-9) {
  groups <- which(spec$grp_b != 0L)
  active <- spec$grp_A
  sizes <- spec$grp_size[groups]
  blocks <- lapply(sort(unique(sizes)), function(size) {
    selected <- groups[sizes == size]
    index <- t(apply(active[selected, , drop = FALSE], 1L, which))
    if (size == 1L) index <- matrix(index, ncol = 1L)
    curvature <- array(0, dim = c(length(selected), size, size))
    for (k in seq_len(size)) {
      for (l in seq_len(size)) {
        curvature[, k, l] <-
          precision[cbind(index[, k], index[, l])] *
          stats$W[cbind(selected, index[, k], index[, l])]
      }
    }
    linear <- matrix(
      stats$C[cbind(rep(selected, size), as.vector(index))],
      nrow = length(selected), ncol = size
    )
    list(groups = selected, index = index, size = size,
         curvature = curvature, linear = linear)
  })

  diagonal <- unlist(lapply(blocks, function(block) {
    unlist(lapply(seq_len(block$size), function(k) {
      block$curvature[, k, k]
    }), use.names = FALSE)
  }), use.names = FALSE)
  positive <- diagonal[is.finite(diagonal) & diagonal > 0]
  rho <- if (length(positive)) stats::median(positive) else 1

  inverse_blocks <- function(value) {
    lapply(blocks, function(block) {
      inverse <- array(0, dim = dim(block$curvature))
      for (i in seq_along(block$groups)) {
        inverse[i, , ] <- solve(
          matrix(block$curvature[i, , ], block$size, block$size) +
            diag(value, block$size)
        )
      }
      inverse
    })
  }

  inverses <- inverse_blocks(rho)
  projected <- .project_mu(start, spec)
  dual <- matrix(0, nrow(projected), ncol(projected))
  primal <- projected
  scale <- sqrt(sum(active[groups, , drop = FALSE]))

  for (iteration in seq_len(max_iterations)) {
    for (j in seq_along(blocks)) {
      block <- blocks[[j]]
      positions <- cbind(rep(block$groups, block$size), as.vector(block$index))
      right <- block$linear + rho * matrix(
        (projected - dual)[positions],
        nrow = length(block$groups), ncol = block$size
      )
      solved <- matrix(0, nrow = length(block$groups), ncol = block$size)
      for (i in seq_along(block$groups)) {
        solved[i, ] <- inverses[[j]][i, , ] %*% right[i, ]
      }
      primal[positions] <- as.vector(solved)
    }

    previous <- projected
    projected <- .project_mu(primal + dual, spec)
    dual <- dual + primal - projected
    primal_residual <- sqrt(sum((primal - projected)^2))
    dual_residual <- rho * sqrt(sum((projected - previous)^2))
    if (primal_residual <= tolerance * scale &&
        dual_residual <= tolerance * scale) break

    if (iteration %% 25L == 0L) {
      if (primal_residual > 10 * dual_residual) {
        rho <- 2 * rho
        dual <- dual / 2
        inverses <- inverse_blocks(rho)
      } else if (dual_residual > 10 * primal_residual) {
        rho <- rho / 2
        dual <- dual * 2
        inverses <- inverse_blocks(rho)
      }
    }
  }
  projected
}

.update_mu <- function(stats, spec, start, precision) {
  if (identical(spec$magnitudes, "shared")) {
    K <- spec$K
    R <- matrix(0, K, K)
    cvec <- colSums(stats$C)
    for (g in which(spec$grp_b != 0L)) R <- R + precision * stats$W[g, , ]
    initial <- vapply(seq_len(K), function(k) {
      values <- start[spec$grp_A[, k], k]
      if (length(values)) mean(values) else spec$a_mu
    }, numeric(1))
    upper <- rep(if (is.null(spec$b_mu)) Inf else spec$b_mu, K)
    shared <- .box_qp(R, cvec, rep(spec$a_mu, K), upper, initial)
    return(spec$grp_A * rep(shared, each = nrow(spec$grp_A)))
  }

  current <- start
  upper <- if (is.null(spec$b_mu)) Inf else spec$b_mu
  for (g in which(spec$grp_b != 0L & stats$Tg > 0)) {
    active <- which(spec$grp_A[g, ])
    block_weight <- matrix(stats$W[g, active, active, drop = TRUE],
                           nrow = length(active), ncol = length(active))
    curvature <- precision[active, active, drop = FALSE] * block_weight
    current[g, active] <- .box_qp(
      curvature, stats$C[g, active], rep(spec$a_mu, length(active)),
      rep(upper, length(active)), start[g, active]
    )
  }
  if (!.ordering_ok(current, spec, tolerance = 0)) {
    current <- .mu_admm(stats, spec, current, precision)
  }
  updated_value <- .q_mu(current, stats, precision)
  initial_value <- .q_mu(start, stats, precision)
  if (!is.finite(updated_value) || updated_value < initial_value) start else current
}

.scatter <- function(es, mu, spec, XtX) {
  H <- spec$pair$H
  group <- spec$pair$group
  means <- H * mu[group, , drop = FALSE]
  selected_means <- matrix(0, nrow = length(es$Sg), ncol = spec$K)
  if (length(es$Sg)) {
    selected_means[cbind(seq_along(es$Sg), es$Sk)] <- mu[cbind(es$Sg, es$Sk)]
  }
  cross_term <- crossprod(es$S, selected_means)
  weighted_means <- means * sqrt(pmax(es$T, 0))
  result <- XtX - cross_term - t(cross_term) + crossprod(weighted_means)
  (result + t(result)) / 2
}

.project_cor <- function(A, kappa, max_iterations = 2000L, tolerance = 1e-13) {
  A <- (as.matrix(A) + t(as.matrix(A))) / 2
  K <- nrow(A)
  Y <- A
  diag(Y) <- 1
  correction <- matrix(0, K, K)
  for (iteration in seq_len(max_iterations)) {
    residual <- Y - correction
    eig <- eigen((residual + t(residual)) / 2, symmetric = TRUE)
    positive <- eig$vectors %*% (pmax(eig$values, kappa) * t(eig$vectors))
    correction <- positive - residual
    next_Y <- (positive + t(positive)) / 2
    diag(next_Y) <- 1
    change <- max(abs(next_Y - Y))
    Y <- next_Y
    if (change <= tolerance) break
  }
  smallest <- min(eigen(Y, symmetric = TRUE, only.values = TRUE)$values)
  if (smallest < kappa) {
    fraction <- (kappa - smallest) / (1 - smallest)
    Y <- (1 - fraction) * Y + fraction * diag(K)
  }
  diag(Y) <- 1
  (Y + t(Y)) / 2
}

.cor_objective <- function(Sigma, scatter) {
  factor <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(factor)) return(Inf)
  precision <- chol2inv(factor)
  2 * sum(log(diag(factor))) + sum(precision * scatter)
}

.update_sigma <- function(scatter, start, kappa, max_iterations = 50L) {
  current <- .project_cor(start, kappa)
  value <- .cor_objective(current, scatter)
  for (iteration in seq_len(max_iterations)) {
    precision <- solve(current)
    gradient <- precision - precision %*% scatter %*% precision
    gradient <- (gradient + t(gradient)) / 2
    diag(gradient) <- 0
    step <- 1
    accepted <- FALSE
    for (trial in seq_len(40L)) {
      proposal <- .project_cor(current - step * gradient, kappa)
      candidate <- .cor_objective(proposal, scatter)
      if (is.finite(candidate) && candidate < value - 1e-4 * step *
          sum(gradient * (current - proposal))) {
        accepted <- TRUE
        break
      }
      step <- step / 2
    }
    if (!accepted) break
    change <- max(abs(proposal - current))
    current <- proposal
    value <- candidate
    if (change <= 1e-9) break
  }
  current
}

.size_means <- function(mu, spec) {
  sizes <- sort(unique(spec$grp_size[spec$grp_b != 0L]))
  result <- matrix(NA_real_, nrow = max(spec$K, max(sizes)), ncol = spec$K)
  for (s in sizes) {
    groups <- which(spec$grp_b != 0L & spec$grp_size == s)
    for (k in seq_len(spec$K)) {
      active <- groups[spec$grp_A[groups, k]]
      if (length(active)) result[s, k] <- mean(mu[active, k])
    }
  }
  result
}

.magnitude_penalty <- function(mu, spec) {
  center <- .size_means(mu, spec)
  groups <- which(spec$grp_b != 0L)
  residual <- mu[groups, , drop = FALSE] - center[spec$grp_size[groups], , drop = FALSE]
  0.5 * sum((residual^2)[spec$grp_A[groups, , drop = FALSE]])
}

.penalize_mu_stats <- function(stats, mu, spec, precision, amount) {
  center <- .size_means(mu, spec)
  for (k in seq_len(spec$K)) {
    groups <- which(spec$grp_b != 0L & spec$grp_A[, k])
    if (!length(groups)) next
    stats$W[groups, k, k] <- stats$W[groups, k, k] + amount / precision[k, k]
    stats$C[groups, k] <- stats$C[groups, k] +
      amount * center[spec$grp_size[groups], k]
  }
  stats
}
