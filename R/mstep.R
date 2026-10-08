# Conditional maximization steps of the constrained ECM algorithm
# (Supplementary Section S2). With responsibilities held fixed, the weight,
# magnitude and correlation blocks are updated in turn; each block does not
# decrease the expected complete-data log-likelihood, so the observed
# log-likelihood does not decrease (Meng and Rubin, 1993).

# ---- Weights: closed form on the closed simplex ----------------------------
.update_weights <- function(Tc, spec, a = 0) {
  cg <- as.vector(rowsum(Tc, spec$pair$grp, reorder = TRUE))
  cg <- pmax(cg, 0)
  cg <- cg / sum(cg)
  if (a > 0) {
    # Within-size pseudo-counts: with W_s the total weight of support size s
    # and n_s its number of support-component pairs, maximize
    #   sum_g c_g log w_g + a * sum_g log(w_g / W_s)
    # (a = alpha / N). The size totals keep their closed-form value W_s = C_s,
    # and each share within a size becomes (c_g + a) / (C_s + n_s a).
    sz <- spec$grp_size
    Cs <- tapply(cg, sz, sum)[as.character(sz)]
    ns <- tapply(rep(1, length(sz)), sz, sum)[as.character(sz)]
    cg <- as.vector(Cs * (cg + a) / (Cs + ns * a))
  }
  cg
}

# Penalty on the weights: -sum_g log(share of g within its support size).
.weight_penalty <- function(varpi, spec) {
  sz <- spec$grp_size
  Ws <- as.vector(tapply(varpi, sz, sum)[as.character(sz)])
  ns <- as.vector(tapply(rep(1, length(sz)), sz, sum)[as.character(sz)])
  # A support size whose total weight is exactly zero (its posterior mass can
  # underflow for large K) has shares 1/n_s, the limit of the weight update
  # (c_g + a) / (C_s + n_s a) as C_s -> 0.
  share <- ifelse(Ws > 0, varpi / Ws, 1 / ns)
  -sum(log(share))
}

# ---- Magnitudes -------------------------------------------------------------
# Sufficient statistics of the magnitude objective
#   Q(mu) = sum_g { c_g' mu_g - 0.5 mu_g' (Omega o W_g) mu_g },
# with W_g = sum_{r in g} T_r h_r h_r' and c_g = sum_{r in g} h_r o (Omega B_r).
# Summing over sign pairs, W_g = sum_pairs T_pair h h' and
# c_g = sum_pairs h o (Omega D_pair), with D_pair = B_+ - B_-.
.mu_stats <- function(es, spec) {
  cp <- spec$pair
  K <- spec$K
  G <- length(spec$grp_b)
  # c_g[k] = sum_{p in g} H_pk (Omega D_p)_k = sum_l S_g[k, l] Omega[l, k]
  cvec <- matrix(0, G, K)
  if (length(es$Sg)) cvec[cbind(es$Sg, es$Sk)] <- rowSums(es$S * es$Omega[es$Sk, , drop = FALSE])
  W <- array(0, c(G, K, K))
  for (k in seq_len(K)) for (l in k:K) {
    w <- as.vector(rowsum(cp$H[, k] * cp$H[, l] * es$T, cp$grp, reorder = TRUE))
    W[, k, l] <- w
    W[, l, k] <- w
  }
  list(c = unname(cvec), W = W, Tg = as.vector(rowsum(es$T, cp$grp, reorder = TRUE)))
}

.Q_mu <- function(MU, st, Omega) {
  K <- ncol(MU)
  val <- sum(st$c * MU)
  for (k in seq_len(K)) for (l in seq_len(K)) {
    val <- val - 0.5 * Omega[k, l] * sum(st$W[, k, l] * MU[, k] * MU[, l])
  }
  val
}

# Minimize 0.5 x'Rx - c'x over lo <= x <= hi (R positive definite) by exact
# cyclic coordinate descent, after trying the unconstrained solution.
.box_qp <- function(R, cvec, lo, hi, x0, tol = 1e-12, max_sweeps = 10000L) {
  x <- tryCatch(solve(R, cvec), error = function(e) NULL)
  if (!is.null(x) && all(is.finite(x)) && all(x >= lo) && all(x <= hi)) return(as.vector(x))
  x <- pmin(hi, pmax(lo, x0))
  d <- diag(R)
  for (it in seq_len(max_sweeps)) {
    delta <- 0
    for (k in seq_along(x)) {
      if (d[k] <= 0) next
      new <- min(hi[k], max(lo[k], x[k] + (cvec[k] - sum(R[k, ] * x)) / d[k]))
      delta <- max(delta, abs(new - x[k]))
      x[k] <- new
    }
    if (delta <= tol * (1 + max(abs(x)))) break
  }
  x
}

# Euclidean projection of the magnitudes onto the box and the ordering
# constraint. For each coordinate k the constraint reads
#   min_m mu[full m, k] >= max_{other g} mu[g, k].
# For a threshold theta in [a_mu, b_mu], the nearest feasible point clips all
# values to the box, raises full-support values to at least theta and lowers
# the others to at most theta. The optimal theta sets the derivative
#   sum_i (theta - u_i)_+ - sum_j (v_j - theta)_+
# to zero, computed with the unclipped values u (full support) and v (others),
# and is then restricted to [a_mu, b_mu].
.threshold <- function(u, v) {
  lo <- min(u)
  hi <- max(v)
  dfun <- function(th) sum(pmax(th - u, 0)) - sum(pmax(v - th, 0))
  for (i in 1:200) {
    mid <- 0.5 * (lo + hi)
    if (dfun(mid) > 0) hi <- mid else lo <- mid
    if (hi - lo <= 1e-15 * (1 + abs(mid))) break
  }
  0.5 * (lo + hi)
}

.project_mu <- function(MU, spec) {
  act <- spec$grp_A
  MU[!act] <- 0
  raw <- MU
  MU[act] <- pmin(spec$b_mu, pmax(spec$a_mu, MU[act]))
  if (length(spec$full) == 0L) return(MU)
  for (k in seq_len(spec$K)) {
    o <- spec$other[[k]]
    if (!length(o)) next
    uc <- MU[spec$full, k]
    vc <- MU[o, k]
    if (min(uc) >= max(vc)) next
    th <- .threshold(raw[spec$full, k], raw[o, k])
    th <- min(spec$b_mu, max(spec$a_mu, th))
    MU[spec$full, k] <- pmax(uc, th)
    MU[o, k] <- pmin(vc, th)
  }
  MU
}

# ADMM for the magnitude QP when the blockwise box solutions violate the
# ordering constraint (Boyd et al., 2011). Splitting x = w with w in the
# feasible set; the x-update solves (R_g + rho I) x_g = c_g + rho (w_g - u_g)
# for all groups at once, vectorized over groups of equal support size.
.mu_admm <- function(st, spec, MU0, Omega, maxit = 3000L, tol = 1e-9) {
  gs <- which(spec$grp_b != 0L)
  A <- spec$grp_A
  sizes <- spec$grp_size[gs]
  blocks <- lapply(sort(unique(sizes)), function(s) {
    g <- gs[sizes == s]
    idx <- t(apply(A[g, , drop = FALSE], 1L, which))
    if (s == 1L) idx <- matrix(idx, ncol = 1L)
    R <- array(0, c(length(g), s, s))
    for (k in seq_len(s)) for (l in seq_len(s)) {
      R[, k, l] <- Omega[cbind(idx[, k], idx[, l])] * st$W[cbind(g, idx[, k], idx[, l])]
    }
    list(g = g, idx = idx, s = s, R = R,
         cvec = matrix(st$c[cbind(rep(g, s), as.vector(idx))], length(g)))
  })
  invert <- function(rho) lapply(blocks, function(b) {
    Inv <- array(0, dim(b$R))
    for (i in seq_along(b$g)) Inv[i, , ] <- solve(matrix(b$R[i, , ], b$s) + diag(rho, b$s))
    Inv
  })
  dg <- unlist(lapply(blocks, function(b) {
    vapply(seq_len(b$s), function(k) b$R[, k, k], numeric(length(b$g)))
  }))
  rho <- stats::median(dg[dg > 0])
  if (!is.finite(rho)) rho <- 1
  Inv <- invert(rho)
  w <- .project_mu(MU0, spec)
  u <- matrix(0, nrow(MU0), ncol(MU0))
  x <- w
  scale <- sqrt(sum(A[gs, ]))
  for (it in seq_len(maxit)) {
    for (j in seq_along(blocks)) {
      b <- blocks[[j]]
      pos <- cbind(rep(b$g, b$s), as.vector(b$idx))
      rhs <- b$cvec + rho * matrix(w[pos] - u[pos], length(b$g))
      xb <- matrix(0, length(b$g), b$s)
      for (k in seq_len(b$s)) for (l in seq_len(b$s)) xb[, k] <- xb[, k] + Inv[[j]][, k, l] * rhs[, l]
      x[pos] <- as.vector(xb)
    }
    w_old <- w
    w <- .project_mu(x + u, spec)
    u <- u + x - w
    r_pri <- sqrt(sum((x - w)^2))
    r_dual <- rho * sqrt(sum((w - w_old)^2))
    if (r_pri <= tol * scale && r_dual <= tol * scale) break
    if (it %% 25L == 0L) {
      if (r_pri > 10 * r_dual) {
        rho <- 2 * rho; u <- u / 2; Inv <- invert(rho)
      } else if (r_dual > 10 * r_pri) {
        rho <- rho / 2; u <- u * 2; Inv <- invert(rho)
      }
    }
  }
  w
}

.update_mu <- function(st, spec, MU0, Omega) {
  G <- nrow(MU0)
  K <- spec$K
  A <- spec$grp_A
  lo <- spec$a_mu
  hi <- spec$b_mu
  nn <- spec$grp_b != 0L
  if (spec$magnitudes == "shared") {
    Wt <- matrix(0, K, K)
    for (k in seq_len(K)) for (l in seq_len(K)) Wt[k, l] <- sum(st$W[nn, k, l])
    mu0 <- vapply(seq_len(K), function(k) max(MU0[nn & A[, k], k]), 0)
    mu <- .box_qp(Omega * Wt, colSums(st$c[nn, , drop = FALSE]), rep(lo, K), rep(hi, K), mu0)
    MU <- A * rep(mu, each = G)
  } else {
    MU <- MU0
    for (g in which(nn & st$Tg > 0)) {
      a <- which(A[g, ])
      Rg <- Omega[a, a, drop = FALSE] * matrix(st$W[g, a, a], length(a))
      MU[g, a] <- .box_qp(Rg, st$c[g, a], rep(lo, length(a)), rep(hi, length(a)), MU0[g, a])
    }
    if (!.ordering_ok(MU, spec, tol = 0)) MU <- .mu_admm(st, spec, MU, Omega)
  }
  q_new <- .Q_mu(MU, st, Omega)
  q_old <- .Q_mu(MU0, st, Omega)
  if (!is.finite(q_new) || q_new < q_old) MU0 else MU
}

# ---- Correlation ------------------------------------------------------------
# Conditional scatter S = N^{-1} sum_{i,r} tau_ir (x_i - m_r)(x_i - m_r)'.
# Over a sign pair, sum_r B_r m_r' = D m_+' and sum_r T_r m_r m_r' = T m_+ m_+'.
.scatter <- function(es, MU, spec, XtX) {
  cp <- spec$pair
  Mean <- cp$H * MU[cp$grp, , drop = FALSE]
  # sum_p D_p Mean_p' = sum over group-sum rows r of S_r MU[g_r, k_r] e_{k_r}'
  Mr <- matrix(0, length(es$Sg), ncol(MU))
  if (length(es$Sg)) Mr[cbind(seq_along(es$Sg), es$Sk)] <- MU[cbind(es$Sg, es$Sk)]
  BM <- crossprod(es$S, Mr)
  S <- XtX - BM - t(BM) + crossprod(Mean * sqrt(es$T))
  (S + t(S)) / 2
}

# Frobenius projection onto {diag(Sigma) = 1, Sigma >= kappa I} by Dykstra's
# alternating projections, followed by an exact feasibility restoration.
.project_cor <- function(A, kappa, maxit = 2000L, tol = 1e-13) {
  K <- nrow(A)
  Y <- (A + t(A)) / 2
  diag(Y) <- 1
  P <- matrix(0, K, K)
  for (it in seq_len(maxit)) {
    Rm <- Y - P
    e <- eigen(Rm, symmetric = TRUE)
    X <- e$vectors %*% (pmax(e$values, kappa) * t(e$vectors))
    P <- X - Rm
    Yn <- (X + t(X)) / 2
    diag(Yn) <- 1
    done <- max(abs(Yn - Y)) < tol
    Y <- Yn
    if (done) break
  }
  lmin <- min(eigen(Y, symmetric = TRUE, only.values = TRUE)$values)
  if (lmin < kappa) {
    beta <- (kappa - lmin) / (1 - lmin)
    Y <- (1 - beta) * Y + beta * diag(K)
  }
  diag(Y) <- 1
  Y
}

.neg_cond_loglik <- function(Sigma, S) {
  U <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(U)) return(Inf)
  2 * sum(log(diag(U))) + sum(chol2inv(U) * S)
}

# Decrease g(Sigma) = log|Sigma| + tr(Sigma^{-1} S) over the feasible set by
# projected gradient steps with Armijo backtracking, starting from the better of
# the current value and the projected normalized scatter.
.update_sigma <- function(S, Sigma0, kappa, maxit = 50L) {
  cur <- Sigma0
  gcur <- .neg_cond_loglik(cur, S)
  d <- sqrt(pmax(diag(S), 1e-12))
  prop <- .project_cor(S / outer(d, d), kappa)
  gp <- .neg_cond_loglik(prop, S)
  if (gp < gcur) {
    cur <- prop
    gcur <- gp
  }
  for (it in seq_len(maxit)) {
    Om <- solve(cur)
    Gr <- Om - Om %*% S %*% Om
    Gr <- (Gr + t(Gr)) / 2
    diag(Gr) <- 0
    eta <- 1
    moved <- FALSE
    while (eta > 1e-10) {
      cand <- .project_cor(cur - eta * Gr, kappa)
      gc <- .neg_cond_loglik(cand, S)
      if (gc <= gcur - 1e-4 * sum(Gr * (cur - cand)) && gc < gcur) {
        moved <- TRUE
        break
      }
      eta <- eta / 2
    }
    if (!moved) break
    change <- max(abs(cand - cur))
    cur <- cand
    gcur <- gc
    if (change < 1e-9) break
  }
  if (.neg_cond_loglik(cur, S) <= .neg_cond_loglik(Sigma0, S)) cur else Sigma0
}

# ---- Magnitude penalty -------------------------------------------------------
# Penalized maximum likelihood shrinks each support's magnitudes toward the
# mean magnitude of the supports of the same size:
#   P(mu) = 0.5 * sum_g sum_{k in A_g} (mu_gk - mbar_{|A_g|, k})^2,
# with mbar the size-specific mean over supports containing k. The objective
# is F(zeta) - (lambda / N) P(mu). P is bounded on the magnitude box, so the
# penalized maximizer is within O(lambda / N) of the maximum likelihood.
.size_means <- function(MU, spec) {
  sizes <- sort(unique(spec$grp_size[spec$grp_b != 0L]))
  m <- matrix(NA_real_, max(spec$K, max(sizes)), spec$K)
  for (s in sizes) for (k in seq_len(spec$K)) {
    g <- which(spec$grp_b != 0L & spec$grp_size == s & spec$grp_A[, k])
    if (length(g)) m[s, k] <- mean(MU[g, k])
  }
  m
}

.mag_penalty <- function(MU, spec) {
  m <- .size_means(MU, spec)
  g <- which(spec$grp_b != 0L)
  A <- spec$grp_A[g, , drop = FALSE]
  target <- m[spec$grp_size[g], , drop = FALSE]
  0.5 * sum(((MU[g, , drop = FALSE] - target)^2)[A])
}

# Minorize-maximize step: with mbar fixed at the current magnitudes,
# -(lambda/N) * 0.5 * sum (mu - mbar)^2 is a lower bound of the penalty term
# that is tight at the current point; adding it to the quadratic magnitude
# objective changes W (diagonal) and c only, so all solvers apply unchanged.
.penalize_mu_stats <- function(st, MU, spec, Omega, eps) {
  m <- .size_means(MU, spec)
  for (k in seq_len(spec$K)) {
    g <- which(spec$grp_b != 0L & spec$grp_A[, k])
    if (!length(g)) next
    st$W[g, k, k] <- st$W[g, k, k] + eps / Omega[k, k]
    st$c[g, k] <- st$c[g, k] + eps * m[spec$grp_size[g], k]
  }
  st
}
