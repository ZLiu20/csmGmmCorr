# Gaussian mixture evaluation: E-step sufficient statistics, log-likelihood,
# and local false discovery rates. Components are handled in sign pairs
# (h, -h), which share a weight and the quadratic term of the log density;
# only their linear terms differ in sign. Features are processed in blocks so
# that each working matrix has at most `block_cells` entries.

.mixture_pieces <- function(spec, p) {
  pr <- spec$pair
  K <- spec$K
  U <- chol(p$Sigma)
  Omega <- chol2inv(U)
  act <- which(p$varpi[pr$grp] > 0)
  Mean <- pr$H[act, , drop = FALSE] * p$mu[pr$grp[act], , drop = FALSE]
  MO <- Mean %*% Omega
  # Group layout of the active pairs for the compiled kernel (0-based).
  ag <- pr$grp[act]
  groups <- unique(ag)
  gsize <- spec$grp_size[groups]
  coords <- lapply(groups, function(g) which(spec$grp_A[g, ]))
  list(
    act = act, Mean = Mean, MO = MO, Omega = Omega, self = pr$self[act],
    colc = log(p$varpi[pr$grp[act]]) - pr$size[act] * log(2) - 0.5 * rowSums(MO * Mean),
    const = -0.5 * K * log(2 * pi) - sum(log(diag(U))),
    gstart = as.integer(match(groups, ag) - 1L), gsize = as.integer(gsize),
    goff = as.integer(c(0L, cumsum(gsize))[seq_along(groups)]),
    gcoord = as.integer(unlist(coords) - 1L),
    gmu = as.numeric(unlist(lapply(seq_along(groups), function(j) p$mu[groups[j], coords[[j]]]))),
    # rows of the group sums S: group (spec index) and coordinate of each row
    srow_g = rep(groups, gsize), srow_k = as.integer(unlist(coords))
  )
}

# Group sums of the D statistics, S[r, ] = sum_{pairs p in group g_r} H[p, k_r] D[p, ],
# for the active coordinates k_r of each group: the only form of D the M-step
# uses (.mu_stats and .scatter). `Da` holds D for the active pairs.
.group_sums <- function(Da, spec, mp) {
  Ha <- spec$pair$H[mp$act, , drop = FALSE]
  nz <- which(Ha != 0, arr.ind = TRUE)
  key <- match(paste(spec$pair$grp[mp$act][nz[, 1]], nz[, 2]), paste(mp$srow_g, mp$srow_k))
  S <- matrix(0, length(mp$srow_g), spec$K)
  if (nrow(nz)) {
    rs <- rowsum(Ha[nz] * Da[nz[, 1], , drop = FALSE], key, reorder = TRUE)
    S[as.integer(rownames(rs)), ] <- rs
  }
  S
}

.block_rows <- function(n_comp, block_cells) max(1L, as.integer(floor(block_cells / max(1L, n_comp))))

# For rows `idx`: scaled exponentials of the log densities (plus log weights)
# of the "+" and "-" members of each pair, their row sum and the row
# log-likelihood. The "-" member of the null pair is absent (weight zero).
.block_logdens <- function(X, idx, mp) {
  Xi <- X[idx, , drop = FALSE]
  n <- length(idx)
  Gm <- tcrossprod(Xi, mp$MO)
  cc <- rep(mp$colc, each = n)
  Lp <- cc + Gm
  Lm <- cc - Gm
  if (any(mp$self)) Lm[, mp$self] <- -Inf
  mx <- pmax(Lp[cbind(seq_len(n), max.col(Lp, ties.method = "first"))],
             Lm[cbind(seq_len(n), max.col(Lm, ties.method = "first"))])
  Ep <- exp(Lp - mx)
  Em <- exp(Lm - mx)
  s <- rowSums(Ep) + rowSums(Em)
  rowc <- -0.5 * rowSums((Xi %*% mp$Omega) * Xi) + mp$const
  list(Ep = Ep, Em = Em, s = s, logf = mx + log(s) + rowc, Xi = Xi)
}

# Returns the average log-likelihood F and, per sign pair, T (sum of the two
# average responsibilities) and D (difference of the two responsibility-
# weighted means, B_+ - B_-); the M-step uses T and the group sums S of D.
.estep_r <- function(X, spec, p, stats = TRUE, block_cells = 5e6) {
  N <- nrow(X)
  K <- spec$K
  P <- nrow(spec$pair$H)
  mp <- .mixture_pieces(spec, p)
  Pa <- length(mp$act)
  nb <- .block_rows(Pa, block_cells)
  Ta <- numeric(Pa)
  Da <- matrix(0, Pa, K)
  ll <- 0
  for (s0 in seq(1L, N, by = nb)) {
    idx <- s0:min(N, s0 + nb - 1L)
    bl <- .block_logdens(X, idx, mp)
    ll <- ll + sum(bl$logf)
    if (stats) {
      Ta <- Ta + as.vector(crossprod(bl$Ep + bl$Em, 1 / bl$s))
      Da <- Da + crossprod(bl$Ep - bl$Em, bl$Xi / bl$s)
    }
  }
  out <- list(F = ll / N, Omega = mp$Omega)
  if (stats) {
    out$T <- numeric(P)
    out$T[mp$act] <- Ta / N
    out$D <- matrix(0, P, K)
    out$D[mp$act, ] <- Da / N
    out$S <- .group_sums(Da / N, spec, mp)
    out$Sg <- mp$srow_g
    out$Sk <- mp$srow_k
  }
  out
}

# ---- Compiled kernel (src/estep.c) -----------------------------------------
# The functions above (.estep_r, .lfdr_core_r) are the plain-R reference
# implementations used in the tests; the fitting code calls the kernel below.
# `cores` threads of the compiled kernel each process a block of features.
.kernel <- function(X, mp, mode, nulls = logical(0), cores = 1L) {
  rowc <- -0.5 * rowSums((X %*% mp$Omega) * X) + mp$const
  .Call(csm_kernel, X, mp$Omega, mp$colc, mp$self, rowc, as.integer(mode), as.logical(nulls),
        mp$gstart, mp$gsize, mp$goff, mp$gcoord, mp$gmu, as.integer(max(1L, cores)))
}

.estep <- function(X, spec, p, stats = TRUE, block_cells = NULL, cores = 1L) {
  N <- nrow(X)
  P <- nrow(spec$pair$H)
  mp <- .mixture_pieces(spec, p)
  r <- .kernel(X, mp, if (stats) 1L else 0L, cores = cores)
  out <- list(F = sum(r[[1]]) / N, Omega = mp$Omega)
  if (stats) {
    out$T <- numeric(P)
    out$T[mp$act] <- r[[2]] / N
    # The kernel returns the group sums S directly (see .group_sums).
    out$S <- r[[3]] / N
    out$Sg <- mp$srow_g
    out$Sk <- mp$srow_k
  }
  out
}

.lfdr_core <- function(Z, spec, p, null_pair, block_cells = NULL, cores = 1L) {
  mp <- .mixture_pieces(spec, p)
  nulls <- null_pair[mp$act]
  if (!any(nulls)) return(numeric(nrow(Z)))
  if (all(nulls)) return(rep(1, nrow(Z)))
  pmin(1, pmax(0, .kernel(Z, mp, 2L, nulls, cores)[[4]]))
}

.check_z <- function(Z, K, name = "Z") {
  if (is.data.frame(Z)) Z <- as.matrix(Z)
  if (is.numeric(Z) && is.null(dim(Z)) && length(Z) == K) Z <- matrix(Z, 1L)
  if (!is.matrix(Z) || !is.numeric(Z) || ncol(Z) != K) {
    stop(sprintf("'%s' must be a numeric matrix with %d columns.", name, K), call. = FALSE)
  }
  if (any(!is.finite(Z))) stop(sprintf("'%s' contains missing or infinite values.", name), call. = FALSE)
  storage.mode(Z) <- "double"
  Z
}

# ---- Null configurations for a target ----
csm_null_components <- function(H, target = c("pleiotropy", "replication", "mediation"),
                                t = 1L, pairs = NULL) {
  target <- match.arg(target)
  H <- as.matrix(H)
  K <- ncol(H)
  if (!is.numeric(t) || length(t) != 1L || t != round(t) || t < 1 || t >= K) {
    stop("'t' must be an integer with 1 <= t < K.", call. = FALSE)
  }
  nz <- rowSums(H != 0)
  switch(target,
    pleiotropy = nz <= t,
    replication = !(nz > t & abs(rowSums(H)) == nz),
    mediation = {
      pairs <- .check_pairs(pairs, K)
      pc <- rowSums((H[, pairs[, 1L], drop = FALSE] != 0) & (H[, pairs[, 2L], drop = FALSE] != 0))
      pc == 0 | nz <= t
    }
  )
}

.check_pairs <- function(pairs, K) {
  if (is.null(pairs)) {
    if (K %% 2L != 0L) stop("Mediation with odd K requires 'pairs'.", call. = FALSE)
    return(matrix(seq_len(K), ncol = 2L, byrow = TRUE))
  }
  pairs <- as.matrix(pairs)
  if (ncol(pairs) != 2L || any(!pairs %in% seq_len(K)) || anyDuplicated(as.vector(pairs))) {
    stop("'pairs' must be a two-column matrix of distinct coordinates.", call. = FALSE)
  }
  pairs
}

.lfdr_core_r <- function(Z, spec, p, null_pair, block_cells = 5e6) {
  mp <- .mixture_pieces(spec, p)
  nulls <- null_pair[mp$act]
  J <- nrow(Z)
  out <- numeric(J)
  if (!any(nulls)) return(out)
  if (all(nulls)) return(rep(1, J))
  nb <- .block_rows(length(mp$act), block_cells)
  for (s0 in seq(1L, J, by = nb)) {
    idx <- s0:min(J, s0 + nb - 1L)
    bl <- .block_logdens(Z, idx, mp)
    out[idx] <- as.vector((bl$Ep + bl$Em) %*% as.numeric(nulls)) / bl$s
  }
  pmin(1, pmax(0, out))
}

# ---- Local false discovery rates ----
csm_lfdr <- function(object, Z, target = c("pleiotropy", "replication", "mediation"),
                     t = 1L, pairs = NULL, spec = NULL, block_cells = 5e6, cores = 1L) {
  target <- match.arg(target)
  if (inherits(object, "csm_fit")) {
    spec <- object$spec
    p <- object$params
  } else {
    p <- object
    if (!inherits(spec, "csm_spec")) stop("'spec' is required when 'object' is a parameter list.", call. = FALSE)
  }
  Z <- .check_z(Z, spec$K)
  # A configuration and its negative are null together for every target.
  null_pair <- csm_null_components(spec$pair$H, target, t, pairs)
  .lfdr_core(Z, spec, p, null_pair, block_cells, cores)
}
