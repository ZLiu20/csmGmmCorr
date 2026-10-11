# Prespecified ranking orders, the monotone adjustment (main-text Definition 1
# and Theorem 4) and incongruity counts.
#
# Every order is represented by a group label and a transformed vector w:
# feature j is at least as strong as feature i (i precedes j) exactly when both
# have the same group label and w_j >= w_i coordinatewise.
#   Pleiotropy and mediation: group = open sign orthant, w = |z|; a vector with
#     a zero coordinate forms its own singleton class.
#   Replication: group = directional branch B+ or B- (Supplementary Section
#     S4), w = z on B+ and w = -z on B-.

.kth_largest <- function(Z, k) apply(Z, 1L, function(v) sort(v, decreasing = TRUE)[k])

.order_repr <- function(Z, target, t) {
  J <- nrow(Z)
  if (target == "replication") {
    k <- t + 1L
    plus <- .kth_largest(Z, k) >= .kth_largest(-Z, k)
    list(group = ifelse(plus, "+", "-"), w = Z * ifelse(plus, 1, -1))
  } else {
    zero <- rowSums(Z == 0) > 0
    g <- apply(Z > 0, 1L, function(v) paste(as.integer(v), collapse = ""))
    g[zero] <- paste0("zero", which(zero))
    list(group = g, w = abs(Z))
  }
}

# v_i = max{ s_j : j in the same group, w_j >= w_i coordinatewise } for i in idx.
.dominance_max <- function(W, s, group, idx = seq_along(s)) {
  out <- s
  for (gname in unique(group[idx])) {
    members <- which(group == gname)
    targets <- idx[group[idx] == gname]
    if (length(members) == 1L) next
    Wm <- W[members, , drop = FALSE]
    sm <- s[members]
    o <- order(Wm[, 1L], decreasing = TRUE)
    Wm <- Wm[o, , drop = FALSE]
    sm <- sm[o]
    first <- Wm[, 1L]
    for (i in targets) {
      n_c <- sum(first >= W[i, 1L])
      if (n_c == 0L) next
      cand <- which(sm[seq_len(n_c)] > out[i])
      if (!length(cand)) next
      Wc <- Wm[cand, , drop = FALSE]
      dom <- rowSums(Wc >= rep(W[i, ], each = length(cand))) == ncol(W)
      if (any(dom)) out[i] <- max(out[i], sm[cand][dom])
    }
  }
  out
}

# ---- Monotone adjustment of local false discovery rates ----
csm_adjust <- function(Z, lfdr, target = c("pleiotropy", "replication", "mediation"), t = 1L,
                       q = 0.1, all = FALSE) {
  target <- match.arg(target)
  if (is.data.frame(Z)) Z <- as.matrix(Z)
  if (!is.matrix(Z) || !is.numeric(Z) || nrow(Z) != length(lfdr) || any(!is.finite(Z))) {
    stop("'Z' must be a finite numeric matrix with one row per lfdr value.", call. = FALSE)
  }
  if (!is.numeric(t) || t < 1 || t >= ncol(Z)) stop("'t' must satisfy 1 <= t < K.", call. = FALSE)
  J <- length(lfdr)
  rep_ <- .order_repr(Z, target, as.integer(t))
  if (all) {
    v <- .dominance_max(rep_$w, lfdr, rep_$group)
    su <- csm_stepup(v, q)
    return(list(adjusted = v, rejected = su$rejected, n_rejected = su$n_rejected,
                exact = rep(TRUE, J)))
  }
  # Decision mode: compute adjusted scores exactly only for the features with the
  # smallest unadjusted scores, and certify that the resulting decisions equal
  # those of the fully adjusted scores. Adjusted scores are lower bounded by the
  # unadjusted ones, so the certificate below is sufficient (see ?csm_adjust).
  o <- order(lfdr, seq_len(J))
  n0 <- min(J, max(100L, 2L * csm_stepup(lfdr, q)$n_rejected))
  repeat {
    C <- o[seq_len(n0)]
    v <- .dominance_max(rep_$w, lfdr, rep_$group, C)
    L <- lfdr
    L[C] <- v[C]
    su <- csm_stepup(L, q)
    exact <- logical(J)
    exact[C] <- TRUE
    if (n0 == J || su$n_rejected == 0L ||
        (all(exact[su$rejected]) && max(L[su$rejected]) < min(lfdr[!exact]))) break
    n0 <- min(J, 2L * n0)
  }
  adjusted <- rep(NA_real_, J)
  adjusted[exact] <- v[exact]
  list(adjusted = adjusted, rejected = su$rejected, n_rejected = su$n_rejected, exact = exact)
}

# ---- Incongruity diagnostic ----
csm_incongruity <- function(Z, scores, target = c("pleiotropy", "replication", "mediation"),
                            t = 1L, max_n = 5000L, tolerance = 0, seed = NULL) {
  target <- match.arg(target)
  if (is.data.frame(Z)) Z <- as.matrix(Z)
  if (!is.matrix(Z) || nrow(Z) != length(scores)) {
    stop("'Z' must have one row per score.", call. = FALSE)
  }
  idx <- seq_len(nrow(Z))
  if (length(idx) > max_n) {
    if (!is.null(seed)) set.seed(seed)
    idx <- sort(sample(idx, max_n))
  }
  Zs <- Z[idx, , drop = FALSE]
  s <- scores[idx]
  rp <- .order_repr(Zs, target, as.integer(t))
  comparable <- 0
  violations <- 0
  for (gname in unique(rp$group)) {
    m <- which(rp$group == gname)
    if (length(m) < 2L) next
    for (i in m) {
      others <- m[m != i]
      dom <- rowSums(rp$w[others, , drop = FALSE] >= rep(rp$w[i, ], each = length(others))) == ncol(Z) &
        rowSums(rp$w[others, , drop = FALSE] != rep(rp$w[i, ], each = length(others))) > 0
      comparable <- comparable + sum(dom)
      v <- sum(s[others[dom]] > s[i] + tolerance)
      violations <- violations + v
    }
  }
  list(n = length(idx), comparable_pairs = comparable, violations = violations,
       violation_rate = if (comparable > 0) violations / comparable else NA_real_)
}
