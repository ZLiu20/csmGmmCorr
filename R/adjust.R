# Prespecified ranking orders, monotone adjustment, and incongruity counts.

.kth_largest <- function(z, k) {
  apply(z, 1L, function(v) sort(v, decreasing = TRUE)[k])
}

.order_representation <- function(z, target, t) {
  if (target == "replication") {
    k <- t + 1L
    positive <- .kth_largest(z, k) >= .kth_largest(-z, k)
    return(list(group = ifelse(positive, "+", "-"),
                w = z * ifelse(positive, 1, -1)))
  }
  has_zero <- rowSums(z == 0) > 0L
  group <- apply(z > 0, 1L, function(v) paste(as.integer(v), collapse = ""))
  group[has_zero] <- paste0("zero", which(has_zero))
  list(group = group, w = abs(z))
}

.dominance_max <- function(w, score, group, idx = seq_along(score)) {
  out <- score
  for (label in unique(group[idx])) {
    members <- which(group == label)
    targets <- idx[group[idx] == label]
    if (length(members) < 2L) next
    wm <- w[members, , drop = FALSE]
    sm <- score[members]
    ordering <- order(wm[, 1L], decreasing = TRUE)
    wm <- wm[ordering, , drop = FALSE]
    sm <- sm[ordering]
    first <- wm[, 1L]
    for (i in targets) {
      n_candidates <- sum(first >= w[i, 1L])
      if (!n_candidates) next
      candidates <- seq_len(n_candidates)
      if (!length(candidates)) next
      candidate_w <- wm[candidates, , drop = FALSE]
      dominates <- rowSums(candidate_w >= rep(w[i, ], each = length(candidates))) == ncol(w)
      if (any(dominates)) out[i] <- max(out[i], sm[candidates][dominates])
    }
  }
  out
}

#' Monotone adjustment of local false discovery rates
#' @export
csm_adjust <- function(z, lfdr,
                       target = c("pleiotropy", "replication", "mediation"),
                       t = 1L, q = 0.1, all = FALSE) {
  target <- match.arg(target)
  if (is.data.frame(z)) z <- as.matrix(z)
  if (!is.matrix(z) || !is.numeric(z) || nrow(z) != length(lfdr) ||
      any(!is.finite(z))) {
    stop("'z' must be a finite numeric matrix with one row per lfdr value.", call. = FALSE)
  }
  if (!is.numeric(t) || length(t) != 1L || t < 1L || t >= ncol(z)) {
    stop("'t' must satisfy 1 <= t < K.", call. = FALSE)
  }
  representation <- .order_representation(z, target, as.integer(t))
  if (all) {
    adjusted <- .dominance_max(representation$w, lfdr, representation$group)
    decision <- csm_stepup(adjusted, q)
    return(list(adjusted = adjusted, rejected = decision$rejected,
                n_rejected = decision$n_rejected, exact = rep(TRUE, length(lfdr))))
  }

  original_order <- order(lfdr, seq_along(lfdr))
  n_current <- min(length(lfdr), max(100L, 2L * csm_stepup(lfdr, q)$n_rejected))
  exact <- logical(length(lfdr))
  repeat {
    candidates <- original_order[seq_len(n_current)]
    adjusted_candidates <- .dominance_max(representation$w, lfdr,
                                           representation$group, candidates)
    local_scores <- lfdr
    local_scores[candidates] <- adjusted_candidates[candidates]
    decision <- csm_stepup(local_scores, q)
    exact[candidates] <- TRUE
    if (n_current == length(lfdr) || decision$n_rejected == 0L ||
        (all(exact[decision$rejected]) &&
         max(local_scores[decision$rejected]) < min(lfdr[!exact]))) break
    n_current <- min(length(lfdr), 2L * n_current)
  }
  adjusted <- rep(NA_real_, length(lfdr))
  adjusted[exact] <- local_scores[exact]
  list(adjusted = adjusted, rejected = decision$rejected,
       n_rejected = decision$n_rejected, exact = exact)
}

#' Count incongruous ordered pairs
#' @export
csm_incongruity <- function(z, scores,
                            target = c("pleiotropy", "replication", "mediation"),
                            t = 1L, max_n = 5000L, tolerance = 0, seed = NULL) {
  target <- match.arg(target)
  if (is.data.frame(z)) z <- as.matrix(z)
  if (!is.matrix(z) || nrow(z) != length(scores)) {
    stop("'z' must have one row per score.", call. = FALSE)
  }
  idx <- seq_len(nrow(z))
  if (length(idx) > max_n) {
    if (!is.null(seed)) set.seed(seed)
    idx <- sort(sample(idx, max_n))
  }
  zs <- z[idx, , drop = FALSE]
  s <- scores[idx]
  representation <- .order_representation(zs, target, as.integer(t))
  comparable <- 0
  violations <- 0
  for (label in unique(representation$group)) {
    members <- which(representation$group == label)
    if (length(members) < 2L) next
    for (i in members) {
      others <- members[members != i]
      dominates <- rowSums(representation$w[others, , drop = FALSE] >=
                             rep(representation$w[i, ], each = length(others))) == ncol(z) &
        rowSums(representation$w[others, , drop = FALSE] !=
                  rep(representation$w[i, ], each = length(others))) > 0L
      comparable <- comparable + sum(dominates)
      violations <- violations + sum(s[others][dominates] > s[i] + tolerance)
    }
  }
  list(n = length(idx), comparable_pairs = comparable, violations = violations,
       violation_rate = if (comparable) violations / comparable else NA_real_)
}
