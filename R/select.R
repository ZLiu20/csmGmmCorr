# ---- Number of magnitude components per support, chosen by BIC ----
# The alternative distribution of each support is a finite Gaussian mixture with
# M_b components. Effect sizes in practice vary continuously, so M_b is an
# approximation; too few components can bias the local false discovery rates
# (a support absorbs the spread of another), while extra components mostly cost
# computation. csm_select fits a small set of prespecified structures (by default
# 1/1, 2/1 and 2/2 components for supports of size one / two or more) and keeps
# the one with the smallest BIC (optionally preferring the richer structure when
# the BIC difference is small).

# M vector (length 2^K, indexed by b + 1) for a structure given by the number of
# components for supports of size 1 and of size >= 2. The null support has one.
csm_structure <- function(K, single = 1L, multiple = 1L) {
  b <- 0:(2L^K - 1L)
  size <- vapply(b, function(x) sum(as.integer(intToBits(x))[seq_len(K)]), 0L)
  M <- ifelse(size == 0L, 1L, ifelse(size == 1L, as.integer(single), as.integer(multiple)))
  as.integer(M)
}

.default_structures <- function() {
  list(A = c(single = 1L, multiple = 1L), B = c(single = 2L, multiple = 1L),
       C = c(single = 2L, multiple = 2L))
}

csm_select <- function(X, K, structures = .default_structures(), bic_tol = 0,
                       starts = c("null", "shared"), control = csm_control(),
                       S = K, a_mu = 0.5, b_mu = NULL, kappa = 0.05, seed = NULL) {
  if (!is.list(structures) || !length(structures) || is.null(names(structures))) {
    stop("'structures' must be a named list of c(single = , multiple = ).", call. = FALSE)
  }
  if (!is.numeric(bic_tol) || length(bic_tol) != 1L || bic_tol < 0) {
    stop("'bic_tol' must be a non-negative number.", call. = FALSE)
  }
  fits <- list()
  tab <- data.frame(structure = names(structures), single = NA_integer_, multiple = NA_integer_,
                    df = NA_real_, loglik = NA_real_, BIC = NA_real_, converged = FALSE,
                    minutes = NA_real_, stringsAsFactors = FALSE)
  for (i in seq_along(structures)) {
    st <- structures[[i]]
    M <- csm_structure(K, st[["single"]], st[["multiple"]])
    sp <- csm_spec(K, M = M, S = S, a_mu = a_mu, b_mu = b_mu, kappa = kappa)
    tm <- system.time(f <- csm_fit(X, sp, starts = starts, control = control, seed = seed))[["elapsed"]]
    fits[[i]] <- f
    tab$single[i] <- st[["single"]]; tab$multiple[i] <- st[["multiple"]]
    tab$converged[i] <- f$converged; tab$minutes[i] <- tm / 60
    if (f$converged) {
      ll <- logLik(f)
      tab$df[i] <- attr(ll, "df"); tab$loglik[i] <- as.numeric(ll)
      tab$BIC[i] <- -2 * as.numeric(ll) + attr(ll, "df") * log(f$n)
    }
  }
  ok <- which(tab$converged & is.finite(tab$BIC))
  if (!length(ok)) stop("No structure converged.", call. = FALSE)
  best <- ok[which.min(tab$BIC[ok])]
  # within bic_tol of the minimum, prefer the structure with the most parameters
  near <- ok[tab$BIC[ok] <= tab$BIC[best] + bic_tol]
  pick <- near[which.max(tab$df[near])]
  tab$delta_BIC <- tab$BIC - tab$BIC[best]
  tab$chosen <- seq_len(nrow(tab)) == pick
  out <- fits[[pick]]
  out$selection <- tab
  out$all_fits <- fits
  out
}
