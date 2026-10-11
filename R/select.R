# Selection of the number of magnitude components by BIC.

.default_structures <- function() {
  list(A = c(single = 1L, multiple = 1L),
       B = c(single = 2L, multiple = 1L),
       C = c(single = 2L, multiple = 2L))
}

#' Select a magnitude-component structure by BIC
#' @export
csm_select <- function(x, K, structures = .default_structures(), bic_tol = 0,
                       starts = c("null", "shared"), control = csm_control(),
                       S = K, a_mu = 0.5, b_mu = NULL, kappa = 0.05,
                       seed = NULL) {
  if (!is.list(structures) || !length(structures) || is.null(names(structures))) {
    stop("'structures' must be a named list of c(single=, multiple=).", call. = FALSE)
  }
  if (!is.numeric(bic_tol) || length(bic_tol) != 1L || bic_tol < 0) {
    stop("'bic_tol' must be nonnegative.", call. = FALSE)
  }
  fits <- vector("list", length(structures))
  table <- data.frame(structure = names(structures), single = NA_integer_,
                      multiple = NA_integer_, df = NA_real_, logLik = NA_real_,
                      BIC = NA_real_, converged = FALSE, minutes = NA_real_)
  for (i in seq_along(structures)) {
    st <- structures[[i]]
    M <- csm_structure(K, st[["single"]], st[["multiple"]])
    spec <- csm_spec(K, M, S, a_mu, b_mu, kappa)
    timed <- system.time(fit <- csm_fit(x, spec, starts = starts,
                                       control = control, seed = seed))["elapsed"]
    fits[[i]] <- fit
    table$single[i] <- st[["single"]]
    table$multiple[i] <- st[["multiple"]]
    table$converged[i] <- isTRUE(fit$converged)
    table$minutes[i] <- as.numeric(timed) / 60
    if (isTRUE(fit$converged)) {
      ll <- logLik(fit)
      table$logLik[i] <- as.numeric(ll)
      table$df[i] <- attr(ll, "df")
      table$BIC[i] <- -2 * as.numeric(ll) + attr(ll, "df") * log(attr(ll, "nobs"))
    }
  }
  valid <- which(table$converged & is.finite(table$BIC))
  if (!length(valid)) stop("No structure converged.", call. = FALSE)
  best <- valid[which.min(table$BIC[valid])]
  close <- valid[table$BIC[valid] <= table$BIC[best] + bic_tol]
  chosen <- close[which.max(table$df[close])]
  table$delta_BIC <- table$BIC - table$BIC[best]
  table$chosen <- seq_len(nrow(table)) == chosen
  out <- fits[[chosen]]
  out$selection <- table
  out$all_fits <- fits
  out
}
