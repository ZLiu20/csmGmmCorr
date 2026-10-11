# Constrained ECM estimation from multiple starting values.

.check_fix_sigma <- function(value) {
  if (is.matrix(value)) {
    if (!is.numeric(value) || nrow(value) != ncol(value) || any(!is.finite(value)) ||
        max(abs(value - t(value))) > 1e-10 || max(abs(diag(value) - 1)) > 1e-10 ||
        min(eigen(value, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      stop("'fix_sigma' must be FALSE, TRUE, or a positive-definite correlation matrix.",
           call. = FALSE)
    }
    return((value + t(value)) / 2)
  }
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("'fix_sigma' must be FALSE, TRUE, or a positive-definite correlation matrix.",
         call. = FALSE)
  }
  value
}

#' Control parameters for correlated csmGmm fitting
#' @export
csm_control <- function(max_iter = 5000L, tol_loglik = NULL, tol_param = Inf,
                        accelerate = TRUE, warm_n = 0L, penalty = 30,
                        weight_penalty = 10, block_cells = 5e6,
                        sigma_steps = 50L, fix_sigma = FALSE, cores = 1L,
                        verbose = FALSE) {
  max_iter <- .check_integer_scalar(max_iter, "max_iter", 1L, .Machine$integer.max)
  if (!is.null(tol_loglik)) tol_loglik <- .check_numeric_scalar(tol_loglik, "tol_loglik", 0, Inf)
  if (!is.numeric(tol_param) || length(tol_param) != 1L || is.na(tol_param) || tol_param < 0) {
    stop("'tol_param' must be a nonnegative number or Inf.", call. = FALSE)
  }
  warm_n <- .check_integer_scalar(warm_n, "warm_n", 0L, .Machine$integer.max)
  penalty <- .check_numeric_scalar(penalty, "penalty", 0, Inf)
  weight_penalty <- .check_numeric_scalar(weight_penalty, "weight_penalty", 0, Inf)
  block_cells <- .check_numeric_scalar(block_cells, "block_cells", 1000, Inf)
  sigma_steps <- .check_integer_scalar(sigma_steps, "sigma_steps", 0L, .Machine$integer.max)
  cores <- .check_integer_scalar(cores, "cores", 1L, .Machine$integer.max)
  list(max_iter = max_iter, tol_loglik = tol_loglik, tol_param = tol_param,
       accelerate = isTRUE(accelerate), warm_n = warm_n, penalty = penalty,
       weight_penalty = weight_penalty, block_cells = block_cells,
       sigma_steps = sigma_steps, fix_sigma = .check_fix_sigma(fix_sigma),
       cores = cores, verbose = isTRUE(verbose))
}

.fixed_sigma <- function(control, K) {
  value <- control$fix_sigma
  if (is.matrix(value)) {
    if (nrow(value) != K) stop("The fixed correlation matrix has the wrong dimension.", call. = FALSE)
    return(value)
  }
  if (isTRUE(value)) diag(K) else NULL
}

.ecm_map <- function(es, params, spec, XtX, control, n) {
  next_params <- params
  next_params$varpi <- .update_weights(es$T, spec, control$weight_penalty / n)
  stats <- .mu_stats(es, spec)
  if (control$penalty > 0) {
    stats <- .penalize_mu_stats(stats, params$mu, spec, es$Omega, control$penalty / n)
  }
  next_params$mu <- .update_mu(stats, spec, params$mu, es$Omega)
  if (is.null(.fixed_sigma(control, spec$K)) && control$sigma_steps > 0L) {
    scatter <- .scatter(es, next_params$mu, spec, XtX)
    next_params$Sigma <- .update_sigma(scatter, params$Sigma, spec$kappa,
                                      control$sigma_steps)
  }
  fixed <- .fixed_sigma(control, spec$K)
  if (!is.null(fixed)) next_params$Sigma <- fixed
  structure(next_params, class = "csm_params")
}

.parameter_distance <- function(a, b) {
  max(abs(c(a$varpi - b$varpi, a$mu - b$mu, a$Sigma - b$Sigma)))
}

.objective <- function(es, params, spec, control, n) {
  es$F - control$penalty * .magnitude_penalty(params$mu, spec) / n -
    control$weight_penalty * .weight_penalty(params$varpi, spec) / n
}

.combine_params <- function(p0, p1, p2, a, b, c, spec) {
  result <- p0
  positive <- p0$varpi > 0 & p1$varpi > 0 & p2$varpi > 0
  log_weights <- rep(-Inf, length(p0$varpi))
  log_weights[positive] <- a * log(p0$varpi[positive]) +
    b * log(p1$varpi[positive]) + c * log(p2$varpi[positive])
  log_weights[positive] <- pmax(log_weights[positive],
                                log(p2$varpi[positive]) - log(1000))
  maximum <- max(log_weights)
  result$varpi <- exp(log_weights - maximum)
  result$varpi <- result$varpi / sum(result$varpi)
  result$mu <- .project_mu(a * p0$mu + b * p1$mu + c * p2$mu, spec)
  result$Sigma <- .project_cor(a * p0$Sigma + b * p1$Sigma + c * p2$Sigma,
                               spec$kappa)
  structure(result, class = "csm_params")
}

.transformed_parameters <- function(params, positive = params$varpi > 0) {
  c(log(params$varpi[positive]), params$mu, params$Sigma)
}

.aitken_gap <- function(values) {
  n <- length(values)
  if (n < 3L) return(Inf)
  first <- values[n] - values[n - 1L]
  previous <- values[n - 1L] - values[n - 2L]
  if (first == 0 && previous == 0) return(0)
  if (first <= 0 || previous <= 0) return(Inf)
  ratio <- first / previous
  if (ratio >= 1) return(Inf)
  first / (1 - ratio)
}

.fit_one <- function(x, spec, initial, control) {
  n <- nrow(x)
  tolerance <- if (is.null(control$tol_loglik)) 0.01 / n else control$tol_loglik
  XtX <- crossprod(x) / n
  evaluate <- function(params) {
    result <- .estep(x, spec, params, block_cells = control$block_cells,
                     cores = control$cores)
    result$J <- .objective(result, params, spec, control, n)
    result
  }
  map <- function(estep, params) {
    .ecm_map(estep, params, spec, XtX, control, n)
  }

  params <- structure(initial, class = "csm_params")
  current <- evaluate(params)
  trace <- current$J
  status <- "max_iter"
  maps <- 0L
  step_limit <- 1
  cycle_values <- current$J
  numerical_tolerance <- function(value) 1e-10 * (1 + abs(value))
  next_report <- 50L

  while (maps < control$max_iter) {
    first_params <- map(current, params)
    first <- evaluate(first_params)
    second_params <- map(first, first_params)
    second <- evaluate(second_params)
    maps <- maps + 2L
    if (!is.finite(first$J) || !is.finite(second$J)) {
      status <- "nonfinite"
      break
    }
    if (first$J < current$J - numerical_tolerance(current$J) ||
        second$J < first$J - numerical_tolerance(first$J)) {
      status <- "likelihood_decrease"
      break
    }
    trace <- c(trace, first$J, second$J)

    gap <- .aitken_gap(c(current$J, first$J, second$J))
    if (control$accelerate) {
      history <- c(cycle_values, second$J)
      if (length(history) > 3L) {
        gap <- max(gap, history[length(history)] - history[length(history) - 3L])
      }
    }
    distance <- .parameter_distance(second_params, first_params)
    if (gap <= tolerance && distance <= control$tol_param) {
      params <- second_params
      current <- second
      status <- "converged"
      break
    }

    accepted_params <- second_params
    accepted <- second

    if (control$accelerate) {
      positive <- params$varpi > 0 & first_params$varpi > 0 & second_params$varpi > 0
      theta0 <- .transformed_parameters(params, positive)
      theta1 <- .transformed_parameters(first_params, positive)
      theta2 <- .transformed_parameters(second_params, positive)
      r <- theta1 - theta0
      v <- theta2 - theta1 - r
      norm_v <- sqrt(sum(v^2))
      if (is.finite(norm_v) && norm_v > 0) {
        alpha <- max(-step_limit, min(-1, -sqrt(sum(r^2)) / norm_v))
        extrapolated <- .combine_params(params, first_params, second_params,
                                        (1 + alpha)^2,
                                        -2 * alpha * (1 + alpha), alpha^2, spec)
        fixed <- .fixed_sigma(control, spec$K)
        if (!is.null(fixed)) extrapolated$Sigma <- fixed
        trial <- tryCatch({
          extrapolated_estep <- evaluate(extrapolated)
          trial_params <- map(extrapolated_estep, extrapolated)
          trial_estep <- evaluate(trial_params)
          list(params = trial_params, estep = trial_estep)
        }, error = function(e) NULL)
        maps <- maps + 1L
        if (!is.null(trial) && is.finite(trial$estep$J) &&
            trial$estep$J >= second$J) {
          accepted_params <- trial$params
          accepted <- trial$estep
          trace <- c(trace, accepted$J)
          if (alpha == -step_limit) step_limit <- 4 * step_limit
        } else {
          step_limit <- max(1, step_limit / 4)
        }
      }
    }

    cycle_values <- c(cycle_values, accepted$J)
    params <- accepted_params
    current <- accepted
    if (control$verbose && maps >= next_report) {
      message(sprintf("ECM maps %d: average log-likelihood %.10f", maps, current$F))
      next_report <- next_report + 50L
    }
  }
  list(params = params, loglik = current$F, objective = current$J,
       iterations = maps, status = status, trace = trace)
}

.warm_rows <- function(n, m) {
  order(((seq_len(n) * 0.6180339887498949) %% 1))[seq_len(m)]
}

.fit_warm <- function(x, spec, initial, control) {
  n <- nrow(x)
  m <- control$warm_n
  if (m < 1L || n <= 3L * m) return(.fit_one(x, spec, initial, control))
  preliminary <- .fit_one(x[.warm_rows(n, m), , drop = FALSE], spec, initial, control)
  final <- .fit_one(x, spec, preliminary$params, control)
  final$warm_iterations <- preliminary$iterations
  final
}

.base_magnitudes <- function(x, spec) {
  out <- apply(abs(x), 2L, function(values) {
    selected <- values[values > 2]
    if (length(selected) >= 5L) mean(selected) else 2.5
  })
  upper <- if (is.null(spec$b_mu)) Inf else spec$b_mu
  pmin(upper, pmax(spec$a_mu, out))
}

.magnitudes_from_base <- function(base, spec) {
  G <- length(spec$grp_b)
  mu <- matrix(0, G, spec$K)
  for (g in which(spec$grp_b != 0L)) {
    components <- spec$M[spec$grp_b[g] + 1L]
    fraction <- if (spec$grp_b[g] == 2^spec$K - 1L) {
      1 + 0.5 * (spec$grp_m[g] - 1L) / max(1L, components - 1L)
    } else {
      0.5 + 0.5 * spec$grp_m[g] / components
    }
    mu[g, spec$grp_A[g, ]] <- base[spec$grp_A[g, ]] * fraction
  }
  .project_mu(mu, spec)
}

.start_params <- function(x, spec, type = c("null", "equal", "shared")) {
  type <- match.arg(type)
  G <- length(spec$grp_b)
  weights <- switch(type,
    equal = rep(1 / G, G),
    null = c(0.5, rep(0.5 / (G - 1L), G - 1L)),
    shared = rep(1 / G, G)
  )
  list(varpi = weights, mu = .magnitudes_from_base(.base_magnitudes(x, spec), spec),
       Sigma = diag(spec$K))
}

.random_start <- function(x, spec, base) {
  G <- length(spec$grp_b)
  weights <- base$varpi * exp(stats::rnorm(G, 0, 0.5))
  weights <- weights / sum(weights)
  mu <- .project_mu(base$mu * exp(matrix(stats::rnorm(G * spec$K, 0, 0.25), G)), spec)
  rho <- stats::runif(1L, 0, 0.3)
  Sigma <- .project_cor((1 - rho) * diag(spec$K) + rho, spec$kappa)
  list(varpi = weights, mu = mu, Sigma = Sigma)
}

.shared_start <- function(x, spec, control) {
  restricted <- .make_spec(spec$K, 1L, spec$S, "shared", spec$a_mu,
                           spec$b_mu, spec$kappa)
  restricted_control <- control
  restricted_control$penalty <- 0
  restricted_control$weight_penalty <- 0
  start <- .start_params(x, restricted, "null")
  fixed <- .fixed_sigma(restricted_control, spec$K)
  if (!is.null(fixed)) start$Sigma <- fixed
  fit <- .fit_warm(x, restricted, start, restricted_control)
  if (fit$status != "converged") stop("The shared-magnitude start did not converge.", call. = FALSE)
  support_weight <- stats::setNames(fit$params$varpi, restricted$grp_b)
  weights <- support_weight[as.character(spec$grp_b)] / spec$M[spec$grp_b + 1L]
  weights <- as.numeric(weights / sum(weights))
  base <- vapply(seq_len(spec$K), function(k) {
    candidates <- fit$params$mu[restricted$grp_A[, k], k]
    if (length(candidates)) max(candidates) else spec$a_mu
  }, numeric(1))
  list(varpi = weights, mu = .magnitudes_from_base(base, spec), Sigma = fit$params$Sigma)
}

#' Fit a correlated csmGmm
#' @export
csm_fit <- function(x, spec, starts = c("null", "equal", "shared"), n_random = 0L,
                    init = NULL, control = csm_control(), seed = NULL) {
  if (!inherits(spec, "csm_spec")) stop("'spec' must come from csm_spec().", call. = FALSE)
  x <- .check_z(x, spec$K, "x")
  if (nrow(x) < 2L) stop("'x' must have at least two rows.", call. = FALSE)
  if (is.null(spec$b_mu)) spec$b_mu <- max(spec$a_mu + 1, ceiling(max(abs(x))) + 1)
  starts <- unique(starts)
  if (!all(starts %in% c("null", "equal", "shared"))) {
    stop("'starts' may contain only 'null', 'equal', and 'shared'.", call. = FALSE)
  }
  if (identical(spec$magnitudes, "shared")) starts <- setdiff(starts, "shared")
  n_random <- .check_integer_scalar(n_random, "n_random", 0L, .Machine$integer.max)
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(init) && inherits(init, "csm_params")) init <- list(init)
  if (is.null(init)) init <- list()

  candidates <- stats::setNames(lapply(starts, function(type) {
    if (type == "shared") {
      tryCatch(.shared_start(x, spec, control), error = function(e) e)
    } else {
      .start_params(x, spec, type)
    }
  }), starts)
  base <- .start_params(x, spec, "null")
  if (n_random) {
    random <- lapply(seq_len(n_random), function(i) .random_start(x, spec, base))
    names(random) <- paste0("random", seq_len(n_random))
    candidates <- c(candidates, random)
  }
  for (i in seq_along(init)) {
    problems <- .feasibility(init[[i]], spec)
    if (length(problems)) {
      stop(sprintf("init[[%d]] is infeasible: %s", i, paste(problems, collapse = "; ")),
           call. = FALSE)
    }
    candidates[[paste0("init", i)]] <- unclass(init[[i]])
  }
  fixed <- .fixed_sigma(control, spec$K)
  if (!is.null(fixed)) {
    if (min(eigen(fixed, symmetric = TRUE, only.values = TRUE)$values) < spec$kappa) {
      stop("The fixed correlation matrix has an eigenvalue below kappa.", call. = FALSE)
    }
    candidates <- lapply(candidates, function(value) {
      if (!inherits(value, "error")) value$Sigma <- fixed
      value
    })
  }
  if (!length(candidates)) stop("No starting values were requested.", call. = FALSE)

  runs <- lapply(candidates, function(candidate) {
    if (inherits(candidate, "error")) {
      return(list(status = "error", message = conditionMessage(candidate),
                  loglik = NA_real_, objective = NA_real_, iterations = 0L,
                  trace = numeric()))
    }
    tryCatch(.fit_warm(x, spec, candidate, control), error = function(e) {
      list(status = "error", message = conditionMessage(e), loglik = NA_real_,
           objective = NA_real_, iterations = 0L, trace = numeric())
    })
  })
  summary <- data.frame(
    start = names(runs),
    status = vapply(runs, `[[`, character(1), "status"),
    iterations = vapply(runs, function(run) as.integer(run$iterations), integer(1)),
    warm_iterations = vapply(runs, function(run) if (is.null(run$warm_iterations)) NA_integer_ else as.integer(run$warm_iterations), integer(1)),
    loglik = vapply(runs, function(run) as.numeric(run$loglik), numeric(1)),
    objective = vapply(runs, function(run) as.numeric(run$objective), numeric(1)),
    stringsAsFactors = FALSE
  )
  eligible <- which(summary$status == "converged")
  if (!length(eligible)) eligible <- which(is.finite(summary$objective))
  selected <- if (length(eligible)) eligible[which.max(summary$objective[eligible])] else NA_integer_
  structure(list(
    params = if (is.na(selected)) NULL else structure(runs[[selected]]$params, class = "csm_params"),
    loglik = if (is.na(selected)) NA_real_ else runs[[selected]]$loglik,
    objective = if (is.na(selected)) NA_real_ else runs[[selected]]$objective,
    selected = if (is.na(selected)) NA_character_ else summary$start[selected],
    converged = !is.na(selected) && summary$status[selected] == "converged",
    runs = summary,
    traces = lapply(runs, `[[`, "trace"),
    spec = spec,
    n = nrow(x),
    control = control
  ), class = "csm_fit")
}

#' @export
print.csm_fit <- function(x, ...) {
  cat("Correlated csmGmm fit\n")
  cat(sprintf("  K = %d, N = %d, largest support size = %d\n", x$spec$K, x$n, x$spec$S))
  if (x$converged) {
    cat(sprintf("  selected start: %s; average log-likelihood: %.6f\n", x$selected, x$loglik))
  } else {
    cat("  no starting value converged; inference defaults to no discoveries\n")
  }
  print(x$runs, row.names = FALSE)
  invisible(x)
}

#' @export
logLik.csm_fit <- function(object, ...) {
  spec <- object$spec
  G <- length(spec$grp_b)
  n_weights <- G - 1L
  n_magnitudes <- if (identical(spec$magnitudes, "shared")) spec$K else sum(spec$grp_A)
  fixed <- .fixed_sigma(object$control, spec$K)
  n_correlation <- if (is.null(fixed)) spec$K * (spec$K - 1L) / 2 else 0
  structure(object$loglik * object$n, df = n_weights + n_magnitudes + n_correlation,
            nobs = object$n, class = "logLik")
}
