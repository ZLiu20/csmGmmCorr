library(csmGmmCorr)

ns <- asNamespace("csmGmmCorr")
assert <- function(condition, label) {
  if (!isTRUE(condition)) stop("FAILED: ", label, call. = FALSE)
  cat("ok -", label, "\n")
}

# Fixed reference example for all scientific targets.
spec4 <- csm_spec(4, M = 1, b_mu = 20)
size_mass <- c(0.60, 0.20, 0.10, 0.06, 0.04)
weights4 <- size_mass[spec4$grp_size + 1L] /
  choose(4, spec4$grp_size)
cor4 <- 0.7 * diag(4) + 0.3
par4 <- csm_params(spec4, varpi = weights4, mu = rep(2, 4), Sigma = cor4)
z4 <- rbind(c(12, 0, 0, 0), c(5, 0, 0, 0),
            c(12, 3, 0, 0), c(5, 3, 0, 0))
assert(max(abs(csm_lfdr(par4, z4, "pleiotropy", 1, spec = spec4) -
                 c(0.00114773, 0.75807354, 0.00376523, 0.08316094))) < 1e-7,
       "paper pleiotropy local-fdr values")
assert(max(abs(csm_lfdr(par4, z4, "mediation", 1, spec = spec4) -
                 c(0.00428877, 0.85201978, 0.04364431, 0.12048256))) < 1e-7,
       "paper mediation local-fdr values")
replication_z <- rbind(c(4, 1.5, 1.5, 1.5), c(12, 1.5, 1.5, 1.5))
assert(max(abs(csm_lfdr(par4, replication_z, "replication", 1, spec = spec4) -
                 c(0.79003448, 0.99829545))) < 1e-7,
       "paper replication local-fdr values")

# Native and R E-step implementations must agree.
set.seed(2)
spec3 <- csm_spec(3, M = 2, b_mu = 10)
G <- length(spec3$grp_b)
weights <- runif(G); weights <- weights / sum(weights)
mu <- matrix(0, G, 3)
for (g in which(spec3$grp_b > 0)) {
  mu[g, spec3$grp_A[g, ]] <- runif(sum(spec3$grp_A[g, ]), 1, 3)
}
mu[spec3$full, ] <- 3 + runif(length(spec3$full) * 3)
cor3 <- matrix(c(1, .4, .2, .4, 1, -.1, .2, -.1, 1), 3)
par3 <- csm_params(spec3, varpi = weights, mu = mu, Sigma = cor3)
z3 <- csm_simulate(3000, par3, spec3)$Z
native <- ns$.estep(z3, spec3, par3)
reference <- ns$.estep_r(z3, spec3, par3)
assert(abs(native$F - reference$F) < 1e-12 &&
         max(abs(native$T - reference$T)) < 1e-12 &&
         max(abs(native$S - reference$S)) < 1e-12,
       "native E-step equals the R reference")

# Magnitude and correlation projections preserve their constraints.
set.seed(3)
es <- ns$.estep(csm_simulate(500, par3, spec3)$Z, spec3, par3)
stats <- ns$.mu_stats(es, spec3)
bad_mu <- mu
bad_mu[spec3$full, ] <- 1.2
start_mu <- ns$.project_mu(bad_mu, spec3)
new_mu <- ns$.update_mu(stats, spec3, start_mu, es$Omega)
assert(ns$.ordering_ok(new_mu, spec3), "magnitude ordering after the M-step")
assert(ns$.q_mu(new_mu, stats, es$Omega) >=
         ns$.q_mu(start_mu, stats, es$Omega) - 1e-8,
       "magnitude M-step is an ascent step")

candidate <- matrix(c(1, .95, .9, .95, 1, .99, .9, .99, 1), 3) + 0.3
projected <- ns$.project_cor(candidate, 0.1)
assert(max(abs(diag(projected) - 1)) < 1e-12 &&
         min(eigen(projected, symmetric = TRUE, only.values = TRUE)$values) >=
           0.1 - 1e-10,
       "correlation projection")

# Full model fitting and conversion.
set.seed(5)
truth_spec <- csm_spec(3, M = 1, b_mu = 8, kappa = 0.05)
truth_weights <- c(.70, .08, .08, .04, .05, .02, .02, .01)
truth_mu <- truth_spec$grp_A * 0
for (g in 2:length(truth_spec$grp_b)) {
  truth_mu[g, truth_spec$grp_A[g, ]] <- 2.5
}
truth_mu[truth_spec$full, ] <- 3.2
truth_cor <- matrix(c(1, .3, .1, .3, 1, .2, .1, .2, 1), 3)
truth <- csm_params(truth_spec, varpi = truth_weights,
                    mu = truth_mu, Sigma = truth_cor)
simulated <- csm_simulate(20000, truth, truth_spec)
fit <- csm_fit(simulated$Z, csm_spec(3, M = 1, b_mu = 8),
               control = csm_control(max_iter = 3000, penalty = 0,
                                     weight_penalty = 0))
assert(fit$converged, "full model converged")
selected_trace <- fit$traces[[match(fit$selected, fit$runs$start)]]
assert(all(diff(selected_trace) >= -1e-10 * (1 + abs(head(selected_trace, -1L)))),
       "fitting objective is nondecreasing")
assert(length(ns$.feasibility(fit$params, fit$spec)) == 0L,
       "fitted parameters are feasible")
assert(max(abs(fit$params$Sigma - truth_cor)) < 0.05,
       "correlation recovery")
assert(max(abs(fit$params$varpi - truth_weights)) < 0.05,
       "mixture-weight recovery")

converted <- as_sun_lists(fit)
assert(abs(sum(unlist(converted$piList)) - 1) < 1e-12,
       "csmGmm block probabilities retain total mass")
assert(max(abs(unlist(converted$piList) - fit$params$varpi)) < 1e-12,
       "csmGmm conversion does not divide block mass over signs")
roundtrip <- csm_params(fit$spec, pi = converted$piList,
                        mu = converted$muList, Sigma = converted$Sigma)
assert(max(abs(roundtrip$varpi - fit$params$varpi)) < 1e-12 &&
         max(abs(roundtrip$mu - fit$params$mu)) == 0,
       "csmGmm list-format round trip")

# Scientific null sets are evaluated directly from signed configurations.
H4 <- rbind(
  c(0, 0, 0, 0), c(1, 0, 0, 0), c(1, 1, 0, 0),
  c(1, -1, 0, 0), c(1, 1, 1, 0), c(1, 0, 1, 0)
)
assert(identical(csm_null_components(H4, "pleiotropy", 1),
                 c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)),
       "pleiotropy null definition")
assert(identical(csm_null_components(H4, "replication", 1),
                 c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE)),
       "same-sign replication null definition")
assert(identical(csm_null_components(H4, "mediation", 1),
                 c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE)),
       "paired-mediation null definition")

# Native thread counts alter scheduling only, not numerical results.
for (K in 2:5) {
  thread_spec <- csm_spec(K, M = 1, b_mu = 10)
  thread_weights <- seq_len(length(thread_spec$grp_b))
  thread_weights <- thread_weights / sum(thread_weights)
  thread_mu <- thread_spec$grp_A * 2.5
  thread_par <- csm_params(thread_spec, varpi = thread_weights,
                           mu = thread_mu, Sigma = diag(K))
  set.seed(100 + K)
  thread_z <- csm_simulate(200, thread_par, thread_spec)$Z
  one <- csm_lfdr(thread_par, thread_z, "pleiotropy", 1,
                  spec = thread_spec, cores = 1)
  three <- csm_lfdr(thread_par, thread_z, "pleiotropy", 1,
                    spec = thread_spec, cores = 3)
  assert(max(abs(one - three)) < 2e-12,
         paste("native thread equivalence K", K))
}

# S3 methods and complete testing output retain their public contracts.
predicted <- predict(fit, simulated$Z[1:25, , drop = FALSE],
                     target = "pleiotropy", t = 1)
likelihood <- logLik(fit)
complete <- csm_test(fit, simulated$Z[1:100, , drop = FALSE],
                     target = "pleiotropy", t = 1, q = 0.1,
                     adjust = TRUE, extrapolation = TRUE)
assert(length(predicted) == 25L && all(is.finite(predicted)),
       "predict method")
assert(inherits(likelihood, "logLik") && attr(likelihood, "nobs") == fit$n,
       "logLik method")
assert(inherits(complete, "csm_test") && length(complete$lfdr) == 100L &&
         complete$n_rejected == sum(complete$rejected),
       "complete test result contract")

# Decision rules and monotone adjustment.
set.seed(7)
for (iteration in seq_len(50)) {
  scores <- round(runif(40)^3, 2)
  result <- csm_stepup(scores, 0.1)
  ordered <- sort(scores)
  expected <- max(c(0, which(cumsum(ordered) / seq_along(ordered) <= 0.1)))
  assert(result$n_rejected == expected, paste("step-up", iteration))
}

set.seed(8)
base_z <- matrix(round(rnorm(900, sd = 2), 1), 300)
scores <- runif(300)^2
for (target in c("pleiotropy", "replication", "mediation")) {
  z <- if (target == "mediation") cbind(base_z, round(rnorm(300, sd = 2), 1)) else base_z
  full <- csm_adjust(z, scores, target, 1, 0.3, all = TRUE)
  decision <- csm_adjust(z, scores, target, 1, 0.3)
  assert(identical(full$rejected, decision$rejected),
         paste("decision adjustment", target))
  assert(csm_incongruity(z, full$adjusted, target, 1)$violations == 0,
         paste("zero incongruity after adjustment", target))
}

# Failed fits default to no discoveries.
failed <- csm_fit(simulated$Z, csm_spec(3, M = 1, b_mu = 8), starts = "equal",
                  control = csm_control(max_iter = 2))
failed_test <- csm_test(failed, simulated$Z, "pleiotropy", 1)
assert(!failed$converged && failed_test$status == "failed" &&
         failed_test$n_rejected == 0L,
       "failed fit makes no discoveries")

cat("All tests passed.\n")
