library(csmGmmCorr)
ns <- asNamespace("csmGmmCorr")
check <- function(cond, what) {
  if (!isTRUE(cond)) stop("FAILED: ", what, call. = FALSE)
  cat("ok -", what, "\n")
}

## 1. Paper example (main-text Sections 3.3 and 4.1), against an independent
##    Python implementation.
sp <- csm_spec(4, M = 1, b_mu = 20)
lam <- c(0.6, 0.2, 0.1, 0.06, 0.04)
varpi <- lam[sp$grp_size + 1] / choose(4, sp$grp_size)
Sig <- 0.7 * diag(4) + 0.3
p <- csm_params(sp, varpi = varpi, mu = rep(2, 4), Sigma = Sig)
pts <- rbind(c(12, 0, 0, 0), c(5, 0, 0, 0), c(12, 3, 0, 0), c(5, 3, 0, 0))
check(max(abs(csm_lfdr(p, pts, "pleiotropy", 1, spec = sp) -
  c(0.00114773, 0.75807354, 0.00376523, 0.08316094))) < 1e-7, "pleiotropy lfdr, paper example")
check(max(abs(csm_lfdr(p, pts, "mediation", 1, spec = sp) -
  c(0.00428877, 0.85201978, 0.04364431, 0.12048256))) < 1e-7, "mediation lfdr, paper example")
check(max(abs(csm_lfdr(p, rbind(c(4, 1.5, 1.5, 1.5), c(12, 1.5, 1.5, 1.5)), "replication", 1, spec = sp) -
  c(0.79003448, 0.99829545))) < 1e-7, "replication lfdr, paper example")
p_ind <- csm_params(sp, varpi = varpi, mu = rep(2, 4), Sigma = diag(4))
check(abs(csm_lfdr(p_ind, pts[1, ], "pleiotropy", 1, spec = sp) - 0.86669956) < 1e-7,
      "pleiotropy lfdr under independence")

## 2. Brute-force two-group posterior for a full model with M_b = 2.
set.seed(2)
K <- 3
sp3 <- csm_spec(K, M = 2, b_mu = 10)
G <- length(sp3$grp_b)
v <- runif(G); v <- v / sum(v)
MU <- matrix(0, G, K)
for (g in which(sp3$grp_b > 0)) MU[g, sp3$grp_A[g, ]] <- runif(sum(sp3$grp_A[g, ]), 1, 3)
MU[sp3$full, ] <- 3 + runif(length(sp3$full) * K)
S3 <- matrix(c(1, .4, .2, .4, 1, -.1, .2, -.1, 1), 3)
p3 <- csm_params(sp3, varpi = v, mu = MU, Sigma = S3)
z <- matrix(rnorm(30, sd = 2), 10)
dmv <- function(x, m, S) exp(-0.5 * sum((x - m) * solve(S, x - m))) / sqrt((2 * pi)^length(x) * det(S))
H <- sp3$comp$H
brute <- apply(z, 1, function(x) {
  d <- vapply(seq_len(nrow(H)), function(r) {
    g <- sp3$comp$grp[r]
    v[g] / 2^sum(H[r, ] != 0) * dmv(x, H[r, ] * MU[g, ], S3)
  }, 0)
  nul <- rowSums(H != 0) <= 1
  sum(d[nul]) / sum(d)
})
check(max(abs(brute - csm_lfdr(p3, z, "pleiotropy", 1, spec = sp3))) < 1e-10,
      "lfdr equals brute-force posterior (K = 3, M_b = 2)")

## 2b. Compiled kernel against the plain-R reference implementation.
zb <- csm_simulate(3000, p3, sp3)$Z
es_c <- ns$.estep(zb, sp3, p3)
es_r <- ns$.estep_r(zb, sp3, p3)
check(abs(es_c$F - es_r$F) < 1e-12 && max(abs(es_c$T - es_r$T)) < 1e-12 &&
      max(abs(es_c$S - es_r$S)) < 1e-12 && identical(es_c$Sg, es_r$Sg),
      "compiled E-step equals the R reference")
## Group sums S reproduce the per-pair forms the M-step used before.
Mean3 <- sp3$pair$H * p3$mu[sp3$pair$grp, , drop = FALSE]
c_old <- rowsum(sp3$pair$H * (es_r$D %*% es_r$Omega), sp3$pair$grp, reorder = TRUE)
check(max(abs(unname(c_old) - ns$.mu_stats(es_c, sp3)$c)) < 1e-12 &&
      max(abs(crossprod(es_r$D, Mean3) - crossprod(es_c$S, {
        Mr <- matrix(0, length(es_c$Sg), sp3$K); Mr[cbind(seq_along(es_c$Sg), es_c$Sk)] <- p3$mu[cbind(es_c$Sg, es_c$Sk)]; Mr }))) < 1e-12,
      "group sums give the same magnitude and scatter statistics")
## Large magnitudes and statistics: the kernel skips groups whose terms are
## below exp(-50) relative to the largest; results agree with the reference.
sp5 <- csm_spec(5, M = 1, b_mu = 12)
p5 <- csm_params(sp5, varpi = ifelse(sp5$grp_b == 0, 0.9, 0.1 / (length(sp5$grp_b) - 1)),
                 mu = sp5$grp_A * (6 + sp5$grp_size), Sigma = 0.3^abs(outer(1:5, 1:5, "-")))
z5 <- csm_simulate(2000, p5, sp5)$Z
e5c <- ns$.estep(z5, sp5, p5); e5r <- ns$.estep_r(z5, sp5, p5)
check(abs(e5c$F - e5r$F) < 1e-12 && max(abs(e5c$T - e5r$T)) < 1e-12 &&
      max(abs(e5c$S - e5r$S)) < 1e-10 * max(1, abs(e5r$S)), "kernel with skipped groups equals the R reference")
if (.Platform$OS.type != "windows") {
  es_2 <- ns$.estep(zb, sp3, p3, cores = 2L)
  check(abs(es_2$F - es_c$F) < 1e-12 && max(abs(es_2$S - es_c$S)) < 1e-12, "two-process E-step equals one process")
}
nm <- csm_null_components(sp3$pair$H, "replication", 1)
check(max(abs(ns$.lfdr_core(zb, sp3, p3, nm) - ns$.lfdr_core_r(zb, sp3, p3, nm))) < 1e-12,
      "compiled lfdr equals the R reference")

## 3. Magnitude update: feasibility, ordering constraint and ascent.
set.seed(3)
es <- ns$.estep(csm_simulate(500, p3, sp3)$Z, sp3, p3)
st <- ns$.mu_stats(es, sp3)
MU_bad <- MU
MU_bad[sp3$full, ] <- 1.2          # violates the ordering constraint
MU_start <- ns$.project_mu(MU_bad, sp3)
check(ns$.ordering_ok(MU_start, sp3), "projection satisfies the ordering constraint")
MU_new <- ns$.update_mu(st, sp3, MU_start, es$Omega)
check(ns$.ordering_ok(MU_new, sp3, 1e-8), "magnitude update keeps the ordering constraint")
check(ns$.Q_mu(MU_new, st, es$Omega) >= ns$.Q_mu(MU_start, st, es$Omega), "magnitude update is an ascent step")
# Projection optimality: no feasible random point is closer than the projection.
Y <- MU_bad + matrix(rnorm(length(MU_bad)), nrow(MU_bad)) * sp3$grp_A
PY <- ns$.project_mu(Y, sp3)
closer <- 0
for (i in 1:200) {
  cand <- ns$.project_mu(PY + 0.05 * matrix(rnorm(length(PY)), nrow(PY)) * sp3$grp_A, sp3)
  if (sum((cand - Y)^2) < sum((PY - Y)^2) - 1e-10) closer <- closer + 1
}
check(closer == 0, "magnitude projection is a nearest feasible point")

## 4. Correlation projection.
A <- matrix(c(1, .95, .9, .95, 1, .99, .9, .99, 1), 3) + 0.3
P <- ns$.project_cor(A, 0.1)
check(max(abs(diag(P) - 1)) < 1e-12 && min(eigen(P, TRUE, TRUE)$values) >= 0.1 - 1e-10,
      "correlation projection is feasible")

## 5. Fitting the full model: monotone likelihood, feasibility, recovery.
set.seed(5)
sp_true <- csm_spec(3, M = 1, b_mu = 8, kappa = 0.05)
Gt <- length(sp_true$grp_b)
vt <- c(0.7, 0.08, 0.08, 0.04, 0.05, 0.02, 0.02, 0.01)
MUt <- sp_true$grp_A * 0
for (g in 2:Gt) MUt[g, sp_true$grp_A[g, ]] <- 2.5
MUt[sp_true$full, ] <- 3.2
St <- matrix(c(1, .3, .1, .3, 1, .2, .1, .2, 1), 3)
pt <- csm_params(sp_true, varpi = vt, mu = MUt, Sigma = St)
sim <- csm_simulate(20000, pt, sp_true)
fit <- csm_fit(sim$Z, csm_spec(3, M = 1, b_mu = 8), control = csm_control(max_iter = 3000, penalty = 0, weight_penalty = 0))
print(fit)
check(fit$converged, "full model converged")
tr <- fit$traces[[fit$selected]]
check(all(diff(tr) >= -1e-10 * (1 + abs(tr[-1]))), "log-likelihood is nondecreasing")
check(length(ns$.feasibility(fit$params, fit$spec)) == 0L, "fitted parameter is feasible")
check(max(abs(fit$params$Sigma - St)) < 0.05, "correlation recovered")
check(max(abs(fit$params$varpi - vt)) < 0.05, "weights recovered")
ll_true <- ns$.estep(sim$Z, fit$spec, unclass(pt), stats = FALSE)$F
check(fit$loglik >= ll_true - 1e-6, "fitted likelihood at least that of the truth")

## 6. Shared-magnitude starting fit (internal) and inference.
sh <- ns$.make_spec(3L, 1L, 3L, "shared", 0.5, 8, 0.05)
fit_sh <- csm_fit(sim$Z, sh, starts = "null", control = csm_control(penalty = 0, weight_penalty = 0))
check(fit_sh$converged && length(ns$.feasibility(fit_sh$params, fit_sh$spec)) == 0L,
      "shared-magnitude starting fit converged and feasible")
check(fit$loglik >= fit_sh$loglik - 1e-6, "full model likelihood at least that of the starting fit")
test <- csm_simulate(5000, pt, sp_true)
for (tg in c("pleiotropy", "replication")) {
  res <- csm_test(fit, test$Z, tg, t = 1, q = 0.1, adjust = TRUE)
  print(res)
  nul <- csm_null_components(test$H, tg, 1)
  fdp <- if (res$n_rejected) mean(nul[res$rejected]) else 0
  check(fdp < 0.2, paste("realized FDP reasonable for", tg))
  check(res$adjusted$n_rejected <= res$n_rejected, paste("adjustment rejects no more for", tg))
}

## 7. Step-up rule against a brute-force search.
set.seed(7)
for (rep in 1:50) {
  s <- round(runif(40)^3, 2)
  su <- csm_stepup(s, 0.1)
  srt <- sort(s)
  r <- max(c(0, which(cumsum(srt) / seq_along(srt) <= 0.1)))
  if (su$n_rejected != r || (r > 0 && sum(su$rejected) != r)) stop("FAILED: step-up")
}
cat("ok - step-up rule\n")

## 8. Monotone adjustment: brute force and decision mode.
set.seed(8)
Zt <- matrix(round(rnorm(900, sd = 2), 1), 300)
s <- runif(300)^2
for (tg in c("pleiotropy", "replication", "mediation")) {
  K <- 3
  if (tg == "mediation") Zt4 <- cbind(Zt, round(rnorm(300, sd = 2), 1)) else Zt4 <- Zt
  rp <- ns$.order_repr(Zt4, tg, 1L)
  brute <- vapply(seq_len(nrow(Zt4)), function(i) {
    j <- which(rp$group == rp$group[i] & rowSums(rp$w >= rep(rp$w[i, ], each = nrow(Zt4))) == ncol(Zt4))
    max(s[j])
  }, 0)
  a_all <- csm_adjust(Zt4, s, tg, 1, 0.3, all = TRUE)
  a_dec <- csm_adjust(Zt4, s, tg, 1, 0.3)
  check(max(abs(a_all$adjusted - brute)) == 0, paste("adjustment equals brute force:", tg))
  check(identical(a_all$rejected, a_dec$rejected), paste("decision mode equals full adjustment:", tg))
  check(csm_incongruity(Zt4, a_all$adjusted, tg, 1)$violations == 0, paste("no incongruity after adjustment:", tg))
}

## 8b. Warm start on a subset and threads give the same maximum.
fit_w <- csm_fit(sim$Z, csm_spec(3, M = 1, b_mu = 8), starts = "null",
                 control = csm_control(warm_n = 2000, cores = 2, penalty = 0, weight_penalty = 0))
fit_c <- csm_fit(sim$Z, csm_spec(3, M = 1, b_mu = 8), starts = "null", control = csm_control(warm_n = 0, penalty = 0, weight_penalty = 0))
check(fit_w$converged && !is.na(fit_w$runs$warm_iterations), "warm start used and converged")
check(abs(fit_w$loglik - fit_c$loglik) < 1e-4, "warm start reaches the same log-likelihood")

## 8c. Penalized fit: objective nondecreasing, magnitudes less dispersed.
fit_p <- csm_fit(sim$Z, csm_spec(3, M = 1, b_mu = 8), starts = "null", control = csm_control(penalty = 50, weight_penalty = 0))
trp <- fit_p$traces[[1]]
check(fit_p$converged && all(diff(trp) >= -1e-10 * (1 + abs(trp[-1]))), "penalized objective is nondecreasing")
check(ns$.mag_penalty(fit_p$params$mu, fit_p$spec) <= ns$.mag_penalty(fit_c$params$mu, fit_c$spec) + 1e-8,
      "penalty reduces magnitude dispersion within support sizes")

## 8d. Default (both penalties): objective nondecreasing, penalties reported, selection by objective.
fit_d <- csm_fit(sim$Z, csm_spec(3, M = 1, b_mu = 8), starts = c("null", "equal"))
trd <- fit_d$traces[[fit_d$selected]]
check(fit_d$converged && all(diff(trd) >= -1e-10 * (1 + abs(trd[-1]))), "default penalized objective is nondecreasing")
ok_d <- fit_d$runs$status == "converged"
check(fit_d$objective == max(fit_d$runs$objective[ok_d]) && fit_d$objective <= fit_d$loglik,
      "selection by the penalized objective")
check(all(fit_d$params$varpi > 0), "weight penalty keeps every weight positive")

## 9. Failed fits make no discoveries.
fit_fail <- csm_fit(sim$Z, csm_spec(3, M = 1, b_mu = 8), starts = "equal", control = csm_control(max_iter = 2))
res <- csm_test(fit_fail, test$Z, "pleiotropy", 1)
check(!fit_fail$converged && res$status == "failed" && res$n_rejected == 0L, "failed fit makes no discoveries")

## 10. Round trip through the csmGmm list format.
L <- as_sun_lists(fit)
p_rt <- csm_params(fit$spec, pi = L$piList, mu = L$muList, Sigma = L$Sigma)
check(max(abs(p_rt$varpi - fit$params$varpi)) < 1e-12 && max(abs(p_rt$mu - fit$params$mu)) == 0,
      "csmGmm list format round trip")
cat("All tests passed.\n")

## M structures and BIC selection
M4 <- csm_structure(4, 2, 1)
check(length(M4) == 16 && M4[1] == 1 && all(M4[c(2, 3, 5, 9)] == 2) && all(M4[-c(1, 2, 3, 5, 9)] == 1),
      "csm_structure: two components for single supports only")
set.seed(21)
sp2 <- csm_spec(2, M = 1, b_mu = 12)
z2 <- rbind(csm_simulate(3000, csm_params(sp2, varpi = c(0.9, 0.04, 0.04, 0.02), mu = sp2$grp_A * 3), sp2)$Z,
            csm_simulate(3000, csm_params(sp2, varpi = c(0.9, 0.04, 0.04, 0.02), mu = sp2$grp_A * 6), sp2)$Z)
fs <- csm_select(z2, 2, structures = list(A = c(single = 1, multiple = 1), B = c(single = 2, multiple = 1)),
                 control = csm_control(max_iter = 500))
check(nrow(fs$selection) == 2 && sum(fs$selection$chosen) == 1 && fs$selection$df[2] > fs$selection$df[1] &&
      isTRUE(all.equal(min(fs$selection$delta_BIC), 0)), "csm_select returns a BIC table and one choice")
check(fs$selection$structure[fs$selection$chosen] == "B",
      "BIC prefers two components when single-support magnitudes are spread (3 and 6)")

# Default structure: two components for size-one supports, one otherwise
spd <- csm_spec(4)
stopifnot(identical(spd$M, csm_structure(4, 2, 1)), sum(spd$grp_size == 1) == 8, sum(spd$grp_size >= 2) == 11)
cat("default structure B: ok\n")

# fix_sigma: the correlation matrix stays at the identity
set.seed(5)
spf <- csm_spec(3, M = 1, b_mu = 8)
pf <- csm_params(spf, mu = spf$grp_A * 3, Sigma = matrix(0.4, 3, 3) + diag(0.6, 3),
                 varpi = c(0.9, rep(0.1 / (length(spf$grp_b) - 1), length(spf$grp_b) - 1)))
simf <- csm_simulate(3000, pf, spf)
ff <- csm_fit(simf$Z, spf, starts = c("null", "shared"), control = csm_control(fix_sigma = TRUE, penalty = 0, weight_penalty = 0))
fc <- csm_fit(simf$Z, spf, starts = c("null", "shared"), control = csm_control(penalty = 0, weight_penalty = 0))
stopifnot(ff$converged, isTRUE(all.equal(ff$params$Sigma, diag(3))), fc$params$Sigma[1, 2] > 0.2,
          attr(logLik(ff), "df") == attr(logLik(fc), "df") - 3, ff$loglik < fc$loglik)
cat("fix_sigma: ok\n")
R0 <- matrix(0.3, 3, 3) + diag(0.7, 3)
fp <- csm_fit(simf$Z, spf, starts = "null", control = csm_control(fix_sigma = R0, penalty = 0, weight_penalty = 0))
stopifnot(fp$converged, isTRUE(all.equal(fp$params$Sigma, R0)), fp$loglik > ff$loglik)
bad <- tryCatch(csm_control(fix_sigma = matrix(2, 3, 3)), error = function(e) "err")
stopifnot(identical(bad, "err"))
cat("fix_sigma (plug-in matrix): ok\n")
