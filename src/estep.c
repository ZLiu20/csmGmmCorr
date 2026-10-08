/* E-step and posterior kernel for csmGmmCorr.
 *
 * Components are handled in sign pairs (h, -h). For feature i and pair p the
 * log densities (plus log weights) are colc[p] + g and colc[p] - g, with
 * g = x_i' Omega m_p; the "-" member of a self-paired (null) component is
 * absent. With W_p = exp(colc[p] - max colc) precomputed, both members need a
 * single exponential: exp(colc +- g - mx) = W_p * exp(+-g) * exp(max colc - mx).
 *
 * Groups: the active pairs are stored group by group. Group j occupies pairs
 * gstart[j], ..., gstart[j] + 2^(s_j - 1) - 1 (one pair if s_j = 0), in the
 * order of the sign patterns of R's expand.grid(c(1, -1), ...): the pattern of
 * local index b has sign -1 on active coordinate i exactly when bit i of b is
 * set. With y = Omega x, the linear term of pair b is
 *   g_b = sum_i h_i mu_i y_{a_i},
 * enumerated in Gray-code order so that each step changes one sign and costs
 * O(1) instead of O(K).
 *
 * Features are split into contiguous blocks processed by POSIX threads, each
 * with private accumulators that are summed at the end. No R API is called
 * from the worker threads.
 *
 * mode 0: row log-likelihoods only
 * mode 1: also T[p] = sum_i (tau_+ + tau_-) and the group sums
 *         S[r, ] = sum_{p in group j} h_pq sum_i (tau_+ - tau_-) x_i,  r = goff[j] + q
 * mode 2: also lfdr[i] = posterior mass of pairs flagged in `nul`
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>

#define GMAX 700.0
#define GINC 300.0
/* A group is skipped for a feature when every member satisfies
 * colc + |g| < mx - PRUNE: each member then contributes less than
 * exp(-PRUNE) relative to the largest term (which contributes 1 to the
 * total), so all skipped members together change the row total by less than
 * 2 P exp(-PRUNE) ~ 1e-17 relative for P = 3e4, below double precision. */
#define PRUNE 50.0

typedef struct {
  /* shared, read-only */
  int n, K, P, NG, mode;
  const double *X, *Om, *colc, *rowc, *W, *gmu, *gcmax;
  const int *self, *nul, *gstart, *gsize, *goff, *gcoord;
  double cmax;
  /* row range for this worker */
  int i0, i1;
  /* outputs: per-row arrays are shared (disjoint rows); T and St are private */
  double *ll, *lf, *T, *St;   /* St: group sums S (ncoord x K, row-major) */
  int ncoord;
  /* private work space */
  double *g, *ep, *em, *x, *y, *cb, *sab, *eu, *ev, *f2, *f2i;
  char *keep;
} work_t;

static void *worker(void *arg) {
  work_t *w = (work_t *) arg;
  const int K = w->K, NG = w->NG, mode = w->mode, n = w->n;
  double *g = w->g, *ep = w->ep, *em = w->em, *x = w->x, *y = w->y;
  double *eu = w->eu, *ev = w->ev, *f2 = w->f2, *f2i = w->f2i, *sab = w->sab;
  char *keep = w->keep;
  double *vq = w->f2i + K;   /* K extra doubles allocated after f2i */
  for (int i = w->i0; i < w->i1; i++) {
    for (int k = 0; k < K; k++) x[k] = w->X[i + (size_t) k * n];
    for (int k = 0; k < K; k++) {
      double s = 0.0;
      for (int l = 0; l < K; l++) s += w->Om[k + (size_t) l * K] * x[l];
      y[k] = s;
    }
    /* Pass 1: per group, c_q = mu_q y_{a_q}, the bound sab = sum |c_q| on |g|,
     * and a lower bound for mx from the pair with all signs aligned. */
    double mlo = -INFINITY;
    for (int j = 0; j < NG; j++) {
      const int s = w->gsize[j], st = w->gstart[j];
      if (s == 0) { sab[j] = 0.0; if (w->colc[st] > mlo) mlo = w->colc[st]; continue; }
      const int *a = w->gcoord + w->goff[j];
      const double *mu = w->gmu + w->goff[j];
      double *c = w->cb + w->goff[j];
      double sa = 0.0;
      int bal = 0;
      const double clast = mu[s - 1] * y[a[s - 1]];
      for (int q = 0; q < s; q++) {
        c[q] = mu[q] * y[a[q]]; sa += fabs(c[q]);
        if (q < s - 1 && (c[q] < 0.0) != (clast < 0.0)) bal |= 1 << q;
      }
      sab[j] = sa;
      const double v = w->colc[st + bal] + sa;
      if (v > mlo) mlo = v;
    }
    /* Pass 2: linear terms and exp(+-g) for the groups that are kept. */
    double mx = -INFINITY;
    for (int j = 0; j < NG; j++) {
      const int s = w->gsize[j], st = w->gstart[j];
      /* lfdr (mode 2): null groups are always kept, so that very small local
       * false discovery rates keep their relative accuracy */
      keep[j] = (w->gcmax[j] + sab[j] >= mlo - PRUNE) || (mode == 2 && w->nul[st]);
      if (!keep[j]) continue;
      if (s == 0) {
        g[st] = 0.0; eu[st] = 1.0; ev[st] = 1.0;
        if (w->colc[st] > mx) mx = w->colc[st];
        continue;
      }
      const double *c = w->cb + w->goff[j];
      double gv = 0.0;
      for (int q = 0; q < s; q++) gv += c[q];
      const int np = 1 << (s - 1);
      /* exp(+-g) along the Gray code: a flip of sign q multiplies exp(g) by
       * exp(-+2 c_q). Used when every |g| <= sab < GINC, so no intermediate
       * value overflows or underflows; otherwise exp is evaluated directly. */
      const int inc = sab[j] < GINC;
      double u = 1.0, v = 1.0;
      if (inc) {
        for (int q = 0; q < s; q++) { f2[q] = exp(2.0 * c[q]); f2i[q] = 1.0 / f2[q]; }
        u = exp(gv); v = 1.0 / u;
      }
      int b = 0;
      for (int t = 0; t < np; t++) {
        if (t > 0) {
          const int q = __builtin_ctz((unsigned) t);   /* the bit that flips */
          b ^= (1 << q);
          if (b & (1 << q)) { gv -= 2.0 * c[q]; if (inc) { u *= f2i[q]; v *= f2[q]; } }
          else              { gv += 2.0 * c[q]; if (inc) { u *= f2[q]; v *= f2i[q]; } }
        }
        const int p = st + b;
        g[p] = gv;
        if (inc) { eu[p] = u; ev[p] = v; }
        else if (fabs(gv) < GMAX) { eu[p] = exp(gv); ev[p] = 1.0 / eu[p]; }
        else eu[p] = -1.0;                              /* direct evaluation below */
        const double a1 = w->colc[p] + fabs(gv);
        if (a1 > mx) mx = a1;
      }
    }
    /* Pass 3: member densities relative to exp(mx). */
    double tot = 0.0, totnull = 0.0;
    const double C = exp(w->cmax - mx);
    for (int j = 0; j < NG; j++) {
      if (!keep[j]) continue;
      const int st = w->gstart[j], np = w->gsize[j] ? 1 << (w->gsize[j] - 1) : 1;
      for (int p = st; p < st + np; p++) {
        double e1, e2 = 0.0;
        if (eu[p] >= 0.0) {
          const double wc = w->W[p] * C;
          e1 = wc * eu[p];
          if (!w->self[p]) e2 = wc * ev[p];
        } else {
          e1 = exp(w->colc[p] + g[p] - mx);
          if (!w->self[p]) e2 = exp(w->colc[p] - g[p] - mx);
        }
        ep[p] = e1;
        em[p] = e2;
        tot += e1 + e2;
        if (mode == 2 && w->nul[p]) totnull += e1 + e2;
      }
    }
    w->ll[i] = mx + log(tot) + w->rowc[i];
    if (mode == 1) {
      /* T per pair; for D only the group sums S_g[q, ] = sum_{p in g} h_pq D_p
       * are needed (M-step), accumulated as v_q x' with v_q = sum_b df_b h_bq
       * and h_bq = -1 exactly when bit q of b is set. */
      const double inv = 1.0 / tot;
      for (int j = 0; j < NG; j++) {
        if (!keep[j]) continue;
        const int s = w->gsize[j], st = w->gstart[j], np = s ? 1 << (s - 1) : 1;
        double sdf = 0.0, sabs = 0.0;
        for (int q = 0; q < s; q++) vq[q] = 0.0;
        for (int b = 0; b < np; b++) {
          const int p = st + b;
          const double sp = ep[p] + em[p];
          if (sp == 0.0) continue;
          w->T[p] += sp * inv;
          const double df = (ep[p] - em[p]) * inv;
          sdf += df; sabs += fabs(df);
          for (unsigned bb = (unsigned) b; bb; bb &= bb - 1) vq[__builtin_ctz(bb)] += df;
        }
        if (s == 0 || sabs <= 1e-17) continue;
        double *S = w->St + (size_t) w->goff[j] * K;
        for (int q = 0; q < s; q++) {
          const double v = sdf - 2.0 * vq[q];
          double *Sq = S + (size_t) q * K;
          for (int k = 0; k < K; k++) Sq[k] += v * x[k];
        }
      }
    } else if (mode == 2) {
      w->lf[i] = totnull / tot;
    }
  }
  return NULL;
}

SEXP csm_kernel(SEXP X_, SEXP Omega_, SEXP colc_, SEXP self_, SEXP rowc_, SEXP mode_, SEXP nul_,
                SEXP gstart_, SEXP gsize_, SEXP goff_, SEXP gcoord_, SEXP gmu_, SEXP threads_) {
  const int n = nrows(X_), K = ncols(X_), P = LENGTH(colc_), mode = asInteger(mode_);
  const int NG = LENGTH(gstart_);
  int nt = asInteger(threads_);
  if (nrows(Omega_) != K || ncols(Omega_) != K || LENGTH(self_) != P || LENGTH(rowc_) != n)
    error("csm_kernel: dimension mismatch");
  if (mode == 2 && LENGTH(nul_) != P) error("csm_kernel: null mask has the wrong length");
  if (nt == NA_INTEGER || nt < 1) nt = 1;
  if (nt > n) nt = n > 0 ? n : 1;

  SEXP ll_ = PROTECT(allocVector(REALSXP, n));
  SEXP T_ = PROTECT(allocVector(REALSXP, mode == 1 ? P : 0));
  const int ncoord = LENGTH(gcoord_);
  SEXP D_ = PROTECT(allocMatrix(REALSXP, mode == 1 ? ncoord : 0, mode == 1 ? K : 0));
  SEXP lf_ = PROTECT(allocVector(REALSXP, mode == 2 ? n : 0));

  const double *colc = REAL(colc_);
  double cmax = -INFINITY;
  for (int p = 0; p < P; p++) if (colc[p] > cmax) cmax = colc[p];
  double *W = (double *) R_alloc(P > 0 ? P : 1, sizeof(double));
  for (int p = 0; p < P; p++) W[p] = exp(colc[p] - cmax);
  /* largest log weight term within each group (for the pruning bound) */
  double *gcmax = (double *) R_alloc(NG > 0 ? NG : 1, sizeof(double));
  {
    const int *gs = INTEGER(gstart_), *gz = INTEGER(gsize_);
    for (int j = 0; j < NG; j++) {
      const int np = gz[j] ? 1 << (gz[j] - 1) : 1;
      double m = -INFINITY;
      for (int p = gs[j]; p < gs[j] + np; p++) if (colc[p] > m) m = colc[p];
      gcmax[j] = m;
    }
  }

  /* allocate all thread work space before starting threads */
  work_t *wk = (work_t *) R_alloc(nt, sizeof(work_t));
  memset(wk, 0, sizeof(work_t) * nt);
  int ok = 1;
  for (int t = 0; t < nt; t++) {
    work_t *w = wk + t;
    w->n = n; w->K = K; w->P = P; w->NG = NG; w->mode = mode;
    w->X = REAL(X_); w->Om = REAL(Omega_); w->colc = colc; w->rowc = REAL(rowc_);
    w->W = W; w->gmu = REAL(gmu_); w->gcmax = gcmax;
    w->self = LOGICAL(self_); w->nul = (mode == 2) ? LOGICAL(nul_) : NULL;
    w->gstart = INTEGER(gstart_); w->gsize = INTEGER(gsize_); w->goff = INTEGER(goff_);
    w->gcoord = INTEGER(gcoord_);
    w->cmax = cmax;
    w->i0 = (int) ((long long) n * t / nt);
    w->i1 = (int) ((long long) n * (t + 1) / nt);
    w->ll = REAL(ll_); w->lf = REAL(lf_);
    w->g = malloc(sizeof(double) * (P > 0 ? P : 1));
    w->ep = malloc(sizeof(double) * (P > 0 ? P : 1));
    w->em = malloc(sizeof(double) * (P > 0 ? P : 1));
    w->x = malloc(sizeof(double) * K);
    w->y = malloc(sizeof(double) * K);
    w->cb = malloc(sizeof(double) * (ncoord > 0 ? ncoord : 1));
    w->sab = malloc(sizeof(double) * (NG > 0 ? NG : 1));
    w->keep = malloc(NG > 0 ? NG : 1);
    w->eu = malloc(sizeof(double) * (P > 0 ? P : 1));
    w->ev = malloc(sizeof(double) * (P > 0 ? P : 1));
    w->f2 = malloc(sizeof(double) * K);
    w->f2i = malloc(sizeof(double) * 2 * K);
    w->T = (mode == 1) ? calloc(P > 0 ? P : 1, sizeof(double)) : NULL;
    w->St = (mode == 1) ? calloc((size_t) (ncoord > 0 ? ncoord : 1) * K, sizeof(double)) : NULL;
    w->ncoord = ncoord;
    if (!w->g || !w->ep || !w->em || !w->x || !w->y || !w->cb || !w->sab || !w->keep || !w->eu || !w->ev || !w->f2 || !w->f2i || (mode == 1 && (!w->T || !w->St))) ok = 0;
  }
  if (ok) {
    if (nt == 1) {
      worker(wk);
    } else {
      pthread_t *th = (pthread_t *) R_alloc(nt, sizeof(pthread_t));
      int *started = (int *) R_alloc(nt, sizeof(int));
      for (int t = 1; t < nt; t++) started[t] = (pthread_create(th + t, NULL, worker, wk + t) == 0);
      worker(wk);
      for (int t = 1; t < nt; t++) {
        if (started[t]) pthread_join(th[t], NULL);
        else worker(wk + t);                    /* thread not available: run inline */
      }
    }
    if (mode == 1) {
      double *T = REAL(T_), *D = REAL(D_);
      memset(T, 0, sizeof(double) * P);
      memset(D, 0, sizeof(double) * (size_t) ncoord * K);
      for (int t = 0; t < nt; t++) {
        for (int p = 0; p < P; p++) T[p] += wk[t].T[p];
        for (int r = 0; r < ncoord; r++)
          for (int k = 0; k < K; k++) D[r + (size_t) k * ncoord] += wk[t].St[(size_t) r * K + k];
      }
    }
  }
  for (int t = 0; t < nt; t++) {
    free(wk[t].g); free(wk[t].ep); free(wk[t].em); free(wk[t].x); free(wk[t].y); free(wk[t].cb); free(wk[t].sab); free(wk[t].keep);
    free(wk[t].eu); free(wk[t].ev); free(wk[t].f2); free(wk[t].f2i);
    free(wk[t].T); free(wk[t].St);
  }
  if (!ok) error("csm_kernel: out of memory");

  SEXP out = PROTECT(allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, ll_);
  SET_VECTOR_ELT(out, 1, T_);
  SET_VECTOR_ELT(out, 2, D_);
  SET_VECTOR_ELT(out, 3, lf_);
  UNPROTECT(5);
  return out;
}

static const R_CallMethodDef callMethods[] = {
  {"csm_kernel", (DL_FUNC) &csm_kernel, 13},
  {NULL, NULL, 0}
};

void R_init_csmGmmCorr(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
