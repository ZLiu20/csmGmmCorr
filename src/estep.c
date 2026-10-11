/* Paired E-step kernel for correlated csmGmm models. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int n, K, P, groups, coordinates, mode, first, last;
    const double *x, *precision, *log_constant, *row_constant, *group_mu;
    const int *self, *null_pair, *group_start, *group_size;
    const int *group_offset, *group_coordinate;
    double *log_density, *local_fdr, *pair_mass, *sufficient;
    int ok;
} worker_data;

static int sign_at(int pattern, int coordinate) {
    return (pattern & (1 << coordinate)) ? -1 : 1;
}

static void *run_worker(void *argument) {
    worker_data *job = (worker_data *) argument;
    const size_t pair_n = (size_t) (job->P > 0 ? job->P : 1);
    const size_t coord_n = (size_t) (job->coordinates > 0 ? job->coordinates : 1);
    double *plus = (double *) malloc(pair_n * sizeof(double));
    double *minus = (double *) malloc(pair_n * sizeof(double));
    double *y = (double *) malloc((size_t) job->K * sizeof(double));
    double *expected_sign = (double *) malloc(coord_n * sizeof(double));

    if (plus == NULL || minus == NULL || y == NULL || expected_sign == NULL) {
        free(plus); free(minus); free(y); free(expected_sign);
        job->ok = 0;
        return NULL;
    }

    for (int i = job->first; i < job->last; ++i) {
        for (int k = 0; k < job->K; ++k) {
            double total = 0.0;
            for (int l = 0; l < job->K; ++l) {
                total += job->precision[k + (size_t) l * job->K] *
                         job->x[i + (size_t) l * job->n];
            }
            y[k] = total;
        }

        double maximum = R_NegInf;
        for (int g = 0; g < job->groups; ++g) {
            const int size = job->group_size[g];
            const int start = job->group_start[g];
            const int offset = job->group_offset[g];
            const int count = size == 0 ? 1 : (1 << (size - 1));
            for (int pattern = 0; pattern < count; ++pattern) {
                const int p = start + pattern;
                double linear = 0.0;
                for (int q = 0; q < size; ++q) {
                    const int coordinate = job->group_coordinate[offset + q];
                    linear += sign_at(pattern, q) * job->group_mu[offset + q] *
                              y[coordinate];
                }
                plus[p] = job->log_constant[p] + linear;
                minus[p] = job->self[p] ? R_NegInf : job->log_constant[p] - linear;
                if (plus[p] > maximum) maximum = plus[p];
                if (minus[p] > maximum) maximum = minus[p];
            }
        }

        double denominator = 0.0;
        double null_mass = 0.0;
        for (int p = 0; p < job->P; ++p) {
            plus[p] = exp(plus[p] - maximum);
            minus[p] = job->self[p] ? 0.0 : exp(minus[p] - maximum);
            const double mass = plus[p] + minus[p];
            denominator += mass;
            if (job->mode == 2 && job->null_pair[p]) null_mass += mass;
        }

        job->log_density[i] = maximum + log(denominator) + job->row_constant[i];
        if (job->mode == 2) job->local_fdr[i] = null_mass / denominator;

        if (job->mode == 1) {
            const double inverse = 1.0 / denominator;
            for (int p = 0; p < job->P; ++p) {
                job->pair_mass[p] += (plus[p] + minus[p]) * inverse;
            }
            memset(expected_sign, 0, coord_n * sizeof(double));
            for (int g = 0; g < job->groups; ++g) {
                const int size = job->group_size[g];
                if (size == 0) continue;
                const int start = job->group_start[g];
                const int offset = job->group_offset[g];
                const int count = 1 << (size - 1);
                for (int pattern = 0; pattern < count; ++pattern) {
                    const int p = start + pattern;
                    const double contrast = (plus[p] - minus[p]) * inverse;
                    for (int q = 0; q < size; ++q) {
                        expected_sign[offset + q] += sign_at(pattern, q) * contrast;
                    }
                }
            }
            for (int r = 0; r < job->coordinates; ++r) {
                for (int k = 0; k < job->K; ++k) {
                    job->sufficient[r + (size_t) k * job->coordinates] +=
                        expected_sign[r] * job->x[i + (size_t) k * job->n];
                }
            }
        }
    }

    free(plus); free(minus); free(y); free(expected_sign);
    return NULL;
}

SEXP csm_kernel(SEXP x_, SEXP precision_, SEXP log_constant_, SEXP self_,
                SEXP row_constant_, SEXP mode_, SEXP null_, SEXP group_start_,
                SEXP group_size_, SEXP group_offset_, SEXP group_coordinate_,
                SEXP group_mu_, SEXP threads_) {
    const int n = nrows(x_);
    const int K = ncols(x_);
    const int P = LENGTH(log_constant_);
    const int groups = LENGTH(group_start_);
    const int coordinates = LENGTH(group_coordinate_);
    const int mode = asInteger(mode_);

    if (!isReal(x_) || !isReal(precision_) || !isReal(log_constant_) ||
        !isLogical(self_) || !isReal(row_constant_) || !isInteger(group_start_) ||
        !isInteger(group_size_) || !isInteger(group_offset_) ||
        !isInteger(group_coordinate_) || !isReal(group_mu_)) {
        error("csm_kernel: incompatible storage mode");
    }
    if (nrows(precision_) != K || ncols(precision_) != K ||
        LENGTH(self_) != P || LENGTH(row_constant_) != n ||
        LENGTH(group_size_) != groups || LENGTH(group_offset_) != groups ||
        LENGTH(group_mu_) != coordinates) {
        error("csm_kernel: incompatible dimensions");
    }
    if (mode < 0 || mode > 2) error("csm_kernel: invalid mode");
    if (mode == 2 && (!isLogical(null_) || LENGTH(null_) != P)) {
        error("csm_kernel: null mask has the wrong length");
    }

    int covered = 0;
    for (int g = 0; g < groups; ++g) {
        const int size = INTEGER(group_size_)[g];
        const int count = size == 0 ? 1 : (1 << (size - 1));
        if (INTEGER(group_start_)[g] != covered ||
            INTEGER(group_offset_)[g] < 0 || size < 0) {
            error("csm_kernel: invalid group layout");
        }
        covered += count;
    }
    if (covered != P) error("csm_kernel: pair groups do not cover the components");

    int thread_count = asInteger(threads_);
    if (thread_count == NA_INTEGER || thread_count < 1) thread_count = 1;
    if (n > 0 && thread_count > n) thread_count = n;
    if (thread_count < 1) thread_count = 1;

    SEXP log_density_ = PROTECT(allocVector(REALSXP, n));
    SEXP pair_mass_ = PROTECT(allocVector(REALSXP, mode == 1 ? P : 0));
    SEXP sufficient_ = PROTECT(allocMatrix(REALSXP,
                                           mode == 1 ? coordinates : 0,
                                           mode == 1 ? K : 0));
    SEXP local_fdr_ = PROTECT(allocVector(REALSXP, mode == 2 ? n : 0));
    if (mode == 1) {
        memset(REAL(pair_mass_), 0, (size_t) P * sizeof(double));
        memset(REAL(sufficient_), 0,
               (size_t) coordinates * K * sizeof(double));
    }

    worker_data *jobs = (worker_data *) calloc((size_t) thread_count,
                                                sizeof(worker_data));
    pthread_t *threads = thread_count > 1 ?
        (pthread_t *) calloc((size_t) thread_count, sizeof(pthread_t)) : NULL;
    int *started = thread_count > 1 ?
        (int *) calloc((size_t) thread_count, sizeof(int)) : NULL;
    if (jobs == NULL || (thread_count > 1 && (threads == NULL || started == NULL))) {
        free(jobs); free(threads); free(started);
        UNPROTECT(4);
        error("csm_kernel: out of memory");
    }

    int allocation_ok = 1;
    for (int t = 0; t < thread_count; ++t) {
        worker_data *job = jobs + t;
        job->n = n; job->K = K; job->P = P; job->groups = groups;
        job->coordinates = coordinates; job->mode = mode;
        job->first = (int) (((long long) n * t) / thread_count);
        job->last = (int) (((long long) n * (t + 1)) / thread_count);
        job->x = REAL(x_); job->precision = REAL(precision_);
        job->log_constant = REAL(log_constant_); job->row_constant = REAL(row_constant_);
        job->group_mu = REAL(group_mu_); job->self = LOGICAL(self_);
        job->null_pair = mode == 2 ? LOGICAL(null_) : NULL;
        job->group_start = INTEGER(group_start_); job->group_size = INTEGER(group_size_);
        job->group_offset = INTEGER(group_offset_);
        job->group_coordinate = INTEGER(group_coordinate_);
        job->log_density = REAL(log_density_);
        job->local_fdr = mode == 2 ? REAL(local_fdr_) : NULL;
        job->pair_mass = mode == 1 ?
            (double *) calloc((size_t) (P > 0 ? P : 1), sizeof(double)) : NULL;
        job->sufficient = mode == 1 ?
            (double *) calloc((size_t) (coordinates > 0 ? coordinates : 1) * K,
                              sizeof(double)) : NULL;
        job->ok = mode != 1 || (job->pair_mass != NULL && job->sufficient != NULL);
        if (!job->ok) allocation_ok = 0;
    }

    if (allocation_ok) {
        if (thread_count == 1) {
            run_worker(jobs);
        } else {
            for (int t = 0; t < thread_count; ++t) {
                started[t] = pthread_create(threads + t, NULL, run_worker, jobs + t) == 0;
                if (!started[t]) run_worker(jobs + t);
            }
            for (int t = 0; t < thread_count; ++t) {
                if (started[t]) pthread_join(threads[t], NULL);
            }
        }
    }

    for (int t = 0; t < thread_count; ++t) {
        if (!jobs[t].ok) allocation_ok = 0;
        if (mode == 1 && jobs[t].pair_mass != NULL && jobs[t].sufficient != NULL) {
            for (int p = 0; p < P; ++p) REAL(pair_mass_)[p] += jobs[t].pair_mass[p];
            for (int index = 0; index < coordinates * K; ++index) {
                REAL(sufficient_)[index] += jobs[t].sufficient[index];
            }
        }
        free(jobs[t].pair_mass);
        free(jobs[t].sufficient);
    }
    free(jobs); free(threads); free(started);
    if (!allocation_ok) {
        UNPROTECT(4);
        error("csm_kernel: out of memory");
    }

    SEXP result = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(result, 0, log_density_);
    SET_VECTOR_ELT(result, 1, pair_mass_);
    SET_VECTOR_ELT(result, 2, sufficient_);
    SET_VECTOR_ELT(result, 3, local_fdr_);
    UNPROTECT(5);
    return result;
}

static const R_CallMethodDef call_methods[] = {
    {"csm_kernel", (DL_FUNC) &csm_kernel, 13},
    {NULL, NULL, 0}
};

void R_init_csmGmmCorr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
