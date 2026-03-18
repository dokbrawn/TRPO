#include "syrk.h"

#include <stddef.h>

int syrk_validate(CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE trans,
                  int n, int k, int lda, int ldc) {
    if ((order != CblasRowMajor && order != CblasColMajor) ||
        (uplo != CblasUpper && uplo != CblasLower) ||
        (trans != CblasNoTrans && trans != CblasTrans) ||
        n < 0 || k < 0 || ldc < n) {
        return 0;
    }

    const int required_lda = (trans == CblasNoTrans) ? k : n;
    return lda >= (required_lda > 0 ? required_lda : 1);
}

static int should_update(CBLAS_UPLO uplo, int row, int col) {
    return (uplo == CblasUpper) ? (row <= col) : (row >= col);
}

static float load_a_f(CBLAS_ORDER order, CBLAS_TRANSPOSE trans,
                      const float *a, int lda, int row, int col) {
    if (order == CblasColMajor) {
        return (trans == CblasNoTrans) ? a[row + col * lda] : a[col + row * lda];
    }
    return (trans == CblasNoTrans) ? a[row * lda + col] : a[col * lda + row];
}

static double load_a_d(CBLAS_ORDER order, CBLAS_TRANSPOSE trans,
                       const double *a, int lda, int row, int col) {
    if (order == CblasColMajor) {
        return (trans == CblasNoTrans) ? a[row + col * lda] : a[col + row * lda];
    }
    return (trans == CblasNoTrans) ? a[row * lda + col] : a[col * lda + row];
}

static float *cell_f(CBLAS_ORDER order, float *c, int ldc, int row, int col) {
    return (order == CblasColMajor) ? &c[row + col * ldc] : &c[row * ldc + col];
}

static double *cell_d(CBLAS_ORDER order, double *c, int ldc, int row, int col) {
    return (order == CblasColMajor) ? &c[row + col * ldc] : &c[row * ldc + col];
}

void syrk_s(CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE trans,
            int n, int k, float alpha, const float *a, int lda,
            float beta, float *c, int ldc) {
    if (!syrk_validate(order, uplo, trans, n, k, lda, ldc) || a == NULL || c == NULL) {
        return;
    }

    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < n; ++col) {
            if (!should_update(uplo, row, col)) {
                continue;
            }

            float sum = 0.0f;
            for (int p = 0; p < k; ++p) {
                sum += load_a_f(order, trans, a, lda, row, p) *
                       load_a_f(order, trans, a, lda, col, p);
            }

            float *dst = cell_f(order, c, ldc, row, col);
            *dst = alpha * sum + beta * (*dst);
        }
    }
}

void syrk_d(CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE trans,
            int n, int k, double alpha, const double *a, int lda,
            double beta, double *c, int ldc) {
    if (!syrk_validate(order, uplo, trans, n, k, lda, ldc) || a == NULL || c == NULL) {
        return;
    }

    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < n; ++col) {
            if (!should_update(uplo, row, col)) {
                continue;
            }

            double sum = 0.0;
            for (int p = 0; p < k; ++p) {
                sum += load_a_d(order, trans, a, lda, row, p) *
                       load_a_d(order, trans, a, lda, col, p);
            }

            double *dst = cell_d(order, c, ldc, row, col);
            *dst = alpha * sum + beta * (*dst);
        }
    }
}
