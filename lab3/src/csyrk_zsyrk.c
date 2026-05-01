#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <string.h>
#include <omp.h>
#include "csyrk_zsyrk.h"

#define SYRC_SUCCESS            0
#define SYRC_ERR_ORDER          1
#define SYRC_ERR_UPLO           2
#define SYRC_ERR_TRANS          3
#define SYRC_ERR_N              4
#define SYRC_ERR_K              5
#define SYRC_ERR_LDA            6
#define SYRC_ERR_LDC            7
#define SYRC_ERR_NULL_PTR       8

static const char* syrc_error_messages[] = {
    "Успешное выполнение",
    "Недопустимое значение Order",
    "Недопустимое значение Uplo",
    "Недопустимое значение Trans",
    "Недопустимое значение N",
    "Недопустимое значение K",
    "Недопустимое значение LDA",
    "Недопустимое значение LDC",
    "Передан NULL указатель"
};

const char* syrc_get_error_message(int code) {
    if (code >= 0 && code <= SYRC_ERR_NULL_PTR) {
        return syrc_error_messages[code];
    }
    return "Неизвестная ошибка";
}

int syrc_validate_params(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo,
                         enum CBLAS_TRANSPOSE Trans, int N, int K,
                         int lda, int ldc) {
    if (Order != CblasRowMajor && Order != CblasColMajor) {
        fprintf(stderr, "ОШИБКА: %s\n", syrc_error_messages[SYRC_ERR_ORDER]);
        return SYRC_ERR_ORDER;
    }
    if (Uplo != CblasUpper && Uplo != CblasLower) {
        fprintf(stderr, "ОШИБКА: %s\n", syrc_error_messages[SYRC_ERR_UPLO]);
        return SYRC_ERR_UPLO;
    }
    if (Trans != CblasNoTrans && Trans != CblasTrans) {
        fprintf(stderr, "ОШИБКА: %s\n", syrc_error_messages[SYRC_ERR_TRANS]);
        return SYRC_ERR_TRANS;
    }
    if (N < 0) {
        fprintf(stderr, "ОШИБКА: %s (N=%d)\n", syrc_error_messages[SYRC_ERR_N], N);
        return SYRC_ERR_N;
    }
    if (K < 0) {
        fprintf(stderr, "ОШИБКА: %s (K=%d)\n", syrc_error_messages[SYRC_ERR_K], K);
        return SYRC_ERR_K;
    }
    
    // Исправленная проверка LDA с учетом порядка матрицы
    int expected_lda;
    if (Order == CblasRowMajor) {
        expected_lda = (Trans == CblasNoTrans) ? K : N;
    } else {
        expected_lda = (Trans == CblasNoTrans) ? N : K;
    }
    
    if (lda < expected_lda) {
        fprintf(stderr, "ОШИБКА: %s (LDA=%d)\n", syrc_error_messages[SYRC_ERR_LDA], lda);
        return SYRC_ERR_LDA;
    }
    if (ldc < N) {
        fprintf(stderr, "ОШИБКА: %s (LDC=%d)\n", syrc_error_messages[SYRC_ERR_LDC], ldc);
        return SYRC_ERR_LDC;
    }
    return SYRC_SUCCESS;
}

int cblas_csyrk_safe(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo,
                     enum CBLAS_TRANSPOSE Trans, int N, int K,
                     const void* alpha, const void* A, int lda,
                     const void* beta, void* C, int ldc) {
    if (!A || !C || !alpha || !beta) {
        fprintf(stderr, "ОШИБКА: %s\n", syrc_error_messages[SYRC_ERR_NULL_PTR]);
        return SYRC_ERR_NULL_PTR;
    }
    
    int err = syrc_validate_params(Order, Uplo, Trans, N, K, lda, ldc);
    if (err != SYRC_SUCCESS) return err;
    
    const cblas_complex_float* Alpha = (const cblas_complex_float*)alpha;
    const cblas_complex_float* Beta = (const cblas_complex_float*)beta;
    const cblas_complex_float* a = (const cblas_complex_float*)A;
    cblas_complex_float* c = (cblas_complex_float*)C;
    
    int i, j, l;
    
    if (Order == CblasRowMajor) {
        if (Trans == CblasNoTrans) {
            #pragma omp parallel for private(i, j, l) schedule(dynamic)
            for (i = 0; i < N; i++) {
                int j_start = (Uplo == CblasUpper) ? i : 0;
                int j_end = (Uplo == CblasUpper) ? N : i + 1;
                for (j = j_start; j < j_end; j++) {
                    cblas_complex_float temp = {0.0f, 0.0f};
                    for (l = 0; l < K; l++) {
                        float ar = a[i * lda + l].real;
                        float ai = a[i * lda + l].imag;
                        float br = a[j * lda + l].real;
                        float bi = a[j * lda + l].imag;
                        temp.real += ar * br - ai * bi;
                        temp.imag += ar * bi + ai * br;
                    }
                    float tr = temp.real, ti = temp.imag;
                    float ar = Alpha->real, ai = Alpha->imag;
                    float br = Beta->real, bi = Beta->imag;
                    float cr = c[i * ldc + j].real, ci = c[i * ldc + j].imag;
                    c[i * ldc + j].real = (ar * tr - ai * ti) + (br * cr - bi * ci);
                    c[i * ldc + j].imag = (ar * ti + ai * tr) + (br * ci + bi * cr);
                    if (i != j) c[j * ldc + i] = c[i * ldc + j];
                }
            }
        } else {
            #pragma omp parallel for private(i, j, l) schedule(dynamic)
            for (i = 0; i < N; i++) {
                int j_start = (Uplo == CblasUpper) ? i : 0;
                int j_end = (Uplo == CblasUpper) ? N : i + 1;
                for (j = j_start; j < j_end; j++) {
                    cblas_complex_float temp = {0.0f, 0.0f};
                    for (l = 0; l < K; l++) {
                        float ar = a[l * lda + i].real;
                        float ai = a[l * lda + i].imag;
                        float br = a[l * lda + j].real;
                        float bi = a[l * lda + j].imag;
                        temp.real += ar * br - ai * bi;
                        temp.imag += ar * bi + ai * br;
                    }
                    float tr = temp.real, ti = temp.imag;
                    float ar = Alpha->real, ai = Alpha->imag;
                    float br = Beta->real, bi = Beta->imag;
                    float cr = c[i * ldc + j].real, ci = c[i * ldc + j].imag;
                    c[i * ldc + j].real = (ar * tr - ai * ti) + (br * cr - bi * ci);
                    c[i * ldc + j].imag = (ar * ti + ai * tr) + (br * ci + bi * cr);
                    if (i != j) c[j * ldc + i] = c[i * ldc + j];
                }
            }
        }
    } else {
        if (Trans == CblasNoTrans) {
            #pragma omp parallel for private(i, j, l) schedule(dynamic)
            for (j = 0; j < N; j++) {
                int i_start = (Uplo == CblasUpper) ? 0 : j;
                int i_end = (Uplo == CblasUpper) ? j + 1 : N;
                for (i = i_start; i < i_end; i++) {
                    cblas_complex_float temp = {0.0f, 0.0f};
                    for (l = 0; l < K; l++) {
                        float ar = a[i + l * lda].real;
                        float ai = a[i + l * lda].imag;
                        float br = a[j + l * lda].real;
                        float bi = a[j + l * lda].imag;
                        temp.real += ar * br - ai * bi;
                        temp.imag += ar * bi + ai * br;
                    }
                    float tr = temp.real, ti = temp.imag;
                    float ar = Alpha->real, ai = Alpha->imag;
                    float br = Beta->real, bi = Beta->imag;
                    float cr = c[i + j * ldc].real, ci = c[i + j * ldc].imag;
                    c[i + j * ldc].real = (ar * tr - ai * ti) + (br * cr - bi * ci);
                    c[i + j * ldc].imag = (ar * ti + ai * tr) + (br * ci + bi * cr);
                    if (i != j) c[j + i * ldc] = c[i + j * ldc];
                }
            }
        } else {
            #pragma omp parallel for private(i, j, l) schedule(dynamic)
            for (j = 0; j < N; j++) {
                int i_start = (Uplo == CblasUpper) ? 0 : j;
                int i_end = (Uplo == CblasUpper) ? j + 1 : N;
                for (i = i_start; i < i_end; i++) {
                    cblas_complex_float temp = {0.0f, 0.0f};
                    for (l = 0; l < K; l++) {
                        float ar = a[l + i * lda].real;
                        float ai = a[l + i * lda].imag;
                        float br = a[l + j * lda].real;
                        float bi = a[l + j * lda].imag;
                        temp.real += ar * br - ai * bi;
                        temp.imag += ar * bi + ai * br;
                    }
                    float tr = temp.real, ti = temp.imag;
                    float ar = Alpha->real, ai = Alpha->imag;
                    float br = Beta->real, bi = Beta->imag;
                    float cr = c[i + j * ldc].real, ci = c[i + j * ldc].imag;
                    c[i + j * ldc].real = (ar * tr - ai * ti) + (br * cr - bi * ci);
                    c[i + j * ldc].imag = (ar * ti + ai * tr) + (br * ci + bi * cr);
                    if (i != j) c[j + i * ldc] = c[i + j * ldc];
                }
            }
        }
    }
    return SYRC_SUCCESS;
}

int cblas_zsyrk_safe(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo,
                     enum CBLAS_TRANSPOSE Trans, int N, int K,
                     const void* alpha, const void* A, int lda,
                     const void* beta, void* C, int ldc) {
    if (!A || !C || !alpha || !beta) {
        fprintf(stderr, "ОШИБКА: %s\n", syrc_error_messages[SYRC_ERR_NULL_PTR]);
        return SYRC_ERR_NULL_PTR;
    }
    
    int err = syrc_validate_params(Order, Uplo, Trans, N, K, lda, ldc);
    if (err != SYRC_SUCCESS) return err;
    
    const cblas_complex_double* Alpha = (const cblas_complex_double*)alpha;
    const cblas_complex_double* Beta = (const cblas_complex_double*)beta;
    const cblas_complex_double* a = (const cblas_complex_double*)A;
    cblas_complex_double* c = (cblas_complex_double*)C;
    
    int i, j, l;
    
    if (Order == CblasRowMajor) {
        if (Trans == CblasNoTrans) {
            #pragma omp parallel for private(i, j, l) schedule(dynamic)
            for (i = 0; i < N; i++) {
                int j_start = (Uplo == CblasUpper) ? i : 0;
                int j_end = (Uplo == CblasUpper) ? N : i + 1;
                for (j = j_start; j < j_end; j++) {
                    cblas_complex_double temp = {0.0, 0.0};
                    for (l = 0; l < K; l++) {
                        double ar = a[i * lda + l].real;
                        double ai = a[i * lda + l].imag;
                        double br = a[j * lda + l].real;
                        double bi = a[j * lda + l].imag;
                        temp.real += ar * br - ai * bi;
                        temp.imag += ar * bi + ai * br;
                    }
                    double tr = temp.real, ti = temp.imag;
                    double ar = Alpha->real, ai = Alpha->imag;
                    double br = Beta->real, bi = Beta->imag;
                    double cr = c[i * ldc + j].real, ci = c[i * ldc + j].imag;
                    c[i * ldc + j].real = (ar * tr - ai * ti) + (br * cr - bi * ci);
                    c[i * ldc + j].imag = (ar * ti + ai * tr) + (br * ci + bi * cr);
                    if (i != j) c[j * ldc + i] = c[i * ldc + j];
                }
            }
        } else {
            #pragma omp parallel for private(i, j, l) schedule(dynamic)
            for (i = 0; i < N; i++) {
                int j_start = (Uplo == CblasUpper) ? i : 0;
                int j_end = (Uplo == CblasUpper) ? N : i + 1;
                for (j = j_start; j < j_end; j++) {
                    cblas_complex_double temp = {0.0, 0.0};
                    for (l = 0; l < K; l++) {
                        double ar = a[l * lda + i].real;
                        double ai = a[l * lda + i].imag;
                        double br = a[l * lda + j].real;
                        double bi = a[l * lda + j].imag;
                        temp.real += ar * br - ai * bi;
                        temp.imag += ar * bi + ai * br;
                    }
                    double tr = temp.real, ti = temp.imag;
                    double ar = Alpha->real, ai = Alpha->imag;
                    double br = Beta->real, bi = Beta->imag;
                    double cr = c[i * ldc + j].real, ci = c[i * ldc + j].imag;
                    c[i * ldc + j].real = (ar * tr - ai * ti) + (br * cr - bi * ci);
                    c[i * ldc + j].imag = (ar * ti + ai * tr) + (br * ci + bi * cr);
                    if (i != j) c[j * ldc + i] = c[i * ldc + j];
                }
            }
        }
    } else {
        if (Trans == CblasNoTrans) {
            #pragma omp parallel for private(i, j, l) schedule(dynamic)
            for (j = 0; j < N; j++) {
                int i_start = (Uplo == CblasUpper) ? 0 : j;
                int i_end = (Uplo == CblasUpper) ? j + 1 : N;
                for (i = i_start; i < i_end; i++) {
                    cblas_complex_double temp = {0.0, 0.0};
                    for (l = 0; l < K; l++) {
                        double ar = a[i + l * lda].real;
                        double ai = a[i + l * lda].imag;
                        double br = a[j + l * lda].real;
                        double bi = a[j + l * lda].imag;
                        temp.real += ar * br - ai * bi;
                        temp.imag += ar * bi + ai * br;
                    }
                    double tr = temp.real, ti = temp.imag;
                    double ar = Alpha->real, ai = Alpha->imag;
                    double br = Beta->real, bi = Beta->imag;
                    double cr = c[i + j * ldc].real, ci = c[i + j * ldc].imag;
                    c[i + j * ldc].real = (ar * tr - ai * ti) + (br * cr - bi * ci);
                    c[i + j * ldc].imag = (ar * ti + ai * tr) + (br * ci + bi * cr);
                    if (i != j) c[j + i * ldc] = c[i + j * ldc];
                }
            }
        } else {
            #pragma omp parallel for private(i, j, l) schedule(dynamic)
            for (j = 0; j < N; j++) {
                int i_start = (Uplo == CblasUpper) ? 0 : j;
                int i_end = (Uplo == CblasUpper) ? j + 1 : N;
                for (i = i_start; i < i_end; i++) {
                    cblas_complex_double temp = {0.0, 0.0};
                    for (l = 0; l < K; l++) {
                        double ar = a[l + i * lda].real;
                        double ai = a[l + i * lda].imag;
                        double br = a[l + j * lda].real;
                        double bi = a[l + j * lda].imag;
                        temp.real += ar * br - ai * bi;
                        temp.imag += ar * bi + ai * br;
                    }
                    double tr = temp.real, ti = temp.imag;
                    double ar = Alpha->real, ai = Alpha->imag;
                    double br = Beta->real, bi = Beta->imag;
                    double cr = c[i + j * ldc].real, ci = c[i + j * ldc].imag;
                    c[i + j * ldc].real = (ar * tr - ai * ti) + (br * cr - bi * ci);
                    c[i + j * ldc].imag = (ar * ti + ai * tr) + (br * ci + bi * cr);
                    if (i != j) c[j + i * ldc] = c[i + j * ldc];
                }
            }
        }
    }
    return SYRC_SUCCESS;
}

void cblas_csyrk(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo,
                 enum CBLAS_TRANSPOSE Trans, int N, int K,
                 const void* alpha, const void* A, int lda,
                 const void* beta, void* C, int ldc) {
    int err = cblas_csyrk_safe(Order, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
    if (err != SYRC_SUCCESS) {
        fprintf(stderr, "КРИТИЧЕСКАЯ ОШИБКА: %s\n", syrc_get_error_message(err));
        exit(EXIT_FAILURE);
    }
}

void cblas_zsyrk(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo,
                 enum CBLAS_TRANSPOSE Trans, int N, int K,
                 const void* alpha, const void* A, int lda,
                 const void* beta, void* C, int ldc) {
    int err = cblas_zsyrk_safe(Order, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
    if (err != SYRC_SUCCESS) {
        fprintf(stderr, "КРИТИЧЕСКАЯ ОШИБКА: %s\n", syrc_get_error_message(err));
        exit(EXIT_FAILURE);
    }
}