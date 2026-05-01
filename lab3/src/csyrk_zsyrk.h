#ifndef CSYRK_ZSYRK_H
#define CSYRK_ZSYRK_H

#include <cblas.h>

typedef struct {
    float  real;
    float  imag;
} cblas_complex_float;

typedef struct {
    double real;
    double imag;
} cblas_complex_double;

#define SYRC_SUCCESS            0
#define SYRC_ERR_ORDER          1
#define SYRC_ERR_UPLO           2
#define SYRC_ERR_TRANS          3
#define SYRC_ERR_N              4
#define SYRC_ERR_K              5
#define SYRC_ERR_LDA            6
#define SYRC_ERR_LDC            7
#define SYRC_ERR_NULL_PTR       8

const char* syrc_get_error_message(int code);

int syrc_validate_params(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo,
                         enum CBLAS_TRANSPOSE Trans, int N, int K,
                         int lda, int ldc);

int cblas_csyrk_safe(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo,
                     enum CBLAS_TRANSPOSE Trans, int N, int K,
                     const void* alpha, const void* A, int lda,
                     const void* beta, void* C, int ldc);

int cblas_zsyrk_safe(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo,
                     enum CBLAS_TRANSPOSE Trans, int N, int K,
                     const void* alpha, const void* A, int lda,
                     const void* beta, void* C, int ldc);

#endif