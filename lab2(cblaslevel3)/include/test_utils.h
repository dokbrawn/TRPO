#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cblas.h>

#ifdef _MSC_VER
    #include <complex.h>
    typedef _Fcomplex cblas_complex_float;
    typedef _Dcomplex cblas_complex_double;
    #define CREAL(c) crealf(c)
    #define CIMAG(c) cimagf(c)
    #define DCREAL(c) creal(c)
    #define DCIMAG(c) cimag(c)
#else
    #include <complex.h>
    typedef float _Complex cblas_complex_float;
    typedef double _Complex cblas_complex_double;
    #define CREAL(c) crealf(c)
    #define CIMAG(c) cimagf(c)
    #define DCREAL(c) creal(c)
    #define DCIMAG(c) cimag(c)
#endif

static inline int float_eq(float a, float b, float eps) {
    return fabsf(a - b) < eps;
}

static inline int double_eq(double a, double b, double eps) {
    return fabs(a - b) < eps;
}

static inline int check_matrix_f(const float *got, const float *exp, int rows, int cols, int ld, float eps) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (!float_eq(got[i + j * ld], exp[i + j * ld], eps)) return 0;
    return 1;
}

static inline int check_matrix_d(const double *got, const double *exp, int rows, int cols, int ld, double eps) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (!double_eq(got[i + j * ld], exp[i + j * ld], eps)) return 0;
    return 1;
}

#define TEST_START(name) printf("▶ %-55s", name); fflush(stdout)
#define TEST_PASS() do { printf("✓ PASS\n"); return 1; } while(0)
#define TEST_FAIL(msg) do { printf("✗ FAIL: %s\n", msg); return 0; } while(0)
#define ASSERT(cond, msg) if (!(cond)) { TEST_FAIL(msg); }

#ifdef USE_THREADS
    #include <openblas_config.h>
    static inline void set_threads(int n) {
    #ifdef OPENBLAS_THREADABLE
        openblas_set_num_threads(n);
        printf(" [threads: %d]", openblas_get_num_threads());
    #endif
    }
#else
    static inline void set_threads(int n) { (void)n; }
#endif

#endif