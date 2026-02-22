#include "test_utils.h"

int test_ssymm_left_upper(void) {
    TEST_START("cblas_ssymm: Left+Upper");
    const int M = 2, N = 2;
    float A[4] = {2,1,1,2}, B[4] = {1,2,3,4}, C[4] = {0};
    float expected[4] = {4,5,10,11};
    cblas_ssymm(CblasColMajor, CblasLeft, CblasUpper, M, N, 1.0f, A, 2, B, 2, 0.0f, C, 2);
    ASSERT(check_matrix_f(C, expected, M, N, 2, 1e-5f), "SSYMM Left+Upper mismatch");
    TEST_PASS();
}

int test_ssymm_left_lower(void) {
    TEST_START("cblas_ssymm: Left+Lower");
    const int M = 2, N = 2;
    float A[4] = {2,1,1,2}, B[4] = {1,2,3,4}, C[4] = {0};
    float expected[4] = {4,5,10,11};
    cblas_ssymm(CblasColMajor, CblasLeft, CblasLower, M, N, 1.0f, A, 2, B, 2, 0.0f, C, 2);
    ASSERT(check_matrix_f(C, expected, M, N, 2, 1e-5f), "SSYMM Left+Lower mismatch");
    TEST_PASS();
}

int test_ssymm_right_upper(void) {
    TEST_START("cblas_ssymm: Right+Upper");
    const int M = 2, N = 2;
    float A[4] = {2,1,1,2}, B[4] = {1,2,3,4}, C[4] = {0};
    cblas_ssymm(CblasColMajor, CblasRight, CblasUpper, M, N, 1.0f, A, 2, B, 2, 0.0f, C, 2);
    ASSERT(C[0] != 0.0f, "SSYMM Right+Upper error");
    TEST_PASS();
}

int test_ssymm_right_lower(void) {
    TEST_START("cblas_ssymm: Right+Lower");
    const int M = 2, N = 2;
    float A[4] = {2,1,1,2}, B[4] = {1,2,3,4}, C[4] = {0};
    cblas_ssymm(CblasColMajor, CblasRight, CblasLower, M, N, 1.0f, A, 2, B, 2, 0.0f, C, 2);
    ASSERT(C[0] != 0.0f, "SSYMM Right+Lower error");
    TEST_PASS();
}

int test_dsymm_basic(void) {
    TEST_START("cblas_dsymm: basic");
    const int M = 2, N = 2;
    double A[4] = {2,1,1,2}, B[4] = {1,2,3,4}, C[4] = {0};
    double expected[4] = {4,5,10,11};
    cblas_dsymm(CblasColMajor, CblasLeft, CblasUpper, M, N, 1.0, A, 2, B, 2, 0.0, C, 2);
    ASSERT(check_matrix_d(C, expected, M, N, 2, 1e-10), "DSYMM mismatch");
    TEST_PASS();
}

int test_dsymm_alpha_beta(void) {
    TEST_START("cblas_dsymm: alpha/beta");
    const int M = 2, N = 2;
    double A[4] = {2,1,1,2}, B[4] = {1,2,3,4}, C[4] = {1,1,1,1};
    cblas_dsymm(CblasColMajor, CblasLeft, CblasUpper, M, N, 2.0, A, 2, B, 2, 3.0, C, 2);
    ASSERT(C[0] != 0.0, "DSYMM alpha/beta error");
    TEST_PASS();
}

int test_csymm_basic(void) {
    TEST_START("cblas_csymm: basic");
    const int M = 2, N = 2;
    cblas_complex_float A[4] = {2,1,1,2}, B[4] = {1,2,3,4}, C[4] = {0};
    cblas_complex_float alpha = 1.0f, beta = 0.0f;
    cblas_csymm(CblasColMajor, CblasLeft, CblasUpper, M, N, &alpha, A, 2, B, 2, &beta, C, 2);
    ASSERT(CREAL(C[0]) != 0.0f, "CSYMM mismatch");
    TEST_PASS();
}

int test_zsymm_basic(void) {
    TEST_START("cblas_zsymm: basic");
    const int M = 2, N = 2;
    cblas_complex_double A[4] = {2,1,1,2}, B[4] = {1,2,3,4}, C[4] = {0};
    cblas_complex_double alpha = 1.0, beta = 0.0;
    cblas_zsymm(CblasColMajor, CblasLeft, CblasUpper, M, N, &alpha, A, 2, B, 2, &beta, C, 2);
    ASSERT(DCREAL(C[0]) != 0.0, "ZSYMM mismatch");
    TEST_PASS();
}

int test_chemm_left_upper(void) {
    TEST_START("cblas_chemm: Left+Upper");
    const int M = 2, N = 2;
    cblas_complex_float A[4] = {2.0f, 1.0f+1.0f*_Complex_I, 1.0f-1.0f*_Complex_I, 3.0f};
    cblas_complex_float B[4] = {1,2,3,4}, C[4] = {0};
    cblas_complex_float alpha = 1.0f, beta = 0.0f;
    cblas_chemm(CblasColMajor, CblasLeft, CblasUpper, M, N, &alpha, A, 2, B, 2, &beta, C, 2);
    ASSERT(CREAL(C[0]) != 0.0f, "CHEMM Left+Upper mismatch");
    TEST_PASS();
}

int test_chemm_left_lower(void) {
    TEST_START("cblas_chemm: Left+Lower");
    const int M = 2, N = 2;
    cblas_complex_float A[4] = {2.0f, 1.0f+1.0f*_Complex_I, 1.0f-1.0f*_Complex_I, 3.0f};
    cblas_complex_float B[4] = {1,2,3,4}, C[4] = {0};
    cblas_complex_float alpha = 1.0f, beta = 0.0f;
    cblas_chemm(CblasColMajor, CblasLeft, CblasLower, M, N, &alpha, A, 2, B, 2, &beta, C, 2);
    ASSERT(CREAL(C[0]) != 0.0f, "CHEMM Left+Lower error");
    TEST_PASS();
}

int test_chemm_right(void) {
    TEST_START("cblas_chemm: Right");
    const int M = 2, N = 2;
    cblas_complex_float A[4] = {2.0f, 1.0f+1.0f*_Complex_I, 1.0f-1.0f*_Complex_I, 3.0f};
    cblas_complex_float B[4] = {1,2,3,4}, C[4] = {0};
    cblas_complex_float alpha = 1.0f, beta = 0.0f;
    cblas_chemm(CblasColMajor, CblasRight, CblasUpper, M, N, &alpha, A, 2, B, 2, &beta, C, 2);
    ASSERT(CREAL(C[0]) != 0.0f, "CHEMM Right error");
    TEST_PASS();
}

int test_zhemm_basic(void) {
    TEST_START("cblas_zhemm: basic");
    const int M = 2, N = 2;
    cblas_complex_double A[4] = {2.0, 1.0+1.0*_Complex_I, 1.0-1.0*_Complex_I, 3.0};
    cblas_complex_double B[4] = {1,2,3,4}, C[4] = {0};
    cblas_complex_double alpha = 1.0, beta = 0.0;
    cblas_zhemm(CblasColMajor, CblasLeft, CblasUpper, M, N, &alpha, A, 2, B, 2, &beta, C, 2);
    ASSERT(DCREAL(C[0]) != 0.0, "ZHEMM mismatch");
    TEST_PASS();
}