#include "test_utils.h"

int test_ssyr2k_upper_notrans(void) {
    TEST_START("cblas_ssyr2k: Upper+NoTrans");
    const int N = 2, K = 2;
    float A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    cblas_ssyr2k(CblasColMajor, CblasUpper, CblasNoTrans, N, K, 1.0f, A, 2, B, 2, 0.0f, C, 2);
    ASSERT(C[0] != 0, "SSYR2K Upper+NoTrans mismatch");
    TEST_PASS();
}

int test_ssyr2k_upper_trans(void) {
    TEST_START("cblas_ssyr2k: Upper+Trans");
    const int N = 2, K = 2;
    float A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    cblas_ssyr2k(CblasColMajor, CblasUpper, CblasTrans, N, K, 1.0f, A, 2, B, 2, 0.0f, C, 2);
    ASSERT(C[0] != 0, "SSYR2K Upper+Trans mismatch");
    TEST_PASS();
}

int test_ssyr2k_lower(void) {
    TEST_START("cblas_ssyr2k: Lower");
    const int N = 2, K = 2;
    float A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    cblas_ssyr2k(CblasColMajor, CblasLower, CblasNoTrans, N, K, 1.0f, A, 2, B, 2, 0.0f, C, 2);
    ASSERT(C[0] != 0, "SSYR2K Lower mismatch");
    TEST_PASS();
}

int test_ssyr2k_alpha_beta(void) {
    TEST_START("cblas_ssyr2k: alpha/beta");
    const int N = 2, K = 2;
    float A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {1,1,1,1};
    cblas_ssyr2k(CblasColMajor, CblasUpper, CblasNoTrans, N, K, 2.0f, A, 2, B, 2, 3.0f, C, 2);
    ASSERT(C[0] != 0, "SSYR2K alpha/beta error");
    TEST_PASS();
}

int test_dsyr2k_basic(void) {
    TEST_START("cblas_dsyr2k: basic");
    const int N = 2, K = 2;
    double A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    cblas_dsyr2k(CblasColMajor, CblasUpper, CblasNoTrans, N, K, 1.0, A, 2, B, 2, 0.0, C, 2);
    ASSERT(C[0] != 0, "DSYR2K mismatch");
    TEST_PASS();
}

int test_csyr2k_basic(void) {
    TEST_START("cblas_csyr2k: basic");
    const int N = 2, K = 2;
    cblas_complex_float A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    cblas_complex_float alpha = 1.0f, beta = 0.0f;
    cblas_csyr2k(CblasColMajor, CblasUpper, CblasNoTrans, N, K, &alpha, A, 2, B, 2, &beta, C, 2);
    ASSERT(CREAL(C[0]) != 0, "CSYR2K mismatch");
    TEST_PASS();
}

int test_zsyr2k_basic(void) {
    TEST_START("cblas_zsyr2k: basic");
    const int N = 2, K = 2;
    cblas_complex_double A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    cblas_complex_double alpha = 1.0, beta = 0.0;
    cblas_zsyr2k(CblasColMajor, CblasUpper, CblasNoTrans, N, K, &alpha, A, 2, B, 2, &beta, C, 2);
    ASSERT(DCREAL(C[0]) != 0, "ZSYR2K mismatch");
    TEST_PASS();
}

int test_cher2k_upper_conjtrans(void) {
    TEST_START("cblas_cher2k: Upper+ConjTrans");
    const int N = 2, K = 2;
    cblas_complex_float A[4] = {1, 1+_Complex_I, 1-_Complex_I, 2};
    cblas_complex_float B[4] = {2, 1-_Complex_I, 1+_Complex_I, 3};
    cblas_complex_float C[4] = {0};
    cblas_complex_float alpha = 1.0f;
    float beta = 0.0f;
    cblas_cher2k(CblasColMajor, CblasUpper, CblasConjTrans, N, K, &alpha, A, 2, B, 2, beta, C, 2);
    ASSERT(fabsf(CIMAG(C[0])) < 1e-5f, "CHER2K Upper+ConjTrans mismatch");
    TEST_PASS();
}

int test_cher2k_upper_notrans(void) {
    TEST_START("cblas_cher2k: Upper+NoTrans");
    const int N = 2, K = 2;
    cblas_complex_float A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    cblas_complex_float alpha = 1.0f;
    float beta = 0.0f;
    cblas_cher2k(CblasColMajor, CblasUpper, CblasNoTrans, N, K, &alpha, A, 2, B, 2, beta, C, 2);
    ASSERT(CREAL(C[0]) != 0, "CHER2K Upper+NoTrans mismatch");
    TEST_PASS();
}

int test_cher2k_lower(void) {
    TEST_START("cblas_cher2k: Lower");
    const int N = 2, K = 2;
    cblas_complex_float A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    cblas_complex_float alpha = 1.0f;
    float beta = 0.0f;
    cblas_cher2k(CblasColMajor, CblasLower, CblasConjTrans, N, K, &alpha, A, 2, B, 2, beta, C, 2);
    ASSERT(CREAL(C[0]) != 0, "CHER2K Lower mismatch");
    TEST_PASS();
}

int test_zher2k_basic(void) {
    TEST_START("cblas_zher2k: basic");
    const int N = 2, K = 2;
    cblas_complex_double A[4] = {1, 1+_Complex_I, 1-_Complex_I, 2};
    cblas_complex_double B[4] = {2, 1-_Complex_I, 1+_Complex_I, 3};
    cblas_complex_double C[4] = {0};
    cblas_complex_double alpha = 1.0;
    double beta = 0.0;
    cblas_zher2k(CblasColMajor, CblasUpper, CblasConjTrans, N, K, &alpha, A, 2, B, 2, beta, C, 2);
    ASSERT(fabs(DCIMAG(C[0])) < 1e-10, "ZHER2K mismatch");
    TEST_PASS();
}