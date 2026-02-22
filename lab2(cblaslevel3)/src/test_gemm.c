#include "test_utils.h"

int test_sgemm_nn(void) {
    TEST_START("cblas_sgemm: NoTrans×NoTrans (ColMajor)");
    const int M = 2, N = 2, K = 2;
    float A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    float expected[4] = {23,34,31,46};
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0f, A, M, B, K, 0.0f, C, M);
    ASSERT(check_matrix_f(C, expected, M, N, M, 1e-5f), "Result mismatch");
    TEST_PASS();
}

int test_sgemm_tn(void) {
    TEST_START("cblas_sgemm: Trans×NoTrans");
    const int M = 2, N = 2, K = 2;
    float A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0f, A, K, B, K, 0.0f, C, M);
    ASSERT(C[0] != 0.0f, "TransA error");
    TEST_PASS();
}

int test_sgemm_nt(void) {
    TEST_START("cblas_sgemm: NoTrans×Trans");
    const int M = 2, N = 2, K = 2;
    float A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0f, A, M, B, N, 0.0f, C, M);
    ASSERT(C[0] != 0.0f, "TransB error");
    TEST_PASS();
}

int test_sgemm_tt(void) {
    TEST_START("cblas_sgemm: Trans×Trans");
    const int M = 2, N = 2, K = 2;
    float A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    cblas_sgemm(CblasColMajor, CblasTrans, CblasTrans, M, N, K, 1.0f, A, K, B, N, 0.0f, C, M);
    ASSERT(C[0] != 0.0f, "TransAB error");
    TEST_PASS();
}

int test_sgemm_alpha_beta(void) {
    TEST_START("cblas_sgemm: alpha/beta (RowMajor)");
    const int M = 2, N = 2, K = 2;
    float A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {1,1,1,1};
    float expected[4] = {41,47,89,103};
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 2.0f, A, K, B, N, 3.0f, C, N);
    ASSERT(check_matrix_f(C, expected, M, N, N, 1e-5f), "Alpha/beta error");
    TEST_PASS();
}

int test_dgemm_nn(void) {
    TEST_START("cblas_dgemm: NoTrans×NoTrans");
    const int M = 2, N = 2, K = 2;
    double A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    double expected[4] = {23,34,31,46};
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0, A, M, B, K, 0.0, C, M);
    ASSERT(check_matrix_d(C, expected, M, N, M, 1e-10), "DGEMM mismatch");
    TEST_PASS();
}

int test_dgemm_transpose(void) {
    TEST_START("cblas_dgemm: Transpose");
    const int M = 2, N = 2, K = 2;
    double A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, A, K, B, K, 0.0, C, M);
    ASSERT(C[0] != 0.0, "DGEMM transpose error");
    TEST_PASS();
}

int test_dgemm_alpha_beta(void) {
    TEST_START("cblas_dgemm: alpha/beta");
    const int M = 2, N = 2, K = 2;
    double A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {1,1,1,1};
    double expected[4] = {41,47,89,103};
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 2.0, A, K, B, N, 3.0, C, N);
    ASSERT(check_matrix_d(C, expected, M, N, N, 1e-10), "DGEMM alpha/beta error");
    TEST_PASS();
}

int test_cgemm_nn(void) {
    TEST_START("cblas_cgemm: NoTrans×NoTrans");
    const int M = 2, N = 2, K = 2;
    cblas_complex_float alpha = 1.0f, beta = 0.0f;
    cblas_complex_float A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    cblas_complex_float expected[4] = {23,34,31,46};
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, K, &alpha, A, M, B, K, &beta, C, M);
    ASSERT(check_matrix_f((float*)C, (float*)expected, M, N, M, 1e-5f), "CGEMM mismatch");
    TEST_PASS();
}

int test_cgemm_conjtrans(void) {
    TEST_START("cblas_cgemm: ConjTrans");
    const int M = 2, N = 2, K = 2;
    cblas_complex_float alpha = 1.0f, beta = 0.0f;
    cblas_complex_float A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    cblas_cgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, M, N, K, &alpha, A, K, B, K, &beta, C, M);
    ASSERT(CREAL(C[0]) != 0.0f, "CGEMM ConjTrans error");
    TEST_PASS();
}

int test_cgemm_alpha_beta(void) {
    TEST_START("cblas_cgemm: alpha/beta");
    const int M = 2, N = 2, K = 2;
    cblas_complex_float alpha = 2.0f, beta = 3.0f;
    cblas_complex_float A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {1,1,1,1};
    cblas_complex_float expected[4] = {41,47,89,103};
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, &alpha, A, K, B, N, &beta, C, N);
    ASSERT(check_matrix_f((float*)C, (float*)expected, M, N, N, 1e-5f), "CGEMM alpha/beta error");
    TEST_PASS();
}

int test_zgemm_nn(void) {
    TEST_START("cblas_zgemm: NoTrans×NoTrans");
    const int M = 2, N = 2, K = 2;
    cblas_complex_double alpha = 1.0, beta = 0.0;
    cblas_complex_double A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    cblas_complex_double expected[4] = {23,34,31,46};
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, K, &alpha, A, M, B, K, &beta, C, M);
    ASSERT(check_matrix_d((double*)C, (double*)expected, M, N, M, 1e-10), "ZGEMM mismatch");
    TEST_PASS();
}

int test_zgemm_conjtrans(void) {
    TEST_START("cblas_zgemm: ConjTrans");
    const int M = 2, N = 2, K = 2;
    cblas_complex_double alpha = 1.0, beta = 0.0;
    cblas_complex_double A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {0};
    cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, M, N, K, &alpha, A, K, B, K, &beta, C, M);
    ASSERT(DCREAL(C[0]) != 0.0, "ZGEMM ConjTrans error");
    TEST_PASS();
}

int test_zgemm_alpha_beta(void) {
    TEST_START("cblas_zgemm: alpha/beta");
    const int M = 2, N = 2, K = 2;
    cblas_complex_double alpha = 2.0, beta = 3.0;
    cblas_complex_double A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4] = {1,1,1,1};
    cblas_complex_double expected[4] = {41,47,89,103};
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, &alpha, A, K, B, N, &beta, C, N);
    ASSERT(check_matrix_d((double*)C, (double*)expected, M, N, N, 1e-10), "ZGEMM alpha/beta error");
    TEST_PASS();
}