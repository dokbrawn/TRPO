#include "test_utils.h"

int test_ssyrk_upper_notrans(void) {
    TEST_START("cblas_ssyrk: Upper+NoTrans");
    const int N = 2, K = 2;
    float A[4] = {1,2,3,4}, C[4] = {0};
    cblas_ssyrk(CblasColMajor, CblasUpper, CblasNoTrans, N, K, 1.0f, A, 2, 0.0f, C, 2);
    ASSERT(C[0] > 0 && C[3] > 0, "SSYRK Upper+NoTrans mismatch");
    TEST_PASS();
}

int test_ssyrk_upper_trans(void) {
    TEST_START("cblas_ssyrk: Upper+Trans");
    const int N = 2, K = 2;
    float A[4] = {1,2,3,4}, C[4] = {0};
    cblas_ssyrk(CblasColMajor, CblasUpper, CblasTrans, N, K, 1.0f, A, 2, 0.0f, C, 2);
    ASSERT(C[0] > 0 && C[3] > 0, "SSYRK Upper+Trans mismatch");
    TEST_PASS();
}

int test_ssyrk_lower_notrans(void) {
    TEST_START("cblas_ssyrk: Lower+NoTrans");
    const int N = 2, K = 2;
    float A[4] = {1,2,3,4}, C[4] = {0};
    cblas_ssyrk(CblasColMajor, CblasLower, CblasNoTrans, N, K, 1.0f, A, 2, 0.0f, C, 2);
    ASSERT(C[0] > 0 && C[3] > 0, "SSYRK Lower+NoTrans mismatch");
    TEST_PASS();
}

int test_ssyrk_alpha_beta(void) {
    TEST_START("cblas_ssyrk: alpha/beta");
    const int N = 2, K = 2;
    float A[4] = {1,2,3,4}, C[4] = {1,1,1,1};
    cblas_ssyrk(CblasColMajor, CblasUpper, CblasNoTrans, N, K, 2.0f, A, 2, 3.0f, C, 2);
    ASSERT(C[0] > 0, "SSYRK alpha/beta error");
    TEST_PASS();
}

int test_dsyrk_basic(void) {
    TEST_START("cblas_dsyrk: basic");
    const int N = 2, K = 2;
    double A[4] = {1,2,3,4}, C[4] = {0};
    cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans, N, K, 1.0, A, 2, 0.0, C, 2);
    ASSERT(C[0] > 0 && C[3] > 0, "DSYRK mismatch");
    TEST_PASS();
}

int test_csyrk_basic(void) {
    TEST_START("cblas_csyrk: basic");
    const int N = 2, K = 2;
    cblas_complex_float A[4] = {1,2,3,4}, C[4] = {0};
    cblas_complex_float alpha = 1.0f, beta = 0.0f;
    cblas_csyrk(CblasColMajor, CblasUpper, CblasNoTrans, N, K, &alpha, A, 2, &beta, C, 2);
    ASSERT(CREAL(C[0]) > 0, "CSYRK mismatch");
    TEST_PASS();
}

int test_zsyrk_basic(void) {
    TEST_START("cblas_zsyrk: basic");
    const int N = 2, K = 2;
    cblas_complex_double A[4] = {1,2,3,4}, C[4] = {0};
    cblas_complex_double alpha = 1.0, beta = 0.0;
    cblas_zsyrk(CblasColMajor, CblasUpper, CblasNoTrans, N, K, &alpha, A, 2, &beta, C, 2);
    ASSERT(DCREAL(C[0]) > 0, "ZSYRK mismatch");
    TEST_PASS();
}

int test_cherk_upper_conjtrans(void) {
    TEST_START("cblas_cherk: Upper+ConjTrans");
    const int N = 2, K = 2;
    cblas_complex_float A[4] = {1, 1+_Complex_I, 1-_Complex_I, 2};
    cblas_complex_float C[4] = {0};
    float alpha = 1.0f, beta = 0.0f;
    cblas_cherk(CblasColMajor, CblasUpper, CblasConjTrans, N, K, alpha, A, 2, beta, C, 2);
    ASSERT(fabsf(CIMAG(C[0])) < 1e-5f && CREAL(C[0]) > 0, "CHERK Upper+ConjTrans mismatch");
    TEST_PASS();
}

int test_cherk_upper_notrans(void) {
    TEST_START("cblas_cherk: Upper+NoTrans");
    const int N = 2, K = 2;
    cblas_complex_float A[4] = {1,2,3,4};
    cblas_complex_float C[4] = {0};
    float alpha = 1.0f, beta = 0.0f;
    cblas_cherk(CblasColMajor, CblasUpper, CblasNoTrans, N, K, alpha, A, 2, beta, C, 2);
    ASSERT(CREAL(C[0]) > 0, "CHERK Upper+NoTrans mismatch");
    TEST_PASS();
}

int test_cherk_lower(void) {
    TEST_START("cblas_cherk: Lower");
    const int N = 2, K = 2;
    cblas_complex_float A[4] = {1,2,3,4};
    cblas_complex_float C[4] = {0};
    float alpha = 1.0f, beta = 0.0f;
    cblas_cherk(CblasColMajor, CblasLower, CblasConjTrans, N, K, alpha, A, 2, beta, C, 2);
    ASSERT(CREAL(C[0]) > 0, "CHERK Lower mismatch");
    TEST_PASS();
}

int test_zherk_basic(void) {
    TEST_START("cblas_zherk: basic");
    const int N = 2, K = 2;
    cblas_complex_double A[4] = {1, 1+_Complex_I, 1-_Complex_I, 2};
    cblas_complex_double C[4] = {0};
    double alpha = 1.0, beta = 0.0;
    cblas_zherk(CblasColMajor, CblasUpper, CblasConjTrans, N, K, alpha, A, 2, beta, C, 2);
    ASSERT(fabs(DCIMAG(C[0])) < 1e-10 && DCREAL(C[0]) > 0, "ZHERK mismatch");
    TEST_PASS();
}