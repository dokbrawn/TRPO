#include "test_utils.h"

int test_strmm_left_upper(void) {
    TEST_START("cblas_strmm: Left+Upper+NoTrans+NonUnit");
    const int M = 2, N = 2;
    float A[4] = {1,0,2,3}, B[4] = {1,2,3,4};
    float expected[4] = {5,6,11,12};
    cblas_strmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, M, N, 1.0f, A, 2, B, 2);
    ASSERT(check_matrix_f(B, expected, M, N, 2, 1e-5f), "STRMM mismatch");
    TEST_PASS();
}

int test_strmm_left_lower(void) {
    TEST_START("cblas_strmm: Left+Lower");
    const int M = 2, N = 2;
    float A[4] = {1,2,0,3}, B[4] = {1,2,3,4};
    cblas_strmm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, M, N, 1.0f, A, 2, B, 2);
    ASSERT(B[0] != 0.0f, "STRMM Lower error");
    TEST_PASS();
}

int test_strmm_right_upper(void) {
    TEST_START("cblas_strmm: Right+Upper");
    const int M = 2, N = 2;
    float A[4] = {1,0,2,3}, B[4] = {1,2,3,4};
    cblas_strmm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, M, N, 1.0f, A, 2, B, 2);
    ASSERT(B[0] != 0.0f, "STRMM Right error");
    TEST_PASS();
}

int test_strmm_trans(void) {
    TEST_START("cblas_strmm: Trans");
    const int M = 2, N = 2;
    float A[4] = {1,2,0,3}, B[4] = {1,2,3,4};
    cblas_strmm(CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasNonUnit, M, N, 1.0f, A, 2, B, 2);
    ASSERT(B[0] != 0.0f, "STRMM Trans error");
    TEST_PASS();
}

int test_strmm_unit_diag(void) {
    TEST_START("cblas_strmm: Unit+Diag");
    const int M = 2, N = 2;
    float A[4] = {1,0,2,1}, B[4] = {1,2,3,4};
    cblas_strmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasUnit, M, N, 1.0f, A, 2, B, 2);
    ASSERT(B[0] != 0.0f, "STRMM Unit error");
    TEST_PASS();
}

int test_strmm_alpha(void) {
    TEST_START("cblas_strmm: alpha");
    const int M = 2, N = 2;
    float A[4] = {1,0,2,3}, B[4] = {1,2,3,4};
    cblas_strmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, M, N, 2.0f, A, 2, B, 2);
    ASSERT(B[0] != 0.0f, "STRMM alpha error");
    TEST_PASS();
}

int test_dtrmm_basic(void) {
    TEST_START("cblas_dtrmm: basic");
    const int M = 2, N = 2;
    double A[4] = {1,0,2,3}, B[4] = {1,2,3,4};
    double expected[4] = {5,6,11,12};
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, M, N, 1.0, A, 2, B, 2);
    ASSERT(check_matrix_d(B, expected, M, N, 2, 1e-10), "DTRMM mismatch");
    TEST_PASS();
}

int test_ctrmm_basic(void) {
    TEST_START("cblas_ctrmm: basic");
    const int M = 2, N = 2;
    cblas_complex_float A[4] = {1,0,2,3}, B[4] = {1,2,3,4};
    cblas_complex_float alpha = 1.0f;
    cblas_ctrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, M, N, &alpha, A, 2, B, 2);
    ASSERT(CREAL(B[0]) != 0.0f, "CTRMM mismatch");
    TEST_PASS();
}

int test_ztrmm_basic(void) {
    TEST_START("cblas_ztrmm: basic");
    const int M = 2, N = 2;
    cblas_complex_double A[4] = {1,0,2,3}, B[4] = {1,2,3,4};
    cblas_complex_double alpha = 1.0;
    cblas_ztrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, M, N, &alpha, A, 2, B, 2);
    ASSERT(DCREAL(B[0]) != 0.0, "ZTRMM mismatch");
    TEST_PASS();
}

int test_strsm_left_upper(void) {
    TEST_START("cblas_strsm: Left+Upper+NoTrans+NonUnit");
    const int M = 2, N = 2;
    float A[4] = {2,0,1,3}, B[4] = {4,6,10,12};
    float expected[4] = {1,2,3,4};
    cblas_strsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, M, N, 1.0f, A, 2, B, 2);
    ASSERT(check_matrix_f(B, expected, M, N, 2, 1e-5f), "STRSM mismatch");
    TEST_PASS();
}

int test_strsm_left_lower(void) {
    TEST_START("cblas_strsm: Left+Lower");
    const int M = 2, N = 2;
    float A[4] = {2,1,0,3}, B[4] = {5,14,4,15};
    cblas_strsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, M, N, 1.0f, A, 2, B, 2);
    ASSERT(B[0] != 0.0f, "STRSM Lower error");
    TEST_PASS();
}

int test_strsm_right(void) {
    TEST_START("cblas_strsm: Right");
    const int M = 2, N = 2;
    float A[4] = {2,0,1,3}, B[4] = {5,8,11,20};
    cblas_strsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, M, N, 1.0f, A, 2, B, 2);
    ASSERT(B[0] != 0.0f, "STRSM Right error");
    TEST_PASS();
}

int test_strsm_trans(void) {
    TEST_START("cblas_strsm: Trans");
    const int M = 2, N = 2;
    float A[4] = {2,1,0,3}, B[4] = {5,11,8,20};
    cblas_strsm(CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasNonUnit, M, N, 1.0f, A, 2, B, 2);
    ASSERT(B[0] != 0.0f, "STRSM Trans error");
    TEST_PASS();
}

int test_strsm_unit_diag(void) {
    TEST_START("cblas_strsm: Unit+Diag");
    const int M = 2, N = 2;
    float A[4] = {1,0,2,1}, B[4] = {5,2,11,4};
    float expected[4] = {1,2,3,4};
    cblas_strsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasUnit, M, N, 1.0f, A, 2, B, 2);
    ASSERT(check_matrix_f(B, expected, M, N, 2, 1e-5f), "STRSM Unit error");
    TEST_PASS();
}

int test_dtrsm_basic(void) {
    TEST_START("cblas_dtrsm: basic");
    const int M = 2, N = 2;
    double A[4] = {2,0,1,3}, B[4] = {4,6,10,12};
    double expected[4] = {1,2,3,4};
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, M, N, 1.0, A, 2, B, 2);
    ASSERT(check_matrix_d(B, expected, M, N, 2, 1e-10), "DTRSM mismatch");
    TEST_PASS();
}

int test_ctrsm_basic(void) {
    TEST_START("cblas_ctrsm: basic");
    const int M = 2, N = 2;
    cblas_complex_float A[4] = {2,0,1,3}, B[4] = {4,6,10,12};
    cblas_complex_float alpha = 1.0f;
    cblas_ctrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, M, N, &alpha, A, 2, B, 2);
    ASSERT(CREAL(B[0]) != 0.0f, "CTRSM mismatch");
    TEST_PASS();
}

int test_ztrsm_basic(void) {
    TEST_START("cblas_ztrsm: basic");
    const int M = 2, N = 2;
    cblas_complex_double A[4] = {2,0,1,3}, B[4] = {4,6,10,12};
    cblas_complex_double alpha = 1.0;
    cblas_ztrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, M, N, &alpha, A, 2, B, 2);
    ASSERT(DCREAL(B[0]) != 0.0, "ZTRSM mismatch");
    TEST_PASS();
}