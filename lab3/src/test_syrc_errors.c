#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include "test_utils.h"
#include "csyrk_zsyrk.h"

int test_csyrk_correct(void) {
    TEST_START("CSYRK: корректный вызов");
    
    const int N = 2, K = 2;
    cblas_complex_float A[4] = {{1.0f, 0.0f}, {2.0f, 0.0f}, {3.0f, 0.0f}, {4.0f, 0.0f}};
    cblas_complex_float C[4] = {{0, 0}};
    cblas_complex_float alpha = {1.0f, 0.0f}, beta = {0.0f, 0.0f};
    
    int err = cblas_csyrk_safe(CblasColMajor, CblasUpper, CblasNoTrans,
                               N, K, &alpha, A, 2, &beta, C, 2);
    
    ASSERT(err == 0, "CSYRK вернул ошибку");
    ASSERT(C[0].real > 0, "CSYRK не вычислил результат");
    
    printf(" [УСПЕХ] C[0]=%.2f%+.2fi\n", C[0].real, C[0].imag);
    TEST_PASS();
}

int test_csyrk_null_ptr(void) {
    TEST_START("CSYRK: NULL указатель");
    
    cblas_complex_float A[4] = {{1, 0}, {2, 0}, {3, 0}, {4, 0}};
    cblas_complex_float C[4] = {{0, 0}};
    cblas_complex_float alpha = {1.0f, 0.0f}, beta = {0.0f, 0.0f};
    
    printf(" [ОЖИДАЕТСЯ ОШИБКА] ");
    int err = cblas_csyrk_safe(CblasColMajor, CblasUpper, CblasNoTrans,
                               2, 2, &alpha, NULL, 2, &beta, C, 2);
    
    ASSERT(err == SYRC_ERR_NULL_PTR, "CSYRK не обнаружил NULL");
    printf(" [УСПЕХ]\n");
    TEST_PASS();
}

int test_csyrk_invalid_n(void) {
    TEST_START("CSYRK: отрицательный N");
    
    cblas_complex_float A[4] = {{1, 0}, {2, 0}, {3, 0}, {4, 0}};
    cblas_complex_float C[4] = {{0, 0}};
    cblas_complex_float alpha = {1.0f, 0.0f}, beta = {0.0f, 0.0f};
    
    printf(" [ОЖИДАЕТСЯ ОШИБКА] ");
    int err = cblas_csyrk_safe(CblasColMajor, CblasUpper, CblasNoTrans,
                               -1, 2, &alpha, A, 2, &beta, C, 2);
    
    ASSERT(err == SYRC_ERR_N, "CSYRK не обнаружил неверный N");
    printf(" [УСПЕХ]\n");
    TEST_PASS();
}

int test_csyrk_invalid_k(void) {
    TEST_START("CSYRK: отрицательный K");
    
    cblas_complex_float A[4] = {{1, 0}, {2, 0}, {3, 0}, {4, 0}};
    cblas_complex_float C[4] = {{0, 0}};
    cblas_complex_float alpha = {1.0f, 0.0f}, beta = {0.0f, 0.0f};
    
    printf(" [ОЖИДАЕТСЯ ОШИБКА] ");
    int err = cblas_csyrk_safe(CblasColMajor, CblasUpper, CblasNoTrans,
                               2, -1, &alpha, A, 2, &beta, C, 2);
    
    ASSERT(err == SYRC_ERR_K, "CSYRK не обнаружил неверный K");
    printf(" [УСПЕХ]\n");
    TEST_PASS();
}

int test_csyrk_invalid_lda(void) {
    TEST_START("CSYRK: неверный LDA");
    
    cblas_complex_float A[4] = {{1, 0}, {2, 0}, {3, 0}, {4, 0}};
    cblas_complex_float C[4] = {{0, 0}};
    cblas_complex_float alpha = {1.0f, 0.0f}, beta = {0.0f, 0.0f};
    
    printf(" [ОЖИДАЕТСЯ ОШИБКА] ");
    int err = cblas_csyrk_safe(CblasColMajor, CblasUpper, CblasNoTrans,
                               2, 2, &alpha, A, 1, &beta, C, 2);
    
    ASSERT(err == SYRC_ERR_LDA, "CSYRK не обнаружил неверный LDA");
    printf(" [УСПЕХ]\n");
    TEST_PASS();
}

int test_zsyrk_correct(void) {
    TEST_START("ZSYRK: корректный вызов");
    
    const int N = 2, K = 2;
    cblas_complex_double A[4] = {{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}};
    cblas_complex_double C[4] = {{0, 0}};
    cblas_complex_double alpha = {1.0, 0.0}, beta = {0.0, 0.0};
    
    int err = cblas_zsyrk_safe(CblasColMajor, CblasUpper, CblasNoTrans,
                               N, K, &alpha, A, 2, &beta, C, 2);
    
    ASSERT(err == 0, "ZSYRK вернул ошибку");
    ASSERT(C[0].real > 0, "ZSYRK не вычислил результат");
    
    printf(" [УСПЕХ] C[0]=%.2f%+.2fi\n", C[0].real, C[0].imag);
    TEST_PASS();
}

int test_zsyrk_invalid_ldc(void) {
    TEST_START("ZSYRK: неверный LDC");
    
    cblas_complex_double A[4] = {{1, 0}, {2, 0}, {3, 0}, {4, 0}};
    cblas_complex_double C[4] = {{0, 0}};
    cblas_complex_double alpha = {1.0, 0.0}, beta = {0.0, 0.0};
    
    printf(" [ОЖИДАЕТСЯ ОШИБКА] ");
    int err = cblas_zsyrk_safe(CblasColMajor, CblasUpper, CblasNoTrans,
                               2, 2, &alpha, A, 2, &beta, C, 1);
    
    ASSERT(err == SYRC_ERR_LDC, "ZSYRK не обнаружил неверный LDC");
    printf(" [УСПЕХ]\n");
    TEST_PASS();
}

int test_csyrk_variations(void) {
    TEST_START("CSYRK: различные конфигурации");
    
    const int N = 3, K = 2;
    cblas_complex_float A[6] = {{1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}};
    cblas_complex_float C[9] = {{0, 0}};
    cblas_complex_float alpha = {2.0f, 0.0f}, beta = {0.5f, 0.0f};
    
    int err = cblas_csyrk_safe(CblasColMajor, CblasUpper, CblasNoTrans,
                               N, K, &alpha, A, N, &beta, C, N);
    ASSERT(err == 0, "CSYRK Upper+NoTrans failed");
    
    err = cblas_csyrk_safe(CblasColMajor, CblasLower, CblasTrans,
                           N, K, &alpha, A, K, &beta, C, N);
    ASSERT(err == 0, "CSYRK Lower+Trans failed");
    
    err = cblas_csyrk_safe(CblasRowMajor, CblasUpper, CblasNoTrans,
                           N, K, &alpha, A, K, &beta, C, N);
    ASSERT(err == 0, "CSYRK RowMajor failed");
    
    printf(" [УСПЕХ]\n");
    TEST_PASS();
}