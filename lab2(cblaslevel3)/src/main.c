#include <stdio.h>
#include <stdlib.h>
#include "test_utils.h"

extern int test_sgemm_nn(void), test_sgemm_tn(void), test_sgemm_nt(void), test_sgemm_tt(void);
extern int test_sgemm_alpha_beta(void);
extern int test_dgemm_nn(void), test_dgemm_transpose(void), test_dgemm_alpha_beta(void);
extern int test_cgemm_nn(void), test_cgemm_conjtrans(void), test_cgemm_alpha_beta(void);
extern int test_zgemm_nn(void), test_zgemm_conjtrans(void), test_zgemm_alpha_beta(void);

extern int test_ssymm_left_upper(void), test_ssymm_left_lower(void);
extern int test_ssymm_right_upper(void), test_ssymm_right_lower(void);
extern int test_dsymm_basic(void), test_dsymm_alpha_beta(void);
extern int test_csymm_basic(void), test_zsymm_basic(void);
extern int test_chemm_left_upper(void), test_chemm_left_lower(void), test_chemm_right(void);
extern int test_zhemm_basic(void);

extern int test_strmm_left_upper(void), test_strmm_left_lower(void), test_strmm_right_upper(void);
extern int test_strmm_trans(void), test_strmm_unit_diag(void), test_strmm_alpha(void);
extern int test_dtrmm_basic(void), test_ctrmm_basic(void), test_ztrmm_basic(void);
extern int test_strsm_left_upper(void), test_strsm_left_lower(void), test_strsm_right(void);
extern int test_strsm_trans(void), test_strsm_unit_diag(void);
extern int test_dtrsm_basic(void), test_ctrsm_basic(void), test_ztrsm_basic(void);

extern int test_ssyrk_upper_notrans(void), test_ssyrk_upper_trans(void);
extern int test_ssyrk_lower_notrans(void), test_ssyrk_alpha_beta(void);
extern int test_dsyrk_basic(void), test_csyrk_basic(void), test_zsyrk_basic(void);
extern int test_cherk_upper_conjtrans(void), test_cherk_upper_notrans(void), test_cherk_lower(void);
extern int test_zherk_basic(void);

extern int test_ssyr2k_upper_notrans(void), test_ssyr2k_upper_trans(void);
extern int test_ssyr2k_lower(void), test_ssyr2k_alpha_beta(void);
extern int test_dsyr2k_basic(void), test_csyr2k_basic(void), test_zsyr2k_basic(void);
extern int test_cher2k_upper_conjtrans(void), test_cher2k_upper_notrans(void), test_cher2k_lower(void);
extern int test_zher2k_basic(void);

typedef int (*test_func_t)(void);
typedef struct { const char *name; test_func_t func; } test_case_t;

int main(void) {
    int passed = 0, total = 0;
    
    printf("\n============================================================\n");
    printf("===       CBLAS Level 3 Interface Tests                  ===\n");
    printf("===       Platform: Windows | OpenBLAS C Interface       ===\n");
    printf("============================================================\n\n");
    
    printf("--- Single-threaded Tests ---\n\n");
    set_threads(1);
    
    test_case_t tests[] = {
        {"cblas_sgemm: NoTrans×NoTrans (ColMajor)", test_sgemm_nn},
        {"cblas_sgemm: Trans×NoTrans", test_sgemm_tn},
        {"cblas_sgemm: NoTrans×Trans", test_sgemm_nt},
        {"cblas_sgemm: Trans×Trans", test_sgemm_tt},
        {"cblas_sgemm: alpha/beta (RowMajor)", test_sgemm_alpha_beta},
        {"cblas_dgemm: NoTrans×NoTrans", test_dgemm_nn},
        {"cblas_dgemm: Transpose", test_dgemm_transpose},
        {"cblas_dgemm: alpha/beta", test_dgemm_alpha_beta},
        {"cblas_cgemm: NoTrans×NoTrans", test_cgemm_nn},
        {"cblas_cgemm: ConjTrans", test_cgemm_conjtrans},
        {"cblas_cgemm: alpha/beta", test_cgemm_alpha_beta},
        {"cblas_zgemm: NoTrans×NoTrans", test_zgemm_nn},
        {"cblas_zgemm: ConjTrans", test_zgemm_conjtrans},
        {"cblas_zgemm: alpha/beta", test_zgemm_alpha_beta},
        {"cblas_ssymm: Left+Upper", test_ssymm_left_upper},
        {"cblas_ssymm: Left+Lower", test_ssymm_left_lower},
        {"cblas_ssymm: Right+Upper", test_ssymm_right_upper},
        {"cblas_ssymm: Right+Lower", test_ssymm_right_lower},
        {"cblas_dsymm: basic", test_dsymm_basic},
        {"cblas_dsymm: alpha/beta", test_dsymm_alpha_beta},
        {"cblas_csymm: basic", test_csymm_basic},
        {"cblas_zsymm: basic", test_zsymm_basic},
        {"cblas_chemm: Left+Upper", test_chemm_left_upper},
        {"cblas_chemm: Left+Lower", test_chemm_left_lower},
        {"cblas_chemm: Right", test_chemm_right},
        {"cblas_zhemm: basic", test_zhemm_basic},
        {"cblas_strmm: Left+Upper+NoTrans+NonUnit", test_strmm_left_upper},
        {"cblas_strmm: Left+Lower", test_strmm_left_lower},
        {"cblas_strmm: Right+Upper", test_strmm_right_upper},
        {"cblas_strmm: Trans", test_strmm_trans},
        {"cblas_strmm: Unit+Diag", test_strmm_unit_diag},
        {"cblas_strmm: alpha", test_strmm_alpha},
        {"cblas_dtrmm: basic", test_dtrmm_basic},
        {"cblas_ctrmm: basic", test_ctrmm_basic},
        {"cblas_ztrmm: basic", test_ztrmm_basic},
        {"cblas_strsm: Left+Upper+NoTrans+NonUnit", test_strsm_left_upper},
        {"cblas_strsm: Left+Lower", test_strsm_left_lower},
        {"cblas_strsm: Right", test_strsm_right},
        {"cblas_strsm: Trans", test_strsm_trans},
        {"cblas_strsm: Unit+Diag", test_strsm_unit_diag},
        {"cblas_dtrsm: basic", test_dtrsm_basic},
        {"cblas_ctrsm: basic", test_ctrsm_basic},
        {"cblas_ztrsm: basic", test_ztrsm_basic},
        {"cblas_ssyrk: Upper+NoTrans", test_ssyrk_upper_notrans},
        {"cblas_ssyrk: Upper+Trans", test_ssyrk_upper_trans},
        {"cblas_ssyrk: Lower+NoTrans", test_ssyrk_lower_notrans},
        {"cblas_ssyrk: alpha/beta", test_ssyrk_alpha_beta},
        {"cblas_dsyrk: basic", test_dsyrk_basic},
        {"cblas_csyrk: basic", test_csyrk_basic},
        {"cblas_zsyrk: basic", test_zsyrk_basic},
        {"cblas_cherk: Upper+ConjTrans", test_cherk_upper_conjtrans},
        {"cblas_cherk: Upper+NoTrans", test_cherk_upper_notrans},
        {"cblas_cherk: Lower", test_cherk_lower},
        {"cblas_zherk: basic", test_zherk_basic},
        {"cblas_ssyr2k: Upper+NoTrans", test_ssyr2k_upper_notrans},
        {"cblas_ssyr2k: Upper+Trans", test_ssyr2k_upper_trans},
        {"cblas_ssyr2k: Lower", test_ssyr2k_lower},
        {"cblas_ssyr2k: alpha/beta", test_ssyr2k_alpha_beta},
        {"cblas_dsyr2k: basic", test_dsyr2k_basic},
        {"cblas_csyr2k: basic", test_csyr2k_basic},
        {"cblas_zsyr2k: basic", test_zsyr2k_basic},
        {"cblas_cher2k: Upper+ConjTrans", test_cher2k_upper_conjtrans},
        {"cblas_cher2k: Upper+NoTrans", test_cher2k_upper_notrans},
        {"cblas_cher2k: Lower", test_cher2k_lower},
        {"cblas_zher2k: basic", test_zher2k_basic},
    };
    
    total = sizeof(tests) / sizeof(tests[0]);
    for (int i = 0; i < total; i++) {
        if (tests[i].func()) passed++;
    }
    
    printf("\n--- Multi-threaded Tests (4 threads) ---\n\n");
    set_threads(4);
    
    int mt_tests = 7;
    total += mt_tests;
    if (test_sgemm_nn()) passed++;
    if (test_dgemm_nn()) passed++;
    if (test_ssymm_left_upper()) passed++;
    if (test_strmm_left_upper()) passed++;
    if (test_ssyrk_upper_notrans()) passed++;
    if (test_ssyr2k_upper_notrans()) passed++;
    if (test_zgemm_nn()) passed++;
    
    printf("\n============================================================\n");
    printf("===                    RESULTS                           ===\n");
    printf("============================================================\n");
    printf("  Total tests:  %d\n", total);
    printf("  Passed:       %d\n", passed);
    printf("  Failed:       %d\n", total - passed);
    printf("  Success rate: %.1f%%\n", (float)passed / total * 100);
    printf("============================================================\n\n");
    
    if (passed == total) {
        printf("✓ All tests PASSED!\n\n");
        return 0;
    } else {
        printf("✗ Some tests FAILED!\n\n");
        return 1;
    }
}