#include <stdio.h>
#include <stdlib.h>
#include "test_utils.h"
#include "csyrk_zsyrk.h"

extern int test_csyrk_correct(void);
extern int test_csyrk_null_ptr(void);
extern int test_csyrk_invalid_n(void);
extern int test_csyrk_invalid_k(void);
extern int test_csyrk_invalid_lda(void);
extern int test_zsyrk_correct(void);
extern int test_zsyrk_invalid_ldc(void);
extern int test_csyrk_variations(void);

typedef int (*test_func_t)(void);

typedef struct {
    const char* name;
    test_func_t func;
} test_case_t;

int main(void) {
    int passed = 0, total = 0;
    
    printf("============================================================\n");
    printf("===       ТЕСТИРОВАНИЕ РЕАЛИЗАЦИИ SYRC                    ===\n");
    printf("===       Complex Symmetric Rank-K Update                  ===\n");
    printf("===       Платформа: Windows + Visual Studio               ===\n");
    printf("============================================================\n\n");
    
    test_case_t tests[] = {
        {"CSYRK: корректный вызов", test_csyrk_correct},
        {"CSYRK: обработка NULL указателя", test_csyrk_null_ptr},
        {"CSYRK: обработка неверного N", test_csyrk_invalid_n},
        {"CSYRK: обработка неверного K", test_csyrk_invalid_k},
        {"CSYRK: обработка неверного LDA", test_csyrk_invalid_lda},
        {"ZSYRK: корректный вызов", test_zsyrk_correct},
        {"ZSYRK: обработка неверного LDC", test_zsyrk_invalid_ldc},
        {"CSYRK: различные конфигурации", test_csyrk_variations},
    };
    
    total = sizeof(tests) / sizeof(tests[0]);
    
    printf("Запуск %d тестов...\n\n", total);
    
    for (int i = 0; i < total; i++) {
        printf("[%d/%d] ", i + 1, total);
        if (tests[i].func()) {
            passed++;
        }
    }
    
    printf("\n============================================================\n");
    printf("                   РЕЗУЛЬТАТЫ ТЕСТОВ\n");
    printf("============================================================\n");
    printf("  Всего тестов:  %d\n", total);
    printf("  Пройдено:      %d\n", passed);
    printf("  Провалено:     %d\n", total - passed);
    printf("  Успешность:    %.1f%%\n", (float)passed / total * 100);
    printf("============================================================\n");
    
    if (passed == total) {
        printf("✓ ВСЕ ТЕСТЫ ПРОЙДЕНЫ!\n");
        return 0;
    } else {
        printf("✗ НЕКОТОРЫЕ ТЕСТЫ ПРОВАЛЕНЫ!\n");
        return 1;
    }
}