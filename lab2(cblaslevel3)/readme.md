# Тестирование интерфейса CBLAS Level 3 (Lab 2)


Отчет по лабораторной работе №2 (ТРПО). Проведено сравнительное тестирование функций уровня 3 библиотеки **CBLAS** на двух реализациях:
1.  **Оригинал (OpenBLAS)** 
2.  **Mock библиотека** — тестовая реализация для проверки обработки ошибок.

---

## 📊 Сравнительная статистика

| Параметр | Оригинальная (OpenBLAS) | Mock библиотека |
| :--- | :---: | :---: |
| **Всего тестов** | 72 | 72 |
| **Успешно ✅** | 72 | 22 |
| **Провалено ❌** | 0 | 50 |
| **Успешность** | **100.0%** | **30.6%** |
| **Статус** | 🟢 Стабильно | 🔴 Критические ошибки |

---

## 1️⃣ Оригинальная библиотека (OpenBLAS)

Оригинальная библиотека показала полную корректность работы во всех режимах (Однопоточный и Многопоточный). Ошибок не выявлено.

### ✅ Ключевые результаты
- **GEMM/SYMM/HEMM:** Все типы данных (`float`, `double`, `complex`) работают корректно.
- **TRMM/TRSM:** Треугольные матрицы обрабатываются верно во всех конфигурациях.
- **SYRK/HERK/SYR2K/HER2K:** Операции ранга 1 и 2 выполнены без ошибок.
- **Многопоточность:** Распараллеливание (4 потока) не нарушает целостность вычислений.

<details>
<summary>📄 <b>Нажмите, чтобы просмотреть лог тестов (Original)</b></summary>

```text
PS C:\Users\tanis\Documents\TRPO\lab2(cblaslevel3)\build> .\cblas_tests.exe

============================================================
===       CBLAS Level 3 Interface Tests                  ===
===       Platform: Windows | OpenBLAS C Interface       ===
============================================================

--- Single-threaded Tests ---

▶ cblas_sgemm: NoTrans×NoTrans (ColMajor)               ✓ PASS
▶ cblas_sgemm: Trans×NoTrans                            ✓ PASS
▶ cblas_sgemm: NoTrans×Trans                            ✓ PASS
▶ cblas_sgemm: Trans×Trans                              ✓ PASS
▶ cblas_sgemm: alpha/beta (RowMajor)                     ✓ PASS
▶ cblas_dgemm: NoTrans×NoTrans                          ✓ PASS
▶ cblas_dgemm: Transpose                                 ✓ PASS
▶ cblas_dgemm: alpha/beta                                ✓ PASS
▶ cblas_cgemm: NoTrans×NoTrans                          ✓ PASS
▶ cblas_cgemm: ConjTrans                                 ✓ PASS
▶ cblas_cgemm: alpha/beta                                ✓ PASS
▶ cblas_zgemm: NoTrans×NoTrans                          ✓ PASS
▶ cblas_zgemm: ConjTrans                                 ✓ PASS
▶ cblas_zgemm: alpha/beta                                ✓ PASS
▶ cblas_ssymm: Left+Upper                                ✓ PASS
▶ cblas_ssymm: Left+Lower                                ✓ PASS
▶ cblas_ssymm: Right+Upper                               ✓ PASS
▶ cblas_ssymm: Right+Lower                               ✓ PASS
▶ cblas_dsymm: basic                                     ✓ PASS
▶ cblas_dsymm: alpha/beta                                ✓ PASS
▶ cblas_csymm: basic                                     ✓ PASS
▶ cblas_zsymm: basic                                     ✓ PASS
▶ cblas_chemm: Left+Upper                                ✓ PASS
▶ cblas_chemm: Left+Lower                                ✓ PASS
▶ cblas_chemm: Right                                     ✓ PASS
▶ cblas_zhemm: basic                                     ✓ PASS
▶ cblas_strmm: Left+Upper+NoTrans+NonUnit                ✓ PASS
▶ cblas_strmm: Left+Lower                                ✓ PASS
▶ cblas_strmm: Right+Upper                               ✓ PASS
▶ cblas_strmm: Trans                                     ✓ PASS
▶ cblas_strmm: Unit+Diag                                 ✓ PASS
▶ cblas_strmm: alpha                                     ✓ PASS
▶ cblas_dtrmm: basic                                     ✓ PASS
▶ cblas_ctrmm: basic                                     ✓ PASS
▶ cblas_ztrmm: basic                                     ✓ PASS
▶ cblas_strsm: Left+Upper+NoTrans+NonUnit                ✓ PASS
▶ cblas_strsm: Left+Lower                                ✓ PASS
▶ cblas_strsm: Right                                     ✓ PASS
▶ cblas_strsm: Trans                                     ✓ PASS
▶ cblas_strsm: Unit+Diag                                 ✓ PASS
▶ cblas_dtrsm: basic                                     ✓ PASS
▶ cblas_ctrsm: basic                                     ✓ PASS
▶ cblas_ztrsm: basic                                     ✓ PASS
▶ cblas_ssyrk: Upper+NoTrans                             ✓ PASS
▶ cblas_ssyrk: Upper+Trans                               ✓ PASS
▶ cblas_ssyrk: Lower+NoTrans                             ✓ PASS
▶ cblas_ssyrk: alpha/beta                                ✓ PASS
▶ cblas_dsyrk: basic                                     ✓ PASS
▶ cblas_csyrk: basic                                     ✓ PASS
▶ cblas_zsyrk: basic                                     ✓ PASS
▶ cblas_cherk: Upper+ConjTrans                           ✓ PASS
▶ cblas_cherk: Upper+NoTrans                             ✓ PASS
▶ cblas_cherk: Lower                                     ✓ PASS
▶ cblas_zherk: basic                                     ✓ PASS
▶ cblas_ssyr2k: Upper+NoTrans                            ✓ PASS
▶ cblas_ssyr2k: Upper+Trans                              ✓ PASS
▶ cblas_ssyr2k: Lower                                    ✓ PASS
▶ cblas_ssyr2k: alpha/beta                               ✓ PASS
▶ cblas_dsyr2k: basic                                    ✓ PASS
▶ cblas_csyr2k: basic                                    ✓ PASS
▶ cblas_zsyr2k: basic                                    ✓ PASS
▶ cblas_cher2k: Upper+ConjTrans                          ✓ PASS
▶ cblas_cher2k: Upper+NoTrans                            ✓ PASS
▶ cblas_cher2k: Lower                                    ✓ PASS
▶ cblas_zher2k: basic                                    ✓ PASS

--- Multi-threaded Tests (4 threads) ---

▶ cblas_sgemm: NoTrans×NoTrans (ColMajor)               ✓ PASS
▶ cblas_dgemm: NoTrans×NoTrans                          ✓ PASS
▶ cblas_ssymm: Left+Upper                                ✓ PASS
▶ cblas_strmm: Left+Upper+NoTrans+NonUnit                ✓ PASS
▶ cblas_ssyrk: Upper+NoTrans                             ✓ PASS
▶ cblas_ssyr2k: Upper+NoTrans                            ✓ PASS
▶ cblas_zgemm: NoTrans×NoTrans                          ✓ PASS

============================================================
===                    RESULTS                           ===
============================================================
  Total tests:  72
  Passed:       72
  Failed:       0
  Success rate: 100.0%
============================================================

✓ All tests PASSED!
```
</details>

---

## 2️⃣ Mock библиотека

Тестовая реализация показала значительное количество ошибок. Основная проблема наблюдается в работе с коэффициентами `alpha/beta`, комплексными числами и симметричными матрицами.

### ❌ Анализ ошибок
| Категория | Статус | Основные проблемы |
| :--- | :---: | :--- |
| **GEMM** | 🔴 Критично | Ошибки `alpha/beta` для `float`/`double`, несоответствия для `complex` |
| **SYMM/HEMM** | 🔴 Критично | Массовые несоответствия (Mismatch) во всех конфигурациях |
| **TRMM/TRSM** | 🟡 Частично | `Single` имеет ошибки, `Complex` проходит лучше |
| **SYRK/HERK** | 🔴 Критично | Ошибки во всех типах данных |
| **SYR2K/HER2K** | 🔴 Критично | Большинство тестов провалено |

### 📉 Детальная статистика Mock
- **Total tests:** 72
- **Passed:** 22
- **Failed:** 50
- **Success rate:** 30.6%

<details>
<summary>📄 <b>Нажмите, чтобы просмотреть лог тестов (Mock)</b></summary>

```text
PS C:\Users\tanis\Documents\TRPO\lab2(cblaslevel3)\build> .\cblas_tests.exe

============================================================
===       CBLAS Level 3 Interface Tests                  ===
===       Platform: Windows | OpenBLAS C Interface       ===
============================================================

--- Single-threaded Tests ---

▶ cblas_sgemm: NoTrans×NoTrans (ColMajor)               ✓ PASS
▶ cblas_sgemm: Trans×NoTrans                            ✓ PASS
▶ cblas_sgemm: NoTrans×Trans                            ✓ PASS
▶ cblas_sgemm: Trans×Trans                              ✓ PASS
▶ cblas_sgemm: alpha/beta (RowMajor)                     ✗ FAIL: Alpha/beta error
▶ cblas_dgemm: NoTrans×NoTrans                          ✓ PASS
▶ cblas_dgemm: Transpose                                 ✓ PASS
▶ cblas_dgemm: alpha/beta                                ✗ FAIL: DGEMM alpha/beta error
▶ cblas_cgemm: NoTrans×NoTrans                          ✗ FAIL: CGEMM mismatch
▶ cblas_cgemm: ConjTrans                                 ✗ FAIL: CGEMM ConjTrans error
▶ cblas_cgemm: alpha/beta                                ✗ FAIL: CGEMM alpha/beta error
▶ cblas_zgemm: NoTrans×NoTrans                          ✗ FAIL: ZGEMM mismatch
▶ cblas_zgemm: ConjTrans                                 ✗ FAIL: ZGEMM ConjTrans error
▶ cblas_zgemm: alpha/beta                                ✗ FAIL: ZGEMM alpha/beta error
▶ cblas_ssymm: Left+Upper                                ✗ FAIL: SSYMM Left+Upper mismatch
▶ cblas_ssymm: Left+Lower                                ✗ FAIL: SSYMM Left+Lower mismatch
▶ cblas_ssymm: Right+Upper                               ✗ FAIL: SSYMM Right+Upper error
▶ cblas_ssymm: Right+Lower                               ✗ FAIL: SSYMM Right+Lower error
▶ cblas_dsymm: basic                                     ✗ FAIL: DSYMM mismatch
▶ cblas_dsymm: alpha/beta                                ✗ FAIL: DSYMM alpha/beta error
▶ cblas_csymm: basic                                     ✗ FAIL: CSYMM mismatch
▶ cblas_zsymm: basic                                     ✗ FAIL: ZSYMM mismatch
▶ cblas_chemm: Left+Upper                                ✗ FAIL: CHEMM Left+Upper mismatch
▶ cblas_chemm: Left+Lower                                ✗ FAIL: CHEMM Left+Lower error
▶ cblas_chemm: Right                                     ✗ FAIL: CHEMM Right error
▶ cblas_zhemm: basic                                     ✗ FAIL: ZHEMM mismatch
▶ cblas_strmm: Left+Upper+NoTrans+NonUnit                ✗ FAIL: STRMM mismatch
▶ cblas_strmm: Left+Lower                                ✓ PASS
▶ cblas_strmm: Right+Upper                               ✓ PASS
▶ cblas_strmm: Trans                                     ✓ PASS
▶ cblas_strmm: Unit+Diag                                 ✓ PASS
▶ cblas_strmm: alpha                                     ✓ PASS
▶ cblas_dtrmm: basic                                     ✗ FAIL: DTRMM mismatch
▶ cblas_ctrmm: basic                                     ✓ PASS
▶ cblas_ztrmm: basic                                     ✓ PASS
▶ cblas_strsm: Left+Upper+NoTrans+NonUnit                ✗ FAIL: STRSM mismatch
▶ cblas_strsm: Left+Lower                                ✓ PASS
▶ cblas_strsm: Right                                     ✓ PASS
▶ cblas_strsm: Trans                                     ✓ PASS
▶ cblas_strsm: Unit+Diag                                 ✗ FAIL: STRSM Unit error
▶ cblas_dtrsm: basic                                     ✗ FAIL: DTRSM mismatch
▶ cblas_ctrsm: basic                                     ✓ PASS
▶ cblas_ztrsm: basic                                     ✓ PASS
▶ cblas_ssyrk: Upper+NoTrans                             ✗ FAIL: SSYRK Upper+NoTrans mismatch
▶ cblas_ssyrk: Upper+Trans                               ✗ FAIL: SSYRK Upper+Trans mismatch
▶ cblas_ssyrk: Lower+NoTrans                             ✗ FAIL: SSYRK Lower+NoTrans mismatch
▶ cblas_ssyrk: alpha/beta                                ✗ FAIL: SSYRK alpha/beta error
▶ cblas_dsyrk: basic                                     ✗ FAIL: DSYRK mismatch
▶ cblas_csyrk: basic                                     ✗ FAIL: CSYRK mismatch
▶ cblas_zsyrk: basic                                     ✗ FAIL: ZSYRK mismatch
▶ cblas_cherk: Upper+ConjTrans                           ✗ FAIL: CHERK Upper+ConjTrans mismatch
▶ cblas_cherk: Upper+NoTrans                             ✗ FAIL: CHERK Upper+NoTrans mismatch
▶ cblas_cherk: Lower                                     ✗ FAIL: CHERK Lower mismatch
▶ cblas_zherk: basic                                     ✗ FAIL: ZHERK mismatch
▶ cblas_ssyr2k: Upper+NoTrans                            ✗ FAIL: SSYR2K Upper+NoTrans mismatch
▶ cblas_ssyr2k: Upper+Trans                              ✗ FAIL: SSYR2K Upper+Trans mismatch
▶ cblas_ssyr2k: Lower                                    ✗ FAIL: SSYR2K Lower mismatch
▶ cblas_ssyr2k: alpha/beta                               ✗ FAIL: SSYR2K alpha/beta error
▶ cblas_dsyr2k: basic                                    ✗ FAIL: DSYR2K mismatch
▶ cblas_csyr2k: basic                                    ✗ FAIL: CSYR2K mismatch
▶ cblas_zsyr2k: basic                                    ✗ FAIL: ZSYR2K mismatch
▶ cblas_cher2k: Upper+ConjTrans                          ✓ PASS
▶ cblas_cher2k: Upper+NoTrans                            ✗ FAIL: CHER2K Upper+NoTrans mismatch
▶ cblas_cher2k: Lower                                    ✗ FAIL: CHER2K Lower mismatch
▶ cblas_zher2k: basic                                    ✓ PASS

--- Multi-threaded Tests (4 threads) ---

▶ cblas_sgemm: NoTrans×NoTrans (ColMajor)               ✓ PASS
▶ cblas_dgemm: NoTrans×NoTrans                          ✓ PASS
▶ cblas_ssymm: Left+Upper                                ✗ FAIL: SSYMM Left+Upper mismatch
▶ cblas_strmm: Left+Upper+NoTrans+NonUnit                ✗ FAIL: STRMM mismatch
▶ cblas_ssyrk: Upper+NoTrans                             ✗ FAIL: SSYRK Upper+NoTrans mismatch
▶ cblas_ssyr2k: Upper+NoTrans                            ✗ FAIL: SSYR2K Upper+NoTrans mismatch
▶ cblas_zgemm: NoTrans×NoTrans                          ✗ FAIL: ZGEMM mismatch

============================================================
===                    RESULTS                           ===
============================================================
  Total tests:  72
  Passed:       22
  Failed:       50
  Success rate: 30.6%
============================================================

✗ Some tests FAILED!
```
</details>

---
