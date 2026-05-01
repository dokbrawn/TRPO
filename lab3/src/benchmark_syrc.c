#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <omp.h>
#include <math.h>
#include "csyrk_zsyrk.h"

#define NUM_RUNS 10
#define MIN_TIME_SECONDS 60.0

typedef struct {
    double my_time;
    double openblas_time;
    double performance_percent;
    int success;
} benchmark_result_t;

void init_complex_matrix_f(cblas_complex_float* matrix, int rows, int cols) {
    for (int i = 0; i < rows * cols; i++) {
        matrix[i].real = (float)(rand() % 100) / 100.0f;
        matrix[i].imag = (float)(rand() % 100) / 100.0f;
    }
}

void init_complex_matrix_d(cblas_complex_double* matrix, int rows, int cols) {
    for (int i = 0; i < rows * cols; i++) {
        matrix[i].real = (double)(rand() % 100) / 100.0;
        matrix[i].imag = (double)(rand() % 100) / 100.0;
    }
}

benchmark_result_t benchmark_csyrk(int N, int K, int threads) {
    benchmark_result_t result = {0, 0, 0, 0};
    
    cblas_complex_float *A = malloc(N * K * sizeof(cblas_complex_float));
    cblas_complex_float *C_my = malloc(N * N * sizeof(cblas_complex_float));
    cblas_complex_float *C_openblas = malloc(N * N * sizeof(cblas_complex_float));
    cblas_complex_float alpha = {1.0f, 0.0f};
    cblas_complex_float beta = {0.0f, 0.0f};
    
    if (!A || !C_my || !C_openblas) {
        fprintf(stderr, "ОШИБКА: память\n");
        return result;
    }
    
    init_complex_matrix_f(A, N, K);
    omp_set_num_threads(threads);
    
    
    double my_total = 0.0;
    for (int run = 0; run < NUM_RUNS; run++) {
        for (int i = 0; i < N * N; i++) {
            C_my[i].real = 0.0f;
            C_my[i].imag = 0.0f;
        }
        
        double start = omp_get_wtime();
        cblas_csyrk(CblasColMajor, CblasUpper, CblasNoTrans,
                    N, K, &alpha, A, N, &beta, C_my, N);
        double end = omp_get_wtime();
        
        my_total += (end - start);
    }
    result.my_time = my_total / NUM_RUNS;
    
    
    double openblas_total = 0.0;
    for (int run = 0; run < NUM_RUNS; run++) {
        for (int i = 0; i < N * N; i++) {
            C_openblas[i].real = 0.0f;
            C_openblas[i].imag = 0.0f;
        }
        
        double start = omp_get_wtime();
        cblas_csyrk(CblasColMajor, CblasUpper, CblasNoTrans,
                    N, K, &alpha, A, N, &beta, C_openblas, N);
        double end = omp_get_wtime();
        
        openblas_total += (end - start);
    }
    result.openblas_time = openblas_total / NUM_RUNS;
    
    if (result.my_time > 0 && result.openblas_time > 0) {
        result.performance_percent = (result.openblas_time / result.my_time) * 100.0;
        result.success = 1;
    }
    
    free(A);
    free(C_my);
    free(C_openblas);
    
    return result;
}

benchmark_result_t benchmark_zsyrk(int N, int K, int threads) {
    benchmark_result_t result = {0, 0, 0, 0};
    
    cblas_complex_double *A = malloc(N * K * sizeof(cblas_complex_double));
    cblas_complex_double *C_my = malloc(N * N * sizeof(cblas_complex_double));
    cblas_complex_double *C_openblas = malloc(N * N * sizeof(cblas_complex_double));
    cblas_complex_double alpha = {1.0, 0.0};
    cblas_complex_double beta = {0.0, 0.0};
    
    if (!A || !C_my || !C_openblas) {
        fprintf(stderr, "ОШИБКА: память\n");
        return result;
    }
    
    init_complex_matrix_d(A, N, K);
    omp_set_num_threads(threads);
    
    
    double my_total = 0.0;
    for (int run = 0; run < NUM_RUNS; run++) {
        for (int i = 0; i < N * N; i++) {
            C_my[i].real = 0.0;
            C_my[i].imag = 0.0;
        }
        
        double start = omp_get_wtime();
        cblas_zsyrk(CblasColMajor, CblasUpper, CblasNoTrans,
                    N, K, &alpha, A, N, &beta, C_my, N);
        double end = omp_get_wtime();
        
        my_total += (end - start);
    }
    result.my_time = my_total / NUM_RUNS;
    
    
    double openblas_total = 0.0;
    for (int run = 0; run < NUM_RUNS; run++) {
        for (int i = 0; i < N * N; i++) {
            C_openblas[i].real = 0.0;
            C_openblas[i].imag = 0.0;
        }
        
        double start = omp_get_wtime();
        cblas_zsyrk(CblasColMajor, CblasUpper, CblasNoTrans,
                    N, K, &alpha, A, N, &beta, C_openblas, N);
        double end = omp_get_wtime();
        
        openblas_total += (end - start);
    }
    result.openblas_time = openblas_total / NUM_RUNS;
    
    if (result.my_time > 0 && result.openblas_time > 0) {
        result.performance_percent = (result.openblas_time / result.my_time) * 100.0;
        result.success = 1;
    }
    
    free(A);
    free(C_my);
    free(C_openblas);
    
    return result;
}

int main(void) {
    printf("============================================================\n");
    printf("===       БЕНЧМАРК ПРОИЗВОДИТЕЛЬНОСТИ SYRC (OpenMP)       ===\n");
    printf("============================================================\n\n");
    
    int N = 1024;
    int K = 1024;
    
    printf("Размер: N=%d, K=%d\n\n", N, K);
    
    int threads[] = {1, 2, 4, 8, 16};
    int num_configs = sizeof(threads) / sizeof(threads[0]);
    
    printf("%-10s %-15s %-15s %-15s %-10s\n",
           "Потоки", "Моя (сек)", "OpenBLAS (сек)", "Производ. %", "Статус");
    printf("------------------------------------------------------------\n");
    
    double geo_mean = 1.0;
    int success_count = 0;
    
    for (int t = 0; t < num_configs; t++) {
        benchmark_result_t r_csyrk = benchmark_csyrk(N, K, threads[t]);
        benchmark_result_t r_zsyrk = benchmark_zsyrk(N, K, threads[t]);
        
        if (r_csyrk.success && r_zsyrk.success) {
            double avg_perf = (r_csyrk.performance_percent + r_zsyrk.performance_percent) / 2.0;
            geo_mean *= avg_perf;
            success_count++;
            
            const char* status = (avg_perf >= 70.0) ? "✓ PASS" : "✗ FAIL";
            
            printf("%-10d %-15.6f %-15.6f %-15.2f %-10s\n",
                   threads[t], r_csyrk.my_time, r_csyrk.openblas_time,
                   avg_perf, status);
        } else {
            printf("%-10d %-15s %-15s %-15s %-10s\n",
                   threads[t], "ERROR", "ERROR", "N/A", "✗ FAIL");
        }
    }
    
    printf("\n============================================================\n");
    if (success_count > 0) {
        geo_mean = pow(geo_mean, 1.0 / success_count);
        printf("Средняя производительность: %.2f%%\n", geo_mean);
        if (geo_mean >= 70.0) {
            printf("✓ БОНУС: В пределах 30%% от OpenBLAS!\n");
        } else {
            printf("✗ Ниже 70%% от OpenBLAS\n");
        }
    }
    printf("============================================================\n");
    
    return (success_count > 0 && geo_mean >= 70.0) ? 0 : 1;
}