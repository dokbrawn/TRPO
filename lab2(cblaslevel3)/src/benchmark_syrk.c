#include "syrk.h"

#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef void (*openblas_set_num_threads_fn)(int);
typedef void (*cblas_ssyrk_fn)(CBLAS_ORDER, CBLAS_UPLO, CBLAS_TRANSPOSE, int, int, float, const float *, int, float, float *, int);
typedef void (*cblas_dsyrk_fn)(CBLAS_ORDER, CBLAS_UPLO, CBLAS_TRANSPOSE, int, int, double, const double *, int, double, double *, int);

typedef struct {
    void *handle;
    openblas_set_num_threads_fn set_threads;
    cblas_ssyrk_fn ssyrk;
    cblas_dsyrk_fn dsyrk;
    const char *library_name;
} openblas_api_t;

typedef struct {
    double min_seconds;
    double avg_seconds;
    double gflops;
} bench_result_t;

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static uint32_t lcg_next(uint32_t *state) {
    *state = (*state * 1664525u) + 1013904223u;
    return *state;
}

static void fill_float(float *data, size_t count) {
    uint32_t state = 0x12345678u;
    for (size_t i = 0; i < count; ++i) {
        data[i] = (float)(lcg_next(&state) % 2001u) / 1000.0f - 1.0f;
    }
}

static void fill_double(double *data, size_t count) {
    uint32_t state = 0x87654321u;
    for (size_t i = 0; i < count; ++i) {
        data[i] = (double)(lcg_next(&state) % 2001u) / 1000.0 - 1.0;
    }
}

#define LOAD_SYMBOL(handle, name, out_fn)            \
    do {                                             \
        void *sym__ = dlsym((handle), (name));       \
        memcpy(&(out_fn), &sym__, sizeof(sym__));    \
    } while (0)

static openblas_api_t load_openblas(void) {
    const char *candidates[] = {
        "libopenblas.so",
        "libopenblas.so.0",
        "libopenblas64_.so",
        NULL
    };

    openblas_api_t api = {0};
    for (int i = 0; candidates[i] != NULL; ++i) {
        void *handle = dlopen(candidates[i], RTLD_LAZY);
        if (handle == NULL) {
            continue;
        }
        api.handle = handle;
        api.library_name = candidates[i];
        LOAD_SYMBOL(handle, "openblas_set_num_threads", api.set_threads);
        LOAD_SYMBOL(handle, "cblas_ssyrk", api.ssyrk);
        LOAD_SYMBOL(handle, "cblas_dsyrk", api.dsyrk);
        if (api.set_threads != NULL && api.ssyrk != NULL && api.dsyrk != NULL) {
            return api;
        }
        dlclose(handle);
        memset(&api, 0, sizeof(api));
    }
    return api;
}

static bench_result_t bench_float_impl(int n, int k, int runs,
                                       void (*fn)(CBLAS_ORDER, CBLAS_UPLO, CBLAS_TRANSPOSE, int, int, float, const float *, int, float, float *, int)) {
    const int lda = k;
    const int ldc = n;
    float *a = (float *)malloc((size_t)n * (size_t)lda * sizeof(float));
    float *c = (float *)malloc((size_t)n * (size_t)ldc * sizeof(float));
    float *seed = (float *)malloc((size_t)n * (size_t)ldc * sizeof(float));
    fill_float(a, (size_t)n * (size_t)lda);
    fill_float(seed, (size_t)n * (size_t)ldc);

    bench_result_t result = {0};
    result.min_seconds = 1e100;
    double total = 0.0;
    for (int run = 0; run < runs; ++run) {
        memcpy(c, seed, (size_t)n * (size_t)ldc * sizeof(float));
        const double t0 = now_seconds();
        fn(CblasColMajor, CblasUpper, CblasNoTrans, n, k, 1.0f, a, lda, 1.0f, c, ldc);
        const double dt = now_seconds() - t0;
        if (dt < result.min_seconds) {
            result.min_seconds = dt;
        }
        total += dt;
    }
    result.avg_seconds = total / (double)runs;
    result.gflops = (double)n * (double)(n + 1) * (double)k / result.min_seconds / 1e9;

    free(seed);
    free(c);
    free(a);
    return result;
}

static bench_result_t bench_double_impl(int n, int k, int runs,
                                        void (*fn)(CBLAS_ORDER, CBLAS_UPLO, CBLAS_TRANSPOSE, int, int, double, const double *, int, double, double *, int)) {
    const int lda = k;
    const int ldc = n;
    double *a = (double *)malloc((size_t)n * (size_t)lda * sizeof(double));
    double *c = (double *)malloc((size_t)n * (size_t)ldc * sizeof(double));
    double *seed = (double *)malloc((size_t)n * (size_t)ldc * sizeof(double));
    fill_double(a, (size_t)n * (size_t)lda);
    fill_double(seed, (size_t)n * (size_t)ldc);

    bench_result_t result = {0};
    result.min_seconds = 1e100;
    double total = 0.0;
    for (int run = 0; run < runs; ++run) {
        memcpy(c, seed, (size_t)n * (size_t)ldc * sizeof(double));
        const double t0 = now_seconds();
        fn(CblasColMajor, CblasUpper, CblasNoTrans, n, k, 1.0, a, lda, 1.0, c, ldc);
        const double dt = now_seconds() - t0;
        if (dt < result.min_seconds) {
            result.min_seconds = dt;
        }
        total += dt;
    }
    result.avg_seconds = total / (double)runs;
    result.gflops = (double)n * (double)(n + 1) * (double)k / result.min_seconds / 1e9;

    free(seed);
    free(c);
    free(a);
    return result;
}

static int parse_arg_or_default(const char *value, int fallback) {
    if (value == NULL || *value == '\0') {
        return fallback;
    }
    char *end = NULL;
    long parsed = strtol(value, &end, 10);
    if (end == value || parsed <= 0L) {
        return fallback;
    }
    return (int)parsed;
}

int main(int argc, char **argv) {
    const int runs = parse_arg_or_default(argc > 1 ? argv[1] : NULL, 10);
    const int n = parse_arg_or_default(argc > 2 ? argv[2] : NULL, 768);
    const int k = parse_arg_or_default(argc > 3 ? argv[3] : NULL, n);
    const int thread_counts[] = {1, 2, 4, 8, 16};

    openblas_api_t openblas = load_openblas();

    printf("SYRK benchmark (n=%d, k=%d, runs=%d)\n", n, k, runs);
    printf("OpenBLAS detected: %s\n", openblas.handle != NULL ? openblas.library_name : "no");

    puts("\nSingle precision:");
    puts("threads,impl,best_sec,avg_sec,gflops,relative_to_openblas_pct");
    for (size_t i = 0; i < sizeof(thread_counts) / sizeof(thread_counts[0]); ++i) {
        const int threads = thread_counts[i];
        bench_result_t custom = bench_float_impl(n, k, runs, syrk_s);
        if (openblas.handle != NULL) {
            openblas.set_threads(threads);
            bench_result_t ref = bench_float_impl(n, k, runs, openblas.ssyrk);
            const double relative = ref.min_seconds > 0.0 ? (ref.min_seconds / custom.min_seconds) * 100.0 : 0.0;
            printf("%d,custom,%0.6f,%0.6f,%0.3f,%0.2f\n", threads, custom.min_seconds, custom.avg_seconds, custom.gflops, relative);
            printf("%d,openblas,%0.6f,%0.6f,%0.3f,100.00\n", threads, ref.min_seconds, ref.avg_seconds, ref.gflops);
        } else {
            printf("%d,custom,%0.6f,%0.6f,%0.3f,n/a\n", threads, custom.min_seconds, custom.avg_seconds, custom.gflops);
        }
    }

    puts("\nDouble precision:");
    puts("threads,impl,best_sec,avg_sec,gflops,relative_to_openblas_pct");
    for (size_t i = 0; i < sizeof(thread_counts) / sizeof(thread_counts[0]); ++i) {
        const int threads = thread_counts[i];
        bench_result_t custom = bench_double_impl(n, k, runs, syrk_d);
        if (openblas.handle != NULL) {
            openblas.set_threads(threads);
            bench_result_t ref = bench_double_impl(n, k, runs, openblas.dsyrk);
            const double relative = ref.min_seconds > 0.0 ? (ref.min_seconds / custom.min_seconds) * 100.0 : 0.0;
            printf("%d,custom,%0.6f,%0.6f,%0.3f,%0.2f\n", threads, custom.min_seconds, custom.avg_seconds, custom.gflops, relative);
            printf("%d,openblas,%0.6f,%0.6f,%0.3f,100.00\n", threads, ref.min_seconds, ref.avg_seconds, ref.gflops);
        } else {
            printf("%d,custom,%0.6f,%0.6f,%0.3f,n/a\n", threads, custom.min_seconds, custom.avg_seconds, custom.gflops);
        }
    }

    if (openblas.handle != NULL) {
        dlclose(openblas.handle);
    }
    return 0;
}
