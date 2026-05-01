#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <omp.h>

#define TEST_START(name) printf("[ТЕСТ] %s... ", name)
#define TEST_PASS() do { printf("✓ PASSED\n"); return 1; } while(0)
#define TEST_FAIL(msg) do { printf("✗ FAILED: %s\n", msg); return 0; } while(0)
#define ASSERT(cond, msg) do { if (!(cond)) { TEST_FAIL(msg); } } while(0)

static inline void set_threads(int n) {
    omp_set_num_threads(n);
}

#endif