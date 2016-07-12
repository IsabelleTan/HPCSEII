//
// Created by Isabelle Tan on 11-07-16.
//

#include <immintrin.h>
#include "complex.h"
#include <iostream>

using namespace std;

void multiply(int N, value_type *x1, value_type *y1, value_type *x2, value_type *y2, value_type *res_r,
              value_type *res_i) {
    // Compute the whole division and the remainder (assuming AVX2 instruction set which supports 256 bit operations)
    int n_vec = 32 / sizeof(value_type);
    int div = N / n_vec;
    int rem = N - div;

    for (int i = 0; i < div; ++i) {
        // Load the data into the registers
        __m256d x1_avx = _mm256_load_pd(x1 + 4 * i);
        __m256d y1_avx = _mm256_load_pd(y1 + 4 * i);
        __m256d x2_avx = _mm256_load_pd(x2 + 4 * i);
        __m256d y2_avx = _mm256_load_pd(y2 + 4 * i);


        // Compute the real part
        __m256d res_xx = _mm256_mul_pd(x1_avx, x2_avx);
        __m256d res_yy = _mm256_mul_pd(y1_avx, y2_avx); // TODO fused multiply-substract possible?
        __m256d real = _mm256_sub_pd(res_xx, res_yy);

        // Compute the imaginary part
        __m256d res_x1y2 = _mm256_mul_pd(x1_avx, y2_avx);
        __m256d res_x2y1 = _mm256_mul_pd(x2_avx, y1_avx);
        __m256d imag = _mm256_add_pd(res_x1y2, res_x2y1);

        // Store the result into res_r and res_i
        _mm256_store_pd(res_r + 4 * i, real);
        _mm256_store_pd(res_i + 4 * i, imag);

        // TODO treat the remaining array entries
    }
}