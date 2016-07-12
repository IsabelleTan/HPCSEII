//
// Created by Isabelle Tan on 11-07-16.
//

#include <immintrin.h>
#include "complexAVX.h"
#include <iostream>

using namespace std;

void multiply(int N, value_type *x1, value_type *y1, value_type *x2, value_type *y2, value_type *res_r,
              value_type *res_i) {
    // Compute the whole division and the remainder (assuming AVX2 instruction set which supports 256 bit operations)
    int n_vec = 32 / sizeof(value_type);
    int div = N / n_vec;
    int rem = N - div * n_vec;

    for (int i = 0; i < div; ++i) {
        // Load the data into the registers
        __m256d x1_avx = _mm256_load_pd(x1 + 4 * i);
        __m256d y1_avx = _mm256_load_pd(y1 + 4 * i);
        __m256d x2_avx = _mm256_load_pd(x2 + 4 * i);
        __m256d y2_avx = _mm256_load_pd(y2 + 4 * i);


        // Compute the real part
        __m256d res_yy = _mm256_mul_pd(y1_avx, y2_avx);
        __m256d real = _mm256_fmsub_pd(x1_avx, x2_avx, res_yy);

        // Compute the imaginary part
        __m256d res_x1y2 = _mm256_mul_pd(x1_avx, y2_avx);
        __m256d imag = _mm256_fmadd_pd(x2_avx, y1_avx, res_x1y2);

        // Store the result into res_r and res_i
        _mm256_store_pd(res_r + 4 * i, real);
        _mm256_store_pd(res_i + 4 * i, imag);

        // Compute the remaining multiplications
        for (int i = div * n_vec; i < N; i++) {
            res_r[i] = x1[i] * x2[i] - y1[i] * y2[i];
            res_i[i] = x1[i] * y2[i] + x2[i] * y1[i];
        }
    }
}