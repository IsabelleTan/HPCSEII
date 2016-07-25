//
// Created by Isabelle Tan on 11-07-16.
//

#include <immintrin.h>
#include "complexAVX.h"
#include <iostream>

using namespace std;


void multiply(__m256d x1, __m256d y1, __m256d x2, __m256d y2, __m256d *res_r,
              __m256d *res_i) {

    // Compute the real part
    __m256d res_yy = _mm256_mul_pd(y1, y2);
    *res_r = _mm256_fmsub_pd(x1, x2, res_yy);

    // Compute the imaginary part
    __m256d res_x1y2 = _mm256_mul_pd(x1, y2);
    *res_i = _mm256_fmadd_pd(x2, y1, res_x1y2);

}

