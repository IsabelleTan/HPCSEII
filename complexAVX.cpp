//
// Created by Isabelle Tan on 11-07-16.
//

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

/*
 * A function to perform the horizontal add of an AVX 256 pack of 4 doubles.
 * Taken from the website https://software.intel.com/en-us/forums/intel-isa-extensions/topic/346695
 */
double HsumAvxDbl(__m256d avx) {
    double sumAVX;

    __m256d hsum = _mm256_add_pd(avx, _mm256_permute2f128_pd(avx, avx, 0x1));
    _mm_store_sd(&sumAVX, _mm_hadd_pd(_mm256_castpd256_pd128(hsum), _mm256_castpd256_pd128(hsum)));

    return sumAVX;
}


