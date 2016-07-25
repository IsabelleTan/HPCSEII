//
// Created by Isabelle Tan on 11-07-16.
//

#ifndef SERIAL_COMPLEXAVX_H
#define SERIAL_COMPLEXAVX_H

// Set the value_type
typedef double value_type;

/*
 * A vectorized multiplication function for complex numbers. Assumes aligned arrays as input.
 */
void multiply(__m256d x1, __m256d y1, __m256d x2, __m256d y2, __m256d *res_r, __m256d *res_i);


#endif //SERIAL_COMPLEXAVX_H
