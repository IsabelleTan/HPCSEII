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
void multiply(int N, value_type *x1, value_type *y1, value_type *x2, value_type *y2, value_type *res_r,
              value_type *res_i);


#endif //SERIAL_COMPLEXAVX_H
