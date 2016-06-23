//
// Created by Isabelle Tan on 04-05-16.
//

#include <algorithm>

#ifndef EXERCISE3_MORTONPARALLEL_H
#define EXERCISE3_MORTONPARALLEL_H

// Set the value_type
typedef double value_type;

// Fill the arrays x and y with N random values between 0 and 1
void initialize(const int N, value_type* x, value_type* y, value_type* mass);

// Compute the size (ext) and location (xmin, ymin) of the root node
void computeExtent(const int N, const value_type* const x, const value_type* const y, value_type& xmin, value_type& ymin, value_type& ext);

// Compute the quadrant of a particle within a window given by (ext, xmin, ymin)
int quadrant(value_type xmin, value_type ymin, value_type ext, value_type x, value_type y);

// Compute the Morton indices for all particles with position (x[],y[]) and store the indices in index
void morton(const int N, const value_type* const x, const value_type* const y, unsigned int* index, int depth);

// Sort the array index and sort the array keys along with the same permutation
void sortit(const int N, uint32_t* index, int* keys);

// Reorder the arrays xsorted and ysorted according to the permutation described by keys
void reorder(const int N, const int* const keys, const value_type* const x, const value_type* const y, const value_type* const mass, value_type *xsorted, value_type *ysorted, value_type* masssorted);



#endif //EXERCISE3_MORTONPARALLEL_H
