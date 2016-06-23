//
// Created by Isabelle Tan on 17-05-16.
//

#include "timer.hpp"
#include <iostream>
#include <fstream>
#include "morton.h"

using namespace std;

// Set the value_type
typedef double value_type;

void timeMorton(int N, int depth){
    // Allocate the timer and output file to measure the computation times of the functions
    timer t;
    ofstream timesFile;

    // Some test particle coordinates
    value_type *x = new value_type[N];
    value_type *y = new value_type[N];
    value_type *mass = new value_type[N];
    // Put some values in the arrays
    initialize(N, x, y, mass);

    // Allocate the index array and compute and store the morton indices in it.
    uint32_t *index = new uint32_t[N];

    t.start();
    morton(N, x, y, index, depth);
    t.stop();
    timesFile <<  t.get_timing() << "\t";

    // Print particles and their morton indices
    /* printf("\n x \t y \t index \n");
    for (int i = 0; i < N; ++i) {
        printf("\n %f \t %f \t %u", x[i], y[i], index[i]);
    }*/

    // Allocation and instantiation of the keys array which we will use to store the Morton index
    // sorting permutation in.
    int *keys = new int[N];
    for (int j = 0; j < N; ++j) {
        keys[j] = j;
    }

    // Sort the array index and store the total permutation in keys
    t.start();
    sortit(N, index, keys);
    t.stop();
    timesFile <<  t.get_timing() << "\t";


    // Allocate the new ordered coordinate arrays
    value_type *xsorted = new value_type[N];
    value_type *ysorted = new value_type[N];
    value_type *masssorted = new value_type[N];

    // Put the coordinates sorted by Morton index in the sorted arrays
    t.start();
    reorder(N, keys, x, y, mass, xsorted, ysorted, masssorted);
    t.stop();
    timesFile <<  t.get_timing() << "\t";

    // Print the sorted coordinate so we can test them
    /* printf("\n \n Sorted particles: \n");
    for (int i = 0; i < N; ++i) {
        printf("%f, \t %f; \n", xsorted[i], ysorted[i]);
    } */

    // Release the memory that we allocated.
    delete[] x;
    delete[] y;
    delete[] mass;
    delete[] xsorted;
    delete[] ysorted;
    delete[] masssorted;
    delete[] index;
    delete[] keys;
}
