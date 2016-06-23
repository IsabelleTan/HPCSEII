//
// Created by Isabelle Tan on 20-06-16.
//

#include "expansion.h"
#include <iostream>
#include <complex.h>

using namespace std;

/*
 * A kernel to compute the particle to expansion computations
 */
void p2e(const value_type* const x, const value_type* const y, const value_type* const q, const int nsources, const int order, const value_type xcom, const value_type ycom, value_type* const rexpansion, value_type* const iexpansion){
    // Make arrays to store the radius of each particle and the values for different orders of r that are used in the multiplication
    complex<value_type>* z = new complex<value_type>[nsources];
    complex<value_type>* r = new complex<value_type>[nsources];
    for (int i = 0; i < nsources; ++i) {
        r[i] = -q[i];
        z[i] = complex<value_type>(x[i],y[i]) - complex<value_type>(xcom, ycom);
    }

    // Compute the alpha's
    for (int j = 0; j < order; ++j) {
        rexpansion[j] = 0;
        iexpansion[j] = 0;
        for (int k = 0; k < nsources; ++k) {
            // Multiply with z to obtain the next power of z
            r[k] *= complex<value_type>(x[k], y[k]);

            // Write real and imaginary part to the arrays containing the expansion alpha's
            rexpansion[j] += r[k].real();
            iexpansion[j] += r[k].imag();
        }

        // Divide the whole sum by j
        rexpansion[j]/=(j+1);
        iexpansion[j]/=(j+1);
    }

    delete[] z;
    delete[] r;
}

/*
 * A kernel to compute the expansion to target particle expansions
 */
value_type e2p(const double xtarget, const double ytarget, const double q, const int order, const double* const rxps, const double* const ixps){
    // Make an array to store the different orders of z_target
    complex<value_type> orders(1,0);
    complex<value_type> z_target(xtarget, ytarget);
    complex<value_type> streamFunction = q * log(z_target);

    for (int i = 0; i < order; ++i) {
        orders *= z_target;
        streamFunction += complex<value_type>(rxps[i], ixps[i]) / orders;
    }

    return streamFunction.real();
}

/*
 * A kernel to compute the interaction between nsources source particles and one target particle
 */
value_type p2p(const value_type* const xsources, const value_type* const ysources, const value_type* const q, const int nsources, const value_type xtarget, const value_type ytarget){
    value_type streamfunction = 0;
    for (int i = 0; i < nsources; ++i) {
        streamfunction += q[i] * log( sqrt( pow(xsources[i] - xtarget,2) + pow(ysources[i] - ytarget,2) ) );
    }
    return streamfunction;
}