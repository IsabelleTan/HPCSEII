//
// Created by Isabelle Tan on 04-05-16.
//
#include <iostream>
#include "test.h"
#include "morton.h"
#include "expansion.h"

using namespace std;

// Set the value_type
typedef double value_type;

/* int main() {
        value_type times[4][20] = {};
        value_type means[4] = {};
        value_type var[4] = {};


        for (int i = 0; i < 20; ++i) {
            cout << "Iteration i: " << i << endl;
            times[0][i] = timeBuild(1000);
            means[0] += times[0][i];
            times[1][i] = timeBuild(10000);
            means[1] += times[1][i];
            times[2][i] = timeBuild(100000);
            means[2] += times[2][i];
            times[3][i] = timeBuild(1000000);
            means[3] += times[3][i];
        }

        for (int j = 0; j < 4; ++j) {
            means[j] /= 20;
            for (int i = 0; i < 20; ++i) {
                var[j] += pow(means[j] - times[j][i], 2);
            }
            var[j] /= 20;

            cout << "\nN = " << pow(10, j) * 1000 << " mean: " << means[j] << " var: " << var[j] << endl;
        }

    cout << "Finished!" << endl;
    return 0;
} */

int main() {
    time(20);

    return 0;
}