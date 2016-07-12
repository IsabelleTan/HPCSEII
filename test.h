//
// Created by Isabelle Tan on 04-05-16.
//

#ifndef EXERCISE3_TEST_H
#define EXERCISE3_TEST_H

#include "quadtree.h"

// Test the array sizes, and the particle domain.
bool testInitialization();

// Test the computation of the size and location of the root node
bool testComputeExtent();

// Test the assignment of the morton methods
bool testMorton();

// Test the assignment of particle groups to children nodes
bool testAssignParticles();

// Test the computation of the center of mass of children nodes
bool testCenterOfMass();

// Test the function split()
bool testBuild();

// Time the function and write the elapsed time to a file named "times.txt"
value_type timeBuild(int N);

// A function to time the p2e kernel
value_type timep2e(int N, int order);

// A function to test the p2e kernel
bool testp2e();

//A function to time the e2p kernel
value_type timee2p(int N, int order);

// A function to test the e2p kernel
bool teste2p();

//A function for the convergence analysis of the p2e+e2p kernels
void convergenceAnalysis(int p_max);

// A function to time the p2p kernel
value_type timep2p(int N);

// A function to test the p2p kernel
bool testp2p();

// A function to test the complex multiplication function.
bool testMultiply();

// A function to time N times that computes the mean and variance
void time(int N);


#endif //EXERCISE3_TEST_H
