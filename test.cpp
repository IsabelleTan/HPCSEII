//
// Created by Isabelle Tan on 04-05-16.
//

#include <iostream>
#include <complex>
#include "morton.h"
#include "quadtree.h"
#include "timer.hpp"
#include "expansion.h"
#include "complexAVX.h"

using namespace std;

// Set the value_type
typedef double value_type;
value_type epsilon = 0.000001;

// Test the initialization function
bool testInitialization(){
    int N = 100;

    // Set return bool to true
    bool result = true;

    // Some test particle coordinates
    value_type *x = new value_type[N];
    value_type *y = new value_type[N];
    value_type *mass = new value_type[N];

    // Put some values in the arrays
    initialize(N, x, y, mass);

    // Check whether the x and y values are within the range (0,1)
    for (int i = 0; i < N; ++i) {
        if(x[i] < 0 || x[i] > 1 || y[i] < 0 || y[i] > 1){
            result = false;
        }
    }

    // Release the allocated memory
    delete[] x;
    delete[] y;

    return result;
}

// Test the computation of the extent and location (xmin,ymin)
bool testComputeExtent(){
    int N = 2;
    // Set the resulting boolean to true
    bool result = true;

    // Initiate arrays
    value_type *x = new value_type[N];
    value_type *y = new value_type[N];
    x[0] = 0.1;
    x[1] = 0.9;
    y[0] = 0.05;
    y[1] = 0.85;

    // Initiate extend and location variables
    value_type xmin;
    value_type ymin;
    value_type ext;

    // Compute the extend and xmin, ymin
    computeExtent(2, x, y, xmin, ymin, ext);

    // Check whether the values are correct
    if(xmin != x[0] || ymin != y[0] || ext - 0.8 > epsilon){
        result = false;
    }

    // Release the allocated memory
    delete[] x;
    delete[] y;

    return result;
}

// Test the morton indices against a test set of 7 particles with known indices.
bool testMorton(){
    int N = 7;
    int depth = 2;
    // Set the resulting boolean to true
    bool result = true;

    // Initiate arrays
    // Create particle coordinate arrays
    value_type x [N] = {0.1, 0.4, 0.8, 0.9, 0.9, 0, 1};
    value_type y [N] = {0.3, 0.6, 0.4, 0.8, 0.1, 0, 1};
    uint32_t *index = new uint32_t[N];
    uint32_t control [N] = {8, 3, 13, 5, 15, 10, 5};

    // Compute indices
    morton(N, x, y, index, depth);

    // test the resulting indices against the control indices
    for (int i = 0; i < N; ++i) {

        if(index[i] != control[i]){
            result = false;
            cout << "Test failed for particle " << i << endl;
        }
    }

    // Print the results
    cout<< "mortonSerial: \t" << " control"<< endl;
    for (int j = 0; j < N; ++j) {
        cout << index[j] << "\t\t\t\t\t" << control[j] << endl;
    }

    // Free the memory
    delete[] index;

    // Print in case of success.
    if (result) {
        cout << "Test succeeded." << endl;
    }

    return result;
}

/*
 * A function to test the assignment of part_start and part_end to child nodes
 */
bool testAssignParticles(){
    // Don't change these variables.
    int N = 7;
    int depth = 2;
    // Set the resulting boolean to true
    bool result = true;

    // Initiate arrays
    // Create index array
    uint32_t index_0[N] = {3, 5, 5, 8, 10, 13, 15};
    uint32_t index_1[N] = {4, 5, 5, 8, 9, 10, 11};
    uint32_t index_2[N] = {3, 13, 13, 13, 14, 15, 15};

    // Make arrays containing the known true values
    int part_start_control_0[4] = {0, 1, 3, 5};
    int part_end_control_0[4] = {0, 2, 4, 6};
    int part_start_control_1[4] = {-1, 0, 3, -1};
    int part_end_control_1[4] = {-1, 2, 6, -1};
    int part_start_control_2[4] = {0, -1, -1, 1};
    int part_end_control_2[4] = {0, -1, -1, 6};

    // Make a root node
    Node root {0, // level
               0, // morton index
               0, // child_id
               0, // part_start
               6, // part_end
               0, // node mass
               0, // x center of mass
               0  // y center of mass
    };


    Node children[4];

    // Compute the indexvalue of this level so we can easily compute the morton-id's of the children
    int indexValue_level = pow(2.0, 2 * (depth - (root.level + 1)));

    // Initialize the children nodes
    Node child_0 = Node {root.level + 1, // level
                         root.morton_id, // morton index
                         -1, // child_id
                         -1, // part_start
                         -1, // part_end
                         1, // node mass
                         1, // x center of mass
                         1  // y center of mass
    };

    Node child_1 = Node {root.level + 1, // level
                         root.morton_id + indexValue_level, // morton index
                         -1, // child_id
                         -1, // part_start
                         -1, // part_end
                         1, // node mass
                         1, // x center of mass
                         1  // y center of mass
    };

    Node child_2 = Node {root.level + 1, // level
                         root.morton_id + 2 * indexValue_level, // morton index
                         -1, // child_id
                         -1, // part_start
                         -1, // part_end
                         1, // node mass
                         1, // x center of mass
                         1  // y center of mass
    };

    Node child_3 = Node {root.level + 1, // level
                         root.morton_id + 3 * indexValue_level, // morton index
                         -1, // child_id
                         -1, // part_start
                         -1, // part_end
                         1, // node mass
                         1, // x center of mass
                         1  // y center of mass
    };

    // Put the children nodes into the array
    children[0] = child_0;
    children[1] = child_1;
    children[2] = child_2;
    children[3] = child_3;


    // Test set 0
    assignParticles(&root, children, depth, index_0);

    // Check the start and end indices of this node
    for (int i = 0; i < 4; ++i) {
        if ((children[i].part_start != part_start_control_0[i]) || (children[i].part_end != part_end_control_0[i])) {
            result = false;
            cout << "Test failed for set 0 and node " << i << endl;
        }
    }

    // Test set 1
    assignParticles(&root, children, depth, index_1);

    // Check the start and end indices of this node
    for (int i = 0; i < 4; ++i) {
        if ((children[i].part_start != part_start_control_1[i]) || (children[i].part_end != part_end_control_1[i])) {
            result = false;
            cout << "Test failed for set 1 and node " << i << endl;
            cout << "children[i].part_start = " << children[i].part_start << ". children[i].part_end = " << children[i].part_end << endl;
            cout << "node " << i << " morton_id = " << children[i].morton_id << " non-inclusive upper bound = " << children[i].morton_id + indexValue_level << endl;
        }
    }

    // Test set 2
    assignParticles(&root, children, depth, index_2);

    // Check the start and end indices of this node
    for (int i = 0; i < 4; ++i) {
        if ((children[i].part_start != part_start_control_2[i]) || (children[i].part_end != part_end_control_2[i])) {
            result = false;
            cout << "Test failed for set 2 and node " << i << endl;
        }
    }

    // Print in case of success.
    if (result) {
        cout << "Test succeeded." << endl;
    }

    return result;
}

// A function to test the function of center of mass
bool testCenterOfMass(){
    bool result = true;
    // Make lists of masses and xy coordinates
    int N = 7;

    value_type xsorted[N] = {0.2, 0.25, 0.4, 0.48, 0.7, 0.72, 0.8};
    value_type ysorted[N] = {0.3, 0.2, 0.2, 0.1, 0.7, 0.4, 0.4};
    value_type masssorted[N] = {0.4, 0.3, 0.3, 0.89, 0.41, 0.1, 0.66};

    // Initiate the children
    Node children[4];

    // Initialize the children nodes
    Node child_0 = Node { 1, // level
                          0, // morton index
                          -1, // child_id
                          0, // part_start
                          2, // part_end
                          1, // node mass
                          1, // x center of mass
                          1  // y center of mass
    };

    Node child_1 = Node { 1, // level
                          0, // morton index
                          -1, // child_id
                          3, // part_start
                          4, // part_end
                          1, // node mass
                          1, // x center of mass
                          1  // y center of mass
    };

    Node child_2 = Node { 1, // level
                          0, // morton index
                          -1, // child_id
                          5, // part_start
                          5, // part_end
                          1, // node mass
                          1, // x center of mass
                          1  // y center of mass
    };

    Node child_3 = Node { 1, // level
                          0, // morton index
                          -1, // child_id
                          6, // part_start
                          6, // part_end
                          1, // node mass
                          1, // x center of mass
                          1  // y center of mass
    };

    // Put the children nodes into the array
    children[0] = child_0;
    children[1] = child_1;
    children[2] = child_2;
    children[3] = child_3;

    centerOfMass(children, xsorted, ysorted, masssorted);

    // Test the values against known true values
    value_type xcom_true[4] = {0.275, 0.5493846, 0.72, 0.8};
    value_type ycom_true[4] = {0.24, 0.2892307, 0.4, 0.4};
    value_type mass_true[4] = {1, 1.3, 0.1, 0.66};

    for (int i = 0; i < 4; ++i) {
        if (children[i].mass - mass_true[i] > 0.000001) {
            cout << "Mass of node " << i << " computed incorrectly." << endl;
            result = false;
        }
        if (children[i].xcom - xcom_true[i] > 0.000001) {
            cout << "xcom of node " << i << " computed incorrectly." << endl;
            result = false;
        }
        if (children[i].ycom - ycom_true[i] > 0.000001) {
            cout << "ycom of node " << i << " computed incorrectly." << endl;
            result = false;
        }
    }

    if (result) {
        cout << "Test succeeded." << endl;
    }
    return result;
}

/*
 * A function to test the split function
 */
bool testBuild(){
    // Test the stopping criterion (#p=7, k=8).
    bool result = true;

    int N = 7;
    int depth = 2;
    int k = 8;

    value_type x[7] = {0.1, 0.15, 0.4, 0.6, 0.7, 0.71, 0.8};
    value_type y[7] = {0.6, 0.85, 0.65, 0.3, 0.6, 0.1, 0.4};
    value_type mass[7] = {1, 2, 3, 4, 5, 6, 7};

    value_type xsorted[7];
    value_type ysorted[7];
    value_type mass_sorted[7];

    // Allocate the tree array containing all the nodes
    int maxNodes = (int) min((float) 8 * N / k, (float) pow(4.0, depth));
    Node* tree = new Node[maxNodes];


    build(x, y, mass, N, k, xsorted, ysorted, mass_sorted, tree, depth);

    // Test the resulting tree
    if (tree[0].child_id != -1) {
        cout << "Test of stopping criterion failed, the root node is not a leaf node." << endl;
        cout << "root.child_id = " << tree[0].child_id << endl;
        result = false;
    } else {
        cout << "Root node is an empty node => The test for the stopping criterion passed." << endl;
    }


    // Test building the tree for k=2
    cout << "Testing the build() function for a tree with N=7, k=2 and depth=2." << endl;
    k = 2;
    build(x, y, mass, N, k, xsorted, ysorted, mass_sorted, tree, depth);

    // Test the resulting tree
    for (int i = 0; i < 13; ++i) {
        printNode(tree[i]);
    }

    if (result) {
        cout << "Check the created nodes to validate the function." << endl;
    }

    delete[] tree;
    return result;
}

/*
 * A function that times the build function and writes the output to a file named "times.txt".
 */
value_type timeBuild(int N){
    // Prepare parameters
    int depth = 16;
    int k = 8;

    value_type *x = new value_type[N];
    value_type *y = new value_type[N];
    value_type *mass = new value_type[N];

    initialize(N, x, y, mass);

    value_type *xsorted = new value_type[N];
    value_type *ysorted = new value_type[N];
    value_type *mass_sorted = new value_type[N];

    // Allocate the tree array containing all the nodes
    int maxNodes = (int) min((float) 8 * N / k, (float) pow(4.0, depth));
    Node* tree = new Node[maxNodes];

    cout << "Building the tree for N = " << N << endl;

    timer t;
    t.start();
    build(x, y, mass, N, k, xsorted, ysorted, mass_sorted, tree, depth);
    t.stop();
    cout << "N = " << N << " \tElapsed time: " << t.get_timing() << endl;


    delete[] tree;
    delete[] x;
    delete[] y;
    delete[] mass;

    delete[] xsorted;
    delete[] ysorted;
    delete[] mass_sorted;


    return t.get_timing();
}


/*
 * A function to test and time the p2e kernel
 */

value_type timep2e(int N, int order){
    timer t;

    // Create arrays with particles
    value_type *x = new value_type[N];
    value_type *y = new value_type[N];
    value_type *mass = new value_type[N];

    value_type xCom = 0.5;
    value_type yCom = 0.5;

    initialize(N, x, y, mass);

    value_type *expansion = new value_type[2*order]{0};

    t.start();
    p2e(x, y, mass, N, order, xCom, yCom, expansion, expansion + order);
    t.stop();

    delete[] x;
    delete[] y;
    delete[] mass;
    delete[] expansion;

    return t.get_timing();
}

/*
 * A function to test the p2e kernel
 */
bool testp2e(){
    bool result = true;
    int order = 2;
    int N = 4;

    value_type x[4] = {1, 2, 3, 4};
    value_type y[4] = {4, 3, 2, 1};
    value_type mass[4] = {5, 6, 7, 8};

    value_type xCom= 2.5;
    value_type yCom = 2.5;

    value_type expansion[4] = {0};

    p2e(x, y, mass, N, order, xCom, yCom, expansion, expansion + order);

    value_type control[4] = {-70, -25, -60, -130};

    for (int i = 0; i < 4; ++i) {
        if (abs(expansion[i] - control[i]) > epsilon){
            cout << "expansion[" << i << "] = " << expansion[i] << " != control[" << i << "] = " << control[i] << endl;
            result = false;
        } else {

            cout << "expansion[" << i << "] = " << expansion[i] << " = control[" << i << "] = " << control[i] << endl;
        }
    }

    if (result){
        cout << "Test succeeded" << endl;
    }

    return result;
}

/*
 * A function to test and time the e2p kernel
 */
value_type timee2p(int N, int order){
    timer t;

    // Create arrays with particles
    value_type *x = new value_type[N];
    value_type *y = new value_type[N];
    value_type *mass = new value_type[N];

    value_type xCom = 0.5;
    value_type yCom = 0.5;

    initialize(N, x, y, mass);

    value_type *expansion = new value_type[2*order]{0};

    p2e(x, y, mass, N, order, xCom, yCom, expansion, expansion + order);

    //Prepare arguments
    value_type xtarget= 3;
    value_type ytarget = 1;
    value_type q = 0;

    // Compute g (sum of masses)
    for (int i = 0; i < N; ++i) {
        q+= mass[i];
    }

    t.start();
    value_type streamfunction = e2p(xtarget, ytarget, q, order, expansion, expansion + order);
    t.stop();


    delete[] x;
    delete[] y;
    delete[] mass;
    delete[] expansion;

    return t.get_timing();
}

/*
 * A function to test the e2p kernel
 */
bool teste2p(){
    bool result = true;
    int order = 2;
    int N = 4;

    value_type x[4] = {1, 2, 3, 4};
    value_type y[4] = {4, 3, 2, 1};
    value_type mass[4] = {5, 6, 7, 8};

    value_type xCom= 2.5;
    value_type yCom = 2.5;

    value_type x_target = 3;
    value_type y_target = 1;

    value_type expansion[2*order] = {0};

    p2e(x, y, mass, N, order, xCom, yCom, expansion, expansion + order);

    value_type streamfunction = e2p(x_target, y_target, 26, order, expansion, expansion + order);

    value_type control = -6.8663938;

    if (abs(streamfunction - control) > epsilon){
        cout << "streamfunction = " << streamfunction << " != control = " << control << endl;
        result = false;
    } else {

        cout << "streamfunction = " << streamfunction << " = control = " << control << endl;
    }


    if (result){
        cout << "Test succeeded" << endl;
    }

    return result;
}

/*
 * A function for the convergence analysis of the p2e+e2p kernels
 */
void convergenceAnalysis(int p_max){
    int N = 100000;

    // Create arrays with particles
    value_type *x = new value_type[N];
    value_type *y = new value_type[N];
    value_type *mass = new value_type[N];

    initialize(N, x, y, mass);

    //Prepare arguments
    value_type xtarget= 3;
    value_type ytarget = 1;
    value_type q = 0;
    value_type xCom = 0;
    value_type yCom = 0;

    // Compute q (sum of masses) and xCom and yCom (center of mass)
    for (int i = 0; i < N; ++i) {
        q+= mass[i];
        xCom+= mass[i]*x[i];
        yCom+= mass[i]*y[i];
    }
    xCom/=q;
    yCom/=q;


    // Compute the streamfunction for different orders and print the resulting value
    for (int p = 3; p <= p_max; ++p) {
        value_type *expansion = new value_type[2*p]{0};

        p2e(x, y, mass, N, p, xCom, yCom, expansion, expansion + p);

        cout << e2p(xtarget, ytarget, q, p, expansion, expansion + p) << endl;

        delete[] expansion;
    }


    delete[] x;
    delete[] y;
    delete[] mass;
}

/*
 * A function to time the p2p kernel
 */
value_type timep2p(int N){
    timer t;

    // Create arrays with particles
    value_type *x = new value_type[N];
    value_type *y = new value_type[N];
    value_type *mass = new value_type[N];

    value_type xtarget = 3;
    value_type ytarget = 1;

    initialize(N, x, y, mass);

    t.start();
    value_type streamfunction = p2p(x, y, mass, N, xtarget, ytarget);
    t.stop();

    delete[] x;
    delete[] y;
    delete[] mass;

    return t.get_timing();

}

/*
 * A function to test the p2p kernel
 */
bool testp2p(){
    bool result = true;
    int N = 4;

    value_type x[4] = {1, 2, 3, 4};
    value_type y[4] = {4, 3, 2, 1};
    value_type mass[4] = {5, 6, 7, 8};

    value_type x_target = 3;
    value_type y_target = 1;


    value_type streamfunction = p2p(x, y, mass, N, x_target, y_target);

    value_type control = 11.240687;

    if (abs(streamfunction - control) > epsilon){
        cout << "streamfunction = " << streamfunction << " != control = " << control << endl;
        result = false;
    } else {

        cout << "streamfunction = " << streamfunction << " = control = " << control << endl;
    }


    if (result){
        cout << "Test succeeded" << endl;
    }

    return result;
}

/*
 * A function to test the complex multiplication function.
 */
bool testMultiply() {
    bool result = true;

    int N = 5;
    void *voidPtr;

    // Allocate aligned data
    posix_memalign(&voidPtr, 64, N * sizeof(value_type));
    value_type *x1 = static_cast<value_type *>(voidPtr);

    posix_memalign(&voidPtr, 64, N * sizeof(value_type));
    value_type *x2 = static_cast<value_type *>(voidPtr);

    posix_memalign(&voidPtr, 64, N * sizeof(value_type));
    value_type *y1 = static_cast<value_type *>(voidPtr);

    posix_memalign(&voidPtr, 64, N * sizeof(value_type));
    value_type *y2 = static_cast<value_type *>(voidPtr);

    // Initialize the data
    for (int i = 0; i < N; ++i) {
        x1[i] = i;
        x2[i] = i + 8;
        y1[i] = 1;
        y2[i] = 2;
    }

    // Allocate space for the resulting values
    posix_memalign(&voidPtr, 64, N * sizeof(value_type));
    value_type *res_r = static_cast<value_type *>(voidPtr);
    posix_memalign(&voidPtr, 64, N * sizeof(value_type));
    value_type *res_i = static_cast<value_type *>(voidPtr);

    int SIMD_width = sizeof(__m256d) / sizeof(value_type);
    int SIMD_blocks = N / SIMD_width;

    for (int k = 0; k < SIMD_blocks; ++k) {
        // Put the numbers in a register
        __m256d x1_avx = _mm256_load_pd(x1 + k * SIMD_width);
        __m256d y1_avx = _mm256_load_pd(y1 + k * SIMD_width);
        __m256d x2_avx = _mm256_load_pd(x2 + k * SIMD_width);
        __m256d y2_avx = _mm256_load_pd(y2 + k * SIMD_width);
        __m256d res_r_avx;
        __m256d res_i_avx;

        // Multiply
        multiply(x1_avx, y1_avx, x2_avx, y2_avx, &res_r_avx, &res_i_avx);

        // Store the result in the res_r and res_i vectors
        _mm256_store_pd(res_r + k * SIMD_width, res_r_avx);
        _mm256_store_pd(res_i + k * SIMD_width, res_i_avx);
    }

    // Compute left-over indices
    for (int l = SIMD_width * SIMD_blocks; l < N; ++l) {
        res_r[l] = x1[l] * x2[l] - y1[l] * y2[l];
        res_i[l] = x1[l] * y2[l] + x2[l] * y1[l];
    }

    double control_r[N] = {-2, 7, 18, 31, 46};
    double control_i[N] = {8, 11, 14, 17, 20};

    for (int j = 0; j < N; ++j) {
        cout << res_r[j] << " " << res_i[j] << endl;
        if (!(res_r[j] == control_r[j] && res_i[j] == control_i[j])) {
            result = false;
        }
    }

    if (result) {
        cout << "Test succeeded!" << endl;
    }

    return result;
}

/*
 * A function to time the multiply() function.
 */
value_type timeMultiply(int N) {
    timer t;

    // INITIALIZATION
    void *voidPtr;

    // Allocate aligned data
    posix_memalign(&voidPtr, 64, N * sizeof(value_type));
    value_type *x1 = static_cast<value_type *>(voidPtr);

    posix_memalign(&voidPtr, 64, N * sizeof(value_type));
    value_type *x2 = static_cast<value_type *>(voidPtr);

    posix_memalign(&voidPtr, 64, N * sizeof(value_type));
    value_type *y1 = static_cast<value_type *>(voidPtr);

    posix_memalign(&voidPtr, 64, N * sizeof(value_type));
    value_type *y2 = static_cast<value_type *>(voidPtr);

    posix_memalign(&voidPtr, 64, N * sizeof(value_type));
    value_type *res_r = static_cast<value_type *>(voidPtr);

    posix_memalign(&voidPtr, 64, N * sizeof(value_type));
    value_type *res_i = static_cast<value_type *>(voidPtr);

    // Fill the arrays with data
    initialize(N, x1, x2, y1);
    initialize(N, y2, res_r, res_i);

    t.start();
    int SIMD_width = sizeof(__m256d) / sizeof(value_type);
    int SIMD_blocks = N / SIMD_width;
    //cout << "SIMD_width = " << SIMD_width << endl;
    //cout << "SIMD_blocks = " << SIMD_blocks << endl;

    for (int i = 0; i < SIMD_blocks; ++i) {
        // Load the values into a register
        __m256d x1_avx = _mm256_load_pd(x1 + i * SIMD_width);
        __m256d y1_avx = _mm256_load_pd(y1 + i * SIMD_width);
        __m256d x2_avx = _mm256_load_pd(x2 + i * SIMD_width);
        __m256d y2_avx = _mm256_load_pd(y2 + i * SIMD_width);
        __m256d res_r_avx;
        __m256d res_i_avx;

        // Perform multiplication
        multiply(x1_avx, y1_avx, x2_avx, y2_avx, &res_r_avx, &res_i_avx);

        // Write the results to an array
        _mm256_store_pd(res_r + i * SIMD_width, res_r_avx);
        _mm256_store_pd(res_i + i * SIMD_width, res_i_avx);
    }

    // Perform multiplications of left-over elements
    for (int j = SIMD_width * SIMD_blocks; j < N; ++j) {
        res_r[j] = x1[j] * x2[j] - y1[j] * y2[j];
        res_i[j] = x1[j] * y2[j] + x2[j] * y1[j];
    }
    t.stop();

    return t.get_timing();
}

/*
 * A function to time the regular multiplication function for complex numbers.
 */
value_type timeMultiplySTD(int N) {
    timer t;

    // Create the arrays
    complex<value_type> *x = new complex<value_type>[N];
    complex<value_type> *y = new complex<value_type>[N];
    complex<value_type> *res = new complex<value_type>[N];

    // Perform the N multiplications
    t.start();
    for (int i = 0; i < N; ++i) {
        res[i] = x[i] * y[i];
    }
    t.stop();

    return t.get_timing();
}

/*
 * A function to time something N times and compute the average and variance
 */
void time(int N){
    value_type *times = new value_type[N];

    // Track the time N times
    for (int i = 0; i < N; ++i) {
        times[i] = timeMultiply(1000000);
    }

    // Compute the average
    value_type mean = 0;
    for (int j = 0; j < N; ++j) {
        mean += times[j];
    }
    mean/=N;

    // Compute the variance
    value_type var = 0;
    for (int k = 0; k < N; ++k) {
        var += pow(times[k] - mean, 2);
    }
    var/=N;

    // Print the results
    cout << "Timed " << N << " times." << endl;
    cout << "Mean = " << mean << endl;
    cout << "Variance = " << var << " is " << var/mean * 100 << "% of the mean" << endl;

    delete[] times;
}