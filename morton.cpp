#include <iostream>
#include <algorithm>

using namespace std;


// Set the value_type
typedef double value_type;

// Put some random values in the arrays
void initialize(const int N, value_type* x, value_type* y, value_type* mass){
    for (int i = 0; i < N; ++i) {
        x[i] = drand48();
        y[i] = drand48();
        mass[i] = drand48();
    }
}

// Compute the size (ext) and location (xmin, ymin) of the root node
void computeExtent(const int N, const value_type* const x, const value_type* const y, value_type& xmin, value_type& ymin, value_type& ext){
    // Prepare the min and max variables
    value_type xmax=x[0];
    value_type ymax=y[0];
    xmin=x[0];
    ymin=y[0];

    // Loop over all particles and compare their coordinates with the min and max variables
    for(int i=1; i<N; i++){
        // Check x coordinate
        if(x[i] > xmax){
            xmax=x[i];
        } else if(x[i]<xmin){
            xmin = x[i];
        }

        // Check y coordinate
        if(y[i] > ymax){
            ymax=y[i];
        } else if(y[i]<ymin){
            ymin = y[i];
        }
    }

    // Compute the extent in the x and y directions
    value_type xext = xmax-xmin;
    value_type yext = ymax-ymin;

    // Compute the largest extent
    ext = max(xext, yext);
}

/*
 * Returns the Morton quadrant that the particle with coordinates x and y lies in
 */
int quadrant(value_type xmin, value_type ymin, value_type ext, value_type x, value_type y){
    // Prepare some variables we will use to distinguish between upper, lower, left and right square.
    value_type halfSize = ext/2;
    value_type centerX = xmin + halfSize;
    value_type centerY = ymin + halfSize;

    // Find the quadrant that this particle lies in
    if(y > centerY){
        if(x < centerX){
            // Top left
            return 0;
        } else{
            // Top right
            return 1;
        }
    } else{
        if(x<centerX){
            // Bottom left
            return 2;
        } else{
            // Bottom right
            return 3;
        }
    }
}

// Provide all particles with a Morton index value
void morton(const int N, const value_type* const x, const value_type* const y, unsigned int* index, int depth){
    // Define the depth (represents the number of subdivisions)
    //int depth = 2;
    value_type xmin, ymin, ext;

    // Compute the coordinates and extent of the root square
    computeExtent(N, x, y, xmin, ymin, ext);
    //cout << "morton2: xmin: "<< xmin << " ymin: " << ymin << " ext: " << ext << endl;

    // Loop over particles to compute their morton index
    {
    for(int i=0; i<N; i++) {
        // Prepare some variables we will use to distinguish between upper, lower, left and right square.
        value_type extTemp = ext;
        value_type xminTemp = xmin;
        value_type yminTemp = ymin;


        // Here indexValue contains the value that is added or not added to the index to distinguish between the upper
        // and lower areas, i.e. the decimal representation of the morton bits.
        uint32_t indexValue = pow(2, ((depth - 1) * 2));
        index[i] = 0;

        // Loop over the number of subdivisions
        for (int j = 1; j <= depth; j++) {
            // Compute the value to be added to the morton index based on this level
            int q = quadrant(xminTemp, yminTemp, extTemp, x[i], y[i]);
            index[i] += q * indexValue;

            // Compute the new xmin, ymin and ext
            extTemp = ext / 2;
            xminTemp += (q % 2) * extTemp;
            yminTemp += (q < 2) ? extTemp : 0;

            // Compute the value for indexValue for the next iteration.
            indexValue = indexValue / 4;
        }
    }
    }
}

// Sort the array index and store the sorting permutation in the array keys
void sortit(const int N, uint32_t* index, int* keys){

    // Sort the array keys according to the values in index.
    sort(keys, keys + N, [&](const int& a, const int& b){
        return (index[a] < index[b]);
    }
    );

    // Sort the array index
    sort(index, index + N);
}

// Reorder the arrays of the particle coordinates according to the Morton Indices of the particles
void reorder(const int N, const int* const keys, const value_type* const x, const value_type* const y, const value_type* const mass, value_type *xsorted, value_type *ysorted, value_type* masssorted){
    // Loop over the array and put the coordinates of the particle that was ranked at that position in the array
    for (int i = 0; i < N; ++i) {
        xsorted[i] = x[keys[i]];
        ysorted[i] = y[keys[i]];
        masssorted[i] = mass[keys[i]];
    }
}