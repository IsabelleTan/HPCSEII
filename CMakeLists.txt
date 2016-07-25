cmake_minimum_required(VERSION 3.5)
project(Serial)

set(CXX g++-5)                                                  # Use g++ version 5 for openMP
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -mavx -mfma -O3")    # Most basic compile options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -ftree-vectorize -fopt-info-vec -O3 -mavx -mfma")   # Use c++11 and openMP and vectorize
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -ftree-vectorize -fopt-info-vec -O -mavx -mfma") # Use AVX intrinsics and FMA

set(SOURCE_FILES
        quadtree.cpp
        morton.cpp
        main.cpp
        test.cpp
        timeMorton.cpp
        expansion.cpp
        complexAVX.cpp)

set(HEADER_FILES
        morton.h
        quadtree.h
        test.h
        timeMorton.h
        timer.hpp
        expansion.h
        complexAVX.h)

add_executable(Serial ${SOURCE_FILES} ${HEADER_FILES})