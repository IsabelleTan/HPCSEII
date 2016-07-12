cmake_minimum_required(VERSION 3.5)
project(Serial)

set(CXX g++-5)                                                  # Use g++ version 5 for openMP
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -ftree-vectorize -fopt-info-vec -O3")   # Use c++11 and openMP
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -ftree-vectorize -fopt-info-vec -O3 -mavx")

set(SOURCE_FILES
        quadtree.cpp
        morton.cpp
        main.cpp
        test.cpp
        timeMorton.cpp
        expansion.cpp
        complex.cpp)

set(HEADER_FILES
        morton.h
        quadtree.h
        test.h
        timeMorton.h
        timer.hpp
        expansion.h
        complex.h)

add_executable(Serial ${SOURCE_FILES} ${HEADER_FILES})