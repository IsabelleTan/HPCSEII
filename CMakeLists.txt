cmake_minimum_required(VERSION 3.5)
project(Serial)

set(CXX g++-5)                                                  # Use g++ version 5 for openMP
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")   # Use c++11 and openMP

set(SOURCE_FILES
        quadtree.cpp
        morton.cpp
        main.cpp
        test.cpp
        timeMorton.cpp
        expansion.cpp
        )

set(HEADER_FILES
        morton.h
        quadtree.h
        test.h
        timeMorton.h
        timer.hpp
        expansion.h
        )

add_executable(Serial ${SOURCE_FILES} ${HEADER_FILES})