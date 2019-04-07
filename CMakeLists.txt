cmake_minimum_required(VERSION 3.10)
project(vaja_III_1)

set(CMAKE_CXX_STANDARD 14)

add_executable(run src/main.cpp src/experiments.hpp src/SimplecticS4.cpp src/SimplecticS4.h)
target_link_libraries(run gflags pthread gsl gslcblas)
