cmake_minimum_required(VERSION 3.3)
project(PolarCode)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

SET(BASEPATH "${CMAKE_SOURCE_DIR}")
INCLUDE_DIRECTORIES("${BASEPATH}")

set(SOURCE_FILES main.cpp PolarCode.cpp PolarCode.h encode.cpp decode.cpp word_error_rate.cpp util.cpp
        build_constraight_matrix.cpp GaloisFieldPolynomial.cpp GaloisFieldPolynomial.h
        )
add_executable(PolarCode ${SOURCE_FILES})