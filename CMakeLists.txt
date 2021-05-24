cmake_minimum_required(VERSION 3.3)
project(PolarCode)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

SET(BASEPATH "${CMAKE_SOURCE_DIR}")
INCLUDE_DIRECTORIES("${BASEPATH}")

set(SOURCE_FILES simulation/main.cpp simulation/PolarCode.cpp simulation/PolarCode.h simulation/encode.cpp simulation/decode.cpp simulation/word_error_rate.cpp simulation/util.cpp
        simulation/build_constraight_matrix.cpp simulation/GaloisFieldPolynomial.cpp simulation/GaloisFieldPolynomial.h
        )
add_executable(PolarCode ${SOURCE_FILES})