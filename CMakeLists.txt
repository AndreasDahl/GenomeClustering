cmake_minimum_required(VERSION 2.8.4)
project(GenomeClustering)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Werror -O4")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/out")

set(SOURCE_FILES
    src/FastaIO.cpp
    src/FastaIO.h
    src/KMeans.cpp
    src/KMeans.h
    src/Main.cpp
    src/GreedyClustering.h
    src/GreedyClustering.cpp
    src/MufDifference.h
    src/MufDifference.cpp
    src/PrintUtils.h
    src/PrintUtils.cpp
    src/FastaComparators.h)

add_executable(GenomeClustering ${SOURCE_FILES})