cmake_minimum_required(VERSION 2.8.4)
project(GenomeClustering)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    src/FastaIO.cpp
    src/FastaIO.h
    src/KMerString.cpp
    src/KMerString.h
    src/KMeans.cpp
    src/KMeans.h
    src/Main.cpp)

add_executable(GenomeClustering ${SOURCE_FILES})