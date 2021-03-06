cmake_minimum_required(VERSION 2.8.4)
project(GenomeClustering)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/out")

set(SOURCE_FILES
    src/FastaIO.cpp
    src/FastaIO.h
    src/KMerHash.cpp
    src/KMerHash.h
    src/Main.cpp
    src/GreedyClustering.h
    src/GreedyClustering.cpp
    src/MufDifference.h
    src/MufDifference.cpp
    src/PrintUtils.h
    src/PrintUtils.cpp
    src/InputParser.h
    src/InputParser.cpp
    src/FastaComparators.h
    src/RecordWriter.h
    src/RecordWriter.cpp)

add_executable(GenomeClustering ${SOURCE_FILES})

add_definitions(-std=c++11 -O3 -Wall -Wextra)
