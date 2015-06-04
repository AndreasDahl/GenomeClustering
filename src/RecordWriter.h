/*
 * @Author: Andreas Dahl
 * @Date: 22/05/15
 */

#ifndef RECORD_WRITER_H
#define RECORD_WRITER_H

#include <string>
#include <iostream>

enum RecordType {
    HIT,
    CENTROID
};

/**
 * Struct containing the resulting data from clustering a single sequence
 */
struct Record {
    RecordType type;    // Record Type.
    unsigned int clusterNumber;     // Cluster number.
    unsigned int sequenceLength;    // Length of the sequence.
    float id;   // Identity to the (as a percent, or * if this is a centroid).
    std::string* query = NULL;  // Fasta label of query sequence.
    std::string* target = NULL; // Fasta label of target centroid.
    /*Record() {
        query = target = NULL;
    }*/
};

std::ostream& operator<< (std::ostream& stream, const Record& record);

#endif

