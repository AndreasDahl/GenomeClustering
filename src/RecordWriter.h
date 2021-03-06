/*
 * Copyright 2015 Andreas Dahl, Christian Muf
 *
 * This file is part of MufDahlClust.
 *
 * MufDahlClust is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MufDahlClust is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MufDahlCLust.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * @Author: Andreas Dahl
 * @Author: Christian Muf
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
};

std::ostream& operator<< (std::ostream& stream, const Record& record);

#endif

