/*
 * @Author: Andreas Dahl
 * @Date: 22/05/15
 */

#include <iostream>
#include <iomanip>
#include <string>
#include "RecordWriter.h"

const char DELIMITER = '\t';

std::ostream& operator<<(std::ostream& out, const Record& record) {
    switch(record.type) {
        case HIT : out << "H"; break;
        case CENTROID : out << "C"; break;
    }
    out << DELIMITER;
    out << record.clusterNumber;
    out << DELIMITER;
    out << record.sequenceLength;
    out << DELIMITER;
    out.precision(3);
    record.type == HIT 
        ? out << record.id 
        : out << '*';
    out << DELIMITER;
    out << record.query;
    out << DELIMITER;
    out << record.target;

    return out;
}


