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

#include <iostream>
#include <iomanip>
#include <string>
#include "RecordWriter.h"

const char DELIMITER = '\t';

std::ostream& operator<<(std::ostream& out, const Record& record) {
    switch (record.type) {
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
    if(record.query) out << *record.query; else out << '*';
    out << DELIMITER;
    if(record.target) out << *record.target; else out << '*';

    return out;
}


