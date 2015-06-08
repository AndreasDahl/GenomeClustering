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
 * @Date: 23/03/15.
 */

#ifndef _GENOMECLUSTERING_PRINTUTILS_H_
#define _GENOMECLUSTERING_PRINTUTILS_H_

#include <string>

typedef unsigned long long timestamp_t;

timestamp_t get_timestamp();
void printProgress(unsigned long current, unsigned long max);
std::string formatDuration(timestamp_t t0, timestamp_t t1);

#endif  // _GENOMECLUSTERING_PRINTUTILS_H_
