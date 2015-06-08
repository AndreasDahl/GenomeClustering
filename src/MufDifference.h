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

/** @file
 * @Author: Christian Muf
 * @Author: Andreas Dahl
 * @Date:   2015-03-17 12:58:09
 */

#ifndef MUF_DIFFERENCE
#define MUF_DIFFERENCE

#include "FastaIO.h"

void initializeMufDifference();
void uninitializeMufDifference();

float mufDifference(FastaContainer& str1, FastaContainer& str2, float threshold);

float distanceLevenshteinFailFast(FastaContainer& kMer1, FastaContainer& kMer2, float threshold);

#endif