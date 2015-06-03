/** @file
* @Author: Christian Muf
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