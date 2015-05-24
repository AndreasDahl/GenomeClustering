/** @file
* @Author: Christian Muf
* @Date:   2015-03-17 12:58:09
*/

#ifndef MUF_DIFFERENCE
#define MUF_DIFFERENCE

#include "FastaIO.h"

float mufDifference(FastaContainer& str1, FastaContainer& str2);
float mufDifference(FastaContainer& str1, FastaContainer& str2, float threshold);

//int createCountedKmer(FastaContainer& str, unsigned int k);
//int createCountedKmer(FastaContainer& str1, FastaContainer& str2, unsigned int k);

float distanceLevenshteinFailFast(FastaContainer& kMer1, FastaContainer& kMer2, float threshold);

#endif