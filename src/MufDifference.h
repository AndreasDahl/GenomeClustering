/** @file
* @Author: Christian Muf
* @Date:   2015-03-17 12:58:09
*/

#ifndef MUF_DIFFERENCE
#define MUF_DIFFERENCE

#include "FastaIO.h"

int createCountedKmer(FastaContainer& str, unsigned int k);
//int createCountedKmer(FastaContainer& str1, FastaContainer& str2, unsigned int k);

float mufDifference(FastaContainer& str1, FastaContainer& str2);
float mufDifferenceFailFast(FastaContainer& str1, FastaContainer& str2);

float kMerDistanceLevenshtein(FastaContainer& kMer1, FastaContainer& kMer2);

#endif