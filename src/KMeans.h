/** @file
* @Author: Andreas Dahl
* @Date:   2015-03-03 01:02:01
*/

#ifndef K_MEANS_H
#define K_MEANS_H

#include <vector>
using std::vector;
#include "KMerString.h"

void kmeans(const vector<KMerString>& data, int k, vector<KMerString>* res);


#endif