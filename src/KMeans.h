/** @file
* @Author: Andreas Dahl
* @Date:   2015-03-03 01:02:01
*/

#ifndef K_MEANS_H
#define K_MEANS_H

#include <vector>

#include "KMerString.h"

void kmeans(const std::vector<KMerString>& data, int k, std::vector<std::vector<KMerString>>& res);


#endif