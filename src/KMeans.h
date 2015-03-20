/** @file
* @Author: Andreas Dahl
* @Date:   2015-03-03 01:02:01
*/

#ifndef K_MEANS_H
#define K_MEANS_H

#include <vector>

#include "FastaIO.h"

void kmeans(std::vector<FastaContainer>& data, int k, std::vector<std::vector<FastaContainer>>& res);


#endif