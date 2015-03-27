/** @file
* @Author: Andreas Dahl
* @Author: Christian Muf
* @Date:   2015-03-10 14:08:00
*/

#ifndef GREEDY_H
#define GREEDY_H

#include <stdio.h>
#include "FastaIO.h"

struct GreedySettings {
    bool greedyPick = true;
    bool lru = true;
    unsigned int cacheSize = 32;
    float similarity;

    GreedySettings(float similarity) {
        this->similarity = similarity;
    }
};

void greedyClustering(FastaIO& dataIO, float (*dist)(FastaContainer &, FastaContainer &), GreedySettings settings, std::ostream* out);

#endif