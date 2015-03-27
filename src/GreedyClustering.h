/** @file
* @Author: Andreas Dahl
* @Author: Christian Muf
* @Date:   2015-03-10 14:08:00
*/

#ifndef SIMPLE_GREEDY_H
#define SIMPLE_GREEDY_H

#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include <list>
#include "PrintUtils.h"

struct GreedySettings {
    bool greedyPick = true;
    bool lru = true;
    unsigned int cacheSize = 32;
    float similarity;

    GreedySettings(float similarity) {
        this->similarity = similarity;
    }
};

template <typename T>
void greedyClustering(FastaIO& dataIO, float (*dist)(T  &, T  &), GreedySettings settings, std::ostream* out) {
    std::vector<int> indexes; // Used for data analysis
    std::list<T*> centroids;
    int c_count = 0;
    int n = 0;
    while(true) {
        FastaContainer* current = new FastaContainer();
        if(dataIO.getNextLine(*current)) {
            delete current;
            break;
        }
        // Output progress FIXME: Move somewhere else.
        if (++n % 100 == 0) {
            printProgress((float)(dataIO.getReadFileRead()) / dataIO.getReadFileLength());
        }
        float bestDist = std::numeric_limits<float>::infinity();
        typename std::list<T*>::iterator it;
        unsigned int i = 0; // Used for data analysis
        unsigned int index = 0; // Used for data analysis
        for (it = centroids.begin(); it != centroids.end(); ++it) {
            ++i;
            float distance = dist(*current, **it);
            if (distance < settings.similarity && distance < bestDist) {
                index = i;
                bestDist = distance;
                if (settings.greedyPick)
                    break;
            }
        }
        if (bestDist > settings.similarity) {
            // Only keep "N" centroids.
            if (centroids.size() >= settings.cacheSize) {
                FastaContainer* temp = centroids.back();
                centroids.pop_back();
                delete temp;
            }
            centroids.push_front(current);
            c_count++;
            indexes.push_back(c_count); // TODO: Work for LRU ?
            if(out) *out << c_count << ' ' << 0 << std::endl;
        } else {
            if (settings.lru) {
                // Move hit to front of cache
                FastaContainer* temp = *it;
                centroids.erase(it);
                centroids.push_front(temp);
            }
            if(out) *out << indexes[c_count - index] << ' ' << index << std::endl;
        }
    }
    std::cout << std::endl << "\r" << "Seq Count:" << n << std::endl;
    std::cout << "\r" << "Cluster Count:" << c_count << std::endl;
}

#endif