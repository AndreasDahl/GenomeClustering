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

struct SimpleGreedySettings {
    bool greedyPick = true;
    bool lru = true;
    unsigned int cacheSize = 32;
    float similarity;

    SimpleGreedySettings(float similarity) {
        this->similarity = similarity;
    }
};

template <typename T>
void simpleGreedyClustering(std::vector<T>& data, float (*dist)(T  &, T  &), SimpleGreedySettings settings, std::ostream& out) {
    std::list<T*> centroids;
    out << data.size() << std::endl;
    int c_count = 0;
    int n = 0;
    for (typename std::vector<T>::iterator current = data.begin(); current != data.end(); ++current) {
        // Output progress FIXME: Move somewhere else.
        if (++n % 1000 == 0) {
            printProgress((float)(n) / data.size());
        }
        float bestDist = std::numeric_limits<float>::infinity();
        typename std::list<T*>::iterator it;
        for (it = centroids.begin(); it != centroids.end(); ++it) {
            float distance = dist(*current, **it);
            if (distance < settings.similarity && distance < bestDist) {
                bestDist = distance;
                if (settings.greedyPick)
                    break;
            }
        }
        if (bestDist > settings.similarity) {
            // Only keep "N" centroids.
            if (centroids.size() >= settings.cacheSize) {
                centroids.pop_back();
            }
            centroids.push_front(&(*current));
            c_count++;
        } else {
            if (settings.lru) {
                // Move hit to front of cache
                centroids.erase(it);
                centroids.push_front(*it);
            }
        }
    }
    out << "\r" << "Cluster Count:" << c_count << std::endl;
}

#endif