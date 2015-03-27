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
    std::vector<int> indexes; // Used for data analysis
    std::list<T*> centroids;
    std::cout << data.size() << std::endl;
    int c_count = 0;
    int n = 0;
    for (typename std::vector<T>::iterator current = data.begin(); current != data.end(); ++current) {
        // Output progress FIXME: Move somewhere else.
        if (++n % 100 == 0) {
            printProgress((float)(n) / data.size());
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
                centroids.pop_back();
            }
            centroids.push_front(&(*current));
            c_count++;
            indexes.push_back(c_count); // TODO: Work for LRU ?
            out << c_count << ' ' << 0 << std::endl;
        } else {
            if (settings.lru) {
                // Move hit to front of cache
                centroids.erase(it);
                centroids.push_front(*it);
            }
            out << indexes[c_count - index] << ' ' << index << std::endl;
        }
    }
    std::cout << "\r" << "Cluster Count:" << c_count << std::endl;
}

#endif