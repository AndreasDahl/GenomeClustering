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

#include "KMerString.h"

template <typename T>
void simpleGreedyClustering(std::vector<T>& data, float (*dist)(T const &, T const &), float similarity, std::vector<std::vector<T>>& res) {
    std::vector<T> centroids;
    int n = 0;
    for (typename std::vector<T>::iterator current = data.begin(); current != data.end(); ++current) {
        if (++n % 100 == 0) {
            std::cout << n << std::endl;
            std::cout.flush();
        }
        std::cout.flush();
        float bestDist = std::numeric_limits<float>::infinity();
        int bestIndex;
        for (int i = 0; i < centroids.size(); ++i) {
            float distance = dist(*current, centroids[i]);
            if (distance < bestDist) {
                bestDist = distance;
                bestIndex = i;
            }
        }
        if (bestDist > similarity) {
            centroids.push_back(*current);
            res.push_back({*current});
        } else {
            res[bestIndex].push_back(*current);
        }
    }
}

#endif