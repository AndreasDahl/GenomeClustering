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

template <typename T>
void simpleGreedyClustering(std::vector<T>& data, float (*dist)(T  &, T  &), float similarity) {
    std::list<T*> centroids;
    std::cout << data.size() << std::endl;
    int c_count = 0;
    int n = 0;
    for (typename std::vector<T>::iterator current = data.begin(); current != data.end(); ++current) {
        // Output progress FIXME: Move somewhere else.
        if (++n % 10000 == 0) {
            std::cout << n << std::endl;
            std::cout.flush();
        }
        float bestDist = std::numeric_limits<float>::infinity();
        T* bestCent;
        for (typename std::list<T*>::iterator it = centroids.begin(); it != centroids.end(); ++it) {
            float distance = dist(*current, **it);
            if (distance < similarity) {
                bestDist = distance;
                bestCent = *it;
                break;
            }
        }
        if (bestDist > similarity) {
            // Only keep "N" centroids.
            if (centroids.size() >= 32) {
                centroids.pop_back();
            }
            centroids.push_front(&(*current));
            c_count++;
        } else {
            // Move hit to front of cache
            centroids.remove(bestCent);
            centroids.push_front(bestCent);
        }
    }
    std::cout << "Cluster Count:" << c_count << std::endl << std::endl;
}

#endif