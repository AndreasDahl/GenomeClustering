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
    std::list<T> centroids;
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
        for (typename std::list<T>::iterator it = centroids.begin(); it != centroids.end(); ++it) {
            float distance = dist(*current, *it);
            if (distance < similarity) {
                bestDist = distance;
                break;
            }
        }
        if (bestDist > similarity) {
            // Only keep "N" centroids.
            if (centroids.size() >= 32) {
                centroids.pop_back();
            }
            centroids.push_front(*current);
            c_count++;
//            std::cout
//                    << "S\t"
//                    << centroids.size() - 1 << "\t"
//                    << "SIZE\t"
//                    << "*\t"
//                    << "*\t"
//                    << "*\t"
//                    << "*\t"
//                    << "*\t"
//                    << "Name_" << n << "\t"
//                    << "*";
        } else {
//            std::cout
//                    << "H\t"
//                    << bestIndex << "\t"
//                    << "SIZE\t"
//                    << bestDist << "\t" // TODO: Percentage value
//                    << "+\t"            // TODO: Determine Strand
//                    << "*\t"            // TODO: 0-based coordinate of alignment start in the query sequence.
//                    << "*\t"            // TODO: 0-based coordinate of alignment start in target (seed) sequence.
//                                        //       If minus strand, Tlo is relative to start of reverse-complemented target.
//                    << "*\t"            // TODO: Alignment?
//                    << "Name_" << n << "\t"
//                    << "Name_" << bestIndex;
        }
//        std::cout << std::endl;
//        std::cout.flush();
    }
    std::cout << "Cluster Count:" << c_count << std::endl << std::endl;
}

#endif