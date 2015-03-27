/*
 * @Author: Andreas Dahl 
 * @Date: 27/03/15
 */

#include "GreedyClustering.h"

#include <vector>
#include <iostream>
#include <list>
#include "PrintUtils.h"

void greedyClustering(FastaIO& dataIO, float (*dist)(FastaContainer &, FastaContainer &), GreedySettings settings, std::ostream* out) {
    std::vector<int> indexes; // Used for data analysis
    std::list<FastaContainer*> centroids;
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
        typename std::list<FastaContainer*>::iterator it;
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