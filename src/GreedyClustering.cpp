/*
 * @Author: Andreas Dahl 
 * @Date: 27/03/15
 */

#include "GreedyClustering.h"

#include <vector>
#include <iostream>
#include <list>
#include "PrintUtils.h"

struct Centroid {
    FastaContainer* fasta;
    unsigned int count = 1;

    Centroid(FastaContainer* fasta) {
        this->fasta=fasta;
    }
};


static bool hasHit(float bestDist, float requirement) {
    return bestDist < requirement;
}

static bool insertIfBig(std::list<Centroid *> &bigCents, Centroid* centroid, unsigned int capacity) {
    // See if fits in bigCache
    bool inserted = false;
    for (std::list<Centroid*>::iterator it = bigCents.begin(); it != bigCents.end(); ++it) {
        if (centroid->count >= (*it)->count) {
            bigCents.insert(it, centroid);
            inserted = true;
            if (bigCents.size() > capacity) {
                Centroid* temp = bigCents.back();
                bigCents.pop_back();
                delete temp->fasta;
                delete temp;
            }
            break;
        }
    }
    if (!inserted && bigCents.size() < capacity) {
        bigCents.push_back(centroid);
        inserted = true;
    }
    return inserted;
}

void greedyClustering(FastaIO& dataIO, float (*dist)(FastaContainer &, FastaContainer &), GreedySettings settings, std::ostream* out) {
    std::vector<int> indexes; // Used for data analysis
    std::list<Centroid *> bigCents;
    std::list<Centroid *> centroids;
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
        bool hitBig = false;
        float bestDist = std::numeric_limits<float>::infinity();
        typename std::list<Centroid *>::iterator it;
        unsigned int i = 0; // Used for data analysis
        unsigned int index = 0; // Used for data analysis
        // Search through 'Centroids' list
        for (it = centroids.begin(); it != centroids.end(); ++it) {
            ++i;
            float distance = dist(*current, *(*it)->fasta);
            if (hasHit(distance, settings.similarity) && distance < bestDist) {
                index = i;
                bestDist = distance;
                if (settings.greedyPick)
                    break;
            }
        }
        // Did not hit in 'Centroids'. Search though bigCents. (Of if running thorough search)
        if (!hasHit(bestDist, settings.similarity) || !settings.greedyPick) {
            for (it = bigCents.begin(); it != bigCents.end(); ++it) {
                ++i;
                float distance = dist(*current, *(*it)->fasta);
                if (hasHit(distance, settings.similarity) && distance < bestDist) {
                    index = i;
                    bestDist = distance;
                    hitBig = true;
                    if (settings.greedyPick)
                        break;
                }
            }
        }

        // Current didn't match anything and is a new cluster
        if (!hasHit(bestDist, settings.similarity)) {
            // Only keep "N" centroids.
            if (centroids.size() >= settings.cacheSize) {
                Centroid* temp = centroids.back();
                centroids.pop_back();

                bool inserted = insertIfBig(bigCents, temp, settings.bigCentCache);
                if (!inserted) {
                    delete temp->fasta;
                    delete temp;
                }
            }
            Centroid* centroid = new Centroid(current);
            centroids.push_front(centroid);

            // Data collection
            c_count++;
            indexes.push_back(c_count); // TODO: Work for LRU ?
            if(out) *out << c_count << ' ' << 0 << std::endl;
        } else { // Current hit a cluster
            if (settings.lru) {
                Centroid *hit = *it;
                hit->count = hit->count + 1;
                // If hit in 'Clusters'
                if (!hitBig) {
                    // Move hit to front of cache
                    centroids.erase(it);
                    centroids.push_front(hit);
                } else {
                    bigCents.erase(it);
                    centroids.push_front(hit);

                    if (centroids.size() > settings.cacheSize) {
                        Centroid* temp = centroids.back();
                        centroids.pop_back();

                        bool inserted = insertIfBig(bigCents, temp, settings.bigCentCache);

                        if (!inserted) {
                            delete temp->fasta;
                            delete temp;
                        }
                    }
                }
            }
            // Data collection
            if(out) *out << indexes[c_count - index] << ' ' << index << std::endl;
        }
    }
    std::cout << std::endl << "\r" << "Seq Count:" << n << std::endl;
    std::cout << "\r" << "Cluster Count:" << c_count << std::endl;
}