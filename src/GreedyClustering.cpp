/*
 * @Author: Andreas Dahl 
 * @Date: 27/03/15
 */

#include "GreedyClustering.h"

#include <vector>
#include <iostream>
#include "PrintUtils.h"


GreedyClustering::GreedyClustering(float similarity) :
    m_greedyPick(true),
    m_lru(true),
    m_cacheSize(32),
    m_bigCentCache(32)
{
    m_similarity = similarity;
}

bool GreedyClustering::hasHit(float bestDist) {
    return bestDist < m_similarity;
}

bool GreedyClustering::insertIfBig(Centroid* centroid) {
    // See if fits in bigCache
    bool inserted = false;
    for (std::list<Centroid*>::iterator it = m_bigCents.begin(); it != m_bigCents.end(); ++it) {
        if (centroid->count >= (*it)->count) {
            m_bigCents.insert(it, centroid);
            inserted = true;
            if (m_bigCents.size() > m_bigCentCache) {
                Centroid* temp = m_bigCents.back();
                m_bigCents.pop_back();
                delete temp->fasta;
                delete temp;
            }
            break;
        }
    }
    if (!inserted && m_bigCents.size() < m_bigCentCache) {
        m_bigCents.push_back(centroid);
        inserted = true;
    }
    return inserted;
}

void GreedyClustering::greedyClustering(FastaIO& dataIO, float (*dist)(FastaContainer &, FastaContainer &), std::ostream* out) {
    std::vector<int> indexes; // Used for data analysis
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
        for (it = m_cache.begin(); it != m_cache.end(); ++it) {
            ++i;
            float distance = dist(*current, *(*it)->fasta);
            if (hasHit(distance) && distance < bestDist) {
                index = i;
                bestDist = distance;
                if (m_greedyPick)
                    break;
            }
        }
        // Did not hit in 'Centroids'. Search though bigCents. (Of if running thorough search)
        if (!hasHit(bestDist) || !m_greedyPick) {
            for (it = m_bigCents.begin(); it != m_bigCents.end(); ++it) {
                ++i;
                float distance = dist(*current, *(*it)->fasta);
                if (hasHit(distance) && distance < bestDist) {
                    index = i;
                    bestDist = distance;
                    hitBig = true;
                    if (m_greedyPick)
                        break;
                }
            }
        }

        // Current didn't match anything and is a new cluster
        if (!hasHit(bestDist)) {
            // Only keep "N" centroids.
            if (m_cache.size() >= m_cacheSize) {
                Centroid* temp = m_cache.back();
                m_cache.pop_back();

                bool inserted = insertIfBig(temp);
                if (!inserted) {
                    delete temp->fasta;
                    delete temp;
                }
            }
            Centroid* centroid = new Centroid(current);
            m_cache.push_front(centroid);

            // Data collection
            c_count++;
            indexes.push_back(c_count); // TODO: Work for LRU ?
            if(out) *out << c_count << ' ' << 0 << std::endl;
        } else { // Current hit a cluster
            if (m_lru) {
                Centroid *hit = *it;
                hit->count = hit->count + 1;
                // If hit in 'Clusters'
                if (!hitBig) {
                    // Move hit to front of cache
                    m_cache.erase(it);
                    m_cache.push_front(hit);
                } else {
                    m_bigCents.erase(it);
                    m_cache.push_front(hit);

                    if (m_cache.size() > m_cacheSize) {
                        Centroid* temp = m_cache.back();
                        m_cache.pop_back();

                        bool inserted = insertIfBig(temp);

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