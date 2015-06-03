/*
 * @Author: Andreas Dahl 
 * @Date: 27/03/15
 */

#include "GreedyClustering.h"

#include <vector>
#include <iostream>
#include <limits>
#include "PrintUtils.h"
#include "RecordWriter.h"

GreedyClustering::GreedyClustering(float similarity) :
    m_greedyPick(true),
    m_lru(true),
    m_cacheSize(32),
    m_longTermCacheSize(32)
{
    m_similarity = similarity;
}

void GreedyClustering::start(
        FastaIO &dataIO,
        float (*dist)(FastaContainer &, FastaContainer &, float),
        std::ostream *out) {
    int n = 0;
    while (true) {
        // Load Sequence
        FastaContainer* current = new FastaContainer();
        if (dataIO.getNextLine(*current)) {
            delete current;
            break;
        }
        // Construct record for capturing results
        Record r;
        r.sequenceLength = (current->sequence).length();
        // Output progress.
        if (++n % 500 == 0) {
            printProgress(dataIO.getReadFileRead(), dataIO.getReadFileLength());
        }
        bool hitBig = false;
        float bestDist = std::numeric_limits<float>::infinity();
        typename std::list<Centroid *>::iterator it;
        unsigned int i = 0;  // Used for data analysis
        // (warning removal) unsigned int index = 0;  // Used for data analysis
        // Search through 'Centroids' list
        for (it = m_cache.begin(); it != m_cache.end(); ++it) {
            ++i;
            float distance = dist(*current, *(*it)->fasta, m_similarity);
            if (isHit(distance) && distance < bestDist) {
                // (warning removal) index = i;
                bestDist = distance;
                if (m_greedyPick)
                    break;
            }
        }
        // Did not hit in 'Centroids'. Search though bigCents. (Or if running thorough search)
        if (!isHit(bestDist) || !m_greedyPick) {
            for (it = m_longTermCache.begin(); it != m_longTermCache.end(); ++it) {
                ++i;
                float distance = dist(*current, *(*it)->fasta, m_similarity);
                if (isHit(distance) && distance < bestDist) {
                    // (warning removal) index = i;
                    bestDist = distance;
                    hitBig = true;
                    if (m_greedyPick)
                        break;
                }
            }
        }

        // Current didn't match anything and is a new cluster
        if (!isHit(bestDist)) {
            Centroid *centroid = new Centroid(current, m_clusterCount);
            pushToCache(centroid);

            r.type = CENTROID;
            r.clusterNumber = centroid->clusterNumber;
            if (out) *out << r << std::endl;

            m_clusterCount++;
        } else {  // Current hit a cluster
            Centroid *hit = *it;
            hit->count = hit->count + 1;
            if (m_lru) {
                // If hit in 'Clusters'
                if (!hitBig) {
                    m_cache.erase(it);
                } else {
                    m_longTermCache.erase(it);
                }
                pushToCache(hit);
            }
            delete current;

            r.type = HIT;
            r.clusterNumber = hit->clusterNumber;
            r.id = 100 - (bestDist * 100);
            if (out) *out << r << std::endl;
        }
    }
    std::cout << std::endl << "\r" << "Seq Count:" << n << std::endl;
    std::cout << "\r" << "Cluster Count:" << m_clusterCount << std::endl;
}

void GreedyClustering::setSimilarity(float similarity) {
    m_similarity = similarity;
}

void GreedyClustering::setCacheSize(unsigned int newCacheSize) {
    m_cacheSize = newCacheSize / 2;
    m_longTermCacheSize = newCacheSize - m_cacheSize;
}

void GreedyClustering::setLRUSize(unsigned int newSize) {
    m_cacheSize = newSize;
}

void GreedyClustering::setLFUSize(unsigned int newSize) {
    m_longTermCacheSize = newSize;
}

unsigned int GreedyClustering::getLRUSize() {
    return m_cacheSize;
}

unsigned int GreedyClustering::getLFUSize() {
    return m_longTermCacheSize;
}

float GreedyClustering::getSimilarity() {
    return m_similarity;
}

unsigned int GreedyClustering::getCacheSize() {
    return m_cacheSize + m_longTermCacheSize;
}

unsigned int GreedyClustering::getClusterCount() {
    return m_clusterCount;
}

// Private Members ------------------------------------------------------------

bool GreedyClustering::isHit(float distance) {
    return distance < m_similarity;
}

void GreedyClustering::pushToCache(Centroid *centroid) {
    m_cache.push_front(centroid);

    // Only keep "m_cacheSize" centroids.
    if (m_cache.size() > m_cacheSize) {
        Centroid* temp = m_cache.back();
        m_cache.pop_back();

        bool inserted = tryInsertIntoLongTerm(temp);
        if (!inserted) {
            delete temp->fasta;
            delete temp;
        }
    }
}

bool GreedyClustering::tryInsertIntoLongTerm(Centroid *centroid) {
    // See if fits in bigCache
    bool inserted = false;
    for (std::list<Centroid*>::iterator it = m_longTermCache.begin(); it != m_longTermCache.end(); ++it) {
        if (centroid->count >= (*it)->count) {
            m_longTermCache.insert(it, centroid);
            inserted = true;
            if (m_longTermCache.size() > m_longTermCacheSize) {
                Centroid* temp = m_longTermCache.back();
                m_longTermCache.pop_back();
                delete temp->fasta;
                delete temp;
            }
            break;
        }
    }
    if (!inserted && m_longTermCache.size() < m_longTermCacheSize) {
        m_longTermCache.push_back(centroid);
        inserted = true;
    }
    return inserted;
}
