/*
 * Copyright 2015 Andreas Dahl, Christian Muf
 *
 * This file is part of MufDahlClust.
 *
 * MufDahlClust is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MufDahlClust is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MufDahlCLust.  If not, see <http://www.gnu.org/licenses/>.
 */

/** @file
* @Author: Andreas Dahl
* @Author: Christian Muf
* @Date:   2015-03-10 14:08:00
*/

#ifndef GREEDY_H
#define GREEDY_H

#include <stdio.h>
#include <list>
#include "FastaIO.h"

struct Centroid {
    FastaContainer* fasta;
    unsigned int count = 1;
    unsigned int clusterNumber;

    Centroid(FastaContainer* fasta, unsigned int clusterNumber) {
        this->fasta=fasta;
        this->clusterNumber = clusterNumber;
    }
};

class GreedyClustering {
    public:
        explicit GreedyClustering(float similarity);

        void start(FastaIO &dataIO, float (*dist)(FastaContainer &, FastaContainer &, float), std::ostream *out);

        void setSimilarity(float similarity);
        void setCacheSize(unsigned int newCacheSize);
        void setLRUSize(unsigned int newSize);
        void setLFUSize(unsigned int newSize);
        void setGreedy(bool greedy);
        void setUsingLRU(bool shouldUseLRU);        

        unsigned int getLRUSize();
        unsigned int getLFUSize();
        float getSimilarity();
        unsigned int getCacheSize();
        unsigned int getClusterCount();
        bool isGreedy();
        bool isUsingLRU();

    private:
        bool m_greedyPick = true;
        bool m_usingLRU = true;
        unsigned int m_cacheSize = 32;
        unsigned int m_longTermCacheSize = 32;
        float m_similarity;

        unsigned int m_clusterCount = 0;

        std::list<Centroid *> m_cache;
        std::list<Centroid *> m_longTermCache;

        /**
         * Return whether the given distance qualifies as a hit.
         *
         * @Param  distance  Distance to test.
         * @Returns  true, if the given distance qualifies as a hit, false otherwise.
         */
        bool isHit(float distance);

        /**
         * Push centroid to cache. If cache is filled it will try to insert the
         * centroid into the long term cache. If not possible the centroid
         * will be deleted.
         *
         * @Param  centroid  Centroid to push to cache.
         */
        void pushToCache(Centroid *centroid);

        /**
         * <p>Insert centroid into long tern cache if it fits.</p>
         *
         * <p>If the given centroid replaces a value into a filled cache, the least
         * desirable centroid will be deleted.</p>
         *
         * @Param  centroid  Centroid to be inserted into long term cache.
         * @Returns  true, if centroid was inserted, false otherwise.
         */
        bool tryInsertIntoLongTerm(Centroid *centroid); // TODO: Better name


};


#endif
