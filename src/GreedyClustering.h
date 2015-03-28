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

    Centroid(FastaContainer* fasta) {
        this->fasta=fasta;
    }
};

class GreedyClustering {
    public:
        GreedyClustering(float similarity);

        void greedyClustering(FastaIO& dataIO, float (*dist)(FastaContainer &, FastaContainer &), std::ostream* out);

    private:
        bool m_greedyPick = true;
        bool m_lru = true;
        unsigned int m_cacheSize = 32;
        unsigned int m_longTermCacheSize = 32;
        float m_similarity;

        std::list<Centroid *> m_cache;
        std::list<Centroid *> m_longTermCache;

        /**
         * Return whether the given distance qualifies as a hit.
         * @Param  distance  Distance to test.
         * @Returns  true, if the given distance qualifies as a hit, false otherwise.
         */
        bool isHit(float distance);

        /**
         * Insert centroid into long tern cache if it fits.
         * @Param  centroid  Centroid to be inserted into long term cache.
         * @Returns  true, if centroid was inserted, false otherwise.
         */
        bool insertIntoLongTerm(Centroid *centroid); // TODO: Better name
};


#endif