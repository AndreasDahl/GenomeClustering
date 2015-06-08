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
 * @Author: Christian Muf
 * @Author: Andreas Dahl
 * @Date:   2015-05-24
 */

#ifndef KMER_HASH_H
#define KMER_HASH_H

#include <stddef.h>

/**
* 
*/
struct KMer
{
    int value;
    int index;
    KMer(int v=0, int i=0) {
        value = v;
        index = i;
    }
};

/**
* Abstract class every kMer data-structure is based on.
*/
template<class T>
class KMerStructure
{
    public:
        virtual void mKerAdd(T) = 0;

        virtual bool iteratorReset() = 0;
        virtual void iteratorIncrease() = 0;
        virtual bool iteratorNotEnded() = 0;
        virtual T* iteratorGet() = 0;
};

/**
* Hash-map structure used for std kMer comparison.
*/
class KMerHashmap : public KMerStructure<KMer>
{
    struct KMerHashmapNode {
        KMerHashmapNode* next;
        KMer element;
        KMerHashmapNode(KMer& elementIn) {
            element = elementIn;
            next = NULL;
        }
    };

    public:
        KMerHashmap();
        virtual ~KMerHashmap();

        KMerHashmap(const KMerHashmap&) = delete;
        KMerHashmap& operator=(const KMerHashmap&) = delete;

        bool isCreated() const;
        void createHashMap(unsigned int numElements, unsigned int kMerBits);
        void mKerAdd(KMer value);

        bool iteratorReset();
        void iteratorIncrease();
        bool iteratorNotEnded();

        KMer* iteratorGet();

    private:
        int m_hashmapSize;
        int m_hashmapMask;
        int m_hashmapShift;
        KMerHashmapNode** m_hashmap;

        int m_iteratorIndex;
        KMerHashmapNode* m_iteratorNode;
};

#endif // KMER_HASH_H