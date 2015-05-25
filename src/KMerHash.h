/** @file
* @Author: Christian Muf
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
        virtual T& iteratorGet() = 0;
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

        KMer& iteratorGet();

    private:
        int m_hashmapSize;
        int m_hashmapMask;
        int m_hashmapShift;
        KMerHashmapNode** m_hashmap;

        int m_iteratorIndex;
        KMerHashmapNode* m_iteratorNode;
};

#endif // KMER_HASH_H