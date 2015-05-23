/** @file
* @Author: Christian Muf
* @Date:   2015-03-17 12:58:21
*/

#include "MufDifference.h"

#include <math.h>
#include <list>

#include <iostream>
using namespace std;


/**
* Creates a integer in the form 2^n that is >= input.
* Examples:
* 00000011 -> 00000100
* 00001111 -> 00010000
* 00001001 -> 00010000
* 00001000 -> 00001000
*/
inline unsigned int bitCap(unsigned int input)
{
    unsigned int store = input;
    unsigned int i = 0;
    while(input != 0) {
        input >>= 1;
        i += 1;
    }
    input = 1 << i;
    return ((input >> 1) == store) ? store : input;
}

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
* Simple array structure used for parallel comparison.
* Template class T must have operators: =(int) and +=(T)
*/
template<class T>
class KMerCountedArray : public KMerStructure<T>
{
    public:
        KMerCountedArray(unsigned int arraySize) :
            m_arraySize(arraySize),
            m_iterator(0)
        {
            m_array = new T[arraySize];
            for(int i = 0; i < arraySize; i++)
                m_array[i] = 0;
        }

        virtual ~KMerCountedArray() {
            delete m_array;
        }

        KMerCountedArray(const KMerCountedArray&) = delete;
        KMerCountedArray& operator=(const KMerCountedArray&) = delete;

        void mKerAdd(T index) { // unsigned int index, 
            m_array[index] += 1;
        }

        bool iteratorReset() {
            m_iterator = 0;
            return true;
        }

        void iteratorIncrease() {
            m_iterator += 1;
        }

        bool iteratorNotEnded() {
            return m_iterator < m_arraySize;
        }

        T& iteratorGet() {
            return m_array[m_iterator];
        }

    private:
        unsigned int m_arraySize;
        T* m_array;

        unsigned int m_iterator;
};

/**
* Hash-map structure used for std kMer comparison.
*/
template<class T>
class KMerHashmap : public KMerStructure<T>
{
    struct KMerHashmapNode {
        KMerHashmapNode* next;
        T element;
        KMerHashmapNode(T& elementIn) {
            element = elementIn;
            next = NULL;
        }
    };

    public:
        KMerHashmap(unsigned int numElements, unsigned int kMerBits) :
            m_hashmapSize(bitCap(numElements / 2))
        {
            m_hashmapMask = m_hashmapSize - 1;

            m_hashmapShift = 0;
            while((1 << (kMerBits - m_hashmapShift)) > m_hashmapSize)
                m_hashmapShift += 1;

            m_hashmap = new KMerHashmapNode*[m_hashmapSize];
            for(int i = 0; i < m_hashmapSize; i++)
                m_hashmap[i] = NULL;
        }

        virtual ~KMerHashmap() {
            KMerHashmapNode* current;
            for(int i = 0; i < m_hashmapSize; i++) {
                current = m_hashmap[i];
                while(current != NULL) {
                    KMerHashmapNode* temp = current;
                    current = current->next;
                    delete temp;
                }
            }
            delete m_hashmap;
        }

        KMerHashmap(const KMerHashmap&) = delete;
        KMerHashmap& operator=(const KMerHashmap&) = delete;

        void mKerAdd(T value) {
            unsigned int index = (value >> m_hashmapShift) & m_hashmapMask;
            KMerHashmapNode* newNode = new KMerHashmapNode(value);
            KMerHashmapNode* currentNode = m_hashmap[index];
            KMerHashmapNode* prevNode = NULL;

            while(currentNode && newNode->element > currentNode->element) {
                prevNode = currentNode;
                currentNode = currentNode->next;
            }
            if(prevNode == NULL) {
                m_hashmap[index] = newNode;
                newNode->next = currentNode;
            }
            else {
                prevNode->next = newNode;
                newNode->next = currentNode;
            }
        }

        bool iteratorReset() {
            m_iteratorIndex = 0;
            m_iteratorNode = NULL;
            while(m_iteratorIndex < m_hashmapSize) {
                m_iteratorNode = m_hashmap[m_iteratorIndex];
                if(m_iteratorNode != NULL)
                    return true;
                m_iteratorIndex += 1;
            }
            return false;
        }

        void iteratorIncrease() {
            if(m_iteratorNode != NULL && m_iteratorNode->next != NULL) {
                m_iteratorNode = m_iteratorNode->next;
                return;
            }
            while(++m_iteratorIndex < m_hashmapSize) {
                m_iteratorNode = m_hashmap[m_iteratorIndex];
                if(m_iteratorNode != NULL)
                    return;
            }
        }

        bool iteratorNotEnded() {
            return m_iteratorIndex < m_hashmapSize;
        }

        T& iteratorGet() {
            return m_iteratorNode->element;
        }

    private:
        unsigned int m_hashmapSize;
        unsigned int m_hashmapMask;
        unsigned int m_hashmapShift;
        KMerHashmapNode** m_hashmap;

        unsigned int m_iteratorIndex;
        KMerHashmapNode* m_iteratorNode;
};

template<class T>
inline int parallelComparison(KMerCountedArray<T>& structure1, KMerCountedArray<T>& structure2)
{
    structure1.iteratorReset();
    structure2.iteratorReset();

    int difference = 0;
    while(structure1.iteratorNotEnded() && structure2.iteratorNotEnded()) {
        difference += abs((int)structure1.iteratorGet() - (int)structure2.iteratorGet());
        structure1.iteratorIncrease();
        structure2.iteratorIncrease();
    }

    return difference;
}

template<class T>
inline int shiftingComparison(KMerHashmap<T>& structure1, KMerHashmap<T>& structure2)
{
    structure1.iteratorReset();
    structure2.iteratorReset();

    int difference = 0;
    while(structure1.iteratorNotEnded() && structure2.iteratorNotEnded()) {
        T t1 = structure1.iteratorGet();
        T t2 = structure2.iteratorGet();
        if(t1 == t2) {
            structure1.iteratorIncrease();
            structure2.iteratorIncrease();
        }
        else {
            difference += 1;
            if(t1 < t2) {
                structure1.iteratorIncrease();
            }
            else { // t1 > t2
                structure2.iteratorIncrease();
            }
        }
    }

    return difference;
}

template<class T>
inline void generateKmer(const std::string& seq, KMerStructure<T>& data, unsigned int k)
{
    unsigned int seqLen = (unsigned int)seq.size()-k+1;
    for(unsigned int i = 0; i < seqLen; i++) {
        T kMer = 0;
        for(unsigned int j = 0; j < k; j++) {
            kMer <<= 2;
            kMer |= ((T)seq[i+j] & 6) >> 1;
        }
        data.mKerAdd(kMer);
    }
}

class TestType
{
    public:
        unsigned int value;
        unsigned int index;

        TestType& operator=(const TestType& other) {
            if(&other != this) {
                value = other.value;
                index = other.index;
            }
            return *this;
        }

        
};

int kMerTest(FastaContainer& str1, FastaContainer& str2)
{
    typedef int UsedType;

    unsigned int k = 5;
    unsigned int arrayLen = 1 << (k << 1); // 4^k == 2^(k*2)

    KMerCountedArray<int> data1(arrayLen);
    KMerCountedArray<int> data2(arrayLen);

    generateKmer<int>(str1.sequence, data1, k);
    generateKmer<int>(str2.sequence, data2, k);

    // shiftingComparison parallelComparison
    int res = parallelComparison<int>(data1, data2);



    KMerHashmap<UsedType> data3(str1.sequence.size(), k<<1);
    KMerHashmap<UsedType> data4(str2.sequence.size(), k<<1);

    generateKmer<UsedType>(str1.sequence, data3, k);
    generateKmer<UsedType>(str2.sequence, data4, k);

    res = shiftingComparison<UsedType>(data3, data4);

    return res;
}


/*
struct FastaContainer
{
    std::string sequence;
    unsigned long lineNumber;
    int* kMer;
    int k; // Size of k
    FastaContainer() {
        sequence = "";
        lineNumber = 0;
        kMer = NULL;
        k = -1;
    }
};
*/

int createCountedKmer(FastaContainer& str, unsigned int k)
{
    // If k has a garbage value, return -1.
    if(k == 0) {
        return -1;
    }

    unsigned int permutations = 1 << (k << 1); // 4^k == 2^(k*2)

    // Delete the counted kMer if it exists.
    str.setKMerLength(permutations);
    str.k = k;

    // Set all elements an the counted kMer to 0.
    for(unsigned int i = 0; i < permutations; i++) {
        str.getKMer()[i] = 0;
    }

    // If the size of the sequence is less than k, no kMer can exist.
    if((unsigned int)str.sequence.size() < k) {
        return -2;
    }

    // Calculate all kMer and increase the values of the counted kMer.
    for(unsigned int i = 0; i < (unsigned int)str.sequence.size()-(k-1); i++) {
        int kMer = 0;
        for(unsigned int j = 0; j < k; j++) {
            kMer <<= 2;
            kMer |= (str.sequence[i+j] & 6) >> 1;
        }
        str.getKMer()[kMer] += 1;
    }

    // Everything went well, return 0.
    return 0;
}

float mufDifference(FastaContainer& str1, FastaContainer& str2)
{
    unsigned int k = 5; // DO THIS BETTER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    unsigned int permutations = 1 << (k << 1); // 4^k == 2^(k*2)

    int strSize1 = (int)str1.sequence.size();
    int strSize2 = (int)str2.sequence.size();

    int sequenceSizeDifference = std::abs(strSize1 - strSize2);
    int difference = -sequenceSizeDifference;

    if(str1.getKMer() == NULL || str1.k != k) {
        createCountedKmer(str1, k);
    }

    if(str2.getKMer() == NULL || str2.k != k) {
        createCountedKmer(str2, k);
    }

    for(unsigned int i = 0; i < permutations; i++) {
        difference += std::abs(str1.getKMer()[i] - str2.getKMer()[i]);
    }

    // Find the real constant instead of 6.0
    return fabs(((float)difference / 6.0f) / (float)((strSize1 < strSize2) ? strSize1 : strSize2));
}

float mufDifferenceFailFast(FastaContainer&, FastaContainer&)
{
    return 0.0f;
}

float mufDifferenceHash(FastaContainer& str1, FastaContainer& str2)
{
    unsigned int k = 9;
    unsigned int skip = 3;

    unsigned int strSize1 = (unsigned int)str1.sequence.size();
    unsigned int strSize2 = (unsigned int)str2.sequence.size();

    unsigned int hashSize1 = strSize1 / skip / 2 + 1;
    unsigned int hashSize2 = strSize2 / skip / 2 + 1;

    unsigned int shift = 0;
    unsigned int temp = hashSize1;
    while(temp > 0) { shift += 1; temp >>= 1; }
    if(hashSize1 > (1 << (shift-1))) hashSize1 = 1 << shift;

    shift = 0;
    temp = hashSize2;
    while(temp > 0) { shift += 1; temp >>= 1; }
    if(hashSize2 > (1 << (shift-1))) hashSize2 = 1 << shift;

    unsigned int hashMask1 = hashSize1 - 1;
    unsigned int hashMask2 = hashSize2 - 1;

    std::list<unsigned int>* hashMap1 = new std::list<unsigned int>[hashSize1];
    std::list<unsigned int>* hashMap2 = new std::list<unsigned int>[hashSize2];

    for(int i = 0; i < hashSize1; i++) hashMap1[i].clear();
    for(int i = 0; i < hashSize2; i++) hashMap2[i].clear();

    int count1 = 0;
    for(unsigned int i = 0; i < strSize1+k; i += skip) {
        const char* ptr = &str1.sequence[i];
        unsigned long kMer = 0;
        for(int j = 0; j < k; j++) {
            kMer <<= 2;
            kMer |= (ptr[j] & 6) >> 1;
        }
        kMer = (kMer * 122949823) >> (k*2);
        hashMap1[kMer & hashMask1].push_back(kMer);
        count1 += 1;
    }

    int count2 = 0;
    for(unsigned int i = 0; i < strSize2+k; i += skip) {
        const char* ptr = &str2.sequence[i];
        unsigned long kMer = 0;
        for(int j = 0; j < k; j++) {
            kMer <<= 2;
            kMer |= (ptr[j] & 6) >> 1;
        }
        kMer = (kMer * 122949823) >> (k*2); 
        hashMap2[kMer & hashMask2].push_back(kMer);
        count2 += 1;
    }

    int countDiff = std::abs(count1 - count2);
    int distance = 0;
    for(int i = 0; i < countDiff; i++) {

    }




    delete hashMap1;
    delete hashMap2;

    return 0.0f;
}

#define MIN2(a, b) (a) < (b) ? (a) : (b)
#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
int kMerDistanceLevenshtein(FastaContainer& kMer1, FastaContainer& kMer2)
{
    unsigned int s1len, s2len, x, y, lastdiag, olddiag, prevdiag;
    char s2c;

    const std::string& s1 = kMer1.sequence;
    const std::string& s2 = kMer2.sequence;

    s1len = s1.size();
    s2len = s2.size();

    unsigned int* column = new unsigned int[s1len+1];

    s2c = s2[0] | 32;
    for(y = 1, lastdiag = 0; y <= s1len; y++)
    {
        olddiag = y;
        column[y] = MIN2(y, lastdiag + ((s1[y-1] | 32) == s2c ? 0 : 1));
        lastdiag = olddiag;
    }

    for(x = 1; x < s2len; x++)
    {
        prevdiag = column[0] = x+1;

        s2c = s2[x] | 32;

        for(y = 1, lastdiag = x; y <= s1len; y++)
        {
            olddiag = column[y];
            prevdiag = column[y] = MIN3(olddiag+1, prevdiag+1, lastdiag + ((s1[y-1] | 32) == s2c ? 0 : 1));
            lastdiag = olddiag;
        }
    }

    unsigned int result = column[s1len];
    delete column;
    return (result - ((s1len - s2len > 0) ? s1len - s2len : s2len - s1len));
    //return (float)(result - ((s1len - s2len > 0) ? s1len - s2len : s2len - s1len)) /
    //    (float)((s1len < s2len) ? s1len : s2len);
}

