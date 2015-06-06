/** @file
 * @Author: Christian Muf
 * @Author: Andreas Dahl
 * @Date:   2015-05-24
 */

#include "KMerHash.h"

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

KMerHashmap::KMerHashmap()
{
    m_hashmapSize = 0;
    m_hashmapMask = 0;
    m_hashmapShift = 0;
    m_hashmap = NULL;

    m_iteratorIndex = 0;
    m_iteratorNode = NULL;
}

KMerHashmap::~KMerHashmap()
{
    if(m_hashmap) {
        KMerHashmapNode* current;
        for(int i = 0; i < m_hashmapSize; i++) {
            current = m_hashmap[i];
            while(current != NULL) {
                KMerHashmapNode* temp = current;
                current = current->next;
                delete temp;
            }
        }
        delete[] m_hashmap;
    }
}

bool KMerHashmap::isCreated() const
{
    return m_hashmap != NULL;
}

void KMerHashmap::createHashMap(unsigned int numElements, unsigned int kMerBits)
{
    m_hashmapSize = (int)bitCap(numElements / 2);
    m_hashmapMask = m_hashmapSize - 1;

    m_hashmapShift = 0;
    while((1 << (kMerBits - m_hashmapShift)) > m_hashmapSize)
        m_hashmapShift += 1;

    if(m_hashmap) delete[] m_hashmap;
    m_hashmap = new KMerHashmapNode*[m_hashmapSize];
    for(int i = 0; i < m_hashmapSize; i++)
        m_hashmap[i] = NULL;

    m_iteratorIndex = 0;
    m_iteratorNode = NULL;
}

void KMerHashmap::mKerAdd(KMer value)
{
    unsigned int index = (value.value >> m_hashmapShift) & m_hashmapMask;
    KMerHashmapNode* newNode = new KMerHashmapNode(value);
    KMerHashmapNode* currentNode = m_hashmap[index];
    KMerHashmapNode* prevNode = NULL;

    while(currentNode && newNode->element.value > currentNode->element.value) {
        prevNode = currentNode;
        currentNode = currentNode->next;
    }
    while(currentNode && newNode->element.value == currentNode->element.value &&
                         newNode->element.index > currentNode->element.index) {
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

bool KMerHashmap::iteratorReset()
{
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

void KMerHashmap::iteratorIncrease()
{
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

bool KMerHashmap::iteratorNotEnded()
{
    return m_iteratorIndex < m_hashmapSize;
}

KMer* KMerHashmap::iteratorGet()
{
    return (m_iteratorIndex < m_hashmapSize) ? &m_iteratorNode->element : NULL;
}