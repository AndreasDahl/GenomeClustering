/** @file
* @Author: Christian Muf
* @Date:   2015-03-17 12:58:21
*/

#include "MufDifference.h"

#include "KMerHash.h"

#include <iostream>
#include <string>
#include <math.h>
#include <vector>
#include <list>

using std::list;
using std::vector;

using namespace std;

struct IndexPair
{
    int i1, i2;
    IndexPair(int I1 = 0, int I2 = 0) {
        i1 = I1; i2 = I2;
    }
};

namespace MufDiffGlobals {
    static int* sortedIndices = NULL;
    static int sortedIndicesSize = 1 << 12; // 4096 seems like a good stating value
}

void initializeMufDifference()
{
    MufDiffGlobals::sortedIndices = new int[MufDiffGlobals::sortedIndicesSize];

    for(int i = 0; i < MufDiffGlobals::sortedIndicesSize; i++) {
        MufDiffGlobals::sortedIndices[i] = -1;
    }
}

void uninitializeMufDifference()
{
    delete[] MufDiffGlobals::sortedIndices;
}

static void reSizeSortedIndices()
{
    delete[] MufDiffGlobals::sortedIndices;

    MufDiffGlobals::sortedIndicesSize <<= 1;
    MufDiffGlobals::sortedIndices = new int[MufDiffGlobals::sortedIndicesSize];

    for(int i = 0; i < MufDiffGlobals::sortedIndicesSize; i++) {
        MufDiffGlobals::sortedIndices[i] = -1;
    }
}

#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
static int levenshteinHelper(const char* str1, int strSize1, const char* str2, int strSize2, int maxErrors)
{
    if(str1 == str2 || strSize1 <= 0 || strSize2 <= 0 || abs(strSize1 - strSize2) >= maxErrors) {
        return abs(strSize1 - strSize2);
    }

    // Perform linear search to find first error
    while(((*str1)|32) == ((*str2)|32)) {
        strSize1 -= 1;
        strSize2 -= 1;
        if(strSize1 <= 0 || strSize2 <= 0) {
            return abs(strSize1 - strSize2);
        }
        str1++;
        str2++;
    }

    int* prevArray = new int[strSize1 + 1];
    int* currArray = new int[strSize1 + 1];

    for(int i = 0; i <= strSize1; i++) {
        prevArray[i] = i;
    }

    int currentMin = 0;
    for(int i = 1; i <= strSize2; i++) {
        currArray[0] = currentMin = i;

        for(int j = 1; j <= strSize1; j++) {
            if((str1[i-1]|32) == (str2[j-1]|32)) {
                currArray[j] = prevArray[j-1];
            } else {
                currArray[j] = MIN3(prevArray[j],   // Deletion
                                    currArray[j-1], // Insertion
                                    prevArray[j-1]) // Substitution
                                    + 1;
            }
            if(currentMin > currArray[j])
                currentMin = currArray[j];
        }

        if(currentMin > maxErrors) { // Fail fast
            delete[] prevArray;
            delete[] currArray;
            return currentMin;
        }

        int* temp = prevArray;
        prevArray = currArray;
        currArray = temp;        
    }

    int result = prevArray[strSize1];

    delete[] prevArray;
    delete[] currArray;

    return result;
}

static void generateHashKmer(const std::string& seq, KMerHashmap& data, unsigned int k)
{
    unsigned int seqLen = (unsigned int)seq.size()-k+1;
    for(unsigned int i = 0; i < seqLen; i++) {
        KMer kMer = KMer(0, i);
        for(unsigned int j = 0; j < k; j++) {
            kMer.value <<= 2;
            kMer.value |= ((int)seq[i+j] & 6) >> 1;
        }
        data.mKerAdd(kMer);
    }
}

static int shiftingComparison(KMerHashmap& bigStruct, KMerHashmap& smallStruct, int* pairList)
{
    bigStruct.iteratorReset();
    smallStruct.iteratorReset();

    int difference = 0;
    while(bigStruct.iteratorNotEnded() && smallStruct.iteratorNotEnded()) {
        KMer* t1 = bigStruct.iteratorGet();
        KMer* t2 = smallStruct.iteratorGet();
        if(t1->value == t2->value) {
            bigStruct.iteratorIncrease();
            smallStruct.iteratorIncrease();
            pairList[t1->index] = t2->index;
        }
        else {
            difference += 1;
            if(t1->value < t2->value) {
                bigStruct.iteratorIncrease();
            }
            else { // t1 > t2
                smallStruct.iteratorIncrease();
            }
        }
    }

    while(bigStruct.iteratorNotEnded()) {
        difference += 1;
        bigStruct.iteratorIncrease();
    }
    while(smallStruct.iteratorNotEnded()) {
        difference += 1;
        smallStruct.iteratorIncrease();
    }

    return difference;
}

// Takes a sorted vector
void removeCrossingIndices(vector<IndexPair>& sortedList)
{
    int lastIndex = -1;
    for(int i = 0; i < (int)sortedList.size(); i++) {
        int currIndex = sortedList[i].i2;
        if(currIndex >= lastIndex) {
            lastIndex = currIndex;
        } else { // Current is less than last, thus crossing
            int c = 1;
            while(true) {
                // Check is deleting current, or the next ones, is a fix 
                if(i+c >= (int)sortedList.size() || sortedList[i+c].i2 >= lastIndex) {
                    sortedList.erase(sortedList.begin()+i, sortedList.begin()+(i+c));
                    i -= 1; // When deleting front set index 1 back
                    break;
                }
                // Check is deleting previous, or ones before that, is a fix 
                else if(i-c-1 < 0 || sortedList[i-c-1].i2 <= currIndex) {
                    sortedList.erase(sortedList.begin()+(i-c), sortedList.begin()+i);
                    lastIndex = (i-c-1 < 0) ? -1 : lastIndex;
                    i -= c+1;
                    break;
                }
                c += 1;
            }
        }
    }
}

// Takes a sorted list
void removeCrossingIndices(list<IndexPair>& sortedList)
{
    int lastIndex = -1;
    list<IndexPair>::iterator it = sortedList.begin();
    while(it != sortedList.end()) {
        int currIndex = it->i2;
        if(currIndex >= lastIndex) {
            lastIndex = currIndex;
            it++;
        } else { // Current is less than last, thus crossing
            list<IndexPair>::iterator back = it; back--;
            list<IndexPair>::iterator front = it;
            if(back == sortedList.begin()) {
                it = sortedList.erase(back);
                lastIndex = -1;
            }
            else while(true) {
                // Check is deleting current, or the next ones, is a fix 
                if(++front == sortedList.end() || front->i2 >= lastIndex) {
                    it = sortedList.erase(it, front);
                    break;
                }
                // Check is deleting previous, or ones before that, is a fix 
                else if(--back == sortedList.begin() || back->i2 <= currIndex) {
                    lastIndex = (back == sortedList.begin()) ? -1 : currIndex;
                    it = sortedList.erase(back, it);
                    break;
                }
            }
        }
    }
}

float mufDifference(FastaContainer& str1, FastaContainer& str2, float threshold)
{
    // Set k for k-mers
    const int k = 8;

    // Return 0% errors if str1 and str2 is the same object.
    if(&str1 == &str2)
        return 0.0f;

    // Find the biggest and smallest of the two inputs.
    FastaContainer *biggest, *smallest;
    int bigSize = (int)str1.sequence.size();
    int smallSize = (int)str2.sequence.size();
    if(bigSize > smallSize) {
        biggest = &str1;
        smallest = &str2;
    } else {
        int temp = bigSize;
        bigSize = smallSize;
        smallSize = temp;
        biggest = &str2;
        smallest = &str1;
    }

    // Build the hashmaps is they don't exists.
    if(!biggest->kMerHash.isCreated()) {
        biggest->kMerHash.createHashMap(bigSize, k<<1);
        generateHashKmer(biggest->sequence, biggest->kMerHash, k);
    }
    if(!smallest->kMerHash.isCreated()) {
        smallest->kMerHash.createHashMap(smallSize, k<<1);
        generateHashKmer(smallest->sequence, smallest->kMerHash, k);
    }

    // Make sure sortedIndices is big enough (used for sorting)
    while(bigSize > MufDiffGlobals::sortedIndicesSize) {
        reSizeSortedIndices();
    }

    // Calculate the allowed number of errors.
    int allowedErrors = (int)ceil((float)smallSize * threshold);

    // Run Manhattan-difference on hashmaps, and save the pair-list
    int manDiff = (shiftingComparison(biggest->kMerHash, smallest->kMerHash, MufDiffGlobals::sortedIndices) -
        (bigSize - smallSize));
    if(manDiff / (2*k) > allowedErrors) {
        for(int i = 0; i < MufDiffGlobals::sortedIndicesSize; i++) {
            MufDiffGlobals::sortedIndices[i] = -1;
        }
        return 1.0f;
    }

    // Sort pairs using counting-sort
    vector<IndexPair> sortedList;
    for(int i = 0; i < MufDiffGlobals::sortedIndicesSize; i++) {
        if(MufDiffGlobals::sortedIndices[i] >= 0) {
            sortedList.push_back(IndexPair(i, MufDiffGlobals::sortedIndices[i]));
            MufDiffGlobals::sortedIndices[i] = -1; // Maintain array
        }
    }

    // Remove crossing values
    removeCrossingIndices(sortedList);

    //cout << biggest->sequence << endl << endl;
    //cout << smallest->sequence << endl << endl;

    /*for(vector<IndexPair>::iterator it = sortedList.begin(); it != sortedList.end(); it++)
        cout << it->i1 << ":" << it->i2 << ", ";
    cout << endl << endl;*/

    //
    int difference = 0;
    int lastHitB = -1;
    int lastHitS = -1;
    int firstHitB = -1;
    int firstHitS = -1;
    for(vector<IndexPair>::iterator it = sortedList.begin(); it != sortedList.end(); it++) {
        if(it->i1-1 != lastHitB && firstHitB >= 0) { // When hit, prev index can't be hit (there must be a hole in hits)
            int missLenB = it->i1 - lastHitB;
            int missLenS = it->i2 - lastHitS;
            difference += levenshteinHelper(&biggest->sequence.data()[lastHitB], missLenB,
                                             &smallest->sequence.data()[lastHitS], missLenS,
                                             allowedErrors - difference);

            if(difference > allowedErrors) {
                return 1.0f;
            }
            /*cout << "B(" << lastHitB << ", " << missLenB << ") S(" << lastHitS << ", " << missLenS << ")\n";
            cout << biggest->sequence.substr(lastHitB, missLenB) << endl;
            cout << smallest->sequence.substr(lastHitS, missLenS) << endl;
            cout << "Diff: " << difference << endl << endl;*/
        }
        if(lastHitB == -1) { // Set first hit
            firstHitB = it->i1;
            firstHitS = it->i2;
        }
        lastHitB = it->i1;
        lastHitS = it->i2;
    }

    if(lastHitB >= 0) {
        int tMin = (firstHitB < firstHitS) ? firstHitB : firstHitS;
        difference += levenshteinHelper(&biggest->sequence.data()[firstHitB-tMin], tMin,
                                        &smallest->sequence.data()[firstHitS-tMin], tMin,
                                        allowedErrors - difference);

        /*cout << "B(" << 0 << ", " << firstHitB << ") S(" << -firstDiff << ", " << firstHitB << ")\n";
        cout << biggest->sequence.substr(0, firstHitB) << endl;
        cout << smallest->sequence.substr(-firstDiff, firstHitB) << endl;
        cout << "Diff: " << difference << endl << endl;*/

        if(difference > allowedErrors) {
            return 1.0f;
        }

        // Handle sides (last)
        int endB = bigSize-lastHitB;
        int endS = smallSize-lastHitS;
        int firstToEnd = (endB > endS) ? endS : endB;
        if(firstToEnd > 0) { // Should always give true, but it's better to be safe.
            difference += levenshteinHelper(&biggest->sequence.data()[lastHitB], firstToEnd,
                                             &smallest->sequence.data()[lastHitS], firstToEnd,
                                             allowedErrors - difference);
            if(difference > allowedErrors) {
                return 1.0f;
            }
        }

        // Add difference if small is more than a sub-string of big
        // mmmhhhhhhhhh
        // hhhhhhhhhmmm
        if((firstHitB < firstHitS && endB > endS) ||
                (firstHitB > firstHitS && endB < endS)) {
            difference += abs(firstHitB - firstHitS) + abs(lastHitB - lastHitS);
        }
    }
    else { // 0% hit -> at least (100-(1/k))% error : (return as error)
        difference += bigSize;
    }

    return (float)(difference) / (float)smallSize;
}

float distanceLevenshteinFailFast(FastaContainer& seq1, FastaContainer& seq2, float threshold)
{
    std::string str1 = seq1.sequence;
    std::string str2 = seq2.sequence;
    // degenerate cases
    if (&str1 == &str2) return 0.0f;
    if (str1.length() == 0) return str2.length();
    if (str2.length() == 0) return str1.length();

    // create two work vectors of integer distances
    int* v0 = new int[str2.length() + 1];
    int* v1 = new int[str2.length() + 1];

    // Get data for calculating relative distance.
    int lengthDiff = (int) (str1.length() < str2.length()
                            ? str2.length() - str1.length()
                            : str1.length() - str2.length());
    int smallestLength = (int) (str1.length() < str2.length() ? str1.length() : str2.length());
    // Maximal errors before fail-fast kicks in
    float maxErrors = threshold * smallestLength + lengthDiff;

    // initialize v0 (the previous row of distances)
    // this row is A[0][i]: edit distance for an empty str1
    // the distance is just the number of characters to delete from str2
    for (int i = 0; i < (int)str2.length() + 1; i++) {
        v0[i] = i;
    }

    for (int i = 0; i < (int)str1.length(); i++)
    {
        // calculate v1 (current row distances) from the previous row v0

        // first element of v1 is A[i+1][0]
        //   edit distance is delete (i+1) chars from s to match empty t
        v1[0] = i + 1;
        int minErrors = v1[0];
        
        // use formula to fill in the rest of the row
        for (int j = 0; j < (int)str2.length(); j++)
        {
            int cost = ((str1[i]|32) == (str2[j]|32)) ? 0 : 1;
            v1[j + 1] = MIN3(v1[j] + 1, v0[j + 1] + 1, v0[j] + cost);
            if (v1[j + 1] < minErrors)
                minErrors = v1[j + 1];
        }
       
        if (minErrors > maxErrors) {
            delete[] v0;
            delete[] v1;
            // Fail Fast
            return 1.0f;
        }

        // Swap pointer v0 and v1
        int* tmp = v0;
        v0 = v1;
        v1 = tmp;
    }

    int result = v0[str2.length()];

    delete[] v0;
    delete[] v1;

    return (float)(result - lengthDiff) / smallestLength;
}

float pureKMer(FastaContainer& kMer1, FastaContainer& kMer2, float threshold)
{
    const int k = 8;

    int smallSize = (kMer1.sequence.size() < kMer2.sequence.size()) ?
                     kMer1.sequence.size() : kMer2.sequence.size();
    int allowedErrors = (int)ceil((float)smallSize * threshold) * k * 2;

    if(!kMer1.kMerHash.isCreated()) {
        kMer1.kMerHash.createHashMap((int)kMer1.sequence.size(), k<<1);
        generateHashKmer(kMer1.sequence, kMer1.kMerHash, k);
    }
    if(!kMer2.kMerHash.isCreated()) {
        kMer2.kMerHash.createHashMap((int)kMer2.sequence.size(), k<<1);
        generateHashKmer(kMer2.sequence, kMer2.kMerHash, k);
    }

    kMer1.kMerHash.iteratorReset();
    kMer2.kMerHash.iteratorReset();

    int difference = 0;
    while(kMer1.kMerHash.iteratorNotEnded() && kMer2.kMerHash.iteratorNotEnded()) {
        KMer* t1 = kMer1.kMerHash.iteratorGet();
        KMer* t2 = kMer2.kMerHash.iteratorGet();
        if(t1->value == t2->value) {
            kMer1.kMerHash.iteratorIncrease();
            kMer2.kMerHash.iteratorIncrease();
        }
        else {
            if(++difference > allowedErrors) return 1.0f;
            if(t1->value < t2->value) {
                kMer1.kMerHash.iteratorIncrease();
            }
            else { // t1 > t2
                kMer2.kMerHash.iteratorIncrease();
            }
        }
    }

    while(kMer1.kMerHash.iteratorNotEnded()) {
        if(++difference > allowedErrors) return 1.0f;
        kMer1.kMerHash.iteratorIncrease();
    }
    while(kMer2.kMerHash.iteratorNotEnded()) {
        if(++difference > allowedErrors) return 1.0f;
        kMer2.kMerHash.iteratorIncrease();
    }

    return (float)difference / (float)smallSize;
}