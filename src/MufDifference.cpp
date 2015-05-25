/** @file
* @Author: Christian Muf
* @Date:   2015-03-17 12:58:21
*/

#include "MufDifference.h"

#include "KMerHash.h"

#include <iostream>
#include <string>
#include <math.h>
#include <list>


using namespace std;

inline void generateHashKmer(const std::string& seq, KMerHashmap& data, unsigned int k)
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

int levenshteinHelper(const std::string& s1, const std::string& s2, float threshold);

inline int shiftingComparison(KMerHashmap& bigStruct, KMerHashmap& smallStruct, short* indexArray)
{
    bigStruct.iteratorReset();
    smallStruct.iteratorReset();

    int difference = 0;
    while(bigStruct.iteratorNotEnded() && smallStruct.iteratorNotEnded()) {
        KMer t1 = bigStruct.iteratorGet();
        KMer t2 = smallStruct.iteratorGet();
        if(t1.value == t2.value) {
            bigStruct.iteratorIncrease();
            smallStruct.iteratorIncrease();
            indexArray[t1.index] = t2.index;
        }
        else {
            difference += 1;
            if(t1.value < t2.value) {
                bigStruct.iteratorIncrease();
            }
            else { // t1 > t2
                smallStruct.iteratorIncrease();
            }
        }
    }

    return difference;
}

float mufDifference(FastaContainer& str1, FastaContainer& str2)
{
    return mufDifference(str1, str2, 1.0f);
}

float mufDifference(FastaContainer& str1, FastaContainer& str2, float threshold)
{
    // Set constant and static variables; and k for k-mers
    const int k = 9;
    const int s_arrayIndexSize = 16384;
    static short s_indexArray[s_arrayIndexSize];
    static bool s_initialized = false;

    // Return 0% errors if str1 and str2 is the same object.
    if(&str1 == &str2)
        return 0.0f;

    // Initialized s_indexArray. s_indexArray is maintained by the main loop.
    if(!s_initialized) {
        s_initialized = true;
        for(int i = 0; i < s_arrayIndexSize; i++)
            s_indexArray[i] = -1;
    }

    // Find the biggest of the two inputs.
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

    // Calculate the allowed number of errors.
    int allowedErrors = (int)ceil((float)smallSize * threshold);

    // Build the hashmaps is they don't exists.
    if(!biggest->kMerHash.isCreated()) {
        biggest->kMerHash.createHashMap(bigSize, k<<1);
        generateHashKmer(biggest->sequence, biggest->kMerHash, k);
    }
    if(!smallest->kMerHash.isCreated()) {
        smallest->kMerHash.createHashMap(smallSize, k<<1);
        generateHashKmer(smallest->sequence, smallest->kMerHash, k);
    }

    
    int kMerManhattanDistance = shiftingComparison(biggest->kMerHash, smallest->kMerHash, s_indexArray) >> 1;
    if(kMerManhattanDistance/k > allowedErrors)
        return (float)(kMerManhattanDistance/k);

    bool debug = false;

    /*if(debug) cout << biggest->sequence << endl;
    if(debug) cout << smallest->sequence << endl << endl;

    if(debug) {
        cout << bigSize << "-" << smallSize << " difference: " << bigSize-smallSize << endl;
        cout << "K-mer distance: " << res << endl;
        for(int i = 0; i < bigSize; i++) {
            if(s_indexArray[i] >= 0) cout << i << ":" << s_indexArray[i] << ", ";
        }
        cout << endl << endl;
    }*/

    int difference = 0;
    int lastHitB = -1;
    int lastHitS = -1;
    int firstHitB = -1;
    int firstHitS = -1;
    for(int i = 0; i < bigSize; i++) {
        if(s_indexArray[i] >= 0 && // Hit
            ((i-1 > 0 && s_indexArray[i-1] >= 0) || (i+1 < bigSize && s_indexArray[i+1] >= 0)) && // Not solo hit: --+---
            s_indexArray[i] > lastHitS && // Not less then previous
            true) // Not heigher than next TODO
        { 
            if(i-1 != lastHitB && firstHitB >= 0) { // When hit, last index can't be hit (there must be a hole in hits)
                int missLenB = i - lastHitB;
                int missLenS = s_indexArray[i] - lastHitS;
                difference += levenshteinHelper(
                    biggest->sequence.substr(lastHitB+k, missLenB-k),
                    smallest->sequence.substr(lastHitS+k, missLenS-k), 1.0f);

                /*if(debug) cout << "B(" << lastHitB << ", " << missLenB << ") S(" << lastHitS << ", " << missLenS << ")\n";
                if(debug) cout << biggest->sequence.substr(lastHitB+k, missLenB-k) << endl;
                if(debug) cout << smallest->sequence.substr(lastHitS+k, missLenS-k) << endl;
                if(debug) cout << "Diff: " << difference << endl;*/
            }
            if(lastHitB == -1) { // Set first hit
                firstHitB = i;
                firstHitS = s_indexArray[i];
            }
            lastHitB = i;
            lastHitS = s_indexArray[i];
            s_indexArray[i] = -1;
        }
    }
    //if(debug) cout << "--Start the end" << endl << flush;
    if(lastHitB >= 0) {
        // Handle sides (first)
        int firstDiff = firstHitB - firstHitS;
        if(firstDiff == 0 && firstHitB > 0) {
            difference += levenshteinHelper(biggest->sequence.substr(0, firstHitB),
                                            smallest->sequence.substr(0, firstHitB), 1.0f);
            if(debug) cout << "B(" << 0 << ", " << firstHitB << ") S(" << 0 << ", " << firstHitB << ")\n";
            if(debug) cout << biggest->sequence.substr(0, firstHitB) << endl;
            if(debug) cout << smallest->sequence.substr(0, firstHitB) << endl;
        }
        else if(firstDiff > 0 && firstHitS > 0) {
            difference += levenshteinHelper(biggest->sequence.substr(firstDiff, firstHitS),
                                            smallest->sequence.substr(0, firstHitS), 1.0f);
            if(debug) cout << "B(" << firstDiff << ", " << firstHitS << ") S(" << 0 << ", " << firstHitS << ")\n";
            if(debug) cout << biggest->sequence.substr(firstDiff, firstHitS) << endl;
            if(debug) cout << smallest->sequence.substr(0, firstHitS) << endl;
        }
        else if(firstDiff < 0 && firstHitB > 0) {
            difference += levenshteinHelper(biggest->sequence.substr(0, firstHitB),
                                            smallest->sequence.substr(-firstDiff, firstHitB), 1.0f);
            if(debug) cout << "B(" << 0 << ", " << firstHitB << ") S(" << -firstDiff << ", " << firstHitB << ")\n";
            if(debug) cout << biggest->sequence.substr(0, firstHitB) << endl;
            if(debug) cout << smallest->sequence.substr(-firstDiff, firstHitB) << endl;
        }
        if(debug) cout << "Diff: " << difference << endl;

        // Handle sides (last)
        int endB = bigSize-lastHitB;
        int endS = smallSize-lastHitS;
        int firstToEnd = (endB > endS) ? endS : endB;
        if(firstToEnd > 0) { // Should always give true, but it's better to be safe.
            difference += levenshteinHelper(biggest->sequence.substr(lastHitB, firstToEnd),
                                            smallest->sequence.substr(lastHitS, firstToEnd), 1.0f);
            if(debug) cout << "B(" << lastHitB << ", " << firstToEnd << ") S(" << lastHitS << ", " << firstToEnd << ")\n";
            if(debug) cout << biggest->sequence.substr(lastHitB, firstToEnd) << endl;
            if(debug) cout << smallest->sequence.substr(lastHitS, firstToEnd) << endl;
        }
        if(debug) cout << "Diff: " << difference << endl << endl;

        // Add difference if small is more than a sub-string of big
        // mmmhhhhhhhhh
        // hhhhhhhhhmmm
        if(debug) cout << firstHitB << ":" << firstHitS << endl << lastHitB << ":" << lastHitS << endl;
        if((firstHitB < firstHitS && endB > endS) ||
                (firstHitB > firstHitS && endB < endS)) {
            difference += abs(firstHitB - firstHitS) + abs(lastHitB - lastHitS);
        }
    }
    else { // 0% hit -> at least (100-(1/k))% error : (return as error)
        difference += bigSize;
    }

    return (float)difference / (float)smallSize;
}

#define MIN2(a, b) (a) < (b) ? (a) : (b)
#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
int levenshteinHelper(const std::string& str1, const std::string& str2, float threshold)
{
    // degenerate cases
    if (&str1 == &str2) return 0;
    if (str1.length() == 0) return str1.length();
    if (str2.length() == 0) return str2.length();

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
            //cout << v1[j+1] << " ";
        }
        //cout << endl;
       
        if (minErrors > maxErrors) {
            // Fail Fast
            return 1.0f;
        }

        // Swap pointer v0 and v1
        int* tmp = v0;
        v0 = v1;
        v1 = tmp;
    }

    return v0[str2.length()];
}

float distanceLevenshteinFailFast(FastaContainer& kMer1, FastaContainer& kMer2, float threshold)
{
    std::string str1 = kMer1.sequence;
    std::string str2 = kMer2.sequence;
    // degenerate cases
    if (&str1 == &str2) return 0;
    if (str1.length() == 0) return str1.length();
    if (str2.length() == 0) return str2.length();

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
            // Fail Fast
            return 1.0f;
        }

        // Swap pointer v0 and v1
        int* tmp = v0;
        v0 = v1;
        v1 = tmp;
    }

    return (float)(v0[str2.length()] - lengthDiff) / smallestLength;
}

