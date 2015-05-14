/** @file
* @Author: Christian Muf
* @Date:   2015-03-17 12:58:21
*/

#include "MufDifference.h"

#include <math.h>
#include <string>
#include <iostream>

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

#define MIN2(a, b) (a) < (b) ? (a) : (b)
#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
float kMerDistanceLevenshtein(FastaContainer& kMer1, FastaContainer& kMer2)
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
    return (float)(result - ((s1len - s2len > 0) ? s1len - s2len : s2len - s1len)) /
        (float)((s1len < s2len) ? s1len : s2len);
}


float distanceLevenshteinFailFast(FastaContainer& kMer1, FastaContainer& kMer2, float threshold)
{
    std::string str1 = kMer1.sequence;
    std::string str2 = kMer2.sequence;
    // degenerate cases
    if (str1 == str2) return 0;
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
    for (int i = 0; i < str2.length() + 1; i++) {
        v0[i] = i;
    }

    for (int i = 0; i < str1.length(); i++)
    {
        // calculate v1 (current row distances) from the previous row v0

        // first element of v1 is A[i+1][0]
        //   edit distance is delete (i+1) chars from s to match empty t
        v1[0] = i + 1;
        int minErrors = v1[0];
        
        // use formula to fill in the rest of the row
        for (int j = 0; j < str2.length(); j++)
        {
            int cost = (str1[i] == str2[j]) ? 0 : 1;
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

    return (float)(v1[str2.length()] - lengthDiff) / smallestLength;
}

