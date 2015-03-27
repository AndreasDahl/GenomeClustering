/** @file
* @Author: Christian Muf
* @Date:   2015-03-17 12:58:21
*/

#include "MufDifference.h"

#include <math.h>

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

