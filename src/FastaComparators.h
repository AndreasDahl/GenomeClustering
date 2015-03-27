/*
 * @Author: Andreas Dahl
 * @Author: Christian Muf
 * @Date: 23/03/15
 */

#ifndef _GENOMECLUSTERING_FASTACOMPARATORS_H_
#define _GENOMECLUSTERING_FASTACOMPARATORS_H_

struct LengthThenKmer
{
    static bool less_than(const FastaContainer& fa1, const FastaContainer& fa2)
    {
        unsigned int permutations = 1 << (fa1.k << 1); // 4^k == 2^(k*2)
        if (fa1.sequence.size() < fa2.sequence.size()) {
            return true;
        } else if (fa1.sequence.size() == fa2.sequence.size()) {
            for (unsigned int i = 0; i < permutations; ++i) {
                if (fa1.getKMer()[i] < fa2.getKMer()[i]) {
                    return true;
                }
            }
        }
        return false;
    }

    static bool greater_than(const FastaContainer& fa1, const FastaContainer& fa2) {
        return LengthThenKmer::less_than(fa2, fa1);
    }
};

struct Kmer
{
    static bool less_than(const FastaContainer& fa1, const FastaContainer& fa2)
    {
        unsigned int permutations = 1 << (fa1.k << 1); // 4^k == 2^(k*2)
        for (unsigned int i = 0; i < permutations; ++i) {
            if (fa1.getKMer()[i] < fa2.getKMer()[i]) {
                return true;
            }
        }
        return false;
    }

    static bool greater_than(const FastaContainer& fa1, const FastaContainer& fa2) {
        return Kmer::less_than(fa2, fa1);
    }
};

#endif //_GENOMECLUSTERING_FASTACOMPARATORS_H_
