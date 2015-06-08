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
        unsigned int permutations = 1 << (fa1.k << 1);  // 4^k == 2^(k*2)
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
        unsigned int permutations = 1 << (fa1.k << 1);  // 4^k == 2^(k*2)
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

#endif  // _GENOMECLUSTERING_FASTACOMPARATORS_H_
