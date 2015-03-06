/** @file
* @Author: Christian Muf
* @Date:   2015-03-04 19:43:58
* @Last Modified time: 2015-03-04 20:19:00
* @Version: 0.0
*/

#ifndef K_MER_STRING_H
#define K_MER_STRING_H

#include <string>

class KMerString
{
	public:
		static const unsigned int c_sortedLength = 64;

	public:
		KMerString();
		KMerString(const std::string& sequence);
		KMerString(const KMerString& kMer);
		~KMerString();

		KMerString& operator=(const KMerString& kMer);

		bool isGenerated() const;
		void gererateKMer();

		std::string& getSequenceRef();
		const std::string& getSequenceRef() const; 

		const int* kMer() const;
		unsigned int kMerLength() const;

		const int* kMerSorted() const;
		unsigned int kMerSortedLength() const;

	private:
		void deleteKMer();

	private:
		std::string m_sequence;

		unsigned int m_kMerLength;
		int* m_kMer;
		int* m_kMerSorted;
};

int kMerDistanceTest(const KMerString& kMer1, const KMerString& kMer2);

float kMerDistanceHellinger(const KMerString& kMer1, const KMerString& kMer2);
int kMerDistanceManhattan(const KMerString& kMer1, const KMerString& kMer2);

int kMerDistanceLevenshtein(const KMerString& kMer1, const KMerString& kMer2);


#endif