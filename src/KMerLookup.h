/** @file
* @Author: Christian Muf
* @Date:   2015-03-04 13:46:30
* @Last Modified time: 2015-03-04 20:49:36
* @Version: 0.0
*/

#ifndef K_MER_LOOKUP_H
#define K_MER_LOOKUP_H

#include <string>

#include "KMerString.h"

class KMerLookup
{
	public:
		KMerLookup(const std::string& sequence = "");
		KMerLookup(const KMerLookup& kMer);
		~KMerLookup();

		KMerLookup& operator=(const KMerLookup& kMer);

		int lookupDistance(const KMerString& other);

	private:
		void generateLookup();
		void deleteLookup();

	private:
		std::string m_sequence;

		int m_tableLength;
		int*** m_lookupTable;
};

#endif