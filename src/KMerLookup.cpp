/** @file
* @Author: Christian Muf
* @Date:   2015-03-04 13:46:36
* @Last Modified time: 2015-03-04 20:12:30
* @Version: 0.0
*/

#include "KMerLookup.h"

KMerLookup::KMerLookup(const std::string& sequence) :
	m_sequence(sequence),
	m_tableLength(0),
	m_lookupTable(NULL)
{
	generateLookup();
}

KMerLookup::KMerLookup(const KMerLookup& kMer) :
	m_sequence(kMer.m_sequence),
	m_tableLength(kMer.m_tableLength),
	m_lookupTable(NULL)
{
	generateLookup();
}

KMerLookup::~KMerLookup()
{
	deleteLookup();
}

KMerLookup& KMerLookup::operator=(const KMerLookup& kMer)
{
	if(this != &kMer) {
		m_sequence = kMer.m_sequence;
		m_tableLength = kMer.m_tableLength;
		m_lookupTable = NULL;

		generateLookup();
	}
	return *this;
}

int KMerLookup::lookupDistance(const KMerString&)
{
	return 0;
}

void KMerLookup::generateLookup()
{
	deleteLookup();

	if(m_tableLength > 0)
	{
		m_lookupTable = new int**[m_tableLength];
		for(int i = 0; i < m_tableLength; i++)
			m_lookupTable[i] = NULL;
	}
}

void KMerLookup::deleteLookup()
{
	int** leftPointer;

	if(m_lookupTable != NULL)
	{
		for(int i = 0; i < m_tableLength; i++)
		{
			leftPointer = m_lookupTable[i];
			if(leftPointer != NULL)
			{
				for(int j = 0; j <= i; j++)
				{
					if(leftPointer[j] != NULL)
						delete leftPointer[j];
				}

				delete leftPointer;
			}
		}

		delete m_lookupTable;
		m_lookupTable = NULL;
	}
}