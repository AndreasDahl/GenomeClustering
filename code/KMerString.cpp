/** @file
* @Author: Christian Muf
* @Date:   2015-03-04 19:44:05
* @Last Modified time: 2015-03-04 20:30:43
* @Version: 0.0
*/

#include "KMerString.h"

KMerString::KMerString() :
	m_sequence(""),
	m_kMerLength(0),
	m_kMer(NULL),
	m_kMerSorted(NULL)
{}

KMerString::KMerString(const std::string& sequence) :
	m_sequence(sequence),
	m_kMerLength(0),
	m_kMer(NULL),
	m_kMerSorted(NULL)
{}

KMerString::KMerString(const KMerString& kMer) :
	m_sequence(kMer.m_sequence),
	m_kMerLength(kMer.m_kMerLength),
	m_kMer(NULL),
	m_kMerSorted(NULL)
{
	if(m_kMerLength > 0)
		gererateKMer();
}

KMerString::~KMerString()
{
	deleteKMer();
}

KMerString& KMerString::operator=(const KMerString& kMer)
{
	if(this != &kMer) {
		m_sequence = kMer.m_sequence;
		m_kMerLength = 0;
		m_kMer = NULL;
		m_kMerSorted = NULL;

		if(kMer.m_kMerLength > 0)
			gererateKMer();
	}
	return *this;
}

void KMerString::gererateKMer()
{
	int a, b, c, current;

	deleteKMer();

	m_kMerLength = m_sequence.size() - 2;

	m_kMer = new int[m_kMerLength];
	m_kMerSorted = new int[c_sortedLength];

	for(unsigned int i = 0; i < c_sortedLength; i++)
		m_kMerSorted[i] = 0;

	a = m_sequence[0] & 6;
	b = m_sequence[1] & 6;
	for(unsigned int i = 0; i < m_kMerLength; i++)
	{
		c = m_sequence[i+2] & 6;

		current = m_kMer[i] = (a << 3) | (b << 1) | (c >> 1);
		m_kMerSorted[current] += 1;

		a = b;
		b = c;
	}
}

std::string& KMerString::getSequenceRef()
{
	return m_sequence;
}

const std::string& KMerString::getSequenceRef() const
{
	return m_sequence;
}

const int* KMerString::kMer() const
{
	return m_kMer;
}

unsigned int KMerString::kMerLength() const
{
	return m_kMerLength;
}

const int* KMerString::kMerSorted() const
{
	return m_kMerSorted;
}

unsigned int KMerString::kMerSortedLength() const
{
	return c_sortedLength;
}

void KMerString::deleteKMer()
{
	m_kMerLength = 0;

	if(m_kMer != NULL)
		delete m_kMer;

	if(m_kMerSorted != NULL)
		delete m_kMerSorted;

	m_kMer = NULL;
	m_kMerSorted = NULL;
}
