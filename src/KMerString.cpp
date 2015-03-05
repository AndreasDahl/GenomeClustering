/** @file
* @Author: Christian Muf
* @Date:   2015-03-04 19:44:05
* @Last Modified time: 2015-03-04 20:30:43
* @Version: 0.0
*/

#include "KMerString.h"

#include <math.h>
#include <stdlib.h>

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
	m_kMerLength(0),
	m_kMer(NULL),
	m_kMerSorted(NULL)
{
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

		gererateKMer();
	}
	return *this;
}

bool KMerString::isGenerated() const
{
	return m_sequence.size()-2 == m_kMerLength;
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

float kMerDistanceHellinger(const KMerString& kMer1, const KMerString& kMer2)
{
	if(!kMer1.isGenerated() || !kMer2.isGenerated())
		throw -1;

	float hellinger = 0;
	for(unsigned int i = 0; i < KMerString::c_sortedLength; i++) {
		float x = (float)kMer1.kMerSorted()[i];
		float y = (float)kMer2.kMerSorted()[i];
		hellinger += x + y - 2.0f * sqrt(x * y);
	}

	return hellinger;
}

int kMerDistanceManhattan(const KMerString& kMer1, const KMerString& kMer2)
{
	if(!kMer1.isGenerated() || !kMer2.isGenerated())
		throw -1;

	int manhattan = 0;
	for(unsigned int i = 0; i < KMerString::c_sortedLength; i++)
		manhattan += abs(kMer1.kMerSorted()[i] - kMer2.kMerSorted()[i]);

	return manhattan;
}

#define MIN2(a, b) (a) < (b) ? (a) : (b)
#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
int kMerDistanceLevenshtein(const KMerString& kMer1, const KMerString& kMer2)
{
	unsigned int s1len, s2len, x, y, lastdiag, olddiag, prevdiag;
	char s2c;

	const std::string& s1 = kMer1.getSequenceRef();
	const std::string& s2 = kMer2.getSequenceRef();

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
	return result;
}

