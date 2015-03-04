/** @file
* @Author: Christian Muf
* @Date:   2015-03-03 00:59:41
* @Last Modified time: 2015-03-04 20:36:47
* @Version: 0.0
*/

#include <stdlib.h>
#include <stdio.h>
#include <string>

#define hypot _hypot
#include <math.h>

#include "FastaIO.h"
#include "KMerString.h"

#include <iostream>

#define MIN2(a, b) (a) < (b) ? (a) : (b)
#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
int levenshtein(const std::string& s1, const std::string& s2)
{
	unsigned int s1len, s2len, x, y, lastdiag, olddiag, prevdiag;
	char s2c;

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

// U == T
// Convert 3 chars (acgt/ACGT) to a 6-bit 3-mer char. 
int char4_3mer(char a, char b, char c)
{
	return ((a & 6) >> 1) | ((b & 6) << 1) | ((c & 6) << 3);
}

void print3mer(std::string& in, int* sum)
{
	for(int i = 0; i < 64; i++)
		sum[i] = 0;

	for(unsigned int i = 0; i < in.size()-2; i++)
		sum[char4_3mer(in[i], in[i+1], in[i+2])] += 1;

	std::cout << in.size() << " :: ";

	for(int i = 0; i < 64; i++)
		std::cout << sum[i] << ", ";

	std::cout << std::endl;
}

int main(int argc, char** argv)
{
	if(argc < 2)
		return -1;

	KMerString test1;
	KMerString test2;

	FastaIO fastaIO;
	fastaIO.openRead(argv[1]);

	fastaIO.getNextLine(test1.getSequenceRef());
	fastaIO.getNextLine(test2.getSequenceRef());

	test1.gererateKMer();
	test2.gererateKMer();

	std::cout << test1.getSequenceRef() << std::endl;
	std::cout << test2.getSequenceRef() << std::endl;

	std::cout << levenshtein(test1.getSequenceRef(), test2.getSequenceRef()) << std::endl;

	float Manhattan = 0;
	float Hellinger = 0;
	for(int i = 0; i < 64; i++) {
		Manhattan += (float)abs(test1.kMerSorted()[i] - test2.kMerSorted()[i]);

		float x = (float)test1.kMerSorted()[i];
		float y = (float)test2.kMerSorted()[i];
		Hellinger += x + y - 2.0f * sqrt(x * y);
	}

	std::cout << "Manhattan: " << Manhattan << std::endl;
	std::cout << "Hellinger: " << Hellinger << std::endl;

	/*
	unsigned int maxLen = 0;
	int maxDif = 0;
	int minDif = 0;
	int count = 0;
	while(!fastaIO.getNextLine(test2)) {
		int dif = levenshtein(test1, test2);
		if(test2.size() > maxLen) maxLen = test2.size();
		if(dif > maxDif) maxDif = dif;
		if(dif < minDif) minDif = dif;
		if(++count >= 10000) break;
	}

	std::cout << "Count: " << count << std::endl;
	std::cout << "Max length: " << maxLen << std::endl;
	std::cout << "Max difference: " << maxDif << std::endl;
	std::cout << "Min difference: " << minDif << std::endl;
	*/

	fastaIO.closeRead();

	return 0;
}