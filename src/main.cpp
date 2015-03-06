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
#include "KMeans.h"

#include <iostream>
#include <vector>
using std::vector;

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

void kmeans_test(char* file_path) {
	FastaIO fastaIO;
	fastaIO.openRead(file_path);

	vector<KMerString> strings(1000);

	for (std::vector<KMerString>::iterator it = strings.begin(); it != strings.end(); ++it) {
		fastaIO.getNextLine(it->getSequenceRef());
		it->gererateKMer();
	}

	int k = 5;
	vector<vector<KMerString>> result(k);

	kmeans(strings, k, result);

	// Print resulting clusters
	for (unsigned int i = 0; i < result.size(); ++i) {
		std::cout << "Cluster nr: " << i << " with " << result[i].size() << " strings" << std::endl;
	}

	fastaIO.closeRead();
}

int main(int argc, char** argv)
{
	if(argc < 2)
		return -1;

	kmeans_test(argv[1]);

	return 0;
}