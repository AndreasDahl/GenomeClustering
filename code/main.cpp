/** @file
* @Author: Christian Muf
* @Date:   2015-03-03 00:59:41
* @Last Modified time: 2015-03-03 00:59:41
* @Version: 0.0
*/

#include <stdio.h>
#include <string>

#include "FastaIO.h"

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

	s2c = s2[0];
	for(y = 1, lastdiag = 0; y <= s1len; y++)
	{
		olddiag = y;
		column[y] = MIN2(y, lastdiag + (s1[y-1] == s2c ? 0 : 1));
		lastdiag = olddiag;
	}

	for(x = 1; x < s2len; x++)
	{
		column[0] = x+1;
		s2c = s2[x];

		for(y = 1, lastdiag = x; y <= s1len; y++)
		{
			olddiag = column[y];
			prevdiag = column[y] = MIN3(olddiag+1, prevdiag+1, lastdiag + (s1[y-1] == s2c ? 0 : 1));
			lastdiag = olddiag;
		}
	}

	unsigned int result = column[s1len];
	delete column;
	return result;
}

int main(int argc, char** argv)
{
	if(argc < 2)
		return -1;

	std::string test1;
	std::string test2;

	FastaIO fastaIO = FastaIO();
	fastaIO.openRead(argv[1]);

	fastaIO.getNextLine(test1);

	int maxLen = 0;
	int maxDif = 0;
	int minDif = 0;
	int count = 0;
	while(!fastaIO.getNextLine(test2)) {
		int dif = levenshtein(test1, test2);
		if(test2.size() > maxLen) maxLen = test2.size();
		if(dif > maxDif) maxDif = dif;
		if(dif < minDif) minDif = dif;
		if(++count >= 100000) break;
	}

	std::cout << "Count: " << count << std::endl;
	std::cout << "Max length: " << maxLen << std::endl;
	std::cout << "Max difference: " << maxDif << std::endl;
	std::cout << "Min difference: " << minDif << std::endl;

	fastaIO.closeRead();

	return 0;
}