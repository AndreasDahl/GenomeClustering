/** @file
* @Author: Christian Muf
* @Author: Andreas Dahl
* @Date:   2015-03-03 00:59:41
*/

#include <stdlib.h>
#include <stdio.h>
#include <string>

#include <math.h>
#include <sys/time.h>

#include "FastaIO.h"
#include "MufDifference.h"
#include "KMeans.h"
#include "SimpleGreedyClustering.h"
#include "FastaComparators.h"

#include <iostream>
#include <vector>

using std::vector;

using std::list;

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp()
{
	struct timeval now;
	gettimeofday (&now, NULL);
	return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

void simpleGreedyClusteringTest(char* file_path) {
	FastaIO fastaIO;
	fastaIO.openRead(file_path);

	vector<FastaContainer> strings;

	timestamp_t t0 = get_timestamp();
	while (true) {
		strings.push_back(FastaContainer());
		if(fastaIO.getNextLine(strings.back())) {
			strings.pop_back();
			break;
		}
		createCountedKmer(strings.back(), 5);
	}
	std::sort(strings.begin(), strings.end(), Kmer::greater_than);

	SimpleGreedySettings settings = SimpleGreedySettings(0.03f);
//	settings.greedyPick = true;
//	settings.cacheSize = 64;
	simpleGreedyClustering<FastaContainer>(strings, mufDifference, settings);
	timestamp_t t1 = get_timestamp();

	std::cout << (t1 - t0) / 1000000.0L << std::endl;

	fastaIO.closeRead();
}

int main(int argc, char** argv)
{
	if(argc < 2)
		return -1;

	simpleGreedyClusteringTest(argv[1]);

	return 0;
}