/** @file
* @Author: Christian Muf
* @Author: Andreas Dahl
* @Date:   2015-03-03 00:59:41
*/

#include <stdlib.h>
#include <stdio.h>
#include <string>

#include <math.h>

#include "FastaIO.h"
#include "MufDifference.h"
#include "KMeans.h"
#include "GreedyClustering.h"
#include "FastaComparators.h"

#include <iostream>
#include <vector>

using std::vector;

using std::list;

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
	}

	greedySettings settings = greedySettings(0.03f);
	std::ofstream nullstream = std::ofstream(0);
	greedyClustering<FastaContainer>(strings, mufDifference, settings, nullstream);
	timestamp_t t1 = get_timestamp();

	std::cout << "Execution took " << formatDuration(t0, t1) << " to complete." << std::endl;

	fastaIO.closeRead();
}

void distance_challenge(char* file_path) {
	FastaIO fastaIO;
	fastaIO.openRead(file_path);

	vector<FastaContainer> strings;

	for (unsigned int i = 0; i < 500; ++i) {
		strings.push_back(FastaContainer());
		if(fastaIO.getNextLine(strings.back())) {
			strings.pop_back();
			break;
		}
	}

	std::ofstream myfile;
	myfile.open ("challenge.csv");
	timestamp_t t0 = get_timestamp();
	unsigned int i = 0;
	for (std::vector<FastaContainer>::iterator it1 = strings.begin(); it1 != strings.end(); ++it1) {
		for (std::vector<FastaContainer>::iterator it2 = strings.begin(); it2 != strings.end(); ++it2) {
			myfile << mufDifference(*it1, *it2);
		}
	}
	timestamp_t t1 = get_timestamp();
	myfile.close();

	std::cout << "Execution took " << formatDuration(t0, t1) << " to complete";
}

int main(int argc, char** argv)
{
	if(argc < 2)
		return -1;

	simpleGreedyClusteringTest(argv[1]);

//	distance_challenge(argv[1]);

	return 0;
}