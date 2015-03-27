/** @file
* @Author: Christian Muf
* @Author: Andreas Dahl
* @Date:   2015-03-03 00:59:41
*/

#include <stdlib.h>
#include <string>

#include "FastaIO.h"
#include "MufDifference.h"
#include "KMeans.h"
#include "GreedyClustering.h"
#include "PrintUtils.h"

#include <iostream>
#include <list>

using std::vector;

using std::list;

void greedyClusteringTest(char* file_path) {
	FastaIO fastaIO;
	fastaIO.openRead(file_path);

	timestamp_t t0 = get_timestamp();

	GreedySettings settings = GreedySettings(0.03f);
	greedyClustering(fastaIO, mufDifference, settings, NULL);
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
	for (std::vector<FastaContainer>::iterator it1 = strings.begin(); it1 != strings.end(); ++it1) {
		for (std::vector<FastaContainer>::iterator it2 = strings.begin(); it2 != strings.end(); ++it2) {
			myfile << mufDifference(*it1, *it2);
		}
	}
	timestamp_t t1 = get_timestamp();
	myfile.close();

	std::cout << "Execution took " << formatDuration(t0, t1) << " to complete";
}

void compareLevenshteinKmer(char* file_path) {
	FastaIO fastaIO;
	fastaIO.openRead(file_path);

	vector<FastaContainer> strings;

	while (true) {
		strings.push_back(FastaContainer());
		if(fastaIO.getNextLine(strings.back())) {
			strings.pop_back();
			break;
		}
	}

	std::cout << "size " << strings.size() << std::flush;
	std::ofstream myfile;
	myfile.open ("compare.csv");
	for (unsigned int i = 0; i < strings.size() - 1; ++i) {
		for (unsigned int j = i + 1; j < strings.size(); ++j) {
			myfile << mufDifference(strings[i], strings[j]);
			myfile << " ";
			myfile << kMerDistanceLevenshtein(strings[i], strings[j]);
			myfile << std::endl;
		}
		printProgress((float) i / strings.size());
	}
	std::cout.flush();
	myfile.close();

}

int main(int argc, char** argv)
{
	if(argc < 2)
		return -1;

	greedyClusteringTest(argv[1]);

//	distance_challenge(argv[1]);

//	compareLevenshteinKmer(argv[1]);

	return 0;
}