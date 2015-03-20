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
#include "SimpleGreedyClustering.h"

#include <iostream>
#include <vector>

using std::vector;
#include <list>
using std::list;

/*void kmeans_test(char* file_path) {
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
}*/

void simpleGreedyClusteringTest(char* file_path) {
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

	simpleGreedyClustering<FastaContainer>(strings, mufDifference, 0.03f);

	fastaIO.closeRead();
}

int main(int argc, char** argv)
{
	if(argc < 2)
		return -1;

//	kmeans_test(argv[1]);
//	kMerTest1();
//	kMerTest2(argv[1]);
//	findBiggest(argv[1]);
	simpleGreedyClusteringTest(argv[1]);


	return 0;
}