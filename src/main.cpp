/** @file
* @Author: Christian Muf
* @Author: Andreas Dahl
* @Date:   2015-03-03 00:59:41
* @Last Modified time: 2015-03-04 20:36:47
* @Version: 0.0
*/

#include <stdlib.h>
#include <stdio.h>
#include <string>

#include <math.h>

#include "FastaIO.h"
#include "KMerString.h"
#include "KMeans.h"
#include "SimpleGreedyClustering.h"

#include <iostream>
#include <vector>

using std::vector;

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

void simpleGreedyClusteringTest(char* file_path) {
	FastaIO fastaIO;
	fastaIO.openRead(file_path);

	vector<KMerString> strings(1000);

	for (std::vector<KMerString>::iterator it = strings.begin(); it != strings.end(); ++it) {
		fastaIO.getNextLine(it->getSequenceRef());
		it->gererateKMer();
	}

	vector<vector<KMerString>> result;

	simpleGreedyClustering<KMerString>(strings, kMerDistanceHellinger, 6.0f, result);

	// Print resulting clusters
	for (unsigned int i = 0; i < result.size(); ++i) {
		std::cout << "Cluster nr: " << i << " with " << result[i].size() << " strings" << std::endl;
	}

	fastaIO.closeRead();
}

void kMerTest1()
{
	std::string test1, test2;
	test1  = "gtatggtgcaagcgttatccggatttactgggtgtaaagggagcgtagacggAAAAGCAAGTC";
	test1 += "TGGAGTGAAAGCCCGGGGCTCAACCCCGGGACTGCTTTGGAAACTGTTATTCTTGAGTGCCGG";
	test1 += "AGAGGTAAGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAAAACCAGTGG";
	test1 += "CGAAGGCGGCTTACTGGACGGTAACTGACGTTGaggctcgaaagcgtggggagcaaacagg";

	test2  =    "tgttgcaagcgttatccggatttactgggtgtaaagggagcgtagacggAAAAGCAAGTC"; //
	test2 += "TGGAGTGAAAGCCCGGGGCTCAACCCCGGGACTGCTTTGGAAACTGTTATTCTTGAGTGCCGG";
	test2 += "AGAGGTAAGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAAAACCAGTGG";
	test2 += "CGAAGGCGGCTTACTGGACGGTAACTGACGTTGaggctcgaaagcgtggggagcaaac"; //

	KMerString mKer1(test1);
	KMerString mKer2(test2);

	mKer1.gererateKMer();
	mKer2.gererateKMer();

	std::cout << "Levenshtein: " << kMerDistanceLevenshtein(mKer1, mKer2) << std::endl;
	std::cout << "Manhattan: " << kMerDistanceManhattan(mKer1, mKer2) << std::endl;
	std::cout << "Hellinger: " << kMerDistanceHellinger(mKer1, mKer2) << std::endl;
	std::cout << "Test: " << kMerDistanceTest(mKer1, mKer2) << std::endl;
}

void kMerTest2(char* filePath)
{
	FastaIO fastaIO;
	fastaIO.openRead(filePath);

	KMerString mKer1, mKer2;

	fastaIO.getNextLine(mKer1.getSequenceRef());
	mKer1.gererateKMer();

	long mean = 0.0;
	long counter = 1;
	int max = 0;
	std::string maxString;
	while(!fastaIO.getNextLine(mKer2.getSequenceRef()) && counter < 100000)
	{
		mKer2.gererateKMer();
		mean += (long)mKer2.getSequenceRef().size();

		int temp = kMerDistanceTest(mKer1, mKer2);

		if(temp > max) {
			max = temp;
			maxString = mKer2.getSequenceRef();
		}

		counter += 1;
	}

	mKer2.getSequenceRef() = maxString;
	mKer2.gererateKMer();

	std::cout << "Count: " << counter << std::endl;
	std::cout << "Mean: " << mean / counter << std::endl;
	std::cout << "Max differnce: " << max << std::endl;
	std::cout << "Max differnce real: " << kMerDistanceLevenshtein(mKer1, mKer2) << std::endl;
	std::cout << std::endl << mKer1.getSequenceRef() << std::endl;
	std::cout << std::endl << mKer2.getSequenceRef() << std::endl;

	fastaIO.closeRead();
}

int main(int argc, char** argv)
{
	if(argc < 2)
		return -1;

	//kmeans_test(argv[1]);

//	kMerTest1();
//	kMerTest2(argv[1]);

	simpleGreedyClusteringTest(argv[1]);

	return 0;
}