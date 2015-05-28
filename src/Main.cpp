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
#include "InputParser.h"

#include <iostream>
#include <list>

using std::vector;

using std::list;

void distance_test(char* file_path) {
    FastaIO fastaIO;
    fastaIO.openRead(file_path);

    const int numTests = 50;
    FastaContainer seq[numTests];
    for(int i = 0; i < numTests; i++) {
        fastaIO.getNextLine(seq[i]);
    }    

    timestamp_t t0 = get_timestamp();
    for(int i = 0; i < numTests; i++)
    for(int j = 0; j < numTests; j++) {
        mufDifference(seq[i], seq[j], 0.05f);
        //std::cout << "MufDiff: " << mufDifference(seq[i], seq[j], 1.0f) << std::endl;
        //std::cout << "Leven: " << distanceLevenshteinFailFast(seq[i], seq[j], 1.0f) << std::endl << std::endl;
    }
    timestamp_t t1 = get_timestamp();
    std::cout << "Execution took " << formatDuration(t0, t1) << " to complete" << std::endl;

    fastaIO.closeRead();
}

void distance_challenge(char* file_path) {
    FastaIO fastaIO;
    fastaIO.openRead(file_path);

    vector<FastaContainer> strings;

    for (unsigned int i = 0; i < 50; ++i) {
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
            myfile << distanceLevenshteinFailFast(*it1, *it2, 1.0f);
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
            myfile << distanceLevenshteinFailFast(strings[i], strings[j], 1.0f);
            myfile << std::endl;
        }
        printProgress(i, strings.size());
    }
    std::cout.flush();
    myfile.close();

}

int main(int argc, char** argv)
{
    initializeMufDifference();
    //return parseInput(argc, argv);

    distance_challenge(argv[1]);

//  compareLevenshteinKmer(argv[1]);
    uninitializeMufDifference();
}
