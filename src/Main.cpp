/** @file
* @Author: Christian Muf
* @Author: Andreas Dahl
* @Date:   2015-03-03 00:59:41
*/

#include <stdlib.h>
#include <string>
#include <iostream>
#include <list>

#include "FastaIO.h"
#include "MufDifference.h"
#include "GreedyClustering.h"
#include "PrintUtils.h"
#include "InputParser.h"

void file_test(char* file_path) {
    FastaIO fastaIO;
    fastaIO.openRead(file_path);

    FastaContainer* fir = new FastaContainer();
    FastaContainer* sec = new FastaContainer();

    fastaIO.getNextLine(*fir);
    while(!fastaIO.getNextLine(*sec)) {
        mufDifference(*fir, *sec, 0.05f);
        delete sec;
        sec = new FastaContainer();
    }

    delete fir;
    delete sec;

    fastaIO.closeRead();
}

void distance_test(char* file_path) {
    FastaIO fastaIO;
    fastaIO.openRead(file_path);

    const int numTests = 500;
    FastaContainer seq[numTests];
    for(int i = 0; i < numTests; i++) {
        fastaIO.getNextLine(seq[i]);
    }    

    timestamp_t t0 = get_timestamp();
    float fail = 0.05f;
    for(int i = 0; i < numTests; i++)
    for(int j = 0; j < numTests; j++) {
        float d1 = mufDifference(seq[i], seq[j], fail);
        float d2 = distanceLevenshteinFailFast(seq[i], seq[j], fail);
        if((d1 > fail && d2 < fail) || (d1 < fail && d2 > fail))
            std::cout << "Error: " << i << ":" << j << " _ " << d1 << " : " << d2 << std::endl;
    }
    timestamp_t t1 = get_timestamp();
    std::cout << "Execution took " << formatDuration(t0, t1) << " to complete" << std::endl;

    fastaIO.closeRead();
}

void compareLevenshteinKmer(char* file_path) {
    const unsigned int s_count = 250;
    FastaIO fastaIO;
    fastaIO.openRead(file_path);


    FastaContainer seq[s_count];
    for(unsigned int i = 0; i < s_count; i++) {
        fastaIO.getNextLine(seq[i]);
    }    


    std::cout << "size " << s_count << std::flush;
    std::ofstream myfile;
    myfile.open ("compare.csv");
    for (unsigned int i = 0; i < s_count - 1; ++i) {
        printProgress(i, s_count);
        for (unsigned int j = i + 1; j < s_count; ++j) {
            timestamp_t t0 = get_timestamp();
            myfile << 1.0f - mufDifference(seq[i], seq[j], 1.0f);
            timestamp_t t1 = get_timestamp();
            myfile << ';' << t1 - t0 << ';';
            timestamp_t t2 = get_timestamp();
            myfile << 1.0f - distanceLevenshteinFailFast(seq[i], seq[j], 1.0f);
            timestamp_t t3 = get_timestamp();
            myfile << ';' << t3 - t2;
            myfile << std::endl;
        }
    }
    std::cout.flush();
    myfile.close();
}


int main(int argc, char** argv)
{
    initializeMufDifference();
    parseInput(argc, argv);

    //file_test(argv[1]);

    //compareLevenshteinKmer(argv[1]);

    uninitializeMufDifference();
    return 0;
}
