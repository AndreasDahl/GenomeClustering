/* 
* @Author: andreasdahl
* @Date:   2015-05-21 18:47:47
*/

#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <set>
#include "InputParser.h"
#include "MufDifference.h"
#include "GreedyClustering.h"
#include "PrintUtils.h"

const char* HELP = {
    "usage: %s <fasta in> <cluster out> <similarity> [<args>]\n"
};

void printHelp(char* programName) {
    printf(HELP, programName);
}

int parseInput(int argc, char** argv) {
    try {
        if (argc < 4)
            throw 1;

        // Open in
        FastaIO fastaIO;
        if (fastaIO.openRead(argv[1])) {
            throw 2;
        }

        // Open out
        std::ofstream out;
        out.open(argv[2]);
        if (!out.is_open()) {
            throw 3;
        }

        GreedyClustering setup(1.0f - atof(argv[3]));
             
        for (int i = 4; i < argc; ++i) {
            std::string argument;
            argument = argv[i];
            if (argument == "--cache_size" || argument == "-c") {
                if (++i < argc) {
                    setup.setCacheSize(atoi(argv[i]));
                } else {
                    throw 5;
                }
            } else if (argument == "--lru_size" || argument == "-r") {
                if (++i < argc) {
                    setup.setLRUSize(atoi(argv[i]));
                } else {
                    throw 5;
                }
            } else if (argument == "--lfu_size" || argument == "-f") {
                if (++i < argc) {
                    setup.setLFUSize(atoi(argv[i]));
                } else {
                    throw 5;
                }
            }
        }
        
        timestamp_t t1 = get_timestamp();
        setup.start(fastaIO, mufDifference, &out);
        timestamp_t t2 = get_timestamp();
        std::cout << "Execution took " << formatDuration(t1, t2) << std::endl;
        
        // Open stats
        //std::ofstream stats;
        //stats.open("stats.csv", std::fstream::out | std::fstream::app);
        //if (!stats.is_open()) {
        //    throw 6;
        //}
        //stats << setup.getLRUSize() << ';' 
        //      << setup.getLFUSize() << ';' 
        //      << (t2 - t1) << ';'
        //      << setup.getClusterCount() << std::endl << std::flush; 

        std::cout << std::endl
            << "Clusters: " << setup.getClusterCount() << std::endl
            << "    Time: " << formatDuration(t1, t2) << std::endl
            << std::flush;

        // Close streams
        fastaIO.closeRead();
        out.close();
        //stats.close();
        return 0;
    } catch (int e) {
        std::cout << e;
        printHelp(argv[0]);
        return e;
    }
    return 0;
}
