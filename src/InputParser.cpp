/* 
* @Author: andreasdahl
* @Date:   2015-05-21 18:47:47
*/

#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <set>
#include <regex>
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
    std::regex stringPattern ("^--(.+)");
    std::regex charsPattern ("^-([^-])$");

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

        GreedyClustering setup(atof(argv[3]));
             
        for (int i = 4; i < argc; ++i) {
            std::string argument;
            std::cmatch matches;
            if (std::regex_match(argv[i], matches, stringPattern)
                    || std::regex_match(argv[i], matches, charsPattern)) {
                argument = matches.str(1);
            } else {
                throw 4;
            }
            if (argument == "cache_size" || argument == "c") {
                if (++i < argc) {
                    setup.setCacheSize(atoi(argv[i]));
                } else {
                    throw 5;
                }

            }
        }
        
        timestamp_t t1 = get_timestamp();
        setup.start(fastaIO, mufDifference, &out);
        timestamp_t t2 = get_timestamp();
        std::cout << "Execution took " << formatDuration(t1, t2) << std::endl;

        // Close streams
        fastaIO.closeRead();
        out.close();
        return 0;
    } catch (int e) {
        printHelp(argv[0]);
        return e;
    }
}
