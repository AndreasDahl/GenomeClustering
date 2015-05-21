/* 
* @Author: andreasdahl
* @Date:   2015-05-21 18:47:47
*/

#include <stdlib.h>
#include <iostream>
#include <string>
#include <set>
#include <regex>
#include "InputParser.h"

void parseInput(int argc, char** argv) {
    std::regex stringPattern ("^--(.+)");
    std::regex charsPattern ("^-([^-]+)");

    std::set<std::string> arguments;
    for (int i = 1; i < argc; ++i) {
        std::cmatch matches;
        if (std::regex_match(argv[i], matches, stringPattern)) {
            arguments.insert(matches.str(1));
        } else if (std::regex_match(argv[i], matches, charsPattern)) {
            for (char& c : matches.str(1)) {
                std::string s (1,c);
                arguments.insert(s);
            }
        } else {
            // TODO: Make a fuss about it (Bad input error)
        }
    }

    std::cout << "myset contains:";
    for (std::set<std::string>::iterator it=arguments.begin(); it!=arguments.end(); ++it)
        std::cout << ' ' << *it;

    std::cout << '\n';
}