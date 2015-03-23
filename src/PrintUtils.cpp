/*
 * @Author: Andreas Dahl
 * @Author: Christian Muf
 * @Date: 23/03/15
 */

#include <iostream>
#include <iomanip>
#include "PrintUtils.h"

void printProgress(float progress) {
    unsigned int barLength = 50;
    unsigned long barsChars = (unsigned long)(progress * barLength);
    std::string bars(barsChars, '=');
    std::string space(barLength - barsChars, ' ');
    std::cout << "\r" << std::fixed << std::setprecision(2) << std::setw(6)
            << progress * 100.0f << "% |" << bars << space << "|    " << std::flush;
}