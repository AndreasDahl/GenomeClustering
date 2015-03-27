/*
 * @Author: Andreas Dahl
 * @Author: Christian Muf
 * @Date: 23/03/15
 */

#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include "PrintUtils.h"

timestamp_t get_timestamp()
{
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

void printProgress(float progress) {
    if (progress > 1.0f)
        progress = 0.0f;
    unsigned int barLength = 50;
    unsigned long barsChars = (unsigned long)(progress * barLength);
    std::string bars(barsChars, '=');
    std::string space(barLength - barsChars, ' ');
    std::cout << "\r" << std::fixed << std::setprecision(2) << std::setw(6)
            << progress * 100.0f << "% |" << bars << space << "|    " << std::flush;
}

std::string formatDuration(timestamp_t t0, timestamp_t t1) {
    long double seconds = (t1 - t0) / 1000000.0L;
    unsigned long minutes = (unsigned long) seconds / 60;
    seconds -= minutes * 60;
    std::stringstream duration;

    if (minutes > 0) {
        duration << minutes << "m ";
    }
    duration << seconds << "s";
    return duration.str();
}