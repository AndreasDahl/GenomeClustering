/*
 * Copyright 2015 Andreas Dahl, Christian Muf
 *
 * This file is part of MufDahlClust.
 *
 * MufDahlClust is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MufDahlClust is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MufDahlCLust.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * @Author: Andreas Dahl
 * @Author: Christian Muf
 * @Date: 23/03/15
 */

#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <sys/time.h>
#include "PrintUtils.h"

timestamp_t get_timestamp()
{
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

void printProgress(unsigned long current, unsigned long max) {
    float progress;
    if (current >= max || max == 0) {
        progress = 1.0f;
    } else {
        progress = (float) current / max;
    }
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
