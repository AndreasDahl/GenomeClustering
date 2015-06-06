/*
 * @Author: Andreas Dahl
 * @Author: Christian Muf
 * @Date: 23/03/15.
 */

#ifndef _GENOMECLUSTERING_PRINTUTILS_H_
#define _GENOMECLUSTERING_PRINTUTILS_H_

#include <string>

typedef unsigned long long timestamp_t;

timestamp_t get_timestamp();
void printProgress(unsigned long current, unsigned long max);
std::string formatDuration(timestamp_t t0, timestamp_t t1);

#endif  // _GENOMECLUSTERING_PRINTUTILS_H_
