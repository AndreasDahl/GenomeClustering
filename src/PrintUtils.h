/*
 * Created by Andreas Dahl on 23/03/15.
 */

#ifndef _GENOMECLUSTERING_PRINTUTILS_H_
#define _GENOMECLUSTERING_PRINTUTILS_H_

#include <string>
#include <sys/time.h>

typedef unsigned long long timestamp_t;

timestamp_t get_timestamp();
void printProgress(float progress);
std::string formatDuration(timestamp_t t0, timestamp_t t1);

#endif //_GENOMECLUSTERING_PRINTUTILS_H_
