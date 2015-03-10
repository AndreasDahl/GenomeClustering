/** @file
* @Author: Christian Muf
* @Date:   2015-03-06 15:46:54
*/

#ifndef WINDOWS_COMPAT_H
#define WINDOWS_COMPAT_H

#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
#define _hypot hypot // Workaround for MingW bug in math.h
#endif

#endif