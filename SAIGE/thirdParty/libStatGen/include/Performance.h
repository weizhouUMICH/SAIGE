/*
 * Copyright (c) 2009 Regents of the University of Michigan
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _PERFORMANCE_H
#define _PERFORMANCE_H

#include <sys/time.h>
#include <time.h>

class Timing
{
    timeval startInterval;
    timeval endInterval;
public:
    Timing()
    {
        start();
    }
    void start();
    void end();
    double interval();
};

inline void Timing::start()
{
    gettimeofday(&startInterval, NULL);
}

inline void Timing::end()
{
    gettimeofday(&endInterval, NULL);
}

///
/// Return time interval between start() and end()
/// @return elapsed time in seconds
///
inline double Timing::interval()
{
    return (endInterval.tv_sec + (endInterval.tv_usec/1000000.0)) -
           (startInterval.tv_sec + (startInterval.tv_usec/1000000.0));
}

#endif
