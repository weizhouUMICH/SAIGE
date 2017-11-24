/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


//////////////////////////////////////////////////////////////////////////////
// This file includes code derived from the original Mersenne Twister Code
// by Makoto Matsumoto and Takuji Nishimura
// and is subject to their original copyright notice copied below:
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// COPYRIGHT NOTICE FOR MERSENNE TWISTER CODE
// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. The names of its contributors may not be used to endorse or promote
// products derived from this software without specific prior written
// permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
///////////////////////////////////////////////////////////////////////////////

#include "Random.h"
#include "MathConstant.h"
#include "Error.h"

#include <math.h>

//Constants used internally by Mersenne random number generator
#define MERSENNE_N 624
#define MERSENNE_M 397

// constant vector a
#define MATRIX_A 0x9908b0dfUL

// most significant w-r bits
#define UPPER_MASK 0x80000000UL

// least significant r bits
#define LOWER_MASK 0x7fffffffUL


// Constants used internally by Park-Miller random generator
#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define RNMX (1.0-EPS)

Random::Random(long s)
{
#ifndef __NO_MERSENNE
    mt = new unsigned long [MERSENNE_N];
    mti = MERSENNE_N + 1;
    mersenneMult = 1.0/4294967296.0;
#else
    shuffler = new long [NTAB];
#endif
    Reset(s);
}

Random::~Random()
{
#ifndef __NO_MERSENNE
    delete [] mt;
#else
    delete [] shuffler;
#endif
}

void Random::Reset(long s)
{
    normSaved = 0;

#ifndef __NO_MERSENNE
    InitMersenne(s);
#else
    // 'Continuous' Random Generator
    if ((seed = s) < 1)
        seed = s == 0 ? 1 : -s;  // seed == 0 would be disastrous

    for (int j=NTAB+7; j>=0; j--)   // Warm up and set shuffle table
    {
        long k = seed / IQ;
        seed = IA * (seed - k * IQ) - IR * k;
        if (seed < 0) seed += IM;
        if (j < NTAB) shuffler[j] = seed;
    }
    last=shuffler[0];
#endif
}

// initializes mt[MERSENNE_N] with a seed
void Random::InitMersenne(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti = 1; mti < MERSENNE_N; mti++)
    {
        mt[mti] = (1812433253UL * (mt[mti-1] ^(mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect */
        /* only MSBs of the array mt[].  */
        /* 2002/01/09 modified by Makoto Matsumoto */

        mt[mti] &= 0xffffffffUL;
    }
}

int Random::Binary()
{
    return Next() > 0.5 ? 1 : 0;
}

#ifndef __NO_MERSENNE

double Random::Next()
{
    unsigned long y;

    // mag01[x] = x * MATRIX_A for x=0,1
    static unsigned long mag01[2]={0x0UL, MATRIX_A};

    if (mti >= MERSENNE_N)
    {
        /* generate MERSENNE_N words at one time */
        int kk;

        // If InitMersenne() has not been called, a default initial seed is used
        if (mti == MERSENNE_N+1)
            InitMersenne(5489UL);

        for (kk=0; kk < MERSENNE_N-MERSENNE_M; kk++)
        {
            y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
            mt[kk] = mt[kk+MERSENNE_M] ^(y >> 1) ^ mag01[y & 0x1UL];
        }

        for (; kk < MERSENNE_N-1; kk++)
        {
            y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
            mt[kk] = mt[kk+(MERSENNE_M - MERSENNE_N)] ^(y >> 1) ^ mag01[y & 0x1UL];
        }

        y = (mt[MERSENNE_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[MERSENNE_N-1] = mt[MERSENNE_M-1] ^(y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (mersenneMult *((double) y + 0.5));
}

// Generates a random number on [0,0xffffffff]-interval

unsigned long Random::NextInt()
{
    unsigned long y;

    // mag01[x] = x * MATRIX_A for x=0,1
    static unsigned long mag01[2]={0x0UL, MATRIX_A};

    if (mti >= MERSENNE_N)
    {
        /* generate MERSENNE_N words at one time */
        int kk;

        // If InitMersenne() has not been called, a default initial seed is used
        if (mti == MERSENNE_N + 1)
            InitMersenne(5489UL);

        for (kk= 0; kk < MERSENNE_N - MERSENNE_M; kk++)
        {
            y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
            mt[kk] = mt[kk+MERSENNE_M] ^(y >> 1) ^ mag01[y & 0x1UL];
        }

        for (; kk< MERSENNE_N-1; kk++)
        {
            y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
            mt[kk] = mt[kk+(MERSENNE_M - MERSENNE_N)] ^(y >> 1) ^ mag01[y & 0x1UL];
        }

        y = (mt[MERSENNE_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[MERSENNE_N-1] = mt[MERSENNE_M-1] ^(y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

#else

double Random::Next()
{
    // Compute seed = (IA * seed) % IM without overflows
    // by Schrage's method
    long k = seed / IQ;
    seed = IA * (seed - k * IQ) - IR * k;
    if (seed < 0) seed += IM;

    // Map to 0..NTAB-1
    int j = last/NDIV;

    // Output value is shuffler[j], which is in turn replaced by seed
    last = shuffler[j];
    shuffler[j] = seed;

    // Map to 0.0 .. 1.0 excluding endpoints
    double temp = AM * last;
    if (temp > RNMX) return RNMX;
    return temp;
}

unsigned long Random::NextInt()
{
    // Compute seed = (IA * seed) % IM without overflows
    // by Schrage's method
    long k = seed / IQ;
    seed = IA * (seed - k * IQ) - IR * k;
    if (seed < 0) seed += IM;

    // Map to 0..NTAB-1
    int j = last/NDIV;

    // Output value is shuffler[j], which is in turn replaced by seed
    last = shuffler[j];
    shuffler[j] = seed;

    return last;
}

#endif

double Random::Normal()
{
    double v1, v2, fac, rsq;

    if (!normSaved)   // Do we need new numbers?
    {
        do
        {
            v1 = 2.0 * Next() - 1.0;  // Pick two coordinates from
            v2 = 2.0 * Next() - 1.0;  // -1 to +1 and check if they
            rsq = v1*v1 + v2*v2;  // are in unit circle...
        }
        while (rsq >= 1.0 || rsq == 0.0);

        fac = sqrt(-2.0 * log(rsq)/rsq);  // Apply the Box-Muller
        normStore = v1 * fac;  // transformation and save
        normSaved = 1;  // one deviate for next time
        return v2 * fac;
    }
    else
    {
        normSaved = 0;
        return normStore;
    }
}

void Random::Choose(int * array, int n, int k)
{
    int choices = 1, others = 0;

    if (k > n / 2)
    {
        choices = 0;
        others = 1;
        k = n - k;
    }

    for (int i = 0; i < n; i++)
        array[i] = others;

    while (k > 0)
    {
        int i = NextInt() % n;

        if (array[i] == choices) continue;

        array[i] = choices;
        k--;
    }
}

void Random::Choose(int * array, float * weights, int n, int k)
{
    int choices = 1, others = 0;

    if (k > n / 2)
    {
        choices = 0;
        others = 1;
        k = n - k;
    }

    // First calculate cumulative sums of weights ...
    float * cumulative = new float [n + 1];

    cumulative[0] = 0;
    for (int i = 1; i <= n; i++)
        cumulative[i] = cumulative[i - 1] + weights[i - 1];

    float & sum = cumulative[n], reject = 0.0;

    for (int i = 0; i < n; i++)
        array[i] = others;

    while (k > 0)
    {
        float weight = Next() * sum;

        int hi = n, lo = 0, i = 0;

        while (hi >= lo)
        {
            i = (hi + lo) / 2;

            if (cumulative[i + 1] <= weight)
                lo = i + 1;
            else if (cumulative[i] >= weight)
                hi = i - 1;
            else break;
        }

        if (array[i] == choices) continue;

        array[i] = choices;
        reject += weights[i];

        // After selecting a substantial number of elements, update the cumulative
        // distribution -- to ensure that at least half of our samples produce a hit
        if (reject > sum * 0.50)
        {
            cumulative[0] = 0;
            for (int i = 1; i <= n; i++)
                if (array[i] != choices)
                    cumulative[i] = cumulative[i - 1] + weights[i - 1];
                else
                    cumulative[i] = cumulative[i - 1];

            reject = 0.0;
	    sum = cumulative[n];
        }

        k--;
    }

    delete [] cumulative;
}

Random globalRandom;


