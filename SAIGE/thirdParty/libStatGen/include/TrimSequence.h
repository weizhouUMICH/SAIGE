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

#ifndef _TRIMSEQUENCE_H
#define _TRIMSEQUENCE_H

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>

#ifndef __WIN32__
#include <unistd.h>
#endif

///
/// TrimSequence is a templated function to find bases
/// which are below a certain moving mean threshold,
/// and can be applied to either end of the sequence string.
///
/// @param sequence is the input sequence
/// @param meanValue is the value below which we wish to trim.
/// @return the iterator of the location at which untrimmed values begin
///
/// Details:
///
/// trimFromLeft is a bool indicating which direction we wish
/// to trim.  true -> left to right, false is right to left.
///
/// The code is convoluted enough here, so for implementation
/// and testing sanity, the following definitions are made:
///
/// When trimFromLeft is true:
///    result == sequence.begin() implies no trimming
///    result == sequence.end() implies all values are trimmed
///
/// When trimFromLeft is false:
///    result == sequence.begin() implies all values are trimmed
///    result == sequence.end() no values are trimmed
///
/// result will always be in the range [sequence.begin() , sequence.end())
/// (begin is inclusive, end is exclusive).
///
/// NOTE: See TrimSequence.h and test/TrimSequence_test.cpp for examples
///
/// THIS CODE IS EXCEPTIONALLY FRAGILE.  DO NOT ATTEMPT TO FIX OR
/// IMPROVE WITHOUT INCLUDING DOCUMENTED, UNDERSTANABLE TEST CASES THAT CLEARLY
/// SHOW WHY OR WHY NOT SOMETHING WORKS.
///
template<typename sequenceType, typename meanValueType>
typename sequenceType::iterator trimSequence(sequenceType &sequence, meanValueType meanValue, const bool trimFromLeft)
{
    const int howManyValues = 4;    // this is used in signed arithmetic below
    int windowThreshold = howManyValues * meanValue;
    int64_t sumOfWindow = 0;
    typename sequenceType::iterator it;

    //
    // Sanity check to weed out what otherwise would be
    // a large number of boundary checks below.  If the input
    // is too small, just punt it back to the caller.  Technically,
    // we can still trim, but we'd just do the simple iteration
    // loop.  Save that for when we care.
    //

    if (sequence.size() < (size_t) howManyValues)
        return trimFromLeft? sequence.begin() : sequence.end();

    typename sequenceType::iterator sequenceBegin;
    typename sequenceType::iterator sequenceEnd;

    // The algorithm should be clear and efficient
    // so it does not bother to write codes for two directions.
    // It that way, we avoid thinking trimming from left and right interchangably.
    if (trimFromLeft)
    {
        // sequenceBegin is inclusive, sequenceEnd is exclusive,
        sequenceBegin = sequence.begin();
        sequenceEnd = sequence.end();

        for (it = sequenceBegin; it < sequenceBegin + howManyValues; it++)
            sumOfWindow += *it;

        for (; it < sequenceEnd; it ++)
        {
            if (sumOfWindow > windowThreshold)
                break;
            sumOfWindow += *it;
            sumOfWindow -= *(it - howManyValues);
        }
        // here it is in the range of [sequenceBegin+howManyValues, sequenceEnd] inclusively
        // the range is also [sequence.begin() + howManyValues, sequence.end()]
        while (*(it-1) >= meanValue && (it-1) >= sequenceBegin)
            it--;
    }
    else
    {
        sequenceBegin = sequence.end() - 1;
        sequenceEnd = sequence.begin() - 1;

        for (it = sequenceBegin; it > sequenceBegin - howManyValues; it--)
            sumOfWindow += *it;

        for (; it > sequenceEnd; it--)
        {
            if (sumOfWindow > windowThreshold)
                break;
            sumOfWindow += *it;
            sumOfWindow -= *(it + howManyValues);
        }

        // here it is in the range of [sequenceEnd, sequenceBegin - howManyValues] inclusively
        // the range is also [sequence.begin() -1, sequence.end() - 1 - howManyValues]
        while (*(it+1) >= meanValue && (it+1) <= sequenceBegin)
            it ++;
        // note, the return value should in the range [sequence.begin(), sequence.end()]
        it += 1;
    }

    // 'it' may be sequence.end() in some cases
    assert(it >= sequence.begin() && it <= sequence.end());

    return it;
}

#endif
