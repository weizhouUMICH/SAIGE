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

#ifndef _GREEDY_TUPLE_H
#define _GREEDY_TUPLE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "Generic.h"
#include "CigarRoller.h"

/*
 *

 TODO:
 1. how to efficiently find insertion?

*/

/**
 * Weight includes various penalties(e.g. gap open) used in local alignment
 */
struct Weight
{
public:
    Weight()
    {
        gapOpen = gapExtend = -1; // here we do not use affine gap penalty for simlicity.
        mismatch = -1;
        match= 2;
    };
    int gapOpen;
    int gapExtend;
    int mismatch;
    int match;
};

//
// tuple number is 3, arbitrary number from my guess!
// another reason
//
template <typename QueryType, typename ReferenceType, typename ReferenceIndex>
class GreedyTupleAligner
{
public:
    GreedyTupleAligner(Weight& wt): weight(wt)
    {/* */}

    /**
     * Match 'query' to the 'reference' from 'searchStartIndex' up
     * to 'searchSize', store matched length to 'matchedLength'
     * and number of mismatch to 'mismatch'
     * @param query            input query
     * @param queryLength       length of query
     * @param reference        reference sequence
     * @param searchStartIndex the positino where search starts
     * @param searchSize       the total length in reference sequence that will be examine
     * @param matchedLength    store how many bases are matched
     * @param mismatch         store how many bases are mismatched
     * @return -1 for unsuccess return
     */
    int MatchTuple(
        const QueryType       query,
        const int             queryLength,
        const ReferenceType   reference,
        const ReferenceIndex  searchStartIndex,
        const int             searchSize,
        int&            matchedLength,
        int&            mismatch)
    {
        // now use naive search,
        // TODO: will incorportate KMP serach later
        // TODO: adjust tolerance of mismatches
        const int MAX_MISMATCH=2;
        int bestPos = 0, bestMismatch = queryLength, bestMatchedLength = 0, bestScore=-1;

#if defined(DEBUG_GREEDY_ALIGNER)
        cout << "searchStartIndex == " << searchStartIndex << ", searchSize == " << searchSize << std::endl;
#endif
        // here i is the matching position (inclusive)
        // j is the matched length
        for (int i = 0; i <= searchSize - tupleSize; i++)
        {
            int j = 0;
            mismatch = 0;
            while (j < queryLength)
            {
                if (searchStartIndex + i + j >= reference.getNumberBases())
                    break;
                if (query[j] != reference[searchStartIndex + i + j])
                {
                    mismatch++;
                    if (mismatch >= MAX_MISMATCH)
                        break;
                }
                j++;
            }

            if (j>0 && (j==queryLength)) j--;

            while (searchStartIndex +i +j < reference.getNumberBases()
                    && ((j+1) > mismatch)
                    && mismatch>0
                    && query[j] != reference[searchStartIndex + i+j])
            {
                // if pattern matching goes beyong the preset mismatch cutoff,
                // we will have to go backwards
                j--;
                mismatch--;
            }

            int score = j - mismatch;

            if (score > bestScore)
            {
                bestPos = i;
                bestScore = score;
                bestMismatch = mismatch;
                bestMatchedLength = j+1;
            }
        }

        if (bestScore > 0)
        {
            mismatch = bestMismatch;
            matchedLength = bestMatchedLength;
            return bestPos;
        }
        return -1;
    }

    /**
     * Core local alignment algorithm
     * @param query             input query
     * @param queryLength       length of query
     * @param reference         reference genome
     * @param searchStartIndex  matching starts here
     * @param searchSize        how far we will search
     * @param cigarRoller       store alignment results here
     * @param matchPosition     store match position
     */
    void Align(
        QueryType       query,
        int             queryLength,
        ReferenceType   reference,
        ReferenceIndex  searchStartIndex,
        int             searchSize,
        CigarRoller&    cigarRoller,
        ReferenceIndex& matchPosition)
    {
        // Algorithm:
        // finished align? (should we try different align position?)
        // if not, try next tuple
        //    is the tuple aligned?
        //    yes, extend to previous, mark unmatched part mismatch or gap
        //         extend to next matched part
        int r1 = 0; // a start index: reference starting from r1 (inclusive) will be used
        int queryMatchCount = 0; // query matched # of bases
        int q1 = 0; // to align
        int pos = -1;
        int lastR1 = -1; // index: record last

        cigarRoller.clear();
        matchPosition = -1;

        while (queryMatchCount < queryLength)
        {
            if (r1 == searchSize - 1)   // touched ref right boundary
            {
                cigarRoller.Add(CigarRoller::softClip, queryLength-queryMatchCount);
                break;
            }
            if (queryLength - q1 < tupleSize)
            {
                // XXX this needs to do something more sane
                // printf("some bases left!\n");
                // a simple fix: treat all left-over bases as mismatches/matches
                cigarRoller.Add(CigarRoller::mismatch, queryLength - queryMatchCount);
                break;
            }
            int mismatch = 0;
            int matchedLen = 0;
            if ((pos = MatchTuple(query+q1, queryLength-q1, reference, searchStartIndex + r1, searchSize - r1, matchedLen, mismatch)) // found match position for tuple

                    >= 0)
            {
                // found match position for tuple

                if (lastR1<0)
                    matchPosition = pos;

                //
                // deal with left
                //
                if (lastR1>=0)   // have previously aligned part of the query to the reference genome yet
                {
                    if (pos > 0)
                    {
                        cigarRoller.Add(CigarRoller::del, pos);
                    }
                }
                else
                {
                    lastR1 = pos;
                }

                r1 += pos;
                r1 += matchedLen;
                q1 += matchedLen;

                //
                // deal with right
                //
                cigarRoller.Add(CigarRoller::match, matchedLen);
                queryMatchCount = q1;
                lastR1 = r1;

                continue;
            } // end if

            //
            // try insertion
            // maximum insert ? say 2
            //
            for (int i = 1; i < queryLength - q1 - tupleSize; i++)
            {
                int mismatch = 0;
                int matchedLen = 0;
                // check reference genome broundary
                if (searchStartIndex + r1 >= reference.getNumberBases())
                    return;
                if ((pos = MatchTuple(query+q1 + i ,
                                      queryLength - q1 -i ,
                                      reference,
                                      searchStartIndex + r1,
                                      searchSize - r1,
                                      matchedLen,
                                      mismatch)) // found match position for tuple
                        >= 0)
                {
                    if (matchPosition < 0)
                        matchPosition = pos + q1 + i ;
                }
                queryMatchCount += i;
                q1 += i;
                cigarRoller.Add(CigarRoller::insert, i);

                lastR1 = r1 + pos;
                r1 += pos + tupleSize;
                q1 += tupleSize;

                // deal with right
                while (searchStartIndex + r1 < reference.getNumberBases()
                        && query[q1]==reference[searchStartIndex + r1]
                        && q1 < queryLength)
                {
                    r1++;
                    q1++;
                }
                if (q1 < queryLength)
                {
                    cigarRoller.Add(CigarRoller::match, q1 - queryMatchCount);
                    queryMatchCount = q1;
                }
                else
                {
                    cigarRoller.Add(CigarRoller::match, queryLength - queryMatchCount);
                    queryMatchCount = queryLength ;
                    break ;
                }
            }

            r1++;
            q1++;

            // try next
        } // end while (queryMatchCount < queryLength)
    }
private:
    static const int tupleSize = 3;
    Weight weight;
};


#endif
