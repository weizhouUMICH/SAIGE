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

#include "assert.h"
#include "ctype.h"
#include "stdio.h"
#include "Error.h"


#include "Generic.h"
#include "ReferenceSequence.h"

#include <algorithm>
#include <istream>
#include <fstream>
#include <sstream>
#include <stdexcept>

//
// Given a buffer with a fasta format contents, count
// the number of chromsomes in it and return that value.
//
bool getFastaStats(const char *fastaData, size_t fastaDataSize, uint32_t &chromosomeCount, uint64_t &baseCount)
{
    chromosomeCount = 0;
    baseCount = 0;
    bool atLineStart = true;

    //
    // loop over the fasta file, essentially matching for the
    // pattern '^>.*$' and counting them.
    //
    for (size_t fastaIndex = 0; fastaIndex < fastaDataSize; fastaIndex++)
    {
        switch (fastaData[fastaIndex])
        {
            case '\n':
            case '\r':
                atLineStart = true;
                break;
            case '>':
            {
                if (!atLineStart) break;
                chromosomeCount++;
                //
                // eat the rest of the line
                //
                while (fastaIndex < fastaDataSize && fastaData[fastaIndex]!='\n' && fastaData[fastaIndex]!='\r')
                {
                    fastaIndex++;
                }
                break;
            }
            default:
                baseCount++;
                atLineStart = false;
                break;
        }

    }
    return false;
}

#if 0
// turn this into a template on read/quality/etc...
int GenomeSequence::debugPrintReadValidation(
    std::string &read,
    std::string &quality,
    char   direction,
    genomeIndex_t   readLocation,
    int sumQuality,
    int mismatchCount,
    bool recurse
) 
{
    int validateSumQ = 0;
    int validateMismatchCount = 0;
    int rc = 0;
    std::string genomeData;

    for (uint32_t i=0; i<read.size(); i++)
    {
        if (tolower(read[i]) != tolower((*this)[readLocation + i]))
        {
            validateSumQ += quality[i] - '!';
            // XXX no longer valid:
            if (direction=='F' ? i<24 : (i >= (read.size() - 24))) validateMismatchCount++;
            genomeData.push_back(tolower((*this)[readLocation + i]));
        }
        else
        {
            genomeData.push_back(toupper((*this)[readLocation + i]));
        }
    }
    assert(validateSumQ>=0);
    if (validateSumQ != sumQuality && validateMismatchCount == mismatchCount)
    {
        printf("SUMQ: Original Genome: %s  test read: %s : actual sumQ = %d, test sumQ = %d\n",
               genomeData.c_str(),
               read.c_str(),
               validateSumQ,
               sumQuality
              );
        rc++;
    }
    else if (validateSumQ == sumQuality && validateMismatchCount != mismatchCount)
    {
        printf("MISM: Original Genome: %s  test read: %s : actual mismatch %d test mismatches %d\n",
               genomeData.c_str(),
               read.c_str(),
               validateMismatchCount,
               mismatchCount
              );
        rc++;
    }
    else if (validateSumQ != sumQuality && validateMismatchCount != mismatchCount)
    {
        printf("BOTH: Original Genome: %s  test read: %s : actual sumQ = %d, test sumQ = %d, actual mismatch %d test mismatches %d\n",
               genomeData.c_str(),
               read.c_str(),
               validateSumQ,
               sumQuality,
               validateMismatchCount,
               mismatchCount
              );
        rc++;
    }

    if (recurse && abs(validateMismatchCount - mismatchCount) > (int) read.size()/2)
    {
        printf("large mismatch difference, trying reverse strand: ");
        std::string reverseRead = read;
        std::string reverseQuality = quality;
        getReverseRead(reverseRead);
        reverse(reverseQuality.begin(), reverseQuality.end());
        rc = debugPrintReadValidation(reverseRead, reverseQuality, readLocation, sumQuality, mismatchCount, false);
    }
    return rc;
}
#endif



