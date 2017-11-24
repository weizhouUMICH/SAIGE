/*
 *  Copyright (C) 2010-2011  Regents of the University of Michigan
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

#include <stdio.h>
#include <stdlib.h>
#include "Cigar.h"
#include "STLUtilities.h"

// Initialize INDEX_NA.
const int32_t Cigar::INDEX_NA = -1;


////////////////////////////////////////////////////////////////////////
//
// Cigar Class
//

//
// Set the passed in string to the string reprentation of the Cigar operations
// in this object.
//
void Cigar::getCigarString(std::string& cigarString) const
{
    using namespace STLUtilities;

    std::vector<CigarOperator>::const_iterator i;

    cigarString.clear();  // clear result string

    // Progressively append the character representations of the operations to
    // the cigar string.
    for (i = cigarOperations.begin(); i != cigarOperations.end(); i++)
    {
        cigarString << (*i).count << (*i).getChar();
    }
}

void Cigar::getCigarString(String& cigarString) const
{
    std::string cigar;

    getCigarString(cigar);

    cigarString = cigar.c_str();

    return;
}

void Cigar::getExpandedString(std::string &s) const
{
    s = "";

    std::vector<CigarOperator>::const_iterator i;

    // Progressively append the character representations of the operations to
    // the string passed in

    for (i = cigarOperations.begin(); i != cigarOperations.end(); i++)
    {
        for (uint32_t j = 0; j<(*i).count; j++) s += (*i).getChar();
    }
    return;
}


bool Cigar::operator == (Cigar &rhs) const
{

    if (this->size() != rhs.size()) return false;

    for (int i = 0; i < this->size(); i++)
    {
        if (cigarOperations[i]!=rhs.cigarOperations[i]) return false;
    }
    return true;
}


// return the length of the read that corresponds to
// the current CIGAR string.
int Cigar::getExpectedQueryBaseCount() const
{
    int matchCount = 0;
    std::vector<CigarOperator>::const_iterator i;
    for (i = cigarOperations.begin(); i != cigarOperations.end(); i++)
    {
        switch (i->operation)
        {
            case match:
            case mismatch:
            case softClip:
            case insert:
                matchCount += i->count;
                break;
            default:
                // we only care about operations that are in the query sequence.
                break;
        }
    }
    return matchCount;
}


// return the number of bases in the reference that
// this read "spans"
int Cigar::getExpectedReferenceBaseCount() const
{
    int matchCount = 0;
    std::vector<CigarOperator>::const_iterator i;
    for (i = cigarOperations.begin(); i != cigarOperations.end(); i++)
    {
        switch (i->operation)
        {
            case match:
            case mismatch:
            case del:
            case skip:
                matchCount += i->count;
                break;
            default:
                // we only care about operations that are in the reference sequence.
                break;
        }
    }
    return matchCount;
}


// Return the number of clips that are at the beginning of the cigar.
int Cigar::getNumBeginClips() const
{
    int numBeginClips = 0;
    for (unsigned int i = 0; i != cigarOperations.size(); i++)
    {
        if ((cigarOperations[i].operation == softClip) ||
                (cigarOperations[i].operation == hardClip))
        {
            // Clipping operator, increment the counter.
            numBeginClips += cigarOperations[i].count;
        }
        else
        {
            // Break out of the loop since a non-clipping operator was found.
            break;
        }
    }
    return(numBeginClips);
}


// Return the number of clips that are at the end of the cigar.
int Cigar::getNumEndClips() const
{
    int numEndClips = 0;
    for (int i = (cigarOperations.size() - 1); i >= 0; i--)
    {
        if ((cigarOperations[i].operation == softClip) ||
                (cigarOperations[i].operation == hardClip))
        {
            // Clipping operator, increment the counter.
            numEndClips += cigarOperations[i].count;
        }
        else
        {
            // Break out of the loop since a non-clipping operator was found.
            break;
        }
    }
    return(numEndClips);
}


int32_t Cigar::getRefOffset(int32_t queryIndex)
{
    // If the vectors aren't set, set them.
    if ((queryToRef.size() == 0) || (refToQuery.size() == 0))
    {
        setQueryAndReferenceIndexes();
    }
    if ((queryIndex < 0) || ((uint32_t)queryIndex >= queryToRef.size()))
    {
        return(INDEX_NA);
    }
    return(queryToRef[queryIndex]);
}


int32_t Cigar::getQueryIndex(int32_t refOffset)
{
    // If the vectors aren't set, set them.
    if ((queryToRef.size() == 0) || (refToQuery.size() == 0))
    {
        setQueryAndReferenceIndexes();
    }
    if ((refOffset < 0) || ((uint32_t)refOffset >= refToQuery.size()))
    {
        return(INDEX_NA);
    }
    return(refToQuery[refOffset]);
}


int32_t Cigar::getRefPosition(int32_t queryIndex, int32_t queryStartPos)
{
    // If the vectors aren't set, set them.
    if ((queryToRef.size() == 0) || (refToQuery.size() == 0))
    {
        setQueryAndReferenceIndexes();
    }
    if ((queryIndex < 0) || ((uint32_t)queryIndex >= queryToRef.size()))
    {
        return(INDEX_NA);
    }

    if (queryToRef[queryIndex] != INDEX_NA)
    {
        return(queryToRef[queryIndex] + queryStartPos);
    }
    return(INDEX_NA);
}


// Return the query index associated with the specified reference position
// when the query starts at the specified reference position based on
// this cigar.
int32_t Cigar::getQueryIndex(int32_t refPosition, int32_t queryStartPos)
{
    // If the vectors aren't set, set them.
    if ((queryToRef.size() == 0) || (refToQuery.size() == 0))
    {
        setQueryAndReferenceIndexes();
    }

    int32_t refOffset = refPosition - queryStartPos;
    if ((refOffset < 0) || ((uint32_t)refOffset >= refToQuery.size()))
    {
        return(INDEX_NA);
    }

    return(refToQuery[refOffset]);
}


int32_t Cigar::getExpandedCigarIndexFromQueryIndex(int32_t queryIndex)
{
    // If the vectors aren't set, set them.
    if ((queryToRef.size() == 0) || (refToQuery.size() == 0))
    {
        setQueryAndReferenceIndexes();
    }
    if ((queryIndex < 0) || ((uint32_t)queryIndex >= queryToCigar.size()))
    {
        return(INDEX_NA);
    }
    return(queryToCigar[queryIndex]);
}


int32_t Cigar::getExpandedCigarIndexFromRefOffset(int32_t refOffset)
{
    // If the vectors aren't set, set them.
    if ((queryToRef.size() == 0) || (refToQuery.size() == 0))
    {
        setQueryAndReferenceIndexes();
    }
    if ((refOffset < 0) || ((uint32_t)refOffset >= refToCigar.size()))
    {
        return(INDEX_NA);
    }
    return(refToCigar[refOffset]);
}


int32_t Cigar::getExpandedCigarIndexFromRefPos(int32_t refPosition, 
                                               int32_t queryStartPos)
{
    return(getExpandedCigarIndexFromRefOffset(refPosition - queryStartPos));
}


char Cigar::getCigarCharOp(int32_t expandedCigarIndex)
{
    // Check if the expanded cigar has been set yet
    if ((queryToRef.size() == 0) || (refToQuery.size() == 0))
    {
        // Set the expanded cigar.
        setQueryAndReferenceIndexes();
    }

    // Check to see if the index is in range.
    if((expandedCigarIndex < 0) || 
       ((uint32_t)expandedCigarIndex >= myExpandedCigar.length()))
    {
        return('?');
    }
    return(myExpandedCigar[expandedCigarIndex]);
}


char Cigar::getCigarCharOpFromQueryIndex(int32_t queryIndex)
{
    return(getCigarCharOp(getExpandedCigarIndexFromQueryIndex(queryIndex)));
}


char Cigar::getCigarCharOpFromRefOffset(int32_t refOffset)
{
    return(getCigarCharOp(getExpandedCigarIndexFromRefOffset(refOffset)));
}


char Cigar::getCigarCharOpFromRefPos(int32_t refPosition, int32_t queryStartPos)
{
    return(getCigarCharOp(getExpandedCigarIndexFromRefPos(refPosition, queryStartPos)));
}


// Return the number of bases that overlap the reference and the
// read associated with this cigar that falls within the specified region.
uint32_t Cigar::getNumOverlaps(int32_t start, int32_t end,
                               int32_t queryStartPos)
{
    // Get the overlap info.
    if ((queryToRef.size() == 0) || (refToQuery.size() == 0))
    {
        setQueryAndReferenceIndexes();
    }

    // Get the start and end offsets.
    int32_t startRefOffset = 0;
    // If the specified start is more than the queryStartPos, set
    // the startRefOffset to the appropriate non-zero value.
    // (if start is <= queryStartPos, than startRefOffset is 0 - it should
    // not be set to a negative value.)
    if (start > queryStartPos)
    {
        startRefOffset = start - queryStartPos;
    }

    int32_t endRefOffset = end - queryStartPos;
    if (end  == -1)
    {
        // -1 means that the region goes to the end of the refrerence.
        // So set endRefOffset to the max refOffset + 1 which is the
        // size of the refToQuery vector.
        endRefOffset = refToQuery.size();
    }


    // if endRefOffset is less than 0, then this read does not fall within
    // the specified region, so return 0.
    if (endRefOffset < 0)
    {
        return(0);
    }

    // Get the overlaps for these offsets.
    // Loop through the read counting positions that match the reference
    // within this region.
    int32_t refOffset = 0;
    int32_t numOverlaps = 0;
    for (unsigned int queryIndex = 0; queryIndex < queryToRef.size();
            queryIndex++)
    {
        refOffset = getRefOffset(queryIndex);
        if (refOffset > endRefOffset)
        {
            // Past the end of the specified region, so stop checking
            // for overlaps since there will be no more.
            break;
        }
        else if ((refOffset >= startRefOffset) && (refOffset < endRefOffset))
        {
            // within the region, increment the counter.
            ++numOverlaps;
        }
    }

    return(numOverlaps);
}


// Return whether or not the cigar has an indel
bool Cigar::hasIndel()
{
    for(unsigned int i = 0; i < cigarOperations.size(); i++)
    {
        if((cigarOperations[i].operation == insert) ||
           (cigarOperations[i].operation == del))
        {
            // Found an indel, so return true.
            return(true);
        }
    }
    // Went through all the operations, and found no indel, so return false.
    return(false);
}


// Clear the query index/reference offset index vectors.
void Cigar::clearQueryAndReferenceIndexes()
{
    queryToRef.clear();
    refToQuery.clear();
    refToCigar.clear();
    queryToCigar.clear();
    myExpandedCigar.clear();
}


///////////////////////////////////////////////////////
// Set the query index/reference offset index vectors.
//
// For Cigar: 3M2I2M1D1M
// That total count of cigar elements is 9 (3+2+2+1+1)
//
// The entries that are valid in the query/reference contain the index/offset
// where they are found in the query/reference.  N/A are marked by 'x':
// query indexes:     0123456x7
//                    ---------
// reference offsets: 012xx3456
//
// This shows what query index is associated with which reference offset and
// vice versa.
// For ones where an x appears, -1 would be returned.
//
void Cigar::setQueryAndReferenceIndexes()
{
    // First ensure that the vectors are clear by clearing them.
    clearQueryAndReferenceIndexes();

    int extPos = 0;

    // Process each cigar index.
    for (uint32_t cigarIndex = 0; cigarIndex < cigarOperations.size(); cigarIndex++)
    {
        // Process the cigar operation.
        switch (cigarOperations[cigarIndex].operation)
        {
            case match:
            case mismatch:
                // For match/mismatch, update the maps between query
                // and reference for the number of matches/mismatches.
                for (uint32_t i = 0; i < cigarOperations[cigarIndex].count; i++)
                {
                    // The associated indexes are the next location in
                    // each array, which is equal to the current size.
                    int32_t queryToRefLen = queryToRef.size();
                    int32_t refToQueryLen = refToQuery.size();
                    queryToRef.push_back(refToQueryLen);
                    refToQuery.push_back(queryToRefLen);
                    refToCigar.push_back(extPos);
                    queryToCigar.push_back(extPos++);
                    myExpandedCigar.push_back(cigarOperations[cigarIndex].getChar());
                }
                break;
            case insert:
            case softClip:
                // Add N/A reference offset for each query index that this
                // insert covers.
                for (uint32_t i = 0; i < cigarOperations[cigarIndex].count; i++)
                {
                    queryToRef.push_back(INDEX_NA);
                    queryToCigar.push_back(extPos++);
                    myExpandedCigar.push_back(cigarOperations[cigarIndex].getChar());
                }
                break;
            case del:
            case skip:
                // Add N/A query index for each reference offset that this
                // deletion/skip covers.
                for (uint32_t i = 0; i < cigarOperations[cigarIndex].count; i++)
                {
                    refToQuery.push_back(INDEX_NA);
                    refToCigar.push_back(extPos++);
                    myExpandedCigar.push_back(cigarOperations[cigarIndex].getChar());
                }
                break;
            case hardClip:
            case pad:
            case none:
                for (uint32_t i = 0; i < cigarOperations[cigarIndex].count; i++)
                {
                    myExpandedCigar.push_back(cigarOperations[cigarIndex].getChar());
                    ++extPos;
                }
                break;
        };
    }
}

