/*
 *  Copyright (C) 2012  Regents of the University of Michigan
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

#include "SamRecordHelper.h"
#include <stdexcept>

int SamRecordHelper::checkSequence(SamRecord& record, int32_t pos0Based, 
                                    const char* sequence)
{
    const char* readSeq = record.getSequence();

    // Get the cigar.
    Cigar* cigar = record.getCigarInfo();

    if(cigar == NULL)
    {
        throw std::runtime_error("Failed to get Cigar.");
    }

    int32_t readStartIndex = 
        cigar->getQueryIndex(pos0Based, record.get0BasedPosition());

    // if the read start is negative, this position was deleted, so 
    // return false, it doesn't match.
    if(readStartIndex == Cigar::INDEX_NA)
    {
        return(false);
    }

    // Increment the readSeq start to where this position is found.
    readSeq += readStartIndex;
    if(strncmp(readSeq, sequence, strlen(sequence)) == 0)
    {
        // Match, so return the readStartIndex (cycle).
        return(readStartIndex);
    }
    // Did not match.
    return(-1);
}


bool SamRecordHelper::genSamTagsString(SamRecord& record,
                                       String& returnString,
                                       char delim)
{
    char tag[3];
    char vtype;
    void* value;

    // Reset the tag iterator to ensure that all the tags are written.
    record.resetTagIter();

    // While there are more tags, write them to the recordString.
    bool firstEntry = true;
    bool returnStatus = true;
    while(record.getNextSamTag(tag, vtype, &value) != false)
    {
        if(!firstEntry)
        {
            returnString += delim;
        }
        else
        {
            firstEntry = false;
        }
        returnStatus &= genSamTagString(tag, vtype, value, returnString);
    }
    return(returnStatus);
}


bool SamRecordHelper::genSamTagString(const char* tag, char vtype, 
                                      void* value, String& returnString)
{
    returnString += tag;
    returnString += ":"; 
    returnString += vtype;
    returnString += ":";
    if(SamRecord::isIntegerType(vtype))
    {
        returnString += (int)*(int*)value;
    }
    else if(SamRecord::isFloatType(vtype))
    {
        returnString.appendFullFloat(*(float*)value);
    }
    else if(SamRecord::isCharType(vtype))
    {
        returnString += (char)*(char*)value;
    }
    else if(SamRecord::isStringType(vtype))
    {
        // String type.
        returnString += (String)*(String*)value;
    }
    else
    {
        // Could not determine the type.
        return(false);
    }
    return(true);
}
