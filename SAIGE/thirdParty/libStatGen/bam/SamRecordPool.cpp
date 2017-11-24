/*
 *  Copyright (C) 2011  Regents of the University of Michigan
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

#include <stdexcept>
#include "SamRecordPool.h"

SamRecordPool::SamRecordPool()
    : myFreeSamRecords(),
      myMaxAllowedRecs(-1),
      myAllocatedRecs(0)
{
}


SamRecordPool::SamRecordPool(int maxNumRecs)
    : myFreeSamRecords(),
      myMaxAllowedRecs(maxNumRecs),
      myAllocatedRecs(0)
{
}


SamRecordPool::~SamRecordPool()
{
    // Loop through the stack deleting the free records.
    while (!myFreeSamRecords.empty())
    {
        delete(myFreeSamRecords.front());
        myFreeSamRecords.pop();
    }
}


SamRecord* SamRecordPool::getRecord()
{
    // Get new samRecord.
    SamRecord* returnSam = NULL;
    if(!myFreeSamRecords.empty())
    {
        // have free already allocated records, so get one of those.
        returnSam = myFreeSamRecords.front();
        myFreeSamRecords.pop();
    }
    else if((myMaxAllowedRecs == -1) || (myAllocatedRecs < myMaxAllowedRecs))
    {
        // There were no free records, but either there is no max or
        // there is still room to allocate more.
        returnSam = new SamRecord();
        ++myAllocatedRecs;
        if(returnSam == NULL)
        {
            // Failed allocation.
            throw(std::runtime_error("Failed to allocate SamRecord"));
        }
    }
    else
    {
        // There are no more free ones and we have already hit the
        // max number allowed to be allocated, so return NULL.
        // The user will have to release some or update the max.
        returnSam = NULL;
    }
    return(returnSam);
}


void SamRecordPool::releaseRecord(SamRecord* record)
{
    if(record == NULL)
    {
        // Nothing to release, so just return.
        return;
    }

    // Release the samRecord to be reused.
    myFreeSamRecords.push(record);
}


void SamRecordPool::setMaxAllocatedRecs(int maxNumRecs)
{
    myMaxAllowedRecs = maxNumRecs;
}
