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

#ifndef __SAM_RECORD_POOL_H__
#define __SAM_RECORD_POOL_H__

#include <queue>
#include "SamRecord.h"


class SamRecordPool
{
public:
    /// Constructor that sets there to be no max number of allocated records.
    SamRecordPool();

    /// Constructor that sets the maximum number of allocated records
    /// \param maxNumRecs maximum number of allocated records (-1 means no max)
    SamRecordPool(int maxNumRecs);

    /// Destructor.  Any records that were allocated without calling "releaseRecord"
    /// will not get cleaned up and the user will need to delete them.
    ~SamRecordPool();

    /// Get a SamRecord.  If records are already allocated and free use those, if not
    /// and there are still more that are allowed to be allocated, allocate a new one.
    /// If no more records are allowed to be allocated, NULL is returned.
    /// NOTE: The user should call releaseRecord when done using the record.
    /// If the user deletes the record instead, it still counts as allocated when
    /// comparing against the maxNumRecs but cannot be reused.
    /// \return pointer to a SamRecord available for use, or NULL if no more records
    /// are allowed to be allocated.
    SamRecord* getRecord();

    /// If record is not NULL, adds it back to the free list.  
    /// If record is NULL, nothing is done.
    /// \param record pointer to a record that is no longer being used 
    /// and is available for reuse.
    void releaseRecord(SamRecord* record);

    /// Set the maximum number of records allowed to be allocated.
    /// If more than the new value have already been allocated,
    /// it does not deallocate any, and will continue to reuse
    /// the already allocated records, but it will not allocate
    /// any additional records.
    void setMaxAllocatedRecs(int maxNumRecs);
    

private:

    std::queue<SamRecord*> myFreeSamRecords;
    int myMaxAllowedRecs;
    int myAllocatedRecs;
};

#endif
