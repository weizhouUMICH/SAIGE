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

#ifndef __SAM_COORD_OUTPUT_H__
#define __SAM_COORD_OUTPUT_H__

#include "SamFile.h"
#include "SamRecordPool.h"

/// Class for buffering up output reads to ensure that it is sorted. 
/// They are added in almost sorted order.
/// Flush writes any records that start at/before the specified position.
class SamCoordOutput
{
public:
    /// Create an output buffer returning any written records to the specified pool.
    /// \param pool pool that any written records should be returned to, a pointer
    /// to this pool is stored, so it should not go out of scope until the output buffer
    /// has written all the records.
    SamCoordOutput(SamRecordPool& pool);
    ~SamCoordOutput();

    /// Set the already opened output file to write to when flushed.
    /// The user should not close/delete the SamFile until this class is done
    /// with it.  This class does NOT close/delete the SamFile.
    /// \param outFile pointer to an already opened (and header written) 
    /// SAM/BAM output file.
    /// \param header pointer to an already written header that should be
    /// used for writing the records.
    void setOutputFile(SamFile* outFile, SamFileHeader* header);

    /// Add the specified record to this read buffer.
    bool add(SamRecord* record);

    /// Flush the entire buffer, writing all records.
    /// If no output buffer is set, the files cannot be written, but the
    /// flushed records are removed/freed.
    bool flushAll();

    /// Flush the buffer based on the specified chromosome id/position, writing
    /// any records that start at/before the specified chromosome id/position.
    /// If no output buffer is set, the files cannot be written, but the
    /// flushed records are removed/freed.
    /// A chromID of -1 will flush everything regardless of pos0Based.
    bool flush(int32_t chromID, int32_t pos0Based);

protected:
    

private:
    // Require a sam record pool, so make the constructor with
    // no parameters private.
    SamCoordOutput();

    SamFile* myOutputFile;
    SamFileHeader* myHeader;
    std::multimap<uint64_t, SamRecord*> myReadBuffer;
    SamRecordPool* myPool;
};


#endif
