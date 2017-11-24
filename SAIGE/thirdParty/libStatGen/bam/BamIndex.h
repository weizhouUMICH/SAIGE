/*
 *  Copyright (C) 2010-2012  Regents of the University of Michigan
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

#ifndef __BAM_INDEX_H__
#define __BAM_INDEX_H__

#include <stdint.h>
#include <vector>
#include <map>
#include <stdlib.h>

#include "IndexBase.h"

#include "InputFile.h"
#include "SamStatus.h"

class BamIndex : public IndexBase
{
public:

    BamIndex();
    virtual ~BamIndex();

    /// Reset the member data for a new index file.
    virtual void resetIndex();

    // Read & parse the specified index file.
    /// \param filename the bam index file to be read.
    /// \return the status of the read.
    SamStatus::Status readIndex(const char* filename);

    /// Get the list of chunks associated with this region.
    /// For an entire reference ID, set start and end to -1.
    /// To start at the beginning of the region, set start to 0/-1.
    /// To go to the end of the region, set end to -1.
    bool getChunksForRegion(int32_t refID, int32_t start, int32_t end, 
                            SortedChunkList& chunkList);

    uint64_t getMaxOffset() const;

    /// Get the minimum and maximum file offsets for the specfied reference ID.
    /// \param refID the reference ID to locate in the file.
    /// \param minOffset returns the min file offset for the specified reference
    /// \param maxOffset returns the max file offset for the specified reference
    /// \return whether or not the reference was found in the file
    bool getReferenceMinMax(int32_t refID, 
                            uint64_t& minOffset, 
                            uint64_t& maxOffset) const;

    /// Get the number of mapped reads for this reference id.  Returns -1 for
    /// out of range refIDs.
    /// \param refID reference ID for which to extract the number of mapped reads.
    /// \return number of mapped reads for the specified reference id.
    int32_t getNumMappedReads(int32_t refID);

    /// Get the number of unmapped reads for this reference id.  Returns -1 for
    /// out of range refIDs.
    /// \param refID reference ID for which to extract the number of unmapped reads.
    /// \return number of unmapped reads for the specified reference id
    int32_t getNumUnMappedReads(int32_t refID);

    /// Print the index information.
    /// \param refID reference ID for which to print info for.  -1 means print for all references.
    /// \param summary whether or not to just print a summary (defaults to false).  The summary just contains summary info for each reference and not every bin/chunk.
    void printIndex(int32_t refID, bool summary = false);

    // Number of reference sequences.
    /// The number used for an unknown number of reads.
    static const int32_t UNKNOWN_NUM_READS = -1;

    /// The number used for the reference id of unmapped reads.
    static const int32_t REF_ID_UNMAPPED = -1;

    /// The number used to indicate that all reference ids should be used.
    static const int32_t REF_ID_ALL = -2;

private:
    uint64_t maxOverallOffset;

    int32_t myUnMappedNumReads;
};


#endif
