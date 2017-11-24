/*
 *  Copyright (C) 2011-2012  Regents of the University of Michigan
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

#ifndef __INDEX_BASE_H__
#define __INDEX_BASE_H__

#include <stdint.h>
#include <vector>
#include <map>
#include <stdlib.h>

#include "InputFile.h"
#include "StatGenStatus.h"


class Chunk
{
public:
    uint64_t chunk_beg; // offset of the start of the chunk
    uint64_t chunk_end; // offset of the end of the chunk
    
    static const uint64_t MAX_CHUNK_VALUE = 0xFFFFFFFFFFFFFFFFULL;

    bool operator< (const Chunk& otherChunk) const
    {
        return(this->chunk_beg < otherChunk.chunk_beg);
    }
};


// This class contains chunks that are sorted by the beginning position.
// This class hides how the chunks are actually stored (map, list ,etc),
// so they can be interchanged.
class SortedChunkList
{
public:
    // Returns the first chunk in the list and  removes it.
    Chunk pop();
    bool insert(const Chunk& chunkToInsert);
    void clear();
    bool empty();
    bool mergeOverlapping();

private:
    std::map<uint64_t, Chunk> chunkList;
};

class IndexBase
{
public:

    IndexBase();
    virtual ~IndexBase();

    /// Reset the member data for a new index file.
    virtual void resetIndex();

    // Read & parse the specified index file.
    /// \param filename the bam index file to be read.
    /// \return the status of the read.
    virtual StatGenStatus::Status readIndex(const char* filename) = 0;

    /// Get the number of references in this index.
    /// \return number of references
    int32_t getNumRefs() const;

    // Returns the minimum offset of records that cross the 16K block that
    // contains the specified position for the given reference id.
    bool getMinOffsetFromLinearIndex(int32_t refID, uint32_t position,
                                     uint64_t& minOffset) const;

protected:
    const static uint32_t MAX_NUM_BINS = 37450; // per specs, at most 37450 bins

    // Maximum allowed position (inclusive 512MB - 1)
    // NOTE: CSI index may not have this same max position.
    const static uint32_t MAX_POSITION = 536870911;

    // Number of bits in 1 linear index - how much to shift a position by
    // to determine which offset into the linear index to look for it.
    const static uint32_t LINEAR_INDEX_SHIFT = 14;

    class Bin
    {
    public:
        Bin(){chunks = NULL; reset();}
        ~Bin() {reset();}
        void reset()
        {
            if(chunks != NULL)
            {
                free(chunks);
                chunks = NULL;
            }
            n_chunk = 0; 
            bin = NOT_USED_BIN;
        }
        uint32_t bin; // The bin id.
        int32_t n_chunk; // The number of chunks.
        Chunk* chunks; // The chunks for this bin.
        static const uint32_t NOT_USED_BIN = 0xFFFFFFFF;
    };

    class Reference
    {
        // Add one to the max since there may now be an extra bin containing
        // the mapped/unmapped counts.
    public:
        static const int32_t UNKNOWN_MAP_INFO = -1;
        Reference(){ioffsets = NULL; reset();}
        ~Reference(){reset();}
        void reset()
        { 
            bins.clear(); 
            if(ioffsets != NULL)
            {
                free(ioffsets);
                ioffsets = NULL;
            }
            n_bin = 0; 
            n_intv = 0;
            minChunkOffset = UNSET_MIN_CHUNK_OFFSET;
            maxChunkOffset = 0;
            n_mapped = UNKNOWN_MAP_INFO;
            n_unmapped = UNKNOWN_MAP_INFO;
        }
        int32_t n_bin; // The number of bins.
        int32_t n_intv; // Number of intervals.
        std::vector<Bin> bins;  // The bins for this reference.
        uint64_t* ioffsets; // Offsets of intervals first alignments
        uint64_t minChunkOffset;
        uint64_t maxChunkOffset;
        int32_t n_mapped; // Number of mapped reads.
        int32_t n_unmapped; // Number of unmapped reads.

        static const uint64_t UNSET_MIN_CHUNK_OFFSET = 0xFFFFFFFFFFFFFFFFULL;
    };

    // Set bins in the region to 1 and all other bins to 0.
    // start is incluive, end is exclusive.
    static void getBinsForRegion(uint32_t start, uint32_t end, bool binMap[MAX_NUM_BINS+1]);

    // Number of reference sequences.
    int32_t n_ref;

    // The references.
    std::vector<Reference> myRefs;
};


#endif
