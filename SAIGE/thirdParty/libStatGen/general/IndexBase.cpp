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

#include "IndexBase.h"
#include <iomanip>

Chunk SortedChunkList::pop()
{
    Chunk newChunk = chunkList.begin()->second;
    chunkList.erase(chunkList.begin());
    return(newChunk);
}


bool SortedChunkList::insert(const Chunk& chunkToInsert)
{
    std::pair<std::map<uint64_t, Chunk>::iterator, bool> insertRes;
    // Insert the passed in chunk.
    insertRes = 
        chunkList.insert(std::pair<uint64_t, Chunk>(chunkToInsert.chunk_beg,
                                                    chunkToInsert));

    if(!insertRes.second)
    {
        // Failed to insert the chunk.
        std::cerr << "Failed to insert into the SortedChunkList.\n";
        std::cerr << "\tpreviously found chunk:\tbeg = " << std::hex
                  << insertRes.first->second.chunk_beg 
                  << "\tend = "
                  << insertRes.first->second.chunk_end
                  << "\nnew chunk:\tbeg = " << std::hex
                  << chunkToInsert.chunk_beg 
                  << "\tend = "
                  << chunkToInsert.chunk_end
                  << std::endl;
    }
    // return the result that comes from insertRes.
    return(insertRes.second);
}

void SortedChunkList::clear()
{
    chunkList.clear();
}

bool SortedChunkList::empty()
{
    return(chunkList.empty());
}


// Merge overlapping chunks found in this list.
bool SortedChunkList::mergeOverlapping()
{
    // Start at the beginning of the list and iterate through.
    std::map<uint64_t, Chunk>::iterator currentPos = chunkList.begin();
    std::map<uint64_t, Chunk>::iterator nextPos = chunkList.begin();
    if(nextPos != chunkList.end())
    {
        ++nextPos;
    }
    
    // Loop until the end is reached.
    while(nextPos != chunkList.end())
    {
        // If the next chunk is completely contained within the current 
        // chunk (its end is less than the current chunk's end), 
        // delete it since its position is already covered.
        if(nextPos->second.chunk_end < currentPos->second.chunk_end)
        {
            chunkList.erase(nextPos);
            nextPos = currentPos;
            ++nextPos;
            continue;
        }
        
        // If the next chunk's start position's BGZF block is less than or
        // equal to the BGZF block of the current chunk's end position,
        // combine the two chunks into the current chunk.
        if((nextPos->second.chunk_beg >> 16) <= 
           (currentPos->second.chunk_end >> 16))
        {
            currentPos->second.chunk_end = nextPos->second.chunk_end;
            // nextPos has now been included in the current pos, so
            // remove it.
            chunkList.erase(nextPos);
            nextPos = currentPos;
            ++nextPos;
            continue;
        }
        else
        {
            // Nothing to combine.  So try combining at the next 
            currentPos = nextPos;
            ++nextPos;
        }
    }
    return(true);
}


IndexBase::IndexBase()
    : n_ref(0)
{
    myRefs.clear();
}



IndexBase::~IndexBase()
{
}


// Reset the member data for a new index file.
void IndexBase::resetIndex()
{
    n_ref = 0;
    // Clear the references.
    myRefs.clear();
}

    
// Get the number of references in this index.
int32_t IndexBase::getNumRefs() const
{
    // Return the number of references.
    return(myRefs.size());
}


// The basic logic is from samtools reg2bins and the samtools format specification pdf.
// Set bins in the region to 1 and all other bins to 0.
void IndexBase::getBinsForRegion(uint32_t start, uint32_t end, bool binMap[MAX_NUM_BINS+1])
{
    for(uint32_t index = 0; index < MAX_NUM_BINS+1; index++)
    {
        binMap[index] = false;
    }

    uint32_t binNum = 0;
    --end;
    
    // Check if beg/end go too high, set to max position.
    if(start > MAX_POSITION)
    {
        start = MAX_POSITION;
    }
    if(end > MAX_POSITION)
    {
        end = MAX_POSITION;
    }
    
    // Turn on bins.
    binMap[binNum] = true;
    for (binNum =    1 + (start>>26); binNum <=    1 + (end>>26); ++binNum) 
        binMap[binNum] = true;
    for (binNum =    9 + (start>>23); binNum <=    9 + (end>>23); ++binNum) 
        binMap[binNum] = true;
    for (binNum =   73 + (start>>20); binNum <=   73 + (end>>20); ++binNum)
        binMap[binNum] = true;
    for (binNum =  585 + (start>>17); binNum <=  585 + (end>>17); ++binNum)
        binMap[binNum] = true;
    for (binNum = 4681 + (start>>14); binNum <= 4681 + (end>>14); ++binNum)
        binMap[binNum] = true;
}


// Returns the minimum offset of records that cross the 16K block that
// contains the specified position for the given reference id.
bool IndexBase::getMinOffsetFromLinearIndex(int32_t refID, uint32_t position,
                                            uint64_t& minOffset) const
{
    int32_t linearIndex = position >> LINEAR_INDEX_SHIFT;

    minOffset = 0;

    if(refID > n_ref)
    {
        // out of range of the references, return false.
        return(false);
    }
    // Check to see if the position is out of range of the linear index.
    int32_t linearOffsetSize = myRefs[refID].n_intv;

    // If there are no entries in the linear index, return false.
    // Or if the linear index is not large enough to include
    // the start block, then there can be no records that cross
    // our region, so return false.
    if((linearOffsetSize == 0) || (linearIndex >= linearOffsetSize))

    {
        return(false);
    }

    // The linear index is specified for this block, so return that value.
    minOffset = myRefs[refID].ioffsets[linearIndex];
    
    // If the offset is 0, go to the previous block that has an offset.
    // This is due to a couple of bugs in older sam tools indexes.
    // 1) they add one to the index location (so when reading those, you
    // may be starting earlier than necessary)
    // 2) (the bigger issue) They did not include bins 4681-37449 in
    // the linear index.
    while((minOffset == 0) && (--linearIndex >= 0))
    {
        minOffset = myRefs[refID].ioffsets[linearIndex]; 
    }


    // If the minOffset is still 0 when moving forward,
    // check later indices to find a non-zero since we don't want to return
    // an offset of 0 since the record can't start at 0 we want to at least
    // return the first record position for this reference.
    linearIndex = 0;
    while((minOffset == 0) && (linearIndex < linearOffsetSize))
    {
         minOffset = myRefs[refID].ioffsets[linearIndex]; 
         linearIndex++;
    }
    if(minOffset == 0)
    {
        // Could not find a valid start position for this reference.
        return(false);
    }
    return(true);
}
