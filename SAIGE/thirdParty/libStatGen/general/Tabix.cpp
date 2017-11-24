/*
 *  Copyright (C) 2012-2013  Regents of the University of Michigan
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
#include "Tabix.h"
#include <stdexcept>
#include "StringBasics.h"

Tabix::Tabix()
    : IndexBase(),
      myChromNamesBuffer(NULL)
{
}


Tabix::~Tabix()
{
    if(myChromNamesBuffer != NULL)
    {
        delete[] myChromNamesBuffer;
        myChromNamesBuffer = NULL;
    }
}


// Reset the member data for a new index file.
void Tabix::resetIndex()
{
    IndexBase::resetIndex();
    if(myChromNamesBuffer != NULL)
    {
        delete[] myChromNamesBuffer;
        myChromNamesBuffer = NULL;
    }
    myChromNamesVector.clear();
}


// Read & parse the specified index file.
StatGenStatus::Status Tabix::readIndex(const char* filename)
{
    // Reset the index from anything that may previously be set.
    resetIndex();

    IFILE indexFile = ifopen(filename, "rb");

    // Failed to open the index file.
    if(indexFile == NULL)
    {
        return(StatGenStatus::FAIL_IO);
    }

    // read the tabix index structure.

    // Read the magic string.
    char magic[4];
    if(ifread(indexFile, magic, 4) != 4)
    {
        // Failed to read the magic
        return(StatGenStatus::FAIL_IO);
    }

    // If this is not an index file, set num references to 0. 
    if (magic[0] != 'T' || magic[1] != 'B' || magic[2] != 'I' || magic[3] != 1)
    {
        // Not a Tabix Index file.
        return(StatGenStatus::FAIL_PARSE);
    }

    // It is a tabix index file.
    // Read the number of reference sequences.
    if(ifread(indexFile, &n_ref, 4) != 4)
    {
        // Failed to read.
        return(StatGenStatus::FAIL_IO);
    }

    // Size the references.
    myRefs.resize(n_ref);

    // Read the Format configuration.
    if(ifread(indexFile, &myFormat, sizeof(myFormat)) != sizeof(myFormat))
    {
        // Failed to read.
        return(StatGenStatus::FAIL_IO);
    }

    // Read the length of the chromosome names.
    uint32_t l_nm;

    if(ifread(indexFile, &l_nm, sizeof(l_nm)) != sizeof(l_nm))
    {
        // Failed to read.
        return(StatGenStatus::FAIL_IO);
    }

    // Read the chromosome names.
    myChromNamesBuffer = new char[l_nm];
    if(ifread(indexFile, myChromNamesBuffer, l_nm) != l_nm)
    {
        return(StatGenStatus::FAIL_IO);
    }
    myChromNamesVector.resize(n_ref);

    // Parse out the chromosome names.
    bool prevNull = true;
    int chromIndex = 0;
    for(uint32_t i = 0; i < l_nm; i++)
    {
        if(chromIndex >= n_ref)
        {
            // already set the pointer for the last chromosome name, 
            // so stop looping.
            break;
        }
        if(prevNull == true)
        {
            myChromNamesVector[chromIndex++] = myChromNamesBuffer + i;
            prevNull = false;
        }
        if(myChromNamesBuffer[i] == '\0')
        {
            prevNull = true;
        }
    }

    for(int refIndex = 0; refIndex < n_ref; refIndex++)
    {
        // Read each reference.
        Reference* ref = &(myRefs[refIndex]);
        
        // Read the number of bins.
        if(ifread(indexFile, &(ref->n_bin), 4) != 4)
        {
            // Failed to read the number of bins.
            // Return failure.
            return(StatGenStatus::FAIL_PARSE);
        }

        // Resize the bins.
        ref->bins.resize(ref->n_bin + 1);
        
        // Read each bin.
        for(int binIndex = 0; binIndex < ref->n_bin; binIndex++)
        {
            uint32_t binNumber;

            // Read in the bin number.
            if(ifread(indexFile, &(binNumber), 4) != 4)
            {
                // Failed to read the bin number.
                // Return failure.
                return(StatGenStatus::FAIL_IO);
            }

            // Add the bin to the reference and get the
            // pointer back so the values can be set in it.
            Bin* binPtr = &(ref->bins[binIndex]);
            binPtr->bin = binNumber;
         
            // Read in the number of chunks.
            if(ifread(indexFile, &(binPtr->n_chunk), 4) != 4)
            {
                // Failed to read number of chunks.
                // Return failure.
                return(StatGenStatus::FAIL_IO);
            }

            // Read in the chunks.
            // Allocate space for the chunks.
            uint32_t sizeOfChunkList = binPtr->n_chunk * sizeof(Chunk);
            binPtr->chunks = (Chunk*)malloc(sizeOfChunkList);
            if(ifread(indexFile, binPtr->chunks, sizeOfChunkList) != sizeOfChunkList)
            {
                // Failed to read the chunks.
                // Return failure.
                return(StatGenStatus::FAIL_IO);
            }
        }

        // Read the number of intervals.
        if(ifread(indexFile, &(ref->n_intv), 4) != 4)
        {
            // Failed to read, set to 0.
            ref->n_intv = 0;
            // Return failure.
            return(StatGenStatus::FAIL_IO);
        }

        // Allocate space for the intervals and read them.
        uint32_t linearIndexSize = ref->n_intv * sizeof(uint64_t);
        ref->ioffsets = (uint64_t*)malloc(linearIndexSize);
        if(ifread(indexFile, ref->ioffsets, linearIndexSize) != linearIndexSize)
        {
            // Failed to read the linear index.
            // Return failure.
            return(StatGenStatus::FAIL_IO);
        }
    }

    // Successfully read teh bam index file.
    return(StatGenStatus::SUCCESS);
}


bool Tabix::getStartPos(const char* refName, int32_t start,
                        uint64_t& fileStartPos) const
{
    // Look for the reference name in the list.
    int refID = 0;
    for(refID = 0; refID < n_ref; refID++)
    {
        if(strcmp(refName, myChromNamesVector[refID]) == 0)
        {
            // found the reference
            break;
        }
    }
    if(refID >= n_ref)
    {
        // Didn't find the refName, so return false.
        return(false);
    }

    // Look up in the linear index.
    if(start < 0)
    {
        // Negative index, so start at 0.
        start = 0;
    }
    return(getMinOffsetFromLinearIndex(refID, start, fileStartPos));
}


const char* Tabix::getRefName(unsigned int indexNum) const
{
    if(indexNum >= myChromNamesVector.size())
    {
        String message = "ERROR: Out of range on Tabix::getRefName(";
        message += indexNum;
        message += ")";
        throw(std::runtime_error(message.c_str()));
        return(NULL);
    }
    return(myChromNamesVector[indexNum]);
}
