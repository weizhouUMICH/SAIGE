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

#include <stdexcept>

#include "PileupElementBaseQual.h"

PileupElementBaseQual::PileupElementBaseQual()
    : PileupElement(),
      myBases(NULL),
      myQualities(NULL),
      myAllocatedSize(0),
      myIndex(-1),
      myAddDelAsBase(false)
{
    myAllocatedSize = 1024;
    myBases = (char*)malloc(myAllocatedSize + 1);
    myQualities = (char*)malloc(myAllocatedSize + 1);
    if((myBases == NULL ) || (myQualities == NULL))
    {
        // TODO, check for malloc failures.
        std::cerr << "Failed Memory Allocation\n";
    }
}

// NOTE that this method does not actually copy, it just resets.
PileupElementBaseQual::PileupElementBaseQual(const PileupElementBaseQual& q)
    : PileupElement(),
      myBases(NULL),
      myQualities(NULL),
      myAllocatedSize(0),
      myIndex(-1)
{
    myAllocatedSize = 1024;
    myBases = (char*)malloc(myAllocatedSize + 1);
    myQualities = (char*)malloc(myAllocatedSize + 1);
    myAddDelAsBase = q.myAddDelAsBase;
    if((myBases == NULL ) || (myQualities == NULL))
    {
        // TODO, check for malloc failures.
        std::cerr << "Failed Memory Allocation\n";
    }
}


PileupElementBaseQual::~PileupElementBaseQual()
{
    if(myBases != NULL)
    {
        free(myBases);
        myBases = NULL;
    }
    if(myQualities != NULL)
    {
        free(myQualities);
        myQualities = NULL;
    }
}


// Add an entry to this pileup element.  
void PileupElementBaseQual::addEntry(SamRecord& record)
{
    // Call the base class:
    PileupElement::addEntry(record);

    // Increment the index
    ++myIndex;
    
    // if the index has gone beyond the allocated space, double the size.
    if(myIndex >= myAllocatedSize)
    {
        char* tempBuffer = (char*)realloc(myBases, myAllocatedSize * 2);
        if(tempBuffer == NULL)
        {
            std::cerr << "Memory Allocation Failure\n";
            // TODO
            return;
        }
        myBases = tempBuffer;
        tempBuffer = (char*)realloc(myQualities, myAllocatedSize * 2);
        if(tempBuffer == NULL)
        {
            std::cerr << "Memory Allocation Failure\n";
            // TODO
            return;
        }
        myQualities = tempBuffer;       
        myAllocatedSize = myAllocatedSize * 2;
    }

    Cigar* cigar = record.getCigarInfo();
    
    if(cigar == NULL)
    {
        throw std::runtime_error("Failed to retrieve cigar info from the record.");
    }


    int32_t readIndex = 
        cigar->getQueryIndex(getRefPosition(), record.get0BasedPosition());

    // If the readPosition is N/A, this is a deletion.
    if(readIndex != CigarRoller::INDEX_NA)
    {
        char base = record.getSequence(readIndex);
        char qual = record.getQuality(readIndex);
        if(qual == UNSET_QUAL)
        {
            qual = ' ';
        }
        myBases[myIndex] = base;
        myQualities[myIndex] = qual;
    }
    else if(myAddDelAsBase)
    {
        // This version adds deletions as bases.
        myBases[myIndex] = '-';
        myQualities[myIndex] = '0';
    }
    else
    {
        // Do not add a deletion.
        // Did not add any entries, so decrement the index counter since the
        // index was not used.
        --myIndex;
    }
}

void PileupElementBaseQual::analyze()
{
    if(getRefPosition() != UNSET_POSITION)
    {
        myBases[myIndex+1] = '\0';
        myQualities[myIndex+1] = '\0';
        std::cout << getChromosome() << "\t" << getRefPosition() << "\tN\t" << myIndex+1 << "\t";
        std::cout << myBases << "\t";
        std::cout << myQualities;
        std::cout << "\n";
    }
}

void PileupElementBaseQual::reset(int refPosition)
{
    // Call the base class.
    PileupElement::reset(refPosition);

    myIndex = -1;
}

