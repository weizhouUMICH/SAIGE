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

#include "SamCoordOutput.h"
#include "SamHelper.h"


SamCoordOutput::SamCoordOutput(SamRecordPool& pool)
    : myOutputFile(NULL),
      myHeader(NULL),
      myPool(&pool)
{
}

SamCoordOutput::~SamCoordOutput()
{
    // Flush the rest of the records.
    flush(-1, -1);

    myOutputFile = NULL;
    myHeader = NULL;
}

void SamCoordOutput::setOutputFile(SamFile* outFile, SamFileHeader* header)
{
    myOutputFile = outFile;
    myHeader = header;
}


bool SamCoordOutput::add(SamRecord* record)
{
    if(record != NULL)
    {
        int32_t chrom = record->getReferenceID();
        uint64_t chromPos = 
            SamHelper::combineChromPos(chrom, record->get0BasedPosition());
        myReadBuffer.insert(std::pair<uint64_t, SamRecord*>(chromPos, record));
        return(true);
    }
    return(false);
}


bool SamCoordOutput::flushAll()
{
    return(flush(-1,-1));
}

bool SamCoordOutput::flush(int32_t chromID, int32_t pos0Based)
{
    static std::multimap<uint64_t, SamRecord*>::iterator iter;

    uint64_t chromPos = SamHelper::combineChromPos(chromID, pos0Based);

    bool returnVal = true;
    iter = myReadBuffer.begin();
    
    if((myOutputFile == NULL) || (myHeader == NULL))
    {
        std::cerr <<
            "SamCoordOutput::flush, no output file/header is set, so records removed without being written\n";
        returnVal = false;
    }

    while((iter != myReadBuffer.end()) &&
          (((*iter).first <= chromPos) || (chromID == -1)))
    {
        if((myOutputFile != NULL) && (myHeader != NULL))
        {
            returnVal &= 
                myOutputFile->WriteRecord(*myHeader, *((*iter).second));
        }
        if(myPool != NULL)
        {
            myPool->releaseRecord((*iter).second);
        }
        else
        {
            delete((*iter).second);
        }
        ++iter;
    }
    // Remove the elements from the begining up to,
    // but not including the current iterator position.
    myReadBuffer.erase(myReadBuffer.begin(), iter);

    return(returnVal);
}
