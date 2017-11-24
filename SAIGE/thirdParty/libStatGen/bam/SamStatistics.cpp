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

#include "SamStatistics.h"
#include <iomanip>
#include "SamFlag.h"

SamStatistics::SamStatistics()
{
    reset();
}

SamStatistics::~SamStatistics()
{
    reset();
}

void SamStatistics::reset()
{
    myReadCount = 0;
    myMappedReadCount = 0;
    myPairedReadCount = 0;
    myProperPairedReadCount = 0;
    myBaseCount = 0;
    myMappedReadBases = 0;
    myDupReadCount = 0;
    myQCFailureReadCount = 0;
}


bool SamStatistics::updateStatistics(SamRecord& samRecord)
{
    // Each record has one read, so update the read count.
    ++myReadCount;

    int32_t readLen = samRecord.getReadLength();

    // Get the flag to determine the type or 
    // read (mapped, paired, proper paired).
    uint16_t flag = samRecord.getFlag();

    // If the read is mapped, update the mapped c
    if(SamFlag::isMapped(flag))
    {
        ++myMappedReadCount;
        myMappedReadBases += readLen;
    }
    if(SamFlag::isPaired(flag))
    {
        ++myPairedReadCount;
        if(SamFlag::isProperPair(flag))
        {
            ++myProperPairedReadCount;
        }
    }
    if(SamFlag::isDuplicate(flag))
    {
        ++myDupReadCount;
    }
    if(SamFlag::isQCFailure(flag))
    {
        ++myQCFailureReadCount;
    }
    
    // Increment the total number of bases.
    myBaseCount += readLen;

    return(true);
}


void SamStatistics::print()
{
    double DIVIDE_UNITS = 1000000;
    std::string units = "(e6)";

    std::cerr << std::fixed << std::setprecision(2);

    // If total reads is less than DIVIDE_UNITS, just show the straight number.
    if(myReadCount < DIVIDE_UNITS)
    {
        DIVIDE_UNITS = 1;
        units.clear();
    }
    
    // Read Counts
    std::cerr << "TotalReads" << units << "\t"
              << myReadCount/DIVIDE_UNITS << std::endl;
    std::cerr << "MappedReads" << units << "\t" 
              << myMappedReadCount/DIVIDE_UNITS << std::endl;
    std::cerr << "PairedReads" << units << "\t" 
              << myPairedReadCount/DIVIDE_UNITS << std::endl;
    std::cerr << "ProperPair" << units << "\t" 
              << myProperPairedReadCount/DIVIDE_UNITS << std::endl;
    std::cerr << "DuplicateReads" << units << "\t" 
              << myDupReadCount/DIVIDE_UNITS << std::endl;
    std::cerr << "QCFailureReads" << units << "\t" 
              << myQCFailureReadCount/DIVIDE_UNITS << std::endl;
    std::cerr << std::endl;

    // Read Percentages
    std::cerr << "MappingRate(%)\t" 
              << 100 * myMappedReadCount/(double)myReadCount << std::endl;
    std::cerr << "PairedReads(%)\t" 
              << 100 * myPairedReadCount/(double)myReadCount << std::endl;
    std::cerr << "ProperPair(%)\t" 
              << 100 * myProperPairedReadCount/(double)myReadCount << std::endl;
    std::cerr << "DupRate(%)\t" 
              << 100 * myDupReadCount/(double)myReadCount << std::endl;
    std::cerr << "QCFailRate(%)\t" 
              << 100 * myQCFailureReadCount/(double)myReadCount << std::endl;
    std::cerr << std::endl;

    // Base Counts
    std::cerr << "TotalBases" << units << "\t"
              << myBaseCount/DIVIDE_UNITS << std::endl;
    std::cerr << "BasesInMappedReads" << units << "\t"
              << myMappedReadBases/DIVIDE_UNITS << std::endl;
}
