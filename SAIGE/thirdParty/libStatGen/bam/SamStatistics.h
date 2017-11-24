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

#ifndef __SAM_STATISTICS_H__
#define __SAM_STATISTICS_H__

#include <stdint.h>
#include "SamRecord.h"

class SamStatistics
{
public:
    SamStatistics();
    ~SamStatistics();

    // Reset the statistics - clear them for processing a new file.
    void reset();

    // Method to update the statistics to include the passed in record.
    bool updateStatistics(SamRecord& samRecord);

    void print();

private:

    ///////////////////////////////////////////////////////
    // Read Counts 

    /// The number of reads (records) that were processed.
    uint64_t myReadCount;
    
    /// The number of mapped reads (records).
    uint64_t myMappedReadCount;

    /// The number of paired reads (records).
    uint64_t myPairedReadCount;

    /// The number of proper paired reads (records).
    uint64_t myProperPairedReadCount;

    /// The number of duplicate reads (based on the flag).
    uint64_t myDupReadCount;

    /// The number of QC failure reads (based on the flag).
    uint64_t myQCFailureReadCount;

    ///////////////////////////////////////////////////////
    // Base Counts

    /// The total number of bases in the reads in the file (sum of read lengths)
    uint64_t myBaseCount;

    /// The total number of bases in mapped reads (sum of read lengths for mapped reads).
    uint64_t myMappedReadBases;
};

#endif
