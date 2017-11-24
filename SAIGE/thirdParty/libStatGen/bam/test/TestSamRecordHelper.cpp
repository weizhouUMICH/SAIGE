/*
 *  Copyright (C) 2012  Regents of the University of Michigan
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

#include "TestSamRecordHelper.h"
#include "TestValidate.h"
#include "SamRecordHelper.h"
#include <assert.h>

void testSamRecordHelper()
{
    // Call generic test.
    SamRecordHelperTest::testSamRecordHelper("testFiles/testSam.sam");
    //    SamRecordHelperTest::testSamRecordHelper("testFiles/testBam.bam");
}


void SamRecordHelperTest::testSamRecordHelper(const char* fileName)
{
    SamFile inSam;
    assert(inSam.OpenForRead(fileName));
    SamFileHeader samHeader;
    assert(inSam.ReadHeader(samHeader));
    validateHeader(samHeader);

    SamRecord samRecord;
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead1(samRecord);
    
    // Validate the entire sequence matches.
    assert(SamRecordHelper::checkSequence(samRecord, 
                                          TestValidate::READ1_POS, 
                                          TestValidate::READ1_SEQ.c_str()) == 0);

    // The read start position is 1010.
    // The sequence is CCGAA.
    assert(SamRecordHelper::checkSequence(samRecord, 1010, "CCGAA") == 0);
  
    // Test not matching.
    assert(SamRecordHelper::checkSequence(samRecord, 1010, "NNNNN") == -1);

    // Test match, but not at the start.
    assert(SamRecordHelper::checkSequence(samRecord, 1011, "CGA") == 1);

    // Test match not at the start, but to the end.
    assert(SamRecordHelper::checkSequence(samRecord, 1011, "CGAA") == 1);
  
    // Test run over the end.
    assert(SamRecordHelper::checkSequence(samRecord, 1011, "CGAAC") == -1);
  
}


