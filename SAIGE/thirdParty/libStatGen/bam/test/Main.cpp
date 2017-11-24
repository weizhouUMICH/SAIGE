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

#include "ReadFiles.h"
#include "WriteFiles.h"
#include "ValidationTest.h"
#include "BamIndexTest.h"
#include "ModifyVar.h"
#include "Modify.h"
#include "SamFileTest.h"
#include "TestEquals.h"
#include "TestFilter.h"
#include "ShiftIndels.h"
#include "TestPileup.h"
#include "TestPosList.h"
#include "TestCigarHelper.h"
#include "TestSamRecordPool.h"
#include "TestSamCoordOutput.h"
#include "TestSamRecordHelper.h"

int main(int argc, char ** argv)
{
    if(argc == 1)
    {
        testReadSam();
#ifdef __ZLIB_AVAILABLE__
        testReadBam();
        testReadBam();
#endif
        testAddHeaderAndTagToFile("testFiles/testSam.sam",
                                  "results/addedTagToSam.bam");
        testAddHeaderAndTagToFile("testFiles/testSam.sam",
                                  "results/addedTagToSam.sam");
// Can't read bams without zlib
#ifdef __ZLIB_AVAILABLE__
        testAddHeaderAndTagToFile("testFiles/testBam.bam",
                                  "results/addedTagToBam.sam");
        testAddHeaderAndTagToFile("testFiles/testBam.bam",
                                  "results/addedTagToBam.bam");
#endif

        testValidateSortedRead();
      
        testWrite();

        testSamQNAME();
        testBamRID();
        testEmptyQual();
      
// Can't read bams without zlib
#ifdef __ZLIB_AVAILABLE__
        testBamIndex();
#endif

        testModifyVar();
        testModify();

        testSamFile();

        testSeqEquals();

        testFilter();

        testShiftIndels();

        testPileup();
        
        testPosList();

        testCigarHelper();
        testSamRecordPool();
        testSamCoordOutput();
        testSamRecordHelper();
    }
    else
    {
        modifyFirstBaseLong();
    }
}

