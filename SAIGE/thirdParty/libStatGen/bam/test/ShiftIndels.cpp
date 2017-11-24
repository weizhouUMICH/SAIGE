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

#include "ShiftIndels.h"
#include "SamFile.h"

void testShiftIndels()
{
    ShiftIndelsTest::testShift("testFiles/testShift.sam", "results/testShift.sam");
#ifdef __ZLIB_AVAILABLE__
    ShiftIndelsTest::testShift("testFiles/testShift.bam", "results/testShift.bam");
    ShiftIndelsTest::testShift("testFiles/testShift.bam", "results/testShiftFromBam.sam");
#endif
    ShiftIndelsTest::testShift("testFiles/testShift.sam", "results/testShiftFromSam.bam");
}

void ShiftIndelsTest::testShift(const char* input, const char* output)
{
    SamFile inSam, outSam;

    assert(inSam.OpenForRead(input));
    assert(outSam.OpenForWrite(output));


    // Read the SAM Header.
    SamFileHeader samHeader;
    assert(inSam.ReadHeader(samHeader));
    assert(outSam.WriteHeader(samHeader));


   SamRecord samRecord;
   int readNum = 1;
   bool shiftResult = true;
   while(inSam.ReadRecord(samHeader, samRecord))
   {
       if((readNum == 3)|| (readNum == 5))
       {
           shiftResult = false;
       }
       else
       {
           shiftResult = true;
       }
       ++readNum;

       assert(samRecord.shiftIndelsLeft() == shiftResult);
       assert(outSam.WriteRecord(samHeader, samRecord));
   }

}
