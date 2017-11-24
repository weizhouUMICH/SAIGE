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

#include "SamFile.h"
#include "SamFlag.h"
#include "Modify.h"

void testModify()
{
    modify modTest;
    modTest.testModify("testFiles/testSam.sam");
#ifdef __ZLIB_AVAILABLE__
    modTest.testModify("testFiles/testBam.bam");
#endif
}


void modify::testModify(const char* filename)
{
    myFilename = filename;

    modifyPosition();
    modifyCigar();
    
    modifyFlag();

    modifyTags();
}

void modify::modifyPosition()
{
    openAndRead1Rec();
   
    // Verify the initial bin.
    assert(samRecord.getBin() == 4681);

    // Change the position and verify that the bin is updated.
    assert(samRecord.set0BasedPosition(33768));

    // Verify the bin was updated.
    assert(samRecord.getBin() == 4683);
    assert(samRecord.get0BasedPosition() == 33768);
}


void modify::modifyCigar()
{
    openAndRead1Rec();
   
    // Verify the initial bin.
    assert(samRecord.getBin() == 4681);

    // Change the Cigar such that it modifies the bin.
    assert(samRecord.setCigar("33768M"));

    // Verify the bin was updated.
    assert(samRecord.getBin() == 585);
}


void modify::modifyFlag()
{
    openAndRead1Rec();
   
    // Verify the initial bin.
    uint16_t flag = 73;
    assert(samRecord.getFlag() == flag);

    SamFlag::setDuplicate(flag);
    assert(flag == 1097);
    assert(samRecord.setFlag(flag));
    assert(samRecord.getFlag() == 1097);

    SamFlag::setNotDuplicate(flag);
    assert(flag == 73);
    assert(samRecord.setFlag(flag));
    assert(samRecord.getFlag() == 73);
}


void modify::openAndRead1Rec()
{
    // Open the file for reading.   
    assert(samIn.OpenForRead(myFilename.c_str()));

    // Read the sam header.
    assert(samIn.ReadHeader(samHeader));
   
    // Read the first record.   
    assert(samIn.ReadRecord(samHeader, samRecord));
}


void modify::modifyTags()
{
    assert(samIn.OpenForRead(myFilename.c_str()));
    // Read the sam header.
    assert(samIn.ReadHeader(samHeader));
   
    SamFile samOut;
    SamFile bamOut;

    std::string inputType = myFilename.substr(myFilename.find_last_of('.'));
    std::string outFileBase = "results/updateTagFrom";
    if(inputType == ".bam")
    {
        outFileBase += "Bam";
    }
    else
    {
        outFileBase += "Sam";
    }

    std::string outFile = outFileBase + ".sam";
    assert(samOut.OpenForWrite(outFile.c_str()));
    outFile = outFileBase + ".bam";
    assert(bamOut.OpenForWrite(outFile.c_str()));
    assert(samOut.WriteHeader(samHeader));
    assert(bamOut.WriteHeader(samHeader));

    int count = 0;
    // Read the records.
    while(samIn.ReadRecord(samHeader, samRecord))
    {
        if(count == 0)
        {
            assert(samRecord.rmTag("MD", 'Z'));
        }
        else if(count == 2)
        {
            assert(samRecord.rmTags("XT:A;MD:Z;AB:c;NM:i"));
        }
        else if(count == 4)
        {
            assert(samRecord.rmTags("MD:Z,AB:c,NM:i"));
        }

        assert(bamOut.WriteRecord(samHeader, samRecord));
        assert(samOut.WriteRecord(samHeader, samRecord));
        ++count;
    }
}
