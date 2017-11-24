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

#include "TestValidate.h"
#include "BaseUtilities.h"

const std::string TestValidate::READ1_CIGAR = "5M2D";
const std::string TestValidate::READ1_SEQ = "CCGAA";
const std::string TestValidate::READ1_QUAL = "6>6+4";

const std::string TestValidate::READ6_CIGAR = "3S2H5M";
const std::string TestValidate::READ6_SEQ = "TGCACGTN";
const std::string TestValidate::READ6_QUAL = "453;>>>>";

const std::string TestValidate::READ7_CIGAR = "3S5M1S3H";
const std::string TestValidate::READ7_SEQ = "TGCACGTNG";
const std::string TestValidate::READ7_QUAL = "453;>>>>5";

void validateRead1(SamRecord& samRecord)
{
    //////////////////////////////////////////
    // Validate Record 1
    // Create record structure for validating.
    int expectedBlockSize = 89;
    const char* expectedReferenceName = "1";
    const char* expectedMateReferenceName = "1";
    const char* expectedMateReferenceNameOrEqual = "=";

    bamRecordStruct* expectedRecordPtr =
        (bamRecordStruct *) malloc(expectedBlockSize + sizeof(int));

    char tag[3];
    char type;
    void* value;
    const bamRecordStruct* bufferPtr;
    unsigned char* varPtr;

    expectedRecordPtr->myBlockSize = expectedBlockSize;
    expectedRecordPtr->myReferenceID = 0;
    expectedRecordPtr->myPosition = TestValidate::READ1_POS;
    expectedRecordPtr->myReadNameLength = 23;
    expectedRecordPtr->myMapQuality = 0;
    expectedRecordPtr->myBin = 4681;
    expectedRecordPtr->myCigarLength = 2;
    expectedRecordPtr->myFlag = 73;
    expectedRecordPtr->myReadLength = 5;
    expectedRecordPtr->myMateReferenceID = 0;
    expectedRecordPtr->myMatePosition = 1010;
    expectedRecordPtr->myInsertSize = 0;
   
    assert(samRecord.getString("MD") == "37");
    assert(samRecord.getString("YZ") == "");
    assert(samRecord.getInteger("YZ") == -1);
    float tmpFloat = -1;
    assert(samRecord.getFloatTag("YZ", tmpFloat) == false);
    // Check the alignment end
    assert(samRecord.get0BasedAlignmentEnd() == TestValidate::READ1_ALIGN_END);
    assert(samRecord.get1BasedAlignmentEnd() == (TestValidate::READ1_ALIGN_END + 1));
    assert(samRecord.getAlignmentLength() == TestValidate::READ1_ALIGN_LEN);
    assert(samRecord.get1BasedUnclippedStart() == (TestValidate::READ1_UNCLIP_START + 1));
    assert(samRecord.get0BasedUnclippedStart() == TestValidate::READ1_UNCLIP_START);
    assert(samRecord.get1BasedUnclippedEnd() == (TestValidate::READ1_UNCLIP_END + 1));
    assert(samRecord.get0BasedUnclippedEnd() == TestValidate::READ1_UNCLIP_END);

    // Check the accessors.
    assert(samRecord.getBlockSize() == expectedRecordPtr->myBlockSize);
    assert(samRecord.getReferenceID() == expectedRecordPtr->myReferenceID);
    assert(strcmp(samRecord.getReferenceName(), expectedReferenceName) == 0);
    assert(samRecord.get1BasedPosition() == expectedRecordPtr->myPosition + 1);
    assert(samRecord.get0BasedPosition() == expectedRecordPtr->myPosition);
    assert(samRecord.getReadNameLength() == 
           expectedRecordPtr->myReadNameLength);
    assert(samRecord.getMapQuality() == expectedRecordPtr->myMapQuality);
    assert(samRecord.getBin() == expectedRecordPtr->myBin);
    assert(samRecord.getCigarLength() == expectedRecordPtr->myCigarLength);
    assert(samRecord.getFlag() == expectedRecordPtr->myFlag);
    assert(samRecord.getReadLength() == expectedRecordPtr->myReadLength);
    assert(samRecord.getMateReferenceID() == 
           expectedRecordPtr->myMateReferenceID);
    assert(strcmp(samRecord.getMateReferenceName(),
                  expectedMateReferenceName) == 0);
    assert(strcmp(samRecord.getMateReferenceNameOrEqual(),
                  expectedMateReferenceNameOrEqual) == 0);
    assert(samRecord.get1BasedMatePosition() ==
           expectedRecordPtr->myMatePosition + 1);
    assert(samRecord.get0BasedMatePosition() == 
           expectedRecordPtr->myMatePosition);
    assert(samRecord.getInsertSize() == expectedRecordPtr->myInsertSize);
    assert(strcmp(samRecord.getReadName(), "1:1011:F:255+17M15D20M") == 0);
    assert(samRecord.getCigar() == TestValidate::READ1_CIGAR);
    assert(samRecord.getSequence() == TestValidate::READ1_SEQ);
    assert(samRecord.getQuality() == TestValidate::READ1_QUAL);

    assert(samRecord.getSequence(0) == 'C');
    assert(samRecord.getQuality(0) == '6');
    assert(samRecord.getSequence(1) == 'C');
    assert(samRecord.getQuality(1) == '>');
    assert(samRecord.getSequence(2) == 'G');
    assert(samRecord.getQuality(2) == '6');
    assert(samRecord.getSequence(3) == 'A');
    assert(samRecord.getQuality(3) == '+');
    assert(samRecord.getSequence(4) == 'A');
    assert(samRecord.getQuality(4) == '4');

    bool caught = false;
    try
    {
        samRecord.getSequence(-1);
    }
    catch (std::exception& e) 
    {
        caught = true;
        assert(strcmp(e.what(), "SamRecord::getSequence(-1) is out of range. Index must be between 0 and 4") == 0);
    }
    assert(caught == true);
    caught = false;
    try
    {
        samRecord.getQuality(-1);
    }
    catch (std::exception& e) 
    {
        caught = true;
        assert(strcmp(e.what(), "SamRecord::getQuality(-1) is out of range. Index must be between 0 and 4") == 0);
    }
    assert(caught == true);
    
    caught = false;
    try
    {
        samRecord.getSequence(5);
    }
    catch (std::exception& e) 
    {
        caught = true;
        assert(strcmp(e.what(), "SamRecord::getSequence(5) is out of range. Index must be between 0 and 4") == 0);
    }
    assert(caught == true);
    caught = false;
    try
    {
        samRecord.getQuality(5);
    }
    catch (std::exception& e) 
    {
        caught = true;
        assert(strcmp(e.what(), "SamRecord::getQuality(5) is out of range. Index must be between 0 and 4") == 0);
    }
    assert(caught == true);
    
    assert(samRecord.getNumOverlaps(1010, 1017) == 5);
    assert(samRecord.getNumOverlaps(1010, 1016) == 5);
    assert(samRecord.getNumOverlaps(1012, 1017) == 3);
    assert(samRecord.getNumOverlaps(1015, 1017) == 0);
    assert(samRecord.getNumOverlaps(1017, 1010) == 0);
    assert(samRecord.getNumOverlaps(1013, 1011) == 0);
    assert(samRecord.getNumOverlaps(-1, 1017) == 5);
    assert(samRecord.getNumOverlaps(1010, -1) == 5);

    // Check the tags.
    assert(samRecord.getNextSamTag(tag, type, &value) == true);
    assert(tag[0] == 'A');
    assert(tag[1] == 'M');
    assert(type == 'i');
    assert(*(char*)value == 0);
    assert(samRecord.getNextSamTag(tag, type, &value) == true);
    assert(tag[0] == 'M');
    assert(tag[1] == 'D');
    assert(type == 'Z');
    assert(*(String*)value == "37");
    assert(samRecord.getNextSamTag(tag, type, &value) == true);
    assert(tag[0] == 'N');
    assert(tag[1] == 'M');
    assert(type == 'i');
    assert(*(char*)value == 0);
    assert(samRecord.getNextSamTag(tag, type, &value) == true);
    assert(tag[0] == 'X');
    assert(tag[1] == 'T');
    assert(type == 'A');
    assert(*(char*)value == 'R');
    // No more tags, should return false.
    assert(samRecord.getNextSamTag(tag, type, &value) == false);
    assert(samRecord.getNextSamTag(tag, type, &value) == false);

    // Get the record ptr.   
    bufferPtr = (const bamRecordStruct*)samRecord.getRecordBuffer();
    // Validate the buffers match.
    assert(bufferPtr->myBlockSize == expectedRecordPtr->myBlockSize);
    assert(bufferPtr->myReferenceID == expectedRecordPtr->myReferenceID);
    assert(bufferPtr->myPosition == expectedRecordPtr->myPosition);
    assert(bufferPtr->myReadNameLength == expectedRecordPtr->myReadNameLength);
    assert(bufferPtr->myMapQuality == expectedRecordPtr->myMapQuality);
    assert(bufferPtr->myBin == expectedRecordPtr->myBin);
    assert(bufferPtr->myCigarLength == expectedRecordPtr->myCigarLength);
    assert(bufferPtr->myFlag == expectedRecordPtr->myFlag);
    assert(bufferPtr->myReadLength == expectedRecordPtr->myReadLength);
    assert(bufferPtr->myMateReferenceID == 
           expectedRecordPtr->myMateReferenceID);
    assert(bufferPtr->myMatePosition == expectedRecordPtr->myMatePosition);
    assert(bufferPtr->myInsertSize == expectedRecordPtr->myInsertSize);

    // Validate the variable length fields in the buffer.
    // Set the pointer to the start of the variable fields.
    varPtr = (unsigned char*)(&(bufferPtr->myData[0]));

    // Validate the readname.
    for(int i = 0; i < expectedRecordPtr->myReadNameLength; i++)
    {
        assert(*varPtr == samRecord.getReadName()[i]);
        varPtr++;
    }

    // Validate the cigar.
    // The First cigar is 5M which is 5 << 4 | 0 = 80
    assert(*(unsigned int*)varPtr == 80);
    // Increment the varptr the size of an int.
    varPtr += 4;
    // The 2nd cigar is 2D which is 2 << 4 | 2 = 34
    assert(*(unsigned int*)varPtr == 34);
    // Increment the varptr the size of an int.
    varPtr += 4;
   
    // Validate the sequence.
    // CC = 0x22
    assert(*varPtr == 0x22);
    varPtr++;
    // GA = 0x41
    assert(*varPtr == 0x41);
    varPtr++;
    // A  = 0x10
    assert(*varPtr == 0x10);
    varPtr++;
  
    // Validate the Quality
    for(int i = 0; i < expectedRecordPtr->myReadLength; i++)
    {
        assert(*varPtr == samRecord.getQuality()[i] - 33);
        varPtr++;
    }

    // Validate the tags.  
    assert(*varPtr == 'A');
    varPtr++;
    assert(*varPtr == 'M');
    varPtr++;
    assert(*varPtr == 'C');
    varPtr++;
    assert(*varPtr == 0);
    varPtr++;
    assert(*varPtr == 'M');
    varPtr++;
    assert(*varPtr == 'D');
    varPtr++;
    assert(*varPtr == 'Z');
    varPtr++;
    assert(*varPtr == '3');
    varPtr++;
    assert(*varPtr == '7');
    varPtr++;
    assert(*varPtr == 0);
    varPtr++;
    assert(*varPtr == 'N');
    varPtr++;
    assert(*varPtr == 'M');
    varPtr++;
    assert(*varPtr == 'C');
    varPtr++;
    assert(*varPtr == 0);
    varPtr++;
    assert(*varPtr == 'X');
    varPtr++;
    assert(*varPtr == 'T');
    varPtr++;
    assert(*varPtr == 'A');
    varPtr++;
    assert(*varPtr == 'R');
    varPtr++;
}


void validateRead2(SamRecord& samRecord)
{
    //////////////////////////////////////////
    // Validate Record 2
    // Create record structure for validating.
    int expectedBlockSize = 61;
    const char* expectedReferenceName = "1";
    const char* expectedMateReferenceName = "1";
    const char* expectedMateReferenceNameOrEqual = "=";

    bamRecordStruct* expectedRecordPtr =
        (bamRecordStruct *) malloc(expectedBlockSize + sizeof(int));

    char tag[3];
    char type;
    void* value;
    bamRecordStruct* bufferPtr;
    unsigned char* varPtr;

    expectedRecordPtr->myBlockSize = expectedBlockSize;
    expectedRecordPtr->myReferenceID = 0;
    expectedRecordPtr->myPosition = TestValidate::READ2_POS;
    expectedRecordPtr->myReadNameLength = 23;
    expectedRecordPtr->myMapQuality = 0;
    expectedRecordPtr->myBin = 4681;
    expectedRecordPtr->myCigarLength = 0;
    expectedRecordPtr->myFlag = 133;
    expectedRecordPtr->myReadLength = 4;
    expectedRecordPtr->myMateReferenceID = 0;
    expectedRecordPtr->myMatePosition = 1010;
    expectedRecordPtr->myInsertSize = 0;
   
    // Check the fields.
    bamRecordStruct retrieveRecord;
    String retrieveReadName;
    String retrieveCigar;
    String retrieveSequence;
    String retrieveQuality;

    assert(samRecord.getFields(retrieveRecord, retrieveReadName, 
                               retrieveCigar, retrieveSequence,
                               retrieveQuality) == true);
    assert(retrieveRecord.myBlockSize == expectedRecordPtr->myBlockSize);
    assert(retrieveRecord.myReferenceID == expectedRecordPtr->myReferenceID);
    assert(retrieveRecord.myPosition == expectedRecordPtr->myPosition);
    assert(retrieveRecord.myReadNameLength == 
           expectedRecordPtr->myReadNameLength);
    assert(retrieveRecord.myMapQuality == expectedRecordPtr->myMapQuality);
    assert(retrieveRecord.myBin == expectedRecordPtr->myBin);
    assert(retrieveRecord.myCigarLength == expectedRecordPtr->myCigarLength);
    assert(retrieveRecord.myFlag == expectedRecordPtr->myFlag);
    assert(retrieveRecord.myReadLength == expectedRecordPtr->myReadLength);
    assert(retrieveRecord.myMateReferenceID == 
           expectedRecordPtr->myMateReferenceID);
    assert(retrieveRecord.myMatePosition == expectedRecordPtr->myMatePosition);
    assert(retrieveRecord.myInsertSize == expectedRecordPtr->myInsertSize);

    // Check the alignment end
    assert(samRecord.getAlignmentLength() == 0);
    assert(samRecord.get0BasedAlignmentEnd() == 1011);
    assert(samRecord.get1BasedAlignmentEnd() == 1012);
    assert(samRecord.get0BasedUnclippedStart() == 1011);
    assert(samRecord.get1BasedUnclippedStart() == 1012);
    assert(samRecord.get0BasedUnclippedEnd() == 1011);
    assert(samRecord.get1BasedUnclippedEnd() == 1012);

    // Check the accessors.
    assert(samRecord.getBlockSize() == expectedRecordPtr->myBlockSize);
    assert(samRecord.getReferenceID() == expectedRecordPtr->myReferenceID);
    assert(strcmp(samRecord.getReferenceName(), expectedReferenceName) == 0);
    assert(samRecord.get1BasedPosition() == expectedRecordPtr->myPosition + 1);
    assert(samRecord.get0BasedPosition() == expectedRecordPtr->myPosition);
    assert(samRecord.getReadNameLength() == 
           expectedRecordPtr->myReadNameLength);
    assert(samRecord.getMapQuality() == expectedRecordPtr->myMapQuality);
    assert(samRecord.getBin() == expectedRecordPtr->myBin);
    assert(samRecord.getCigarLength() == expectedRecordPtr->myCigarLength);
    assert(samRecord.getFlag() == expectedRecordPtr->myFlag);
    assert(samRecord.getReadLength() == expectedRecordPtr->myReadLength);
    assert(samRecord.getMateReferenceID() == 
           expectedRecordPtr->myMateReferenceID);
    assert(strcmp(samRecord.getMateReferenceName(), 
                  expectedMateReferenceName) == 0);
    assert(strcmp(samRecord.getMateReferenceNameOrEqual(), 
                  expectedMateReferenceNameOrEqual) == 0);
    assert(samRecord.get1BasedMatePosition() ==
           expectedRecordPtr->myMatePosition + 1);
    assert(samRecord.get0BasedMatePosition() == 
           expectedRecordPtr->myMatePosition);
    assert(samRecord.getInsertSize() == expectedRecordPtr->myInsertSize);
    assert(strcmp(samRecord.getReadName(), "1:1011:F:255+17M15D20M") == 0);
    assert(strcmp(samRecord.getCigar(), "*") == 0);
    assert(strcmp(samRecord.getSequence(), "CTGT") == 0);
    assert(strcmp(samRecord.getQuality(), ">>9>") == 0);

    assert(samRecord.getSequence(0) == 'C');
    assert(samRecord.getQuality(0) == '>');
    assert(samRecord.getSequence(1) == 'T');
    assert(samRecord.getQuality(1) == '>');
    assert(samRecord.getSequence(2) == 'G');
    assert(samRecord.getQuality(2) == '9');
    assert(samRecord.getSequence(3) == 'T');
    assert(samRecord.getQuality(3) == '>');
    bool caught = false;
    try
    {
        samRecord.getSequence(-1);
    }
    catch (std::exception& e) 
    {
        caught = true;
        assert(strcmp(e.what(), "SamRecord::getSequence(-1) is out of range. Index must be between 0 and 3") == 0);
    }
    assert(caught == true);
    caught = false;
    try
    {
        samRecord.getQuality(-1);
    }
    catch (std::exception& e) 
    {
        caught = true;
        assert(strcmp(e.what(), "SamRecord::getQuality(-1) is out of range. Index must be between 0 and 3") == 0);
    }
    assert(caught == true);
    
    caught = false;
    try
    {
        samRecord.getSequence(4);
    }
    catch (std::exception& e) 
    {
        caught = true;
        assert(strcmp(e.what(), "SamRecord::getSequence(4) is out of range. Index must be between 0 and 3") == 0);
    }
    assert(caught == true);
    caught = false;
    try
    {
        samRecord.getQuality(4);
    }
    catch (std::exception& e) 
    {
        caught = true;
        assert(strcmp(e.what(), "SamRecord::getQuality(4) is out of range. Index must be between 0 and 3") == 0);
    }
    assert(caught == true);
    
    assert(samRecord.getNumOverlaps(1011, 1017) == 0);
    assert(samRecord.getNumOverlaps(0, 1116) == 0);

    // No Tags to check, should return false.
    assert(samRecord.getNextSamTag(tag, type, &value) == false);
    assert(samRecord.getNextSamTag(tag, type, &value) == false);

    // Get the record ptr.   
    bufferPtr = (bamRecordStruct*)samRecord.getRecordBuffer();
    // Validate the buffers match.
    assert(bufferPtr->myBlockSize == expectedRecordPtr->myBlockSize);
    assert(bufferPtr->myReferenceID == expectedRecordPtr->myReferenceID);
    assert(bufferPtr->myPosition == expectedRecordPtr->myPosition);
    assert(bufferPtr->myReadNameLength == expectedRecordPtr->myReadNameLength);
    assert(bufferPtr->myMapQuality == expectedRecordPtr->myMapQuality);
    assert(bufferPtr->myBin == expectedRecordPtr->myBin);
    assert(bufferPtr->myCigarLength == expectedRecordPtr->myCigarLength);
    assert(bufferPtr->myFlag == expectedRecordPtr->myFlag);
    assert(bufferPtr->myReadLength == expectedRecordPtr->myReadLength);
    assert(bufferPtr->myMateReferenceID == 
           expectedRecordPtr->myMateReferenceID);
    assert(bufferPtr->myMatePosition == expectedRecordPtr->myMatePosition);
    assert(bufferPtr->myInsertSize == expectedRecordPtr->myInsertSize);

    // Validate the variable length fields in the buffer.
    // Set the pointer to the start of the variable fields.
    varPtr = (unsigned char*)(&(bufferPtr->myData[0]));

    // Validate the readname.
    for(int i = 0; i < expectedRecordPtr->myReadNameLength; i++)
    {
        assert(*varPtr == samRecord.getReadName()[i]);
        varPtr++;
    }

    // No cigar to validate. 
    // Validate the sequence.
    // CT = 0x28
    assert(*varPtr == 0x28);
    varPtr++;
    // GT = 0x48
    assert(*varPtr == 0x48);
    varPtr++;

    // Validate the Quality
    for(int i = 0; i < expectedRecordPtr->myReadLength; i++)
    {
        assert(*varPtr == samRecord.getQuality()[i] - 33);
        varPtr++;
    }

    // No tags.  
}


void validateRead3(SamRecord& samRecord)
{
    //////////////////////////////////////////
    // Validate Record 3
    // Create record structure for validating.
    int expectedBlockSize = 87;
    const char* expectedReferenceName = "1";
    const char* expectedMateReferenceName = "18";
    const char* expectedMateReferenceNameOrEqual = "18";

    bamRecordStruct* expectedRecordPtr =
        (bamRecordStruct *) malloc(expectedBlockSize + sizeof(int));

    char tag[3];
    char type;
    void* value;
    bamRecordStruct* bufferPtr;
    unsigned char* varPtr;

    expectedRecordPtr->myBlockSize = expectedBlockSize;
    expectedRecordPtr->myReferenceID = 0;
    expectedRecordPtr->myPosition = 74;
    expectedRecordPtr->myReadNameLength = 21;
    expectedRecordPtr->myMapQuality = 0;
    expectedRecordPtr->myBin = 4681;
    expectedRecordPtr->myCigarLength = 1;
    expectedRecordPtr->myFlag = 97;
    expectedRecordPtr->myReadLength = 5;
    expectedRecordPtr->myMateReferenceID = 17;
    expectedRecordPtr->myMatePosition = 756;
    expectedRecordPtr->myInsertSize = 0;
   
    // Check the accessors.
    assert(samRecord.getBlockSize() == expectedRecordPtr->myBlockSize);
    assert(samRecord.getReferenceID() == expectedRecordPtr->myReferenceID);
    assert(strcmp(samRecord.getReferenceName(), expectedReferenceName) == 0);
    assert(samRecord.get1BasedPosition() == expectedRecordPtr->myPosition + 1);
    assert(samRecord.get0BasedPosition() == expectedRecordPtr->myPosition);
    assert(samRecord.getReadNameLength() == 
           expectedRecordPtr->myReadNameLength);
    assert(samRecord.getMapQuality() == expectedRecordPtr->myMapQuality);
    assert(samRecord.getBin() == expectedRecordPtr->myBin);
    assert(samRecord.getCigarLength() == expectedRecordPtr->myCigarLength);
    assert(samRecord.getFlag() == expectedRecordPtr->myFlag);
    assert(samRecord.getReadLength() == expectedRecordPtr->myReadLength);
    assert(samRecord.getMateReferenceID() == 
           expectedRecordPtr->myMateReferenceID);
    assert(strcmp(samRecord.getMateReferenceName(), 
                  expectedMateReferenceName) == 0);
    assert(strcmp(samRecord.getMateReferenceNameOrEqual(),
                  expectedMateReferenceNameOrEqual) == 0);
    assert(samRecord.get1BasedMatePosition() == 
           expectedRecordPtr->myMatePosition + 1);
    assert(samRecord.get0BasedMatePosition() ==
           expectedRecordPtr->myMatePosition);
    assert(samRecord.getInsertSize() == expectedRecordPtr->myInsertSize);
    assert(strcmp(samRecord.getReadName(), "18:462+29M5I3M:F:295") == 0);
    assert(strcmp(samRecord.getCigar(), "5M") == 0);
    assert(strcmp(samRecord.getSequence(), "ACGTN") == 0);
    assert(strcmp(samRecord.getQuality(), ";>>>>") == 0);
    assert(samRecord.getNumOverlaps(74, 79) == 5);
    assert(samRecord.getNumOverlaps(74, 78) == 4);
    assert(samRecord.getNumOverlaps(73, 79) == 5);
    assert(samRecord.getNumOverlaps(75, 79) == 4);
    assert(samRecord.getNumOverlaps(0, 179) == 5);
    assert(samRecord.getNumOverlaps(0, 19) == 0);

    // Check the alignment end
    assert(samRecord.get0BasedAlignmentEnd() == 78);
    assert(samRecord.get1BasedAlignmentEnd() == 79);
    assert(samRecord.getAlignmentLength() == 5);
    assert(samRecord.get0BasedUnclippedStart() == 74);
    assert(samRecord.get1BasedUnclippedStart() == 75);
    assert(samRecord.get0BasedUnclippedEnd() == 78);
    assert(samRecord.get1BasedUnclippedEnd() == 79);


    // Check the tags.
    assert(samRecord.getNextSamTag(tag, type, &value) == true);
    assert(tag[0] == 'A');
    assert(tag[1] == 'M');
    assert(type == 'i');
    assert(*(char*)value == 0);
    assert(samRecord.getNextSamTag(tag, type, &value) == true);
    assert(tag[0] == 'M');
    assert(tag[1] == 'D');
    assert(type == 'Z');
    assert(*(String*)value == "30A0C5");
    assert(samRecord.getNextSamTag(tag, type, &value) == true);
    assert(tag[0] == 'N');
    assert(tag[1] == 'M');
    assert(type == 'i');
    assert(*(char*)value == 2);
    assert(samRecord.getNextSamTag(tag, type, &value) == true);
    assert(tag[0] == 'X');
    assert(tag[1] == 'T');
    assert(type == 'A');
    assert(*(char*)value == 'R');
    // No more tags, should return false.
    assert(samRecord.getNextSamTag(tag, type, &value) == false);
    assert(samRecord.getNextSamTag(tag, type, &value) == false);

    // Get the record ptr.   
    bufferPtr = (bamRecordStruct*)samRecord.getRecordBuffer();
    // Validate the buffers match.
    assert(bufferPtr->myBlockSize == expectedRecordPtr->myBlockSize);
    assert(bufferPtr->myReferenceID == expectedRecordPtr->myReferenceID);
    assert(bufferPtr->myPosition == expectedRecordPtr->myPosition);
    assert(bufferPtr->myReadNameLength == expectedRecordPtr->myReadNameLength);
    assert(bufferPtr->myMapQuality == expectedRecordPtr->myMapQuality);
    assert(bufferPtr->myBin == expectedRecordPtr->myBin);
    assert(bufferPtr->myCigarLength == expectedRecordPtr->myCigarLength);
    assert(bufferPtr->myFlag == expectedRecordPtr->myFlag);
    assert(bufferPtr->myReadLength == expectedRecordPtr->myReadLength);
    assert(bufferPtr->myMateReferenceID == 
           expectedRecordPtr->myMateReferenceID);
    assert(bufferPtr->myMatePosition == expectedRecordPtr->myMatePosition);
    assert(bufferPtr->myInsertSize == expectedRecordPtr->myInsertSize);

    // Validate the variable length fields in the buffer.
    // Set the pointer to the start of the variable fields.
    varPtr = (unsigned char*)(&(bufferPtr->myData[0]));

    // Validate the readname.
    for(int i = 0; i < expectedRecordPtr->myReadNameLength; i++)
    {
        assert(*varPtr == samRecord.getReadName()[i]);
        varPtr++;
    }

    // Validate the cigar.
    // The cigar is 5M which is 5 << 4 | 0 = 80
    assert(*(unsigned int*)varPtr == 80);
    // Increment the varptr the size of an int.
    varPtr += 4;
   
    // Validate the sequence.
    // AC = 0x12
    assert(*varPtr == 0x12);
    varPtr++;
    // GT = 0x48
    assert(*varPtr == 0x48);
    varPtr++;
    // N  = 0xF0
    assert(*varPtr == 0xF0);
    varPtr++;

    // Validate the Quality
    for(int i = 0; i < expectedRecordPtr->myReadLength; i++)
    {
        assert(*varPtr == samRecord.getQuality()[i] - 33);
        varPtr++;
    }

    // Validate the tags.  
    assert(*varPtr == 'A');
    varPtr++;
    assert(*varPtr == 'M');
    varPtr++;
    assert(*varPtr == 'C');
    varPtr++;
    assert(*varPtr == 0);
    varPtr++;
    assert(*varPtr == 'M');
    varPtr++;
    assert(*varPtr == 'D');
    varPtr++;
    assert(*varPtr == 'Z');
    varPtr++;
    assert(*varPtr == '3');
    varPtr++;
    assert(*varPtr == '0');
    varPtr++;
    assert(*varPtr == 'A');
    varPtr++;
    assert(*varPtr == '0');
    varPtr++;
    assert(*varPtr == 'C');
    varPtr++;
    assert(*varPtr == '5');
    varPtr++;
    assert(*varPtr == 0);
    varPtr++;
    assert(*varPtr == 'N');
    varPtr++;
    assert(*varPtr == 'M');
    varPtr++;
    assert(*varPtr == 'C');
    varPtr++;
    assert(*varPtr == 2);
    varPtr++;
    assert(*varPtr == 'X');
    varPtr++;
    assert(*varPtr == 'T');
    varPtr++;
    assert(*varPtr == 'A');
    varPtr++;
    assert(*varPtr == 'R');
    varPtr++;
}


void validateRead4(SamRecord& samRecord)
{
    //////////////////////////////////////////
    // Validate Record 4
    // Create record structure for validating.
    int expectedBlockSize = 57;
    const char* expectedReferenceName = "1";
    const char* expectedMateReferenceName = "18";
    const char* expectedMateReferenceNameOrEqual = "18";

    bamRecordStruct* expectedRecordPtr =
        (bamRecordStruct *) malloc(expectedBlockSize + sizeof(int));

    char tag[3];
    char type;
    void* value;
    bamRecordStruct* bufferPtr;
    unsigned char* varPtr;

    expectedRecordPtr->myBlockSize = expectedBlockSize;
    expectedRecordPtr->myReferenceID = 0;
    expectedRecordPtr->myPosition = 74;
    expectedRecordPtr->myReadNameLength = 21;
    expectedRecordPtr->myMapQuality = 0;
    expectedRecordPtr->myBin = 4681;
    expectedRecordPtr->myCigarLength = 0;
    expectedRecordPtr->myFlag = 97;
    expectedRecordPtr->myReadLength = 0;
    expectedRecordPtr->myMateReferenceID = 17;
    expectedRecordPtr->myMatePosition = 756;
    expectedRecordPtr->myInsertSize = 0;
   
    // Check the alignment end
    assert(samRecord.get1BasedUnclippedEnd() == 75);
    assert(samRecord.get0BasedUnclippedEnd() == 74);
    assert(samRecord.get0BasedUnclippedStart() == 74);
    assert(samRecord.get1BasedUnclippedStart() == 75);
    assert(samRecord.get1BasedAlignmentEnd() == 75);
    assert(samRecord.get0BasedAlignmentEnd() == 74);
    assert(samRecord.getAlignmentLength() == 0);

    // Check the accessors.
    assert(samRecord.getBlockSize() == expectedRecordPtr->myBlockSize);
    assert(samRecord.getReferenceID() == expectedRecordPtr->myReferenceID);
    assert(strcmp(samRecord.getReferenceName(), expectedReferenceName) == 0);
    assert(samRecord.get1BasedPosition() == expectedRecordPtr->myPosition + 1);
    assert(samRecord.get0BasedPosition() == expectedRecordPtr->myPosition);
    assert(samRecord.getReadNameLength() ==
           expectedRecordPtr->myReadNameLength);
    assert(samRecord.getMapQuality() == expectedRecordPtr->myMapQuality);
    assert(samRecord.getBin() == expectedRecordPtr->myBin);
    assert(samRecord.getCigarLength() == expectedRecordPtr->myCigarLength);
    assert(samRecord.getFlag() == expectedRecordPtr->myFlag);
    assert(samRecord.getReadLength() == expectedRecordPtr->myReadLength);
    assert(samRecord.getMateReferenceID() ==
           expectedRecordPtr->myMateReferenceID);
    assert(strcmp(samRecord.getMateReferenceName(),
                  expectedMateReferenceName) == 0);
    assert(strcmp(samRecord.getMateReferenceNameOrEqual(),
                  expectedMateReferenceNameOrEqual) == 0);
    assert(samRecord.get1BasedMatePosition() ==
           expectedRecordPtr->myMatePosition + 1);
    assert(samRecord.get0BasedMatePosition() ==
           expectedRecordPtr->myMatePosition);
    assert(samRecord.getInsertSize() == expectedRecordPtr->myInsertSize);
    assert(strcmp(samRecord.getReadName(), "18:462+29M5I3M:F:295") == 0);
    assert(strcmp(samRecord.getCigar(), "*") == 0);
    assert(strcmp(samRecord.getSequence(), "*") == 0);
    assert(strcmp(samRecord.getQuality(), "*") == 0);

    bool caught = false;
    try
    {
        samRecord.getSequence(0);
    }
    catch (std::exception& e) 
    {
        caught = true;
        assert(strcmp(e.what(), "SamRecord::getSequence(0) is not allowed since sequence = '*'") == 0);
    }
    assert(caught == true);
    caught = false;
    try
    {
        assert(samRecord.getQuality(0) == BaseUtilities::UNKNOWN_QUALITY_CHAR);
    }
    catch (std::exception& e) 
    {
        caught = true;
    }
    assert(caught == false);
    try
    {
        samRecord.getSequence(-1);
    }
    catch (std::exception& e) 
    {
        caught = true;
        assert(strcmp(e.what(), "SamRecord::getSequence(-1) is not allowed since sequence = '*'") == 0);
    }
    assert(caught == true);
    caught = false;
    try
    {
        assert(samRecord.getQuality(-1) == BaseUtilities::UNKNOWN_QUALITY_CHAR);
    }
    catch (std::exception& e) 
    {
        caught = true;
    }
    assert(caught == false);
    
    caught = false;
    try
    {
        samRecord.getSequence(5);
    }
    catch (std::exception& e) 
    {
        caught = true;
        assert(strcmp(e.what(), "SamRecord::getSequence(5) is not allowed since sequence = '*'") == 0);
    }
    assert(caught == true);
    caught = false;
    try
    {
        assert(samRecord.getQuality(5) == BaseUtilities::UNKNOWN_QUALITY_CHAR);
    }
    catch (std::exception& e) 
    {
        caught = true;
    }
    assert(caught == false);

    assert(samRecord.getNumOverlaps(74, 79) == 0);
    assert(samRecord.getNumOverlaps(74, 78) == 0);
    assert(samRecord.getNumOverlaps(73, 79) == 0);
    assert(samRecord.getNumOverlaps(75, 79) == 0);
    assert(samRecord.getNumOverlaps(0, 179) == 0);
    assert(samRecord.getNumOverlaps(0, 19) == 0);

    // Check the tag.
    assert(samRecord.getNextSamTag(tag, type, &value) == true);
    assert(tag[0] == 'A');
    assert(tag[1] == 'M');
    assert(type == 'i');
    assert(*(char*)value == 0);
    // No more Tags to check, should return false.
    assert(samRecord.getNextSamTag(tag, type, &value) == false);
    assert(samRecord.getNextSamTag(tag, type, &value) == false);

    // Get the record ptr.   
    bufferPtr = (bamRecordStruct*)samRecord.getRecordBuffer();
    // Validate the buffers match.
    assert(bufferPtr->myBlockSize == expectedRecordPtr->myBlockSize);
    assert(bufferPtr->myReferenceID == expectedRecordPtr->myReferenceID);
    assert(bufferPtr->myPosition == expectedRecordPtr->myPosition);
    assert(bufferPtr->myReadNameLength == expectedRecordPtr->myReadNameLength);
    assert(bufferPtr->myMapQuality == expectedRecordPtr->myMapQuality);
    assert(bufferPtr->myBin == expectedRecordPtr->myBin);
    assert(bufferPtr->myCigarLength == expectedRecordPtr->myCigarLength);
    assert(bufferPtr->myFlag == expectedRecordPtr->myFlag);
    assert(bufferPtr->myReadLength == expectedRecordPtr->myReadLength);
    assert(bufferPtr->myMateReferenceID ==
           expectedRecordPtr->myMateReferenceID);
    assert(bufferPtr->myMatePosition == expectedRecordPtr->myMatePosition);
    assert(bufferPtr->myInsertSize == expectedRecordPtr->myInsertSize);

    // Validate the variable length fields in the buffer.
    // Set the pointer to the start of the variable fields.
    varPtr = (unsigned char*)(&(bufferPtr->myData[0]));

    // Validate the readname.
    for(int i = 0; i < expectedRecordPtr->myReadNameLength; i++)
    {
        assert(*varPtr == samRecord.getReadName()[i]);
        varPtr++;
    }

    // No cigar to validate. 
    // Validate the sequence.
    // No sequence.
    // No Quality.

    // Validate the tags.  
    assert(*varPtr == 'A');
    varPtr++;
    assert(*varPtr == 'M');
    varPtr++;
    assert(*varPtr == 'C');
    varPtr++;
    assert(*varPtr == 0);
    varPtr++;
}


void validateRead5(SamRecord& samRecord)
{
    //////////////////////////////////////////
    // Validate Record 5
    int expectedBlockSize = 87;
    const char* expectedReferenceName = "2";
    const char* expectedMateReferenceName = "18";
    const char* expectedMateReferenceNameOrEqual = "18";

    bamRecordStruct* expectedRecordPtr =
        (bamRecordStruct *) malloc(expectedBlockSize + sizeof(int));

    char tag[3];
    char type;
    void* value;
    bamRecordStruct* bufferPtr;
    unsigned char* varPtr;

    expectedRecordPtr->myBlockSize = expectedBlockSize;
    expectedRecordPtr->myReferenceID = 1;
    expectedRecordPtr->myPosition = 74;
    expectedRecordPtr->myReadNameLength = 21;
    expectedRecordPtr->myMapQuality = 0;
    expectedRecordPtr->myBin = 4681;
    expectedRecordPtr->myCigarLength = 1;
    expectedRecordPtr->myFlag = 97;
    expectedRecordPtr->myReadLength = 5;
    expectedRecordPtr->myMateReferenceID = 17;
    expectedRecordPtr->myMatePosition = 756;
    expectedRecordPtr->myInsertSize = 0;
   
    // Check the accessors.
    assert(samRecord.getBlockSize() == expectedRecordPtr->myBlockSize);
    assert(samRecord.getReferenceID() == expectedRecordPtr->myReferenceID);
    assert(strcmp(samRecord.getReferenceName(), expectedReferenceName) == 0);
    assert(samRecord.get1BasedPosition() == expectedRecordPtr->myPosition + 1);
    assert(samRecord.get0BasedPosition() == expectedRecordPtr->myPosition);
    assert(samRecord.getReadNameLength() ==
           expectedRecordPtr->myReadNameLength);
    assert(samRecord.getMapQuality() == expectedRecordPtr->myMapQuality);
    assert(samRecord.getBin() == expectedRecordPtr->myBin);
    assert(samRecord.getCigarLength() == expectedRecordPtr->myCigarLength);
    assert(samRecord.getFlag() == expectedRecordPtr->myFlag);
    assert(samRecord.getReadLength() == expectedRecordPtr->myReadLength);
    assert(samRecord.getMateReferenceID() ==
           expectedRecordPtr->myMateReferenceID);
    assert(strcmp(samRecord.getMateReferenceName(),
                  expectedMateReferenceName) == 0);
    assert(strcmp(samRecord.getMateReferenceNameOrEqual(),
                  expectedMateReferenceNameOrEqual) == 0);
    assert(samRecord.get1BasedMatePosition() == 
           expectedRecordPtr->myMatePosition + 1);
    assert(samRecord.get0BasedMatePosition() ==
           expectedRecordPtr->myMatePosition);
    assert(samRecord.getInsertSize() == expectedRecordPtr->myInsertSize);
    assert(strcmp(samRecord.getReadName(), "18:462+29M5I3M:F:295") == 0);
    assert(strcmp(samRecord.getCigar(), "5M") == 0);
    assert(strcmp(samRecord.getSequence(), "ACGTN") == 0);
    assert(strcmp(samRecord.getQuality(), "*") == 0);
    assert(samRecord.getNumOverlaps(74, 79) == 5);
    assert(samRecord.getNumOverlaps(74, 78) == 4);
    assert(samRecord.getNumOverlaps(73, 79) == 5);
    assert(samRecord.getNumOverlaps(75, 79) == 4);
    assert(samRecord.getNumOverlaps(0, 179) == 5);
    assert(samRecord.getNumOverlaps(0, 19) == 0);

    assert(samRecord.getSequence(0) == 'A');
    char expChar = BaseUtilities::UNKNOWN_QUALITY_CHAR;
    assert(samRecord.getQuality(0) == expChar);
    assert(samRecord.getSequence(1) == 'C');
    assert(samRecord.getQuality(1) == expChar);
    assert(samRecord.getSequence(2) == 'G');
    assert(samRecord.getQuality(2) == expChar);
    assert(samRecord.getSequence(3) == 'T');
    assert(samRecord.getQuality(3) == expChar);
    assert(samRecord.getSequence(4) == 'N');
    assert(samRecord.getQuality(4) == expChar);

    bool caught = false;
    try
    {
        samRecord.getSequence(-1);
    }
    catch (std::exception& e) 
    {
        caught = true;
        assert(strcmp(e.what(), "SamRecord::getSequence(-1) is out of range. Index must be between 0 and 4") == 0);
    }
    assert(caught == true);
    caught = false;
    try
    {
        samRecord.getQuality(-1);
    }
    catch (std::exception& e) 
    {
        caught = true;
        assert(strcmp(e.what(), "SamRecord::getQuality(-1) is out of range. Index must be between 0 and 4") == 0);
    }
    assert(caught == true);
    
    caught = false;
    try
    {
        samRecord.getSequence(5);
    }
    catch (std::exception& e) 
    {
        caught = true;
        assert(strcmp(e.what(), "SamRecord::getSequence(5) is out of range. Index must be between 0 and 4") == 0);
    }
    assert(caught == true);
    caught = false;
    try
    {
        samRecord.getQuality(5);
    }
    catch (std::exception& e) 
    {
        caught = true;
        assert(strcmp(e.what(), "SamRecord::getQuality(5) is out of range. Index must be between 0 and 4") == 0);
    }
    assert(caught == true);

    // Check the tags.
    assert(samRecord.getNextSamTag(tag, type, &value) == true);
    assert(tag[0] == 'A');
    assert(tag[1] == 'M');
    assert(type == 'i');
    assert(*(char*)value == 0);
    assert(samRecord.getNextSamTag(tag, type, &value) == true);
    assert(tag[0] == 'M');
    assert(tag[1] == 'D');
    assert(type == 'Z');
    assert(*(String*)value == "30A0C5");
    assert(samRecord.getNextSamTag(tag, type, &value) == true);
    assert(tag[0] == 'N');
    assert(tag[1] == 'M');
    assert(type == 'i');
    assert(*(char*)value == 2);
    assert(samRecord.getNextSamTag(tag, type, &value) == true);
    assert(tag[0] == 'X');
    assert(tag[1] == 'T');
    assert(type == 'A');
    assert(*(char*)value == 'R');
    // No more tags, should return false.
    assert(samRecord.getNextSamTag(tag, type, &value) == false);
    assert(samRecord.getNextSamTag(tag, type, &value) == false);

    // Get the record ptr.   
    bufferPtr = (bamRecordStruct*)samRecord.getRecordBuffer();
    // Validate the buffers match.
    assert(bufferPtr->myBlockSize == expectedRecordPtr->myBlockSize);
    assert(bufferPtr->myReferenceID == expectedRecordPtr->myReferenceID);
    assert(bufferPtr->myPosition == expectedRecordPtr->myPosition);
    assert(bufferPtr->myReadNameLength == expectedRecordPtr->myReadNameLength);
    assert(bufferPtr->myMapQuality == expectedRecordPtr->myMapQuality);
    assert(bufferPtr->myBin == expectedRecordPtr->myBin);
    assert(bufferPtr->myCigarLength == expectedRecordPtr->myCigarLength);
    assert(bufferPtr->myFlag == expectedRecordPtr->myFlag);
    assert(bufferPtr->myReadLength == expectedRecordPtr->myReadLength);
    assert(bufferPtr->myMateReferenceID ==
           expectedRecordPtr->myMateReferenceID);
    assert(bufferPtr->myMatePosition == expectedRecordPtr->myMatePosition);
    assert(bufferPtr->myInsertSize == expectedRecordPtr->myInsertSize);

    // Validate the variable length fields in the buffer.
    // Set the pointer to the start of the variable fields.
    varPtr = (unsigned char*)(&(bufferPtr->myData[0]));

    // Validate the readname.
    for(int i = 0; i < expectedRecordPtr->myReadNameLength; i++)
    {
        assert(*varPtr == samRecord.getReadName()[i]);
        varPtr++;
    }

    // Validate the cigar.
    // The cigar is 5M which is 5 << 4 | 0 = 80
    assert(*(unsigned int*)varPtr == 80);
    // Increment the varptr the size of an int.
    varPtr += 4;
   
    // Validate the sequence.
    // AC = 0x12
    assert(*varPtr == 0x12);
    varPtr++;
    // GT = 0x48
    assert(*varPtr == 0x48);
    varPtr++;
    // N  = 0xF0
    assert(*varPtr == 0xF0);
    varPtr++;

    // Validate the Quality
    for(int i = 0; i < expectedRecordPtr->myReadLength; i++)
    {
        assert(*varPtr == 0xFF);
        varPtr++;
    }

    // Validate the tags.  
    assert(*varPtr == 'A');
    varPtr++;
    assert(*varPtr == 'M');
    varPtr++;
    assert(*varPtr == 'C');
    varPtr++;
    assert(*varPtr == 0);
    varPtr++;
    assert(*varPtr == 'M');
    varPtr++;
    assert(*varPtr == 'D');
    varPtr++;
    assert(*varPtr == 'Z');
    varPtr++;
    assert(*varPtr == '3');
    varPtr++;
    assert(*varPtr == '0');
    varPtr++;
    assert(*varPtr == 'A');
    varPtr++;
    assert(*varPtr == '0');
    varPtr++;
    assert(*varPtr == 'C');
    varPtr++;
    assert(*varPtr == '5');
    varPtr++;
    assert(*varPtr == 0);
    varPtr++;
    assert(*varPtr == 'N');
    varPtr++;
    assert(*varPtr == 'M');
    varPtr++;
    assert(*varPtr == 'C');
    varPtr++;
    assert(*varPtr == 2);
    varPtr++;
    assert(*varPtr == 'X');
    varPtr++;
    assert(*varPtr == 'T');
    varPtr++;
    assert(*varPtr == 'A');
    varPtr++;
    assert(*varPtr == 'R');
    varPtr++;
}


void validateRead6(SamRecord& samRecord)
{
    //////////////////////////////////////////
    // Validate Record 6
    // Create record structure for validating.
    int expectedBlockSize = 77;
    const char* expectedReferenceName = "1";
    const char* expectedMateReferenceName = "18";
    const char* expectedMateReferenceNameOrEqual = "18";

    bamRecordStruct* expectedRecordPtr =
        (bamRecordStruct *) malloc(expectedBlockSize + sizeof(int));

    char tag[3];
    char type;
    void* value;
    bamRecordStruct* bufferPtr;
    unsigned char* varPtr;

    expectedRecordPtr->myBlockSize = expectedBlockSize;
    expectedRecordPtr->myReferenceID = 0;
    expectedRecordPtr->myPosition = TestValidate::READ6_POS;
    expectedRecordPtr->myReadNameLength = 21;
    expectedRecordPtr->myMapQuality = 0;
    expectedRecordPtr->myBin = 4681;
    expectedRecordPtr->myCigarLength = 3;
    expectedRecordPtr->myFlag = 97;
    expectedRecordPtr->myReadLength = 8;
    expectedRecordPtr->myMateReferenceID = 17;
    expectedRecordPtr->myMatePosition = 756;
    expectedRecordPtr->myInsertSize = 0;
   
    // Check the accessors.
    assert(samRecord.getBlockSize() == expectedRecordPtr->myBlockSize);
    assert(samRecord.getReferenceID() == expectedRecordPtr->myReferenceID);
    assert(strcmp(samRecord.getReferenceName(), expectedReferenceName) == 0);
    assert(samRecord.get1BasedPosition() == expectedRecordPtr->myPosition + 1);
    assert(samRecord.get0BasedPosition() == expectedRecordPtr->myPosition);
    assert(samRecord.getReadNameLength() == 
           expectedRecordPtr->myReadNameLength);
    assert(samRecord.getMapQuality() == expectedRecordPtr->myMapQuality);
    assert(samRecord.getBin() == expectedRecordPtr->myBin);
    assert(samRecord.getCigarLength() == expectedRecordPtr->myCigarLength);
    assert(samRecord.getFlag() == expectedRecordPtr->myFlag);
    assert(samRecord.getReadLength() == expectedRecordPtr->myReadLength);
    assert(samRecord.getMateReferenceID() ==
           expectedRecordPtr->myMateReferenceID);
    assert(strcmp(samRecord.getMateReferenceName(),
                  expectedMateReferenceName) == 0);
    assert(strcmp(samRecord.getMateReferenceNameOrEqual(),
                  expectedMateReferenceNameOrEqual) == 0);
    assert(samRecord.get1BasedMatePosition() ==
           expectedRecordPtr->myMatePosition + 1);
    assert(samRecord.get0BasedMatePosition() ==
           expectedRecordPtr->myMatePosition);
    assert(samRecord.getInsertSize() == expectedRecordPtr->myInsertSize);
    assert(strcmp(samRecord.getReadName(), "18:462+29M5I3M:F:296") == 0);
    assert(samRecord.getCigar() == TestValidate::READ6_CIGAR);
    assert(samRecord.getSequence() == TestValidate::READ6_SEQ);
    assert(samRecord.getQuality() == TestValidate::READ6_QUAL);
    assert(samRecord.getNumOverlaps(1750, 1755) == 5);
    assert(samRecord.getNumOverlaps(1750, 1754) == 4);
    assert(samRecord.getNumOverlaps(0, 2000) == 5);
    assert(samRecord.getNumOverlaps(1749, 1755) == 5);
    assert(samRecord.getNumOverlaps(1751, 1755) == 4);
    assert(samRecord.getNumOverlaps(0, 1752) == 2);
    assert(samRecord.getNumOverlaps(0, 19) == 0);

    // Check the alignment end
    assert(samRecord.get0BasedAlignmentEnd() == TestValidate::READ6_ALIGN_END);
    assert(samRecord.get1BasedAlignmentEnd() == (TestValidate::READ6_ALIGN_END + 1));
    assert(samRecord.getAlignmentLength() == TestValidate::READ6_ALIGN_LEN);
    assert(samRecord.get0BasedUnclippedStart() == TestValidate::READ6_UNCLIP_START);
    assert(samRecord.get1BasedUnclippedStart() == (TestValidate::READ6_UNCLIP_START + 1));
    assert(samRecord.get0BasedUnclippedEnd() == TestValidate::READ6_UNCLIP_END);
    assert(samRecord.get1BasedUnclippedEnd() == (TestValidate::READ6_UNCLIP_END + 1));

    // No tags.
    assert(samRecord.getNextSamTag(tag, type, &value) == false);

    // Get the record ptr.   
    bufferPtr = (bamRecordStruct*)samRecord.getRecordBuffer();
    // Validate the buffers match.
    assert(bufferPtr->myBlockSize == expectedRecordPtr->myBlockSize);
    assert(bufferPtr->myReferenceID == expectedRecordPtr->myReferenceID);
    assert(bufferPtr->myPosition == expectedRecordPtr->myPosition);
    assert(bufferPtr->myReadNameLength == expectedRecordPtr->myReadNameLength);
    assert(bufferPtr->myMapQuality == expectedRecordPtr->myMapQuality);
    assert(bufferPtr->myBin == expectedRecordPtr->myBin);
    assert(bufferPtr->myCigarLength == expectedRecordPtr->myCigarLength);
    assert(bufferPtr->myFlag == expectedRecordPtr->myFlag);
    assert(bufferPtr->myReadLength == expectedRecordPtr->myReadLength);
    assert(bufferPtr->myMateReferenceID ==
           expectedRecordPtr->myMateReferenceID);
    assert(bufferPtr->myMatePosition == expectedRecordPtr->myMatePosition);
    assert(bufferPtr->myInsertSize == expectedRecordPtr->myInsertSize);

    // Validate the variable length fields in the buffer.
    // Set the pointer to the start of the variable fields.
    varPtr = (unsigned char*)(&(bufferPtr->myData[0]));

    // Validate the readname.
    for(int i = 0; i < expectedRecordPtr->myReadNameLength; i++)
    {
        assert(*varPtr == samRecord.getReadName()[i]);
        varPtr++;
    }

    // Validate the cigar.
    // The cigar is 3S2H5M which is:
    // 3S: 3 << 4 | 4 = 0x34
    assert(*(unsigned int*)varPtr == 0x34);
    // Increment the varptr the size of an int.
    varPtr += 4;
    // 2H: 2 << 4 | 5 = 0x25
    assert(*(unsigned int*)varPtr == 0x25);
    // Increment the varptr the size of an int.
    varPtr += 4;
    // 5M: 5 << 4 | 0 = 0x50
    assert(*(unsigned int*)varPtr == 0x50);
    // Increment the varptr the size of an int.
    varPtr += 4;
   
    // Validate the sequence.
    // TG = 0x84
    assert(*varPtr == 0x84);
    varPtr++;
    // CA = 0x21
    assert(*varPtr == 0x21);
    varPtr++;
    // CG = 0x24
    assert(*varPtr == 0x24);
    varPtr++;
    // TN = 0x8F
    assert(*varPtr == 0x8F);
    varPtr++;

    // Validate the Quality
    for(int i = 0; i < expectedRecordPtr->myReadLength; i++)
    {
        assert(*varPtr == samRecord.getQuality()[i] - 33);
        varPtr++;
    }
}


void validateRead7(SamRecord& samRecord)
{
    //////////////////////////////////////////
    // Validate Record 7
    // Create record structure for validating.
    int expectedBlockSize = 83;
    const char* expectedReferenceName = "2";
    const char* expectedMateReferenceName = "18";
    const char* expectedMateReferenceNameOrEqual = "18";

    bamRecordStruct* expectedRecordPtr =
        (bamRecordStruct *) malloc(expectedBlockSize + sizeof(int));

    char tag[3];
    char type;
    void* value;
    bamRecordStruct* bufferPtr;
    unsigned char* varPtr;

    expectedRecordPtr->myBlockSize = expectedBlockSize;
    expectedRecordPtr->myReferenceID = 1;
    expectedRecordPtr->myPosition = TestValidate::READ7_POS;
    expectedRecordPtr->myReadNameLength = 21;
    expectedRecordPtr->myMapQuality = 0;
    expectedRecordPtr->myBin = 4681;
    expectedRecordPtr->myCigarLength = 4;
    expectedRecordPtr->myFlag = 97;
    expectedRecordPtr->myReadLength = 9;
    expectedRecordPtr->myMateReferenceID = 17;
    expectedRecordPtr->myMatePosition = 756;
    expectedRecordPtr->myInsertSize = 0;
   
    // Check the accessors.
    assert(samRecord.getBlockSize() == expectedRecordPtr->myBlockSize);
    assert(samRecord.getReferenceID() == expectedRecordPtr->myReferenceID);
    assert(strcmp(samRecord.getReferenceName(), expectedReferenceName) == 0);
    assert(samRecord.get1BasedPosition() == expectedRecordPtr->myPosition + 1);
    assert(samRecord.get0BasedPosition() == expectedRecordPtr->myPosition);
    assert(samRecord.getReadNameLength() ==
           expectedRecordPtr->myReadNameLength);
    assert(samRecord.getMapQuality() == expectedRecordPtr->myMapQuality);
    assert(samRecord.getBin() == expectedRecordPtr->myBin);
    assert(samRecord.getCigarLength() == expectedRecordPtr->myCigarLength);
    assert(samRecord.getFlag() == expectedRecordPtr->myFlag);
    assert(samRecord.getReadLength() == expectedRecordPtr->myReadLength);
    assert(samRecord.getMateReferenceID() ==
           expectedRecordPtr->myMateReferenceID);
    assert(strcmp(samRecord.getMateReferenceName(),
                  expectedMateReferenceName) == 0);
    assert(strcmp(samRecord.getMateReferenceNameOrEqual(),
                  expectedMateReferenceNameOrEqual) == 0);
    assert(samRecord.get1BasedMatePosition() ==
           expectedRecordPtr->myMatePosition + 1);
    assert(samRecord.get0BasedMatePosition() ==
           expectedRecordPtr->myMatePosition);
    assert(samRecord.getInsertSize() == expectedRecordPtr->myInsertSize);
    assert(strcmp(samRecord.getReadName(), "18:462+29M5I3M:F:297") == 0);
    assert(samRecord.getCigar() == TestValidate::READ7_CIGAR);
    assert(samRecord.getSequence() == TestValidate::READ7_SEQ);
    assert(samRecord.getQuality() == TestValidate::READ7_QUAL);
    assert(samRecord.getNumOverlaps(1750, 1755) == 5);
    assert(samRecord.getNumOverlaps(1750, 1754) == 4);
    assert(samRecord.getNumOverlaps(0, 2000) == 5);
    assert(samRecord.getNumOverlaps(1749, 1755) == 5);
    assert(samRecord.getNumOverlaps(1751, 1755) == 4);
    assert(samRecord.getNumOverlaps(0, 1752) == 2);
    assert(samRecord.getNumOverlaps(0, 19) == 0);

    // Check the alignment end
    assert(samRecord.get0BasedAlignmentEnd() == TestValidate::READ7_ALIGN_END);
    assert(samRecord.get1BasedAlignmentEnd() == (TestValidate::READ7_ALIGN_END + 1));
    assert(samRecord.getAlignmentLength() == TestValidate::READ7_ALIGN_LEN);
    assert(samRecord.get0BasedUnclippedStart() == TestValidate::READ7_UNCLIP_START);
    assert(samRecord.get1BasedUnclippedStart() == (TestValidate::READ7_UNCLIP_START + 1));
    assert(samRecord.get0BasedUnclippedEnd() == TestValidate::READ7_UNCLIP_END);
    assert(samRecord.get1BasedUnclippedEnd() == (TestValidate::READ7_UNCLIP_END + 1));

    // No tags.
    assert(samRecord.getNextSamTag(tag, type, &value) == false);

    // Get the record ptr.   
    bufferPtr = (bamRecordStruct*)samRecord.getRecordBuffer();
    // Validate the buffers match.
    assert(bufferPtr->myBlockSize == expectedRecordPtr->myBlockSize);
    assert(bufferPtr->myReferenceID == expectedRecordPtr->myReferenceID);
    assert(bufferPtr->myPosition == expectedRecordPtr->myPosition);
    assert(bufferPtr->myReadNameLength == expectedRecordPtr->myReadNameLength);
    assert(bufferPtr->myMapQuality == expectedRecordPtr->myMapQuality);
    assert(bufferPtr->myBin == expectedRecordPtr->myBin);
    assert(bufferPtr->myCigarLength == expectedRecordPtr->myCigarLength);
    assert(bufferPtr->myFlag == expectedRecordPtr->myFlag);
    assert(bufferPtr->myReadLength == expectedRecordPtr->myReadLength);
    assert(bufferPtr->myMateReferenceID ==
           expectedRecordPtr->myMateReferenceID);
    assert(bufferPtr->myMatePosition == expectedRecordPtr->myMatePosition);
    assert(bufferPtr->myInsertSize == expectedRecordPtr->myInsertSize);

    // Validate the variable length fields in the buffer.
    // Set the pointer to the start of the variable fields.
    varPtr = (unsigned char*)(&(bufferPtr->myData[0]));

    // Validate the readname.
    for(int i = 0; i < expectedRecordPtr->myReadNameLength; i++)
    {
        assert(*varPtr == samRecord.getReadName()[i]);
        varPtr++;
    }

    // Validate the cigar.
    // The cigar is 3S5M1S3H which is:
    // 3S: 3 << 4 | 4 = 0x34
    assert(*(unsigned int*)varPtr == 0x34);
    // Increment the varptr the size of an int.
    varPtr += 4;
    // 5M: 5 << 4 | 0 = 0x50
    assert(*(unsigned int*)varPtr == 0x50);
    // Increment the varptr the size of an int.
    varPtr += 4;
    // 1S: 1 << 4 | 4 = 0x14
    assert(*(unsigned int*)varPtr == 0x14);
    // Increment the varptr the size of an int.
    varPtr += 4;
    // 3H: 3 << 4 | 5 = 0x35
    assert(*(unsigned int*)varPtr == 0x35);
    // Increment the varptr the size of an int.
    varPtr += 4;
   
    // Validate the sequence.
    // TG = 0x84
    assert(*varPtr == 0x84);
    varPtr++;
    // CA = 0x21
    assert(*varPtr == 0x21);
    varPtr++;
    // CG = 0x24
    assert(*varPtr == 0x24);
    varPtr++;
    // TN = 0x8F
    assert(*varPtr == 0x8F);
    varPtr++;
    // G  = 0x40
    assert(*varPtr == 0x40);
    varPtr++;

    // Validate the Quality
    for(int i = 0; i < expectedRecordPtr->myReadLength; i++)
    {
        assert(*varPtr == samRecord.getQuality()[i] - 33);
        varPtr++;
    }
}


void validateRead8(SamRecord& samRecord)
{
    //////////////////////////////////////////
    // Validate Record 8
    // Create record structure for validating.
    int expectedBlockSize = 65;
    const char* expectedReferenceName = "*";
    const char* expectedMateReferenceName = "*";
    const char* expectedMateReferenceNameOrEqual = "*";

    bamRecordStruct* expectedRecordPtr =
        (bamRecordStruct *) malloc(expectedBlockSize + sizeof(int));

    char tag[3];
    char type;
    void* value;
    bamRecordStruct* bufferPtr;
    unsigned char* varPtr;

    expectedRecordPtr->myBlockSize = expectedBlockSize;
    expectedRecordPtr->myReferenceID = -1;
    expectedRecordPtr->myPosition = -1;
    expectedRecordPtr->myReadNameLength = 27;
    expectedRecordPtr->myMapQuality = 0;
    expectedRecordPtr->myBin = 4680;
    expectedRecordPtr->myCigarLength = 0;
    expectedRecordPtr->myFlag = 141;
    expectedRecordPtr->myReadLength = 4;
    expectedRecordPtr->myMateReferenceID = -1;
    expectedRecordPtr->myMatePosition = -1;
    expectedRecordPtr->myInsertSize = 0;
   
    // Check the alignment end
    assert(samRecord.get0BasedAlignmentEnd() == -1);
    assert(samRecord.get1BasedAlignmentEnd() == 0);
    assert(samRecord.getAlignmentLength() == 0);
    assert(samRecord.get0BasedUnclippedStart() == -1);
    assert(samRecord.get1BasedUnclippedStart() == 0);
    assert(samRecord.get0BasedUnclippedEnd() == -1);
    assert(samRecord.get1BasedUnclippedEnd() == 0);

    // Check the accessors.
    assert(samRecord.getBlockSize() == expectedRecordPtr->myBlockSize);
    assert(samRecord.getReferenceID() == expectedRecordPtr->myReferenceID);
    assert(strcmp(samRecord.getReferenceName(), expectedReferenceName) == 0);
    assert(samRecord.get1BasedPosition() == expectedRecordPtr->myPosition + 1);
    assert(samRecord.get0BasedPosition() == expectedRecordPtr->myPosition);
    assert(samRecord.getReadNameLength() ==
           expectedRecordPtr->myReadNameLength);
    assert(samRecord.getMapQuality() == expectedRecordPtr->myMapQuality);
    assert(samRecord.getBin() == expectedRecordPtr->myBin);
    assert(samRecord.getCigarLength() == expectedRecordPtr->myCigarLength);
    assert(samRecord.getFlag() == expectedRecordPtr->myFlag);
    assert(samRecord.getReadLength() == expectedRecordPtr->myReadLength);
    assert(samRecord.getMateReferenceID() ==
           expectedRecordPtr->myMateReferenceID);
    assert(strcmp(samRecord.getMateReferenceName(),
                  expectedMateReferenceName) == 0);
    assert(strcmp(samRecord.getMateReferenceNameOrEqual(),
                  expectedMateReferenceNameOrEqual) == 0);
    assert(samRecord.get1BasedMatePosition() == 
           expectedRecordPtr->myMatePosition + 1);
    assert(samRecord.get0BasedMatePosition() ==
           expectedRecordPtr->myMatePosition);
    assert(samRecord.getInsertSize() == expectedRecordPtr->myInsertSize);
    assert(strcmp(samRecord.getReadName(), "Y:16597235+13M13I11M:F:181") == 0);
    assert(strcmp(samRecord.getCigar(), "*") == 0);
    assert(strcmp(samRecord.getSequence(), "AACT") == 0);
    assert(strcmp(samRecord.getQuality(), "==;;") == 0);
    assert(samRecord.getNumOverlaps(1750, 1755) == 0);
    assert(samRecord.getNumOverlaps(1750, 1754) == 0);
    assert(samRecord.getNumOverlaps(0, 2000) == 0);
    assert(samRecord.getNumOverlaps(1749, 1755) == 0);
    assert(samRecord.getNumOverlaps(1751, 1755) == 0);
    assert(samRecord.getNumOverlaps(0, 1752) == 0);
    assert(samRecord.getNumOverlaps(0, 19) == 0);
    assert(samRecord.getNumOverlaps(-1, 4) == 0);

    // No Tags to check, should return false.
    assert(samRecord.getNextSamTag(tag, type, &value) == false);
    assert(samRecord.getNextSamTag(tag, type, &value) == false);

    // Get the record ptr.   
    bufferPtr = (bamRecordStruct*)samRecord.getRecordBuffer();
    // Validate the buffers match.
    assert(bufferPtr->myBlockSize == expectedRecordPtr->myBlockSize);
    assert(bufferPtr->myReferenceID == expectedRecordPtr->myReferenceID);
    assert(bufferPtr->myPosition == expectedRecordPtr->myPosition);
    assert(bufferPtr->myReadNameLength == expectedRecordPtr->myReadNameLength);
    assert(bufferPtr->myMapQuality == expectedRecordPtr->myMapQuality);
    assert(bufferPtr->myBin == expectedRecordPtr->myBin);
    assert(bufferPtr->myCigarLength == expectedRecordPtr->myCigarLength);
    assert(bufferPtr->myFlag == expectedRecordPtr->myFlag);
    assert(bufferPtr->myReadLength == expectedRecordPtr->myReadLength);
    assert(bufferPtr->myMateReferenceID ==
           expectedRecordPtr->myMateReferenceID);
    assert(bufferPtr->myMatePosition == expectedRecordPtr->myMatePosition);
    assert(bufferPtr->myInsertSize == expectedRecordPtr->myInsertSize);

    // Validate the variable length fields in the buffer.
    // Set the pointer to the start of the variable fields.
    varPtr = (unsigned char*)(&(bufferPtr->myData[0]));

    // Validate the readname.
    for(int i = 0; i < expectedRecordPtr->myReadNameLength; i++)
    {
        assert(*varPtr == samRecord.getReadName()[i]);
        varPtr++;
    }

    // No cigar to validate. 
    // Validate the sequence.
    // AA = 0x11
    assert(*varPtr == 0x11);
    varPtr++;
    // CT = 0x28
    assert(*varPtr == 0x28);
    varPtr++;

    // Validate the Quality
    for(int i = 0; i < expectedRecordPtr->myReadLength; i++)
    {
        assert(*varPtr == samRecord.getQuality()[i] - 33);
        varPtr++;
    }

    // No tags.  
}


void validateRead9(SamRecord& samRecord)
{
    //////////////////////////////////////////
    // Validate Record 9
    // Create record structure for validating.
    int expectedBlockSize = 77;
    const char* expectedReferenceName = "3";
    const char* expectedMateReferenceName = "18";
    const char* expectedMateReferenceNameOrEqual = "18";

    bamRecordStruct* expectedRecordPtr =
        (bamRecordStruct *) malloc(expectedBlockSize + sizeof(int));

    char tag[3];
    char type;
    void* value;
    bamRecordStruct* bufferPtr;
    unsigned char* varPtr;

    expectedRecordPtr->myBlockSize = expectedBlockSize;
    expectedRecordPtr->myReferenceID = 2;
    expectedRecordPtr->myPosition = 74;
    expectedRecordPtr->myReadNameLength = 21;
    expectedRecordPtr->myMapQuality = 0;
    expectedRecordPtr->myBin = 4681;
    expectedRecordPtr->myCigarLength = 3;
    expectedRecordPtr->myFlag = 97;
    expectedRecordPtr->myReadLength = 8;
    expectedRecordPtr->myMateReferenceID = 17;
    expectedRecordPtr->myMatePosition = 756;
    expectedRecordPtr->myInsertSize = 0;
   
    // Check the accessors.
    assert(samRecord.getBlockSize() == expectedRecordPtr->myBlockSize);
    assert(samRecord.getStatus() == SamStatus::SUCCESS);
    assert(samRecord.getReferenceID() == expectedRecordPtr->myReferenceID);
    assert(strcmp(samRecord.getReferenceName(), expectedReferenceName) == 0);
    assert(samRecord.get1BasedPosition() == expectedRecordPtr->myPosition + 1);
    assert(samRecord.get0BasedPosition() == expectedRecordPtr->myPosition);
    assert(samRecord.getReadNameLength() ==
           expectedRecordPtr->myReadNameLength);
    assert(samRecord.getMapQuality() == expectedRecordPtr->myMapQuality);
    assert(samRecord.getBin() == expectedRecordPtr->myBin);
    assert(samRecord.getCigarLength() == expectedRecordPtr->myCigarLength);
    assert(samRecord.getFlag() == expectedRecordPtr->myFlag);
    assert(samRecord.getReadLength() == expectedRecordPtr->myReadLength);
    assert(samRecord.getMateReferenceID() ==
           expectedRecordPtr->myMateReferenceID);
    assert(strcmp(samRecord.getMateReferenceName(),
                  expectedMateReferenceName) == 0);
    assert(strcmp(samRecord.getMateReferenceNameOrEqual(),
                  expectedMateReferenceNameOrEqual) == 0);
    assert(samRecord.get1BasedMatePosition() ==
           expectedRecordPtr->myMatePosition + 1);
    assert(samRecord.get0BasedMatePosition() ==
           expectedRecordPtr->myMatePosition);
    assert(samRecord.getInsertSize() == expectedRecordPtr->myInsertSize);
    assert(strcmp(samRecord.getReadName(), "18:462+29M5I3M:F:298") == 0);
    assert(strcmp(samRecord.getCigar(), "3S5M4H") == 0);
    assert((strcmp(samRecord.getSequence(), "TGCACGTN") == 0) || 
           (strcmp(samRecord.getSequence(), "tgcacgtn") == 0));
    assert(strcmp(samRecord.getQuality(), "453;>>>>") == 0);
    assert(samRecord.getNumOverlaps(74, 79) == 5);
    assert(samRecord.getNumOverlaps(73, 79) == 5);
    assert(samRecord.getNumOverlaps(75, 78) == 3);
    assert(samRecord.getNumOverlaps(0, 1017) == 5);
    assert(samRecord.getNumOverlaps(79, 85) == 0);
    assert(samRecord.getNumOverlaps(78, 85) == 1);
    assert(samRecord.getNumOverlaps(-1, 1017) == 5);

    // Check the alignment end
    assert(samRecord.get0BasedAlignmentEnd() == 78);
    assert(samRecord.get1BasedAlignmentEnd() == 79);
    assert(samRecord.getAlignmentLength() == 5);
    assert(samRecord.get0BasedUnclippedStart() == 71);
    assert(samRecord.get1BasedUnclippedStart() == 72);
    assert(samRecord.get0BasedUnclippedEnd() == 82);
    assert(samRecord.get1BasedUnclippedEnd() == 83);

    // No tags.
    assert(samRecord.getNextSamTag(tag, type, &value) == false);

    // Get the record ptr.   
    bufferPtr = (bamRecordStruct*)samRecord.getRecordBuffer();
    assert(bufferPtr != NULL);
    // Validate the buffers match.
    assert(bufferPtr->myBlockSize == expectedRecordPtr->myBlockSize);
    assert(bufferPtr->myReferenceID == expectedRecordPtr->myReferenceID);
    assert(bufferPtr->myPosition == expectedRecordPtr->myPosition);
    assert(bufferPtr->myReadNameLength == expectedRecordPtr->myReadNameLength);
    assert(bufferPtr->myMapQuality == expectedRecordPtr->myMapQuality);
    assert(bufferPtr->myBin == expectedRecordPtr->myBin);
    assert(bufferPtr->myCigarLength == expectedRecordPtr->myCigarLength);
    assert(bufferPtr->myFlag == expectedRecordPtr->myFlag);
    assert(bufferPtr->myReadLength == expectedRecordPtr->myReadLength);
    assert(bufferPtr->myMateReferenceID == 
           expectedRecordPtr->myMateReferenceID);
    assert(bufferPtr->myMatePosition == expectedRecordPtr->myMatePosition);
    assert(bufferPtr->myInsertSize == expectedRecordPtr->myInsertSize);

    // Validate the variable length fields in the buffer.
    // Set the pointer to the start of the variable fields.
    varPtr = (unsigned char*)(&(bufferPtr->myData[0]));

    // Validate the readname.
    for(int i = 0; i < expectedRecordPtr->myReadNameLength; i++)
    {
        assert(*varPtr == samRecord.getReadName()[i]);
        varPtr++;
    }

    // Validate the cigar.
    // The cigar is 3S5M1S3H which is:
    // 3S: 3 << 4 | 4 = 0x34
    assert(*(unsigned int*)varPtr == 0x34);
    // Increment the varptr the size of an int.
    varPtr += 4;
    // 5M: 5 << 4 | 0 = 0x50
    assert(*(unsigned int*)varPtr == 0x50);
    // Increment the varptr the size of an int.
    varPtr += 4;
    // 4H: 4 << 4 | 5 = 0x45
    assert(*(unsigned int*)varPtr == 0x45);
    // Increment the varptr the size of an int.
    varPtr += 4;
   
    // Validate the sequence.
    // TG = 0x84
    assert(*varPtr == 0x84);
    varPtr++;
    // CA = 0x21
    assert(*varPtr == 0x21);
    varPtr++;
    // CG = 0x24
    assert(*varPtr == 0x24);
    varPtr++;
    // TN = 0x8F
    assert(*varPtr == 0x8F);
    varPtr++;

    // Validate the Quality
    for(int i = 0; i < expectedRecordPtr->myReadLength; i++)
    {
        assert(*varPtr == samRecord.getQuality()[i] - 33);
        varPtr++;
    }
}


void validateRead10(SamRecord& samRecord)
{
    //////////////////////////////////////////
    // Validate Record 10
    // Create record structure for validating.
    int expectedBlockSize = 59;
    const char* expectedReferenceName = "*";
    const char* expectedMateReferenceName = "*";
    const char* expectedMateReferenceNameOrEqual = "*";

    bamRecordStruct* expectedRecordPtr =
        (bamRecordStruct *) malloc(expectedBlockSize + sizeof(int));

    char tag[3];
    char type;
    void* value;
    bamRecordStruct* bufferPtr;
    unsigned char* varPtr;

    expectedRecordPtr->myBlockSize = expectedBlockSize;
    expectedRecordPtr->myReferenceID = -1;
    expectedRecordPtr->myPosition = -1;
    expectedRecordPtr->myReadNameLength = 27;
    expectedRecordPtr->myMapQuality = 0;
    expectedRecordPtr->myBin = 4680;
    expectedRecordPtr->myCigarLength = 0;
    expectedRecordPtr->myFlag = 141;
    expectedRecordPtr->myReadLength = 0;
    expectedRecordPtr->myMateReferenceID = -1;
    expectedRecordPtr->myMatePosition = -1;
    expectedRecordPtr->myInsertSize = 0;
   
    // Check the alignment end
    assert(samRecord.get0BasedUnclippedStart() == -1);
    assert(samRecord.get1BasedUnclippedStart() == 0);
    assert(samRecord.get0BasedUnclippedEnd() == -1);
    assert(samRecord.get1BasedUnclippedEnd() == 0);
    assert(samRecord.get1BasedAlignmentEnd() == 0);
    assert(samRecord.get0BasedAlignmentEnd() == -1);
    assert(samRecord.getAlignmentLength() == 0);

    // Check the accessors.
    assert(samRecord.getBlockSize() == expectedRecordPtr->myBlockSize);
    assert(samRecord.getReferenceID() == expectedRecordPtr->myReferenceID);
    assert(strcmp(samRecord.getReferenceName(), expectedReferenceName) == 0);
    assert(samRecord.get1BasedPosition() == expectedRecordPtr->myPosition + 1);
    assert(samRecord.get0BasedPosition() == expectedRecordPtr->myPosition);
    assert(samRecord.getReadNameLength() == 
           expectedRecordPtr->myReadNameLength);
    assert(samRecord.getMapQuality() == expectedRecordPtr->myMapQuality);
    assert(samRecord.getBin() == expectedRecordPtr->myBin);
    assert(samRecord.getCigarLength() == expectedRecordPtr->myCigarLength);
    assert(samRecord.getFlag() == expectedRecordPtr->myFlag);
    assert(samRecord.getReadLength() == expectedRecordPtr->myReadLength);
    assert(samRecord.getMateReferenceID() == 
           expectedRecordPtr->myMateReferenceID);
    assert(strcmp(samRecord.getMateReferenceName(), 
                  expectedMateReferenceName) == 0);
    assert(strcmp(samRecord.getMateReferenceNameOrEqual(),
                  expectedMateReferenceNameOrEqual) == 0);
    assert(samRecord.get1BasedMatePosition() == 
           expectedRecordPtr->myMatePosition + 1);
    assert(samRecord.get0BasedMatePosition() == 
           expectedRecordPtr->myMatePosition);
    assert(samRecord.getInsertSize() == expectedRecordPtr->myInsertSize);
    assert(strcmp(samRecord.getReadName(), "Y:16597235+13M13I11M:F:181") == 0);
    assert(strcmp(samRecord.getCigar(), "*") == 0);
    assert(strcmp(samRecord.getSequence(), "*") == 0);
    assert(strcmp(samRecord.getQuality(), "*") == 0);
    assert(samRecord.getNumOverlaps(74, 79) == 0);
    assert(samRecord.getNumOverlaps(73, 79) == 0);
    assert(samRecord.getNumOverlaps(75, 78) == 0);
    assert(samRecord.getNumOverlaps(0, 1017) == 0);
    assert(samRecord.getNumOverlaps(79, 85) == 0);
    assert(samRecord.getNumOverlaps(78, 85) == 0);
    assert(samRecord.getNumOverlaps(-1, 1017) == 0);

    // No Tags to check, should return false.
    assert(samRecord.getNextSamTag(tag, type, &value) == false);
    assert(samRecord.getNextSamTag(tag, type, &value) == false);

    // Get the record ptr.   
    bufferPtr = (bamRecordStruct*)samRecord.getRecordBuffer();
    // Validate the buffers match.
    assert(bufferPtr->myBlockSize == expectedRecordPtr->myBlockSize);
    assert(bufferPtr->myReferenceID == expectedRecordPtr->myReferenceID);
    assert(bufferPtr->myPosition == expectedRecordPtr->myPosition);
    assert(bufferPtr->myReadNameLength == expectedRecordPtr->myReadNameLength);
    assert(bufferPtr->myMapQuality == expectedRecordPtr->myMapQuality);
    assert(bufferPtr->myBin == expectedRecordPtr->myBin);
    assert(bufferPtr->myCigarLength == expectedRecordPtr->myCigarLength);
    assert(bufferPtr->myFlag == expectedRecordPtr->myFlag);
    assert(bufferPtr->myReadLength == expectedRecordPtr->myReadLength);
    assert(bufferPtr->myMateReferenceID == 
           expectedRecordPtr->myMateReferenceID);
    assert(bufferPtr->myMatePosition == expectedRecordPtr->myMatePosition);
    assert(bufferPtr->myInsertSize == expectedRecordPtr->myInsertSize);

    // Validate the variable length fields in the buffer.
    // Set the pointer to the start of the variable fields.
    varPtr = (unsigned char*)(&(bufferPtr->myData[0]));

    // Validate the readname.
    for(int i = 0; i < expectedRecordPtr->myReadNameLength; i++)
    {
        assert(*varPtr == samRecord.getReadName()[i]);
        varPtr++;
    }

    // No cigar to validate. 
    // No sequence.
    // No Quality.
    // No Tags.
}


void validateHeader(SamFileHeader& samHeader)
{
    validateHeaderFields(samHeader);
    validateHeaderString(samHeader);
}


void validateHeaderFields(SamFileHeader& samHeader)
{
    const char* value;

    ////////////////////////////////////////////////////////
    // Test getting a specific HD Tag value from the header
    // that does not exist.
    value = samHeader.getHDTagValue("GO");
    assert(strcmp(value, "") == 0);

    ////////////////////////////////////////////////////////
    // Test getting a specific PG Tag value from the header
    // that does not exist.
    value = samHeader.getPGTagValue("CL", "1");
    assert(strcmp(value, "") == 0);

    ////////////////////////////////////////////////////////
    // Test getting a specific SQ Tag value from the header
    value = samHeader.getSQTagValue("LN", "1");
    assert(value != NULL);
    assert(strcmp(value, "247249719") == 0);
    value = samHeader.getSQTagValue("LN", "22");
    assert(value != NULL);
    assert(strcmp(value, "49691432") == 0);

    ////////////////////////////////////////////////////////
    // Test getting a specific SQ Tag value from the header
    // that does not exist.
    value = samHeader.getSQTagValue("LN", "1000");
    assert(strcmp(value, "") == 0);

    ////////////////////////////////////////////////////////
    // Test getting a specific SQ Tag value from the header
    // that does not exist - sq exists, but not with that tag.
    value = samHeader.getSQTagValue("AS", "1");
    assert(strcmp(value, "") == 0);

    ////////////////////////////////////////////////////////
    // Test getting a specific RG Tag value from the header
    value = samHeader.getRGTagValue("LB", "myID2");
    assert(value != NULL);
    assert(strcmp(value, "library2") == 0);
    value = samHeader.getRGTagValue("LB", "myID");
    assert(value != NULL);
    assert(strcmp(value, "library") == 0);

    ////////////////////////////////////////////////////////
    // Test getting a specific SQ from the header
    // Then pulling the tags out of it.
    SamHeaderSQ* sq = samHeader.getSQ("10");
    assert(strcmp(sq->getTagValue("SN"), "10") == 0);
    assert(strcmp(sq->getTagValue("LN"), "135374737") == 0);
   
    // Test pulling a tag that does not exist.
    assert(strcmp(sq->getTagValue("DD"), "") == 0);
   

    ////////////////////////////////////////////////////////
    // Test getting a specific RG from the header
    // Then pulling the tags out of it.
    const SamHeaderRG* rg = samHeader.getRG("myID");
    assert(strcmp(rg->getTagValue("ID"), "myID") == 0);
    assert(strcmp(rg->getTagValue("SM"), "sample") == 0);
    assert(strcmp(rg->getTagValue("LB"), "library") == 0);
   
    // Test pulling a tag that does not exist.
    assert(strcmp(rg->getTagValue("DD"), "") == 0);
   
    ////////////////////////////////////////////////////////
    // Test getting a specific RG from the header that does not exist.
    rg = samHeader.getRG("noExist");
    assert(rg == NULL);

    ////////////////////////////////////////////////////////
    // Test getting a specific SQ from the header that does not exist.
    sq = samHeader.getSQ("noExist");
    assert(sq == NULL);

    ////////////////////////////////////////////////////////
    // Test getting the reference ID.
    assert(samHeader.getReferenceID("2") == 1);
    std::string refIDStdString = "X";
    assert(samHeader.getReferenceID(refIDStdString.c_str()) == 22);
    String refIDString = "22";
    assert(samHeader.getReferenceID(refIDString) == 21);
    assert(samHeader.getReferenceID(refIDString.c_str()) == 21);
    assert(samHeader.getReferenceID("Z") == SamReferenceInfo::NO_REF_ID);
    assert(samHeader.getReferenceID("Z", true) == 23);
    assert(samHeader.getReferenceID("*") == -1);
    refIDString = "*";
    assert(samHeader.getReferenceID(refIDString) == -1);
    assert(samHeader.getReferenceID(refIDString.c_str()) == -1);
}

void validateHeaderString(SamFileHeader& samHeader)
{
    // Check the header line.
    std::string headerString = "";
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tSN:1\tLN:247249719\n@SQ\tSN:2\tLN:242951149\n@SQ\tSN:3\tLN:199501827\n@SQ\tSN:4\tLN:191273063\n@SQ\tSN:5\tLN:180857866\n@SQ\tSN:6\tLN:170899992\n@SQ\tSN:7\tLN:158821424\n@SQ\tSN:8\tLN:146274826\n@SQ\tSN:9\tLN:140273252\n@SQ\tSN:10\tLN:135374737\n@SQ\tSN:11\tLN:134452384\n@SQ\tSN:12\tLN:132349534\n@SQ\tSN:13\tLN:114142980\n@SQ\tSN:14\tLN:106368585\n@SQ\tSN:15\tLN:100338915\n@SQ\tSN:16\tLN:88827254\n@SQ\tSN:17\tLN:78774742\n@SQ\tSN:18\tLN:76117153\n@SQ\tSN:19\tLN:63811651\n@SQ\tSN:20\tLN:62435964\n@SQ\tSN:21\tLN:46944323\n@SQ\tSN:22\tLN:49691432\n@SQ\tSN:X\tLN:154913754\n@RG\tID:myID\tLB:library\tSM:sample\n@RG\tID:myID2\tSM:sample2\tLB:library2\n@CO\tComment 1\n@CO\tComment 2\n");
}
