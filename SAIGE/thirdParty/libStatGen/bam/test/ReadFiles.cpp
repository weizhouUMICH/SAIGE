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
#include "TestValidate.h"
#include "SamTags.h"
#include <assert.h>

void testReadSam()
{
    SamFileHeader samHdr;
    SamFile inSam1("testFiles/testSam1.sam", SamFile::READ, &samHdr);
    assert(!inSam1.IsEOF());
    SamRecord rec;
    assert(inSam1.ReadRecord(samHdr, rec) == true);
    assert(inSam1.IsEOF());

    SamFile inSam;
    assert(inSam.OpenForRead("testFiles/testSam.sam"));

    // Call generic test which since the sam and bam are identical, should
    // contain the same results.
    testRead(inSam);

    inSam.Close();

    testFlagRead("testFiles/testSam.sam");
}

void testReadBam()
{
    SamFileHeader samHdr;
    SamFile inSam1("testFiles/testBam1.bam", SamFile::READ, &samHdr);
    assert(!inSam1.IsEOF());
    SamRecord rec;
    assert(inSam1.ReadRecord(samHdr, rec) == true);
    assert(inSam1.IsEOF());

    SamFile inSam;
    assert(inSam.OpenForRead("testFiles/testBam.bam"));

    // Call generic test which since the sam and bam are identical, should
    // contain the same results.
    testRead(inSam);

    inSam.Close();

    testFlagRead("testFiles/testBam.bam");
}

void testRead(SamFile &inSam)
{
    // Read the SAM Header.
    SamFileHeader samHeader;
    assert(inSam.ReadHeader(samHeader));
    assert(!inSam.IsEOF());

    validateHeader(samHeader);

    testCopyHeader(samHeader);    

    testModHeader(samHeader);

    SamRecord samRecord;
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead1(samRecord);

    // Set a new quality and get the buffer.
    samRecord.setQuality("ABCDE");
    validateRead1ModQuality(samRecord);
    //   void* buffer = samRecord.getRecordBuffer();

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead2(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead3(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead4(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead5(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead6(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead7(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead8(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead9(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead10(samRecord);
}


void testAddHeaderAndTagToFile(const char* inputName, const char* outputName)
{
    SamFile inSam, outSam;
    assert(inSam.OpenForRead(inputName));
    assert(outSam.OpenForWrite(outputName));

    // Read the SAM Header.
    SamFileHeader samHeader;
    assert(inSam.ReadHeader(samHeader));

    // Add a header line.
    assert(samHeader.addHeaderLine("@RG\tID:myID\tSM:mySM") == false);
    assert(samHeader.addHeaderLine("@RG\tID:myID3\tSM:mySM") == true);

    // Write Header
    assert(outSam.WriteHeader(samHeader));

    SamRecord samRecord;
    assert(inSam.ReadRecord(samHeader, samRecord));
    //   validateRead1(samRecord);
    // Add two tags.
    assert(samRecord.addIntTag("XA", 123));
    assert(samRecord.addIntTag("XA", 456));
    assert(samRecord.addTag("RR", 'Z', "myID1"));
    assert(samRecord.addTag("RR", 'Z', "myID2"));

    // Write as Sam.
    assert(outSam.WriteRecord(samHeader, samRecord));

    // TODO, add test to verify it was written correctly.

    // Read a couple of records to make sure it properly can read them even
    // if they are bigger than the original.
    assert(inSam.ReadRecord(samHeader, samRecord));
    assert(inSam.ReadRecord(samHeader, samRecord));

    //  Check the MD tag, which requires the reference.
    GenomeSequence reference("testFiles/chr1_partial.fa");
    assert(SamTags::isMDTagCorrect(samRecord, reference) == false);
    String newMDTag;
    SamTags::createMDTag(newMDTag, samRecord, reference);
    assert(newMDTag == "2T1N0");
    assert(SamTags::updateMDTag(samRecord, reference));
    // Write as Sam.
    assert(outSam.WriteRecord(samHeader, samRecord));
}


// Test reading a file, validating it is sorted.
void testValidateSortedRead()
{
    // Open a file for reading.
    SamFile inSam(ErrorHandler::RETURN);
    assert(inSam.OpenForRead("testFiles/testSam.sam"));

    // Set the validation to COORDINATE.
    inSam.setSortedValidation(SamFile::COORDINATE);

    // Read the SAM Header.
    SamFileHeader samHeader;
    assert(inSam.ReadHeader(samHeader));
    
    SamRecord samRecord;
    // Succeed, first record.
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead1(samRecord);

    // Succeed, higher coordinate.
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead2(samRecord);

    // Failed sort order - due to coord.
    assert(inSam.ReadRecord(samHeader, samRecord) == false);
    validateRead3(samRecord);

    // Failed sort order - due to coord.
    assert(inSam.ReadRecord(samHeader, samRecord) == false);
    validateRead4(samRecord);

    // Succeed, new reference id
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead5(samRecord);

    // Fail, previous reference id.
    assert(inSam.ReadRecord(samHeader, samRecord) == false);
    validateRead6(samRecord);

    // Succeed, same reference id, higher coord.
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead7(samRecord);

    // Succeed, *, new reference id.
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead8(samRecord);

    // Fail, reference id is not *
    assert(inSam.ReadRecord(samHeader, samRecord) == false);
    validateRead9(samRecord);

    // Succeed, valid reference id, and no coordinate.
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead10(samRecord);


    ////////////////////////////////////////////
    // Reopen the file for reading
    assert(inSam.OpenForRead("testFiles/testSam.sam"));

    // Set the validation to QUERY_NAME.
    inSam.setSortedValidation(SamFile::QUERY_NAME);

    // Read the SAM Header.
    assert(inSam.ReadHeader(samHeader));
   
    // Succeed, first record.
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead1(samRecord);

    // Succeed, same name.
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead2(samRecord);

    // Succeeds - numeric sort
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead3(samRecord);

    // Succeeds - numeric sort
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead4(samRecord);

    // Succeeds - numeric sort
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead5(samRecord);

    // Succeeds - numeric sort
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead6(samRecord);

    // Succeeds - numeric sort
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead7(samRecord);

    // Succeed - std sort
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead8(samRecord);

    // Succeed - numeric sort (Y<18)
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead9(samRecord);

    // Succeed - std sort
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead10(samRecord);

    ////////////////////////////////////////////
    // Reopen the file for reading
    assert(inSam.OpenForRead("testFiles/testSam.sam"));

    // Set the validation to the SO Flag.  Not set, so it is UNSORTED, so 
    // all reads should pass.
    inSam.setSortedValidation(SamFile::FLAG);

    // Read the SAM Header.
    assert(inSam.ReadHeader(samHeader));
   
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead1(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead2(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead3(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead4(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead5(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead6(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead7(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead8(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead9(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead10(samRecord);

    ////////////////////////////////////////////
    // Reopen for reading SO FLAG set to coordinate.
    assert(inSam.OpenForRead("testFiles/testSamSOcoord.sam"));

    // Set the validation to SO FLAG which is set to coordinate.
    inSam.setSortedValidation(SamFile::FLAG);

    // Read the SAM Header.
    assert(inSam.ReadHeader(samHeader));

    // Succeed, first record.
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead1(samRecord);

    // Succeed, higher coordinate.
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead2(samRecord);

    // Failed sort order - due to coord.
    assert(inSam.ReadRecord(samHeader, samRecord) == false);
    validateRead3(samRecord);

    // Failed sort order - due to coord.
    assert(inSam.ReadRecord(samHeader, samRecord) == false);
    validateRead4(samRecord);

    // Succeed, new reference id
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead5(samRecord);

    // Fail, previous reference id.
    assert(inSam.ReadRecord(samHeader, samRecord) == false);
    validateRead6(samRecord);

    // Succeed, same reference id, higher coord.
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead7(samRecord);

    // Succeed, *, new reference id.
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead8(samRecord);

    // Fail, reference id is not *
    assert(inSam.ReadRecord(samHeader, samRecord) == false);
    validateRead9(samRecord);

    // Succeed, valid reference id, and no coordinate.
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead10(samRecord);


    ////////////////////////////////////////////
    // Reopen the file for reading
    assert(inSam.OpenForRead("testFiles/testSamSOquery.sam"));

    // Set the validation to FLAG, SO set to queryname.
    inSam.setSortedValidation(SamFile::FLAG);

    // Read the SAM Header.
    assert(inSam.ReadHeader(samHeader));
   
    // Succeed, first record.
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead1(samRecord);

    // Succeed, same name.
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead2(samRecord);

    // Succeeds - numeric sort
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead3(samRecord);

    // Succeeds - numeric sort
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead4(samRecord);

    // Succeeds - numeric sort
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead5(samRecord);

    // Succeeds - numeric sort
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead6(samRecord);

    // Succeeds - numeric sort
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead7(samRecord);

    // Succeed - std sort
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead8(samRecord);

    // Succeed - numeric sort (Y<18)
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead9(samRecord);

    // Succeed - std sort
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead10(samRecord);

    ////////////////////////////////////////////
    // Reopen the file for reading, SO flag set to junk.
    assert(inSam.OpenForRead("testFiles/testSamSOinvalid.sam"));

    // Set the validation to the SO Flag.  Not set to anything valid,
    // so it is considered UNSORTED, so all reads should pass.
    inSam.setSortedValidation(SamFile::FLAG);

    // Read the SAM Header.
    assert(inSam.ReadHeader(samHeader));
   
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead1(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead2(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead3(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead4(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead5(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead6(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead7(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead8(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead9(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead10(samRecord);
}



void validateRead1ModQuality(SamRecord& samRecord)
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
    bamRecordStruct* bufferPtr;
    unsigned char* varPtr;

    expectedRecordPtr->myBlockSize = expectedBlockSize;
    expectedRecordPtr->myReferenceID = 0;
    expectedRecordPtr->myPosition = 1010;
    expectedRecordPtr->myReadNameLength = 23;
    expectedRecordPtr->myMapQuality = 0;
    expectedRecordPtr->myBin = 4681;
    expectedRecordPtr->myCigarLength = 2;
    expectedRecordPtr->myFlag = 73;
    expectedRecordPtr->myReadLength = 5;
    expectedRecordPtr->myMateReferenceID = 0;
    expectedRecordPtr->myMatePosition = 1010;
    expectedRecordPtr->myInsertSize = 0;
   
    // Check the alignment end
    assert(samRecord.get0BasedAlignmentEnd() == 1016);
    assert(samRecord.get1BasedAlignmentEnd() == 1017);
    assert(samRecord.getAlignmentLength() == 7);
    assert(samRecord.get0BasedUnclippedStart() == 1010);
    assert(samRecord.get1BasedUnclippedStart() == 1011);
    assert(samRecord.get0BasedUnclippedEnd() == 1016);
    assert(samRecord.get1BasedUnclippedEnd() == 1017);

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
    assert(strcmp(samRecord.getCigar(), "5M2D") == 0);
    assert(strcmp(samRecord.getSequence(), "CCGAA") == 0);
    assert(strcmp(samRecord.getQuality(), "ABCDE") == 0);
    assert(samRecord.getNumOverlaps(1010, 1017) == 5);
    assert(samRecord.getNumOverlaps(1010, 1016) == 5);
    assert(samRecord.getNumOverlaps(1012, 1017) == 3);
    assert(samRecord.getNumOverlaps(1015, 1017) == 0);
    assert(samRecord.getNumOverlaps(1017, 1010) == 0);
    assert(samRecord.getNumOverlaps(1013, 1011) == 0);
    assert(samRecord.getNumOverlaps(-1, 1017) == 5);

    // Reset the tag iter, since the tags have already been read.
    samRecord.resetTagIter();

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


void testModHeader(SamFileHeader& samHeader)
{
    // Check the header line.
    std::string headerString = "";
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tSN:1\tLN:247249719\n@SQ\tSN:2\tLN:242951149\n@SQ\tSN:3\tLN:199501827\n@SQ\tSN:4\tLN:191273063\n@SQ\tSN:5\tLN:180857866\n@SQ\tSN:6\tLN:170899992\n@SQ\tSN:7\tLN:158821424\n@SQ\tSN:8\tLN:146274826\n@SQ\tSN:9\tLN:140273252\n@SQ\tSN:10\tLN:135374737\n@SQ\tSN:11\tLN:134452384\n@SQ\tSN:12\tLN:132349534\n@SQ\tSN:13\tLN:114142980\n@SQ\tSN:14\tLN:106368585\n@SQ\tSN:15\tLN:100338915\n@SQ\tSN:16\tLN:88827254\n@SQ\tSN:17\tLN:78774742\n@SQ\tSN:18\tLN:76117153\n@SQ\tSN:19\tLN:63811651\n@SQ\tSN:20\tLN:62435964\n@SQ\tSN:21\tLN:46944323\n@SQ\tSN:22\tLN:49691432\n@SQ\tSN:X\tLN:154913754\n@RG\tID:myID\tLB:library\tSM:sample\n@RG\tID:myID2\tSM:sample2\tLB:library2\n@CO\tComment 1\n@CO\tComment 2\n");

    // Remove a tag - by setting it to "".
    assert(samHeader.setRGTag("LB", "", "myID2") == true);


    // Check the header line.
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tSN:1\tLN:247249719\n@SQ\tSN:2\tLN:242951149\n@SQ\tSN:3\tLN:199501827\n@SQ\tSN:4\tLN:191273063\n@SQ\tSN:5\tLN:180857866\n@SQ\tSN:6\tLN:170899992\n@SQ\tSN:7\tLN:158821424\n@SQ\tSN:8\tLN:146274826\n@SQ\tSN:9\tLN:140273252\n@SQ\tSN:10\tLN:135374737\n@SQ\tSN:11\tLN:134452384\n@SQ\tSN:12\tLN:132349534\n@SQ\tSN:13\tLN:114142980\n@SQ\tSN:14\tLN:106368585\n@SQ\tSN:15\tLN:100338915\n@SQ\tSN:16\tLN:88827254\n@SQ\tSN:17\tLN:78774742\n@SQ\tSN:18\tLN:76117153\n@SQ\tSN:19\tLN:63811651\n@SQ\tSN:20\tLN:62435964\n@SQ\tSN:21\tLN:46944323\n@SQ\tSN:22\tLN:49691432\n@SQ\tSN:X\tLN:154913754\n@RG\tID:myID\tLB:library\tSM:sample\n@RG\tID:myID2\tSM:sample2\n@CO\tComment 1\n@CO\tComment 2\n");

    //  Add an HD tag.
    SamHeaderHD* hd = new SamHeaderHD();
    assert(hd->setTag("VN", "1.3") == true);
    assert(samHeader.addHD(hd) == true);
    assert(strcmp(samHeader.getHDTagValue("VN"), "1.3") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tSN:1\tLN:247249719\n@SQ\tSN:2\tLN:242951149\n@SQ\tSN:3\tLN:199501827\n@SQ\tSN:4\tLN:191273063\n@SQ\tSN:5\tLN:180857866\n@SQ\tSN:6\tLN:170899992\n@SQ\tSN:7\tLN:158821424\n@SQ\tSN:8\tLN:146274826\n@SQ\tSN:9\tLN:140273252\n@SQ\tSN:10\tLN:135374737\n@SQ\tSN:11\tLN:134452384\n@SQ\tSN:12\tLN:132349534\n@SQ\tSN:13\tLN:114142980\n@SQ\tSN:14\tLN:106368585\n@SQ\tSN:15\tLN:100338915\n@SQ\tSN:16\tLN:88827254\n@SQ\tSN:17\tLN:78774742\n@SQ\tSN:18\tLN:76117153\n@SQ\tSN:19\tLN:63811651\n@SQ\tSN:20\tLN:62435964\n@SQ\tSN:21\tLN:46944323\n@SQ\tSN:22\tLN:49691432\n@SQ\tSN:X\tLN:154913754\n@RG\tID:myID\tLB:library\tSM:sample\n@RG\tID:myID2\tSM:sample2\n@HD\tVN:1.3\n@CO\tComment 1\n@CO\tComment 2\n");

    // Try adding another HD tag.
    SamHeaderHD* hd2 = new SamHeaderHD();
    assert(hd2->setTag("VN", "1.4") == true);
    assert(samHeader.addHD(hd2) == false);
    assert(strcmp(samHeader.getHDTagValue("VN"), "1.4") != 0);
    assert(strcmp(samHeader.getHDTagValue("VN"), "1.3") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tSN:1\tLN:247249719\n@SQ\tSN:2\tLN:242951149\n@SQ\tSN:3\tLN:199501827\n@SQ\tSN:4\tLN:191273063\n@SQ\tSN:5\tLN:180857866\n@SQ\tSN:6\tLN:170899992\n@SQ\tSN:7\tLN:158821424\n@SQ\tSN:8\tLN:146274826\n@SQ\tSN:9\tLN:140273252\n@SQ\tSN:10\tLN:135374737\n@SQ\tSN:11\tLN:134452384\n@SQ\tSN:12\tLN:132349534\n@SQ\tSN:13\tLN:114142980\n@SQ\tSN:14\tLN:106368585\n@SQ\tSN:15\tLN:100338915\n@SQ\tSN:16\tLN:88827254\n@SQ\tSN:17\tLN:78774742\n@SQ\tSN:18\tLN:76117153\n@SQ\tSN:19\tLN:63811651\n@SQ\tSN:20\tLN:62435964\n@SQ\tSN:21\tLN:46944323\n@SQ\tSN:22\tLN:49691432\n@SQ\tSN:X\tLN:154913754\n@RG\tID:myID\tLB:library\tSM:sample\n@RG\tID:myID2\tSM:sample2\n@HD\tVN:1.3\n@CO\tComment 1\n@CO\tComment 2\n");

    // Remove the entire HD Tag.
    assert(samHeader.removeHD() == true);
    assert(strcmp(samHeader.getHDTagValue("VN"), "") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tSN:1\tLN:247249719\n@SQ\tSN:2\tLN:242951149\n@SQ\tSN:3\tLN:199501827\n@SQ\tSN:4\tLN:191273063\n@SQ\tSN:5\tLN:180857866\n@SQ\tSN:6\tLN:170899992\n@SQ\tSN:7\tLN:158821424\n@SQ\tSN:8\tLN:146274826\n@SQ\tSN:9\tLN:140273252\n@SQ\tSN:10\tLN:135374737\n@SQ\tSN:11\tLN:134452384\n@SQ\tSN:12\tLN:132349534\n@SQ\tSN:13\tLN:114142980\n@SQ\tSN:14\tLN:106368585\n@SQ\tSN:15\tLN:100338915\n@SQ\tSN:16\tLN:88827254\n@SQ\tSN:17\tLN:78774742\n@SQ\tSN:18\tLN:76117153\n@SQ\tSN:19\tLN:63811651\n@SQ\tSN:20\tLN:62435964\n@SQ\tSN:21\tLN:46944323\n@SQ\tSN:22\tLN:49691432\n@SQ\tSN:X\tLN:154913754\n@RG\tID:myID\tLB:library\tSM:sample\n@RG\tID:myID2\tSM:sample2\n@CO\tComment 1\n@CO\tComment 2\n");

    // Remove an entire SQ Tag.
    assert(strcmp(samHeader.getSQTagValue("LN", "11"), "134452384") == 0);
    assert(samHeader.removeSQ("11") == true);
    assert(strcmp(samHeader.getSQTagValue("LN", "11"), "") == 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tSN:1\tLN:247249719\n@SQ\tSN:2\tLN:242951149\n@SQ\tSN:3\tLN:199501827\n@SQ\tSN:4\tLN:191273063\n@SQ\tSN:5\tLN:180857866\n@SQ\tSN:6\tLN:170899992\n@SQ\tSN:7\tLN:158821424\n@SQ\tSN:8\tLN:146274826\n@SQ\tSN:9\tLN:140273252\n@SQ\tSN:10\tLN:135374737\n@SQ\tSN:12\tLN:132349534\n@SQ\tSN:13\tLN:114142980\n@SQ\tSN:14\tLN:106368585\n@SQ\tSN:15\tLN:100338915\n@SQ\tSN:16\tLN:88827254\n@SQ\tSN:17\tLN:78774742\n@SQ\tSN:18\tLN:76117153\n@SQ\tSN:19\tLN:63811651\n@SQ\tSN:20\tLN:62435964\n@SQ\tSN:21\tLN:46944323\n@SQ\tSN:22\tLN:49691432\n@SQ\tSN:X\tLN:154913754\n@RG\tID:myID\tLB:library\tSM:sample\n@RG\tID:myID2\tSM:sample2\n@CO\tComment 1\n@CO\tComment 2\n");

    // Try adding a null HD tag.
    hd = NULL;
    assert(samHeader.addHD(hd) == false);
    assert(strcmp(samHeader.getHDTagValue("VN"), "") == 0);
    assert(strcmp(samHeader.getHDTagValue("VN"), "1.4") != 0);
    assert(strcmp(samHeader.getHDTagValue("VN"), "1.3") != 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tSN:1\tLN:247249719\n@SQ\tSN:2\tLN:242951149\n@SQ\tSN:3\tLN:199501827\n@SQ\tSN:4\tLN:191273063\n@SQ\tSN:5\tLN:180857866\n@SQ\tSN:6\tLN:170899992\n@SQ\tSN:7\tLN:158821424\n@SQ\tSN:8\tLN:146274826\n@SQ\tSN:9\tLN:140273252\n@SQ\tSN:10\tLN:135374737\n@SQ\tSN:12\tLN:132349534\n@SQ\tSN:13\tLN:114142980\n@SQ\tSN:14\tLN:106368585\n@SQ\tSN:15\tLN:100338915\n@SQ\tSN:16\tLN:88827254\n@SQ\tSN:17\tLN:78774742\n@SQ\tSN:18\tLN:76117153\n@SQ\tSN:19\tLN:63811651\n@SQ\tSN:20\tLN:62435964\n@SQ\tSN:21\tLN:46944323\n@SQ\tSN:22\tLN:49691432\n@SQ\tSN:X\tLN:154913754\n@RG\tID:myID\tLB:library\tSM:sample\n@RG\tID:myID2\tSM:sample2\n@CO\tComment 1\n@CO\tComment 2\n");

    // Try adding a null SQ tag.
    SamHeaderSQ* sq = NULL;
    assert(samHeader.addSQ(sq) == false);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tSN:1\tLN:247249719\n@SQ\tSN:2\tLN:242951149\n@SQ\tSN:3\tLN:199501827\n@SQ\tSN:4\tLN:191273063\n@SQ\tSN:5\tLN:180857866\n@SQ\tSN:6\tLN:170899992\n@SQ\tSN:7\tLN:158821424\n@SQ\tSN:8\tLN:146274826\n@SQ\tSN:9\tLN:140273252\n@SQ\tSN:10\tLN:135374737\n@SQ\tSN:12\tLN:132349534\n@SQ\tSN:13\tLN:114142980\n@SQ\tSN:14\tLN:106368585\n@SQ\tSN:15\tLN:100338915\n@SQ\tSN:16\tLN:88827254\n@SQ\tSN:17\tLN:78774742\n@SQ\tSN:18\tLN:76117153\n@SQ\tSN:19\tLN:63811651\n@SQ\tSN:20\tLN:62435964\n@SQ\tSN:21\tLN:46944323\n@SQ\tSN:22\tLN:49691432\n@SQ\tSN:X\tLN:154913754\n@RG\tID:myID\tLB:library\tSM:sample\n@RG\tID:myID2\tSM:sample2\n@CO\tComment 1\n@CO\tComment 2\n");

    // Try adding an HD tag again.
    assert(samHeader.addHD(hd2) == true);
    assert(strcmp(samHeader.getHDTagValue("VN"), "1.4") == 0);
    assert(strcmp(samHeader.getHDTagValue("VN"), "1.3") != 0);
    assert(samHeader.getHeaderString(headerString) == true);
    assert(headerString == "@SQ\tSN:1\tLN:247249719\n@SQ\tSN:2\tLN:242951149\n@SQ\tSN:3\tLN:199501827\n@SQ\tSN:4\tLN:191273063\n@SQ\tSN:5\tLN:180857866\n@SQ\tSN:6\tLN:170899992\n@SQ\tSN:7\tLN:158821424\n@SQ\tSN:8\tLN:146274826\n@SQ\tSN:9\tLN:140273252\n@SQ\tSN:10\tLN:135374737\n@SQ\tSN:12\tLN:132349534\n@SQ\tSN:13\tLN:114142980\n@SQ\tSN:14\tLN:106368585\n@SQ\tSN:15\tLN:100338915\n@SQ\tSN:16\tLN:88827254\n@SQ\tSN:17\tLN:78774742\n@SQ\tSN:18\tLN:76117153\n@SQ\tSN:19\tLN:63811651\n@SQ\tSN:20\tLN:62435964\n@SQ\tSN:21\tLN:46944323\n@SQ\tSN:22\tLN:49691432\n@SQ\tSN:X\tLN:154913754\n@RG\tID:myID\tLB:library\tSM:sample\n@RG\tID:myID2\tSM:sample2\n@HD\tVN:1.4\n@CO\tComment 1\n@CO\tComment 2\n");


    // TODO Get the comments.

}



void testFlagRead(const char* fileName)
{
    SamFile inSam;
    SamFileHeader samHeader;
    SamRecord samRecord;

    ////////////////////////////////////////////////////////////
    // Required flag 0x48  (only flag 73 matches)
    // Exclude nothing
    assert(inSam.OpenForRead(fileName));
    assert(inSam.ReadHeader(samHeader));
    validateHeader(samHeader);
    inSam.SetReadFlags(0x48, 0x0);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead1(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == false);

    inSam.Close();

    ////////////////////////////////////////////////////////////
    // No required flags.
    // Exclude 0x48.  This leaves just the one read with flag 133.
    assert(inSam.OpenForRead(fileName));
    assert(inSam.ReadHeader(samHeader));
    validateHeader(samHeader);
    inSam.SetReadFlags(0x0, 0x48);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead2(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == false);

    inSam.Close();

    ////////////////////////////////////////////////////////////
    // Required flag 0x40 
    // Exclude 0x48.
    // This will not find any records since the exclude and required conflict.
    assert(inSam.OpenForRead(fileName));
    assert(inSam.ReadHeader(samHeader));
    validateHeader(samHeader);
    inSam.SetReadFlags(0x40, 0x48);

    assert(inSam.ReadRecord(samHeader, samRecord) == false);

    inSam.Close();

    ////////////////////////////////////////////////////////////
    // Required flag 0x4
    // Exclude 0x8.
    // Only finds flag 133.
    assert(inSam.OpenForRead(fileName));
    assert(inSam.ReadHeader(samHeader));
    validateHeader(samHeader);
    inSam.SetReadFlags(0x4, 0x8);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead2(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == false);

    inSam.Close();

     ////////////////////////////////////////////////////////////
    // Required flag 0x4
    // Exclude nothing
    // Finds flags 133 & 141.
    assert(inSam.OpenForRead(fileName));
    assert(inSam.ReadHeader(samHeader));
    validateHeader(samHeader);
    inSam.SetReadFlags(0x4, 0x0);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead2(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead8(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead10(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == false);

    inSam.Close();
}


void testCopyHeader(SamFileHeader& samHeader)
{
    // Copy the header.
    SamFileHeader samHeader2;
    
    SamHeaderRecord* recPtr = samHeader.getNextHeaderRecord();
    while(recPtr != NULL)
    {
        samHeader2.addRecordCopy(*recPtr);
        recPtr = samHeader.getNextHeaderRecord();
    }
    // Add the comments.
    std::string nextComment = samHeader.getNextComment();
    while(nextComment != SamFileHeader::EMPTY_RETURN)
    {
        samHeader2.addComment(nextComment.c_str());
        nextComment = samHeader.getNextComment();
    }
    // Validate the header.
    validateHeader(samHeader2);
}
