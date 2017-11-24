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
#include "ModifyVar.h"
#include "SamValidation.h"

void modifyFirstBase()
{
    SamFile samIn;
    // Open the file for reading.   
    assert(samIn.OpenForRead("testFiles/testVar.bam"));

    SamFile samOut;
    // Open the file for writing.   
    assert(samOut.OpenForWrite("results/updateVar.bam"));

    // Read the sam header.
    SamFileHeader samHeader;
    assert(samIn.ReadHeader(samHeader));

    assert(samOut.WriteHeader(samHeader));

    // Read the sam records.
    SamRecord samRecord;

    // Keep reading records until the end of the file is reached.
    while(samIn.ReadRecord(samHeader, samRecord))
    {
        // Successfully read a record from the file, so check to see
        // if it is valid.
        SamValidationErrors samValidationErrors;
        assert(SamValidator::isValid(samHeader, samRecord, 
                                     samValidationErrors));
      
        // Get the sequence.
        const char* sequence = samRecord.getSequence();
        assert(strcmp(sequence, "") != 0);
        std::string upSeq = sequence;
        upSeq[0] = 'N';
        assert(samRecord.setSequence(upSeq.c_str()));
      
        // write the sequence.
        assert(samOut.WriteRecord(samHeader, samRecord));
    }

    // Should have exited only when done reading.
    assert(samIn.GetStatus() == SamStatus::NO_MORE_RECS);
}


void modifyFirstBaseLong()
{
    SamFile samIn;
    // Open the file for reading.
    assert(samIn.OpenForRead("statusTests/HalfG.bam"));

    SamFile samOut;
    // Open the file for writing.   
    assert(samOut.OpenForWrite("results/updateSeq.bam"));

    // Read the sam header.
    SamFileHeader samHeader;
    assert(samIn.ReadHeader(samHeader));

    assert(samOut.WriteHeader(samHeader));

    // Read the sam records.
    SamRecord samRecord;

    // Keep reading records until the end of the file is reached.
    while(samIn.ReadRecord(samHeader, samRecord))
    {
        // Successfully read a record from the file, so check to see
        // if it is valid.
        SamValidationErrors samValidationErrors;
        assert(SamValidator::isValid(samHeader, samRecord, 
                                     samValidationErrors));
      
        // Get the sequence.
        const char* sequence = samRecord.getSequence();
        assert(strcmp(sequence, "") != 0);
        std::string upSeq = sequence;
        upSeq[0] = 'N';
        assert(samRecord.setSequence(upSeq.c_str()));
      
        // write the sequence.
        assert(samOut.WriteRecord(samHeader, samRecord));
    }

    // Should have exited only when done reading.
    assert(samIn.GetStatus() == SamStatus::NO_MORE_RECS);
}


void testModifyVar()
{
#ifdef __ZLIB_AVAILABLE__
    modifyFirstBase();
#endif
    modifyVar modTest;
    modTest.testModifyVar("testFiles/testSam.sam", true);
    modTest.testModifyVar("testFiles/testSam.sam", false);
#ifdef __ZLIB_AVAILABLE__
    modTest.testModifyVar("testFiles/testBam.bam", true);
    modTest.testModifyVar("testFiles/testBam.bam", false);
#endif
}


void modifyVar::testModifyVar(const char* filename, bool valBufFirst)
{
    myFilename = filename;
    myValBufFirst = valBufFirst;

    testModifyReadNameOnlySameLength();
    testModifyCigarOnlySameLength();
    testModifySequenceOnlySameLength();
    testModifyQualityOnlySameLength();
    testRemoveQuality();
    testShortenQuality();
    testLengthenQuality();

    testShortenReadName();
    testShortenCigar();
    testShortenSequence();

    testLengthenReadName();
    testLengthenCigar();
    testLengthenSequence();
   
    testRemoveCigar();
    testRemoveSequence();
    testLengthenSequenceAndQuality();
}


void modifyVar::testModifyReadNameOnlySameLength()
{
    resetExpected();
    openAndRead1Rec();

    // Set the Read Name - same length, just different name.
    expectedReadNameString = "1:1011:G:255+17M15D20M";
    samRecord.setReadName(expectedReadNameString.c_str());

    validate();
}


void modifyVar::testModifyCigarOnlySameLength()
{
    resetExpected();
    openAndRead1Rec();

    // Set the Cigar - same length, just different values.
    expectedCigarString = "3M2I";
    samRecord.setCigar(expectedCigarString.c_str());

    // The new Cigar for record 1 is 3M2I
    // 3M = 3 << 4 | 0 = 0x30
    // 2I = 2 << 4 | 1 = 0x21
    expectedCigarBufLen = 2;
    expectedCigarBuffer[0] = 0x30;
    expectedCigarBuffer[1] = 0x21;
   
    validate();
}


void modifyVar::testModifySequenceOnlySameLength()
{
    resetExpected();
    openAndRead1Rec();

    // Set the Sequence - same length, just different values.
    expectedSequenceString = "NCGAN";
    samRecord.setSequence(expectedSequenceString.c_str());

    // NCGAN = NC GA N = 0xF2 0x41 0xF0
    expectedSequenceBuffer[0] = 0xF2;
    expectedSequenceBuffer[1] = 0x41;
    expectedSequenceBuffer[2] = 0xF0;

    validate();
}


void modifyVar::testModifyQualityOnlySameLength()
{
    resetExpected();
    openAndRead1Rec();
   
    // Set the Quality - same length, just different values.
    expectedQualityString = "!>6+!";
    samRecord.setQuality(expectedQualityString.c_str());

    validate();
}


void modifyVar::testRemoveQuality()
{
    resetExpected();
    openAndRead1Rec();
   
    // Set the Quality - to "*" - does not affect the length since the 
    // sequence field drives the length.
    expectedQualityString = "*";
    samRecord.setQuality(expectedQualityString.c_str());

    validate();
}


void modifyVar::testShortenQuality()
{
    resetExpected();
    openAndRead1Rec();

    // Set the Quality - shorten, but doesn't affect the length since
    // the sequence field drives the length.
    expectedQualityString = "!!";
    samRecord.setQuality(expectedQualityString.c_str());

    validate();
}


void modifyVar::testLengthenQuality()
{
    resetExpected();
    openAndRead1Rec();

    // Set the Quality - lengthen, but doesn't affect the length since
    // the sequence field drives the length.
    expectedQualityString = "!!@@##";
    samRecord.setQuality(expectedQualityString.c_str());

    validate();
}


void modifyVar::testShortenReadName()
{
    resetExpected();
    openAndRead1Rec();

    // Set the Read Name - shorter length
    expectedReadNameString = "1:1011:G:255";
    samRecord.setReadName(expectedReadNameString.c_str());

    validate();
}


void modifyVar::testShortenCigar()
{
    resetExpected();
    openAndRead1Rec();

    // Set the Cigar - shorter length
    expectedCigarString = "5M";
    samRecord.setCigar(expectedCigarString.c_str());

    // The new Cigar for record 1 is 5M
    // 5M = 5 << 4 | 0 = 0x50
    expectedCigarBufLen = 1;
    expectedCigarBuffer[0] = 0x50;
   
    validate();
}


void modifyVar::testShortenSequence()
{
    resetExpected();
    openAndRead1Rec();

    // Set the Sequence - shorter length.
    expectedSequenceString = "CCGA";
    samRecord.setSequence(expectedSequenceString.c_str());

    // CCGA = CC GA = 0x22 0x41
    expectedSequenceBuffer[0] = 0x22;
    expectedSequenceBuffer[1] = 0x41;

    validate();
}


void modifyVar::testLengthenReadName()
{
    resetExpected();
    openAndRead1Rec();

    // Set the Read Name - longer.
    expectedReadNameString = "1:1011:G:255+17M15D20M:1111111";
    samRecord.setReadName(expectedReadNameString.c_str());

    validate();
}


void modifyVar::testLengthenCigar()
{
    resetExpected();
    openAndRead1Rec();

    // Set the Cigar - longer length.
    expectedCigarString = "3M2D2I";
    samRecord.setCigar(expectedCigarString.c_str());

    // The new Cigar for record 1 is 3M2I
    // 3M = 3 << 4 | 0 = 0x30
    // 2D = 2 << 2 | 1 = 0x22
    // 2I = 2 << 4 | 1 = 0x21
    expectedCigarBufLen = 3;
    expectedCigarBuffer[0] = 0x30;
    expectedCigarBuffer[1] = 0x22;
    expectedCigarBuffer[2] = 0x21;
   
    validate();
}


void modifyVar::testLengthenSequence()
{
    resetExpected();
    openAndRead1Rec();

    // Set the Sequence - longer length.
    expectedSequenceString = "CCGAATT";
    samRecord.setSequence(expectedSequenceString.c_str());

    // CCGAATT = CC GA AT T = 0x22 0x41 0x18 0x80
    expectedSequenceBuffer[0] = 0x22;
    expectedSequenceBuffer[1] = 0x41;
    expectedSequenceBuffer[2] = 0x18;
    expectedSequenceBuffer[3] = 0x80;

    validate();
}


void modifyVar::testRemoveCigar()
{
    resetExpected();
    openAndRead1Rec();

    // Set the Cigar - same length, just different values.
    expectedCigarString = "*";
    expectedCigarBufLen = 0;
    samRecord.setCigar(expectedCigarString.c_str());

    validate();
}


void modifyVar::testRemoveSequence()
{
    resetExpected();
    openAndRead1Rec();

    // Set the Sequence - shorter length.
    expectedSequenceString = "*";
    samRecord.setSequence(expectedSequenceString.c_str());

    validate();
}


void modifyVar::testLengthenSequenceAndQuality()
{
    resetExpected();
    openAndRead1Rec();

    // Set the Sequence & quality - longer.
    expectedSequenceString = "CCGAATT";
    expectedQualityString = "!@#$%^&";
    samRecord.setSequence(expectedSequenceString.c_str());
    samRecord.setQuality(expectedQualityString.c_str());

    // CCGAATT = CC GA AT T = 0x22 0x41 0x18 0x80
    expectedSequenceBuffer[0] = 0x22;
    expectedSequenceBuffer[1] = 0x41;
    expectedSequenceBuffer[2] = 0x18;
    expectedSequenceBuffer[3] = 0x80;

    validate();
}


void modifyVar::validate()
{
    if(myValBufFirst)
    {
        // get the record.
        const bamRecordStruct* recordBuffer = 
            (const bamRecordStruct*)samRecord.getRecordBuffer();
      
        // Validate the buffer.
        validateReadName(recordBuffer);
        validateCigar(recordBuffer);
        validateSequence(recordBuffer);
        validateQuality(recordBuffer);
        validateTags(recordBuffer);
      
        // Validate the strings.
        validateReadNameString();
        validateCigarString();
        validateSequenceString();
        validateQualityString();
        validateTagsString();
    }
    else
    {
        // get the record.
        const bamRecordStruct* recordBuffer = 
            (const bamRecordStruct*)samRecord.getRecordBuffer();
      
        // Validate the buffer.
        validateReadName(recordBuffer);
        validateCigar(recordBuffer);
        validateSequence(recordBuffer);
        validateQuality(recordBuffer);
        validateTags(recordBuffer);
      
        // Validate the strings.
        validateReadNameString();
        validateCigarString();
        validateSequenceString();
        validateQualityString();
        validateTagsString();
    }
}

void modifyVar::validateReadName(const bamRecordStruct* recordBuffer)
{
    const char* varPtr = (const char*)&(recordBuffer->myData);

    unsigned int len = expectedReadNameString.length();
    for(unsigned int i = 0; i < len; i++)
    {
        assert(varPtr[i] == expectedReadNameString[i]);
    }

    // Verify ending null.
    assert(varPtr[len] == 0);
   
    // verify the length - add one for the terminating null.
    assert(recordBuffer->myReadNameLength == 
           expectedReadNameString.length() + 1);
}


void modifyVar::validateCigar(const bamRecordStruct* recordBuffer)
{
    const unsigned char* cigarStart = 
        (const unsigned char*)&(recordBuffer->myData) 
        + recordBuffer->myReadNameLength;

    unsigned int* varPtr = (unsigned int*)cigarStart;
   
    for(int i = 0; i < expectedCigarBufLen; i++)
    {
        assert(varPtr[i] == expectedCigarBuffer[i]);
    }
   
    assert(recordBuffer->myCigarLength == expectedCigarBufLen);
}


void modifyVar::validateSequence(const bamRecordStruct* recordBuffer)
{
    // Calculate the sequence length.
    int expectedReadLen = expectedSequenceString.length();
    int seqLen = (expectedReadLen + 1)/2;
    if(expectedSequenceString == "*")
    {
        expectedReadLen = 0;
        seqLen = 0;
    }
   
    const unsigned char* sequencePtr = 
        (const unsigned char*)&(recordBuffer->myData) 
        + recordBuffer->myReadNameLength + (recordBuffer->myCigarLength * 4);

    for(int i = 0; i < seqLen; i++)
    {
        assert(sequencePtr[i] == expectedSequenceBuffer[i]);
    }

    assert(recordBuffer->myReadLength == expectedReadLen);
}


void modifyVar::validateQuality(const bamRecordStruct* recordBuffer)
{
    int expectedReadLen = expectedSequenceString.length();
    int seqLen = (expectedReadLen + 1)/2;
    if(expectedSequenceString == "*")
    {
        expectedReadLen = 0;
        seqLen = 0;
    }

    const uint8_t* qualityPtr = 
        (const unsigned char*)&(recordBuffer->myData) 
        + recordBuffer->myReadNameLength + (recordBuffer->myCigarLength * 4)
        + seqLen;
 
    int qualityLen = expectedQualityString.length();

    for(int i = 0; i < expectedReadLen; i++)
    {
        if(expectedQualityString == "*")
        {
            // no quality, so check for 0xFF.
            assert(qualityPtr[i] == 0xFF);
        }
        else if(i >= qualityLen)
        {
            // Quality is shorter than the sequence, so should be padded with
            // 0xFF.
            assert(qualityPtr[i] == 0xFF);
        }
        else
        {
            assert(qualityPtr[i] == (expectedQualityString[i] - 33));
        }
    }

    assert(recordBuffer->myReadLength == expectedReadLen);
}


void modifyVar::validateTags(const bamRecordStruct* recordBuffer)
{
    const unsigned char* tagsPtr = 
        (const unsigned char*)&(recordBuffer->myData) 
        + recordBuffer->myReadNameLength + (recordBuffer->myCigarLength * 4)
        + (recordBuffer->myReadLength + 1)/2 + recordBuffer->myReadLength;
 
    for(int i = 0; i < expectedTagsLen; i++)
    {
        assert(tagsPtr[i] == expectedTagsBuffer[i]);
    }

    // Calculate expected block size - from the start of the buffer to the
    // start of the tags plus the tags length - minus the size of the blocksize
    // field.
    int32_t expectedBlockSize = tagsPtr - (const unsigned char*)(recordBuffer) 
        + expectedTagsLen - 4;
    assert(recordBuffer->myBlockSize == expectedBlockSize);
}


void modifyVar::validateTagsString()
{
    char tag[3];
    char type;
    void* value;
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
}

void modifyVar::validateReadNameString()
{
    assert(samRecord.getReadName() == expectedReadNameString);
}

void modifyVar::validateCigarString()
{
    assert(samRecord.getCigar() == expectedCigarString);
}

void modifyVar::validateSequenceString()
{
    assert(samRecord.getSequence() == expectedSequenceString);
}

void modifyVar::validateQualityString()
{
    assert(samRecord.getQuality() == expectedQualityString);
}


void modifyVar::resetExpected()
{
    expectedReadNameString = "1:1011:F:255+17M15D20M";
    expectedCigarString = "5M2D";
    expectedSequenceString = "CCGAA";
    expectedQualityString = "6>6+4";

    // The default Cigar for record 1 is 5M2D
    // 5M = 5 << 4 | 0 = 0x50
    // 2D = 2 << 4 | 2 = 0x22
    expectedCigarBufLen = 2;
    expectedCigarBuffer[0] = 0x50;
    expectedCigarBuffer[1] = 0x22;

    // CCGAA = CC GA A = 0x22 0x41 0x10
    expectedSequenceBuffer[0] = 0x22;
    expectedSequenceBuffer[1] = 0x41;
    expectedSequenceBuffer[2] = 0x10;

    expectedTagsLen = 18;
    expectedTagsBuffer[0] = 'A';
    expectedTagsBuffer[1] =  'M';
    expectedTagsBuffer[2] =  'C';
    expectedTagsBuffer[3] =  0;
    expectedTagsBuffer[4] =  'M';
    expectedTagsBuffer[5] =  'D';
    expectedTagsBuffer[6] =  'Z';
    expectedTagsBuffer[7] =  '3';
    expectedTagsBuffer[8] =  '7';
    expectedTagsBuffer[9] =  0;
    expectedTagsBuffer[10] =  'N';
    expectedTagsBuffer[11] =  'M';
    expectedTagsBuffer[12] =  'C';
    expectedTagsBuffer[13] =  0;
    expectedTagsBuffer[14] =  'X';
    expectedTagsBuffer[15] =  'T';
    expectedTagsBuffer[16] =  'A';
    expectedTagsBuffer[17] =  'R';
}


void modifyVar::openAndRead1Rec()
{
    // Open the file for reading.   
    assert(samIn.OpenForRead(myFilename));

    // Read the sam header.
    assert(samIn.ReadHeader(samHeader));
   
    // Read the first record.   
    assert(samIn.ReadRecord(samHeader, samRecord));
}
