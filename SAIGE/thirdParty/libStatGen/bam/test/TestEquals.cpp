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

#include "TestEquals.h"
#include "SamQuerySeqWithRefHelper.h"
#include <assert.h>



void testSeqEquals()
{
    // Call generic test which since the sam and bam are identical, should
    // contain the same results.
    EqualsTest::testEq(EqualsTest::SAM);
#ifdef __ZLIB_AVAILABLE__
    EqualsTest::testEq(EqualsTest::BAM);
#endif
}


const char* EqualsTest::READ_NAMES[] = 
    {"01:====", "02:===X", "03:==X=", "04:==XX", "05:=X==", "06:=X=X",
     "07:=XX=", "08:=XXX", "09:X===", "10:X==X", "11:X=X=", "12:X=XX",
     "13:XX==", "14:XX=X", "15:XXX=", "16:XXXX", "Read:GGCCTA;Ref:CCTA",
     "Read:CCTA;Ref:CCTA", "Read:CCGTxxxC;Ref:CCxTAACC", 
     "Read:CCxxAC;Ref:CCTAACC", "chromNotInRef", "chromNotInRef1"};

const char* EqualsTest::READ_SEQS_BASES[] = 
    {"CCTA", "CCTT", "CCAA", "CCAT", "CTTA", "CTTT", "CTAA", "CTAT", "TCTA", "TCTT", "TCAA", "TCAT", "TTTA", "TTTT", "TTAA", "TTAT", "GGCCTA",
     "CCTA", "CCGTC", "CCAC", "CCTA", "CC=A"};

const char* EqualsTest::READ_SEQS_EQUALS[] = 
    {"====", "===T", "==A=", "==AT", "=T==", "=T=T", "=TA=", "=TAT", "T===", "T==T", "T=A=", "T=AT", "TT==", "TT=T", "TTA=", "TTAT", "GG====",
     "====", "==G==", "====", "CCTA", "CC=A"};

const char* EqualsTest::READ_SEQS_MIXED[] = 
    {"C===", "=C=T", "==AA", "==AT", "=TTA", "CT=T", "=TAA", "=TAT", "T=TA", "TC=T", "TCA=", "TCAT", "TT=A", "TT=T", "TTA=", "TTAT", "GGC=T=",
     "C=T=", "C=GT=", "C=A=", "CCTA", "CC=A"};

const char* EqualsTest::expectedReferenceName;
const char* EqualsTest::expectedMateReferenceName;
const char* EqualsTest::expectedMateReferenceNameOrEqual;
const char* EqualsTest::expectedCigar;
const char* EqualsTest::expectedQuality;

std::vector<unsigned int> EqualsTest::expectedCigarHex;

int EqualsTest::expected0BasedAlignmentEnd;
int EqualsTest::expected1BasedAlignmentEnd;
int EqualsTest::expectedAlignmentLength;
int EqualsTest::expected0BasedUnclippedStart;
int EqualsTest::expected1BasedUnclippedStart;
int EqualsTest::expected0BasedUnclippedEnd;
int EqualsTest::expected1BasedUnclippedEnd;
bamRecordStruct EqualsTest::expectedRecord;

void EqualsTest::testEq(FileType inputType)
{
    reset();
    SamFile inSam;

    std::string outputBase = "results/out";

    if(inputType == SAM)
    {
        assert(inSam.OpenForRead("testFiles/testEq.sam"));
        outputBase += "SamEq";
    }
    else
    {
        assert(inSam.OpenForRead("testFiles/testEq.bam"));
        outputBase += "BamEq";
    }

   // Read the SAM Header.
    SamFileHeader samHeader;
    assert(inSam.ReadHeader(samHeader));

    std::string outputName = outputBase + "Bases.sam";
    SamFile outBasesSam( outputName.c_str(), SamFile::WRITE);
    outputName = outputBase + "Equals.sam";
    SamFile outEqualsSam(outputName.c_str(), SamFile::WRITE);
    outputName = outputBase + "Orig.sam";
    SamFile outOrigSam(  outputName.c_str(), SamFile::WRITE);
    outputName = outputBase + "Bases.bam";
    SamFile outBasesBam( outputName.c_str(), SamFile::WRITE);
    outputName = outputBase + "Equals.bam";
    SamFile outEqualsBam(outputName.c_str(), SamFile::WRITE);
    outputName = outputBase + "Orig.bam";
    SamFile outOrigBam(  outputName.c_str(), SamFile::WRITE);
    assert(outBasesSam.WriteHeader(samHeader));
    assert(outEqualsSam.WriteHeader(samHeader));
    assert(outOrigSam.WriteHeader(samHeader));
    assert(outBasesBam.WriteHeader(samHeader));
    assert(outEqualsBam.WriteHeader(samHeader));
    assert(outOrigBam.WriteHeader(samHeader));

    outBasesSam.SetWriteSequenceTranslation(SamRecord::BASES);
    outEqualsSam.SetWriteSequenceTranslation(SamRecord::EQUAL);
    outOrigSam.SetWriteSequenceTranslation(SamRecord::NONE);
    outBasesBam.SetWriteSequenceTranslation(SamRecord::BASES);
    outEqualsBam.SetWriteSequenceTranslation(SamRecord::EQUAL);
    outOrigBam.SetWriteSequenceTranslation(SamRecord::NONE);

    GenomeSequence reference("testFiles/chr1_partial.fa");

    inSam.SetReference(&reference);

    SamRecord samRecord;

    // The set of 16 variations are repeated 3 times: once with all charcters
    // 1) Matches have the actual bases in them.
    // 2) Matches have '='
    // 3) Matches are mixed between bases and '='
    // Since Sequences are 4 characters long, there are 16 variations
    // of match/mismatch.
    for(int j = 0; j < 16; j++)
    {
        assert(inSam.ReadRecord(samHeader, samRecord) == true);
        validateEqRead(samRecord, j, READ_SEQS_BASES[j]);
        assert(outBasesSam.WriteRecord(samHeader, samRecord));
        assert(outEqualsSam.WriteRecord(samHeader, samRecord));
        assert(outOrigSam.WriteRecord(samHeader, samRecord));
        assert(outBasesBam.WriteRecord(samHeader, samRecord));
        assert(outEqualsBam.WriteRecord(samHeader, samRecord));
        assert(outOrigBam.WriteRecord(samHeader, samRecord));
    }
    for(int j = 0; j < 16; j++)
    {
        assert(inSam.ReadRecord(samHeader, samRecord) == true);
        validateEqRead(samRecord, j, READ_SEQS_EQUALS[j]);
        assert(outBasesSam.WriteRecord(samHeader, samRecord));
        assert(outEqualsSam.WriteRecord(samHeader, samRecord));
        assert(outOrigSam.WriteRecord(samHeader, samRecord));
        assert(outBasesBam.WriteRecord(samHeader, samRecord));
        assert(outEqualsBam.WriteRecord(samHeader, samRecord));
        assert(outOrigBam.WriteRecord(samHeader, samRecord));
    }
    for(int j = 0; j < 16; j++)
    {
        assert(inSam.ReadRecord(samHeader, samRecord) == true);
        validateEqRead(samRecord, j, READ_SEQS_MIXED[j]);
        assert(outBasesSam.WriteRecord(samHeader, samRecord));
        assert(outEqualsSam.WriteRecord(samHeader, samRecord));
        assert(outOrigSam.WriteRecord(samHeader, samRecord));
        assert(outBasesBam.WriteRecord(samHeader, samRecord));
        assert(outEqualsBam.WriteRecord(samHeader, samRecord));
        assert(outOrigBam.WriteRecord(samHeader, samRecord));
    }

    expectedCigar = "2S4M";
    expectedCigarHex.clear();
    expectedCigarHex.push_back(0x24);
    expectedCigarHex.push_back(0x40);
    expected0BasedUnclippedStart = expectedRecord.myPosition-2;
    expected1BasedUnclippedStart = expected0BasedUnclippedStart + 1;
    expectedRecord.myBlockSize = 70;
    expectedRecord.myReadNameLength = 21;
    expectedRecord.myCigarLength = 2;
    expectedRecord.myReadLength = 6;
    expectedQuality = "??I00?";
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateEqRead(samRecord, 16, READ_SEQS_MIXED[16]);
    assert(outBasesSam.WriteRecord(samHeader, samRecord));
    assert(outEqualsSam.WriteRecord(samHeader, samRecord));
    assert(outOrigSam.WriteRecord(samHeader, samRecord));
    assert(outBasesBam.WriteRecord(samHeader, samRecord));
    assert(outEqualsBam.WriteRecord(samHeader, samRecord));
    assert(outOrigBam.WriteRecord(samHeader, samRecord));
    
    expectedCigar = "4M4H";
    expectedCigarHex.clear();
    expectedCigarHex.push_back(0x40);
    expectedCigarHex.push_back(0x45);
    expected0BasedUnclippedStart = expectedRecord.myPosition;
    expected1BasedUnclippedStart = expected0BasedUnclippedStart + 1;
    expected0BasedUnclippedEnd = expectedRecord.myPosition + 7;
    expected1BasedUnclippedEnd = expected0BasedUnclippedEnd + 1;
    expectedRecord.myBlockSize = 65;
    expectedRecord.myReadNameLength = 19;
    expectedRecord.myCigarLength = 2;
    expectedRecord.myReadLength = 4;
    expectedQuality = "I00?";
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateEqRead(samRecord, 17, READ_SEQS_MIXED[17]);
    assert(outBasesSam.WriteRecord(samHeader, samRecord));
    assert(outEqualsSam.WriteRecord(samHeader, samRecord));
    assert(outOrigSam.WriteRecord(samHeader, samRecord));
    assert(outBasesBam.WriteRecord(samHeader, samRecord));
    assert(outEqualsBam.WriteRecord(samHeader, samRecord));
    assert(outOrigBam.WriteRecord(samHeader, samRecord));
    
    expectedCigar = "1M1P1M1I1M3D1M";
    expectedCigarHex.clear();
    expectedCigarHex.push_back(0x10);
    expectedCigarHex.push_back(0x16);
    expectedCigarHex.push_back(0x10);
    expectedCigarHex.push_back(0x11);
    expectedCigarHex.push_back(0x10);
    expectedCigarHex.push_back(0x32);
    expectedCigarHex.push_back(0x10);
    expected0BasedAlignmentEnd = expectedRecord.myPosition + 6;
    expected1BasedAlignmentEnd = expected0BasedAlignmentEnd + 1;
    expectedAlignmentLength = 7;
    expected0BasedUnclippedStart = expectedRecord.myPosition;
    expected1BasedUnclippedStart = expected0BasedUnclippedStart + 1;
    expected0BasedUnclippedEnd = expected0BasedAlignmentEnd;
    expected1BasedUnclippedEnd = expected0BasedUnclippedEnd + 1;
    expectedRecord.myBlockSize = 95;
    expectedRecord.myReadNameLength = 27;
    expectedRecord.myCigarLength = 7;
    expectedRecord.myReadLength = 5;
    expectedQuality = "I00??";
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateEqRead(samRecord, 18, READ_SEQS_MIXED[18]);
    assert(outBasesSam.WriteRecord(samHeader, samRecord));
    assert(outEqualsSam.WriteRecord(samHeader, samRecord));
    assert(outOrigSam.WriteRecord(samHeader, samRecord));
    assert(outBasesBam.WriteRecord(samHeader, samRecord));
    assert(outEqualsBam.WriteRecord(samHeader, samRecord));
    assert(outOrigBam.WriteRecord(samHeader, samRecord));

    expectedCigar = "2M2N2M";
    expectedCigarHex.clear();
    expectedCigarHex.push_back(0x20);
    expectedCigarHex.push_back(0x23);
    expectedCigarHex.push_back(0x20);
    expected0BasedAlignmentEnd = expectedRecord.myPosition + 5;
    expected1BasedAlignmentEnd = expected0BasedAlignmentEnd + 1;
    expectedAlignmentLength = 6;
    expected0BasedUnclippedStart = expectedRecord.myPosition;
    expected1BasedUnclippedStart = expected0BasedUnclippedStart + 1;
    expected0BasedUnclippedEnd = expected0BasedAlignmentEnd;
    expected1BasedUnclippedEnd = expected0BasedUnclippedEnd + 1;
    expectedRecord.myBlockSize = 74;
    expectedRecord.myReadNameLength = 24;
    expectedRecord.myCigarLength = 3;
    expectedRecord.myReadLength = 4;
    expectedQuality = "I00?";
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateEqRead(samRecord, 19, READ_SEQS_MIXED[19]);
    assert(outBasesSam.WriteRecord(samHeader, samRecord));
    assert(outEqualsSam.WriteRecord(samHeader, samRecord));
    assert(outOrigSam.WriteRecord(samHeader, samRecord));
    assert(outBasesBam.WriteRecord(samHeader, samRecord));
    assert(outEqualsBam.WriteRecord(samHeader, samRecord));
    assert(outOrigBam.WriteRecord(samHeader, samRecord));

    // Test getNextMatchMismatch.
    SamSingleBaseMatchInfo matchTest;
    SamQuerySeqWithRefIter queryIter(samRecord, reference, true);
    assert(queryIter.getNextMatchMismatch(matchTest) == true);
    assert(matchTest.getType() == SamSingleBaseMatchInfo::MATCH);
    assert(matchTest.getQueryIndex() == 0);
    assert(queryIter.getNextMatchMismatch(matchTest) == true);
    assert(matchTest.getType() == SamSingleBaseMatchInfo::MATCH);
    assert(matchTest.getQueryIndex() == 1);
    assert(queryIter.getNextMatchMismatch(matchTest) == true);
    assert(matchTest.getType() == SamSingleBaseMatchInfo::MATCH);
    assert(matchTest.getQueryIndex() == 2);
    assert(queryIter.getNextMatchMismatch(matchTest) == true);
    assert(matchTest.getType() == SamSingleBaseMatchInfo::MATCH);
    assert(matchTest.getQueryIndex() == 3);
    assert(queryIter.getNextMatchMismatch(matchTest) == false);

    // Check the read that is on a different chormosome not
    // found in the reference.
    reset();
    expectedRecord.myBlockSize = 56;
    expectedRecord.myReadNameLength = 14;
    expectedRecord.myReferenceID = 1;
    expectedReferenceName = "2";
    expectedRecord.myMateReferenceID = 1;
    expectedMateReferenceName = "2";

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateEqRead(samRecord, 20, READ_SEQS_MIXED[20]);
    assert(outBasesSam.WriteRecord(samHeader, samRecord));
    assert(outEqualsSam.WriteRecord(samHeader, samRecord));
    assert(outOrigSam.WriteRecord(samHeader, samRecord));
    assert(outBasesBam.WriteRecord(samHeader, samRecord));
    assert(outEqualsBam.WriteRecord(samHeader, samRecord));
    assert(outOrigBam.WriteRecord(samHeader, samRecord));   

    // Check the read that is on a different chormosome and
    // has '=', but the chromosome is not found in the reference.
    reset();
    expectedRecord.myBlockSize = 57;
    expectedRecord.myReadNameLength = 15;
    expectedRecord.myReferenceID = 1;
    expectedReferenceName = "2";
    expectedRecord.myMateReferenceID = 1;
    expectedMateReferenceName = "2";

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateEqRead(samRecord, 21, READ_SEQS_MIXED[21]);
    assert(outBasesSam.WriteRecord(samHeader, samRecord));
    assert(outEqualsSam.WriteRecord(samHeader, samRecord));
    assert(outOrigSam.WriteRecord(samHeader, samRecord));
    assert(outBasesBam.WriteRecord(samHeader, samRecord));
    assert(outEqualsBam.WriteRecord(samHeader, samRecord));
    assert(outOrigBam.WriteRecord(samHeader, samRecord));   

    SamQuerySeqWithRefIter queryIter2(samRecord, reference, true);
    assert(queryIter2.getNextMatchMismatch(matchTest) == false);
}


void EqualsTest::reset()
{
    expectedReferenceName = "1";
    expectedMateReferenceName = "1";
    expectedMateReferenceNameOrEqual = "=";
    expectedCigar = "4M";
    expectedQuality = "I00?";

    // The First cigar is 4M which is 4 << 4 | 0 = 0x40 = 64
    expectedCigarHex.clear();
    expectedCigarHex.push_back(0x40);

    expectedRecord.myBlockSize = 50;
    expectedRecord.myReferenceID = 0;
    expectedRecord.myPosition = 10010;
    expectedRecord.myReadNameLength = 8;
    expectedRecord.myMapQuality = 0;
    expectedRecord.myBin = 4681;
    expectedRecord.myCigarLength = 1;
    expectedRecord.myFlag = 73;
    expectedRecord.myReadLength = 4;
    expectedRecord.myMateReferenceID = 0;
    expectedRecord.myMatePosition = 10008;
    expectedRecord.myInsertSize = 0;

    expected0BasedAlignmentEnd = 10013;
    expected1BasedAlignmentEnd = expected0BasedAlignmentEnd + 1;
    expectedAlignmentLength = 4;
    expected0BasedUnclippedStart = expectedRecord.myPosition;
    expected1BasedUnclippedStart = expected0BasedUnclippedStart + 1;
    expected0BasedUnclippedEnd = expected0BasedAlignmentEnd;
    expected1BasedUnclippedEnd = expected1BasedAlignmentEnd;
}

void EqualsTest::validateEqRead(SamRecord& samRecord, 
                                int readIndex,
                                const char* actualExpectedSequence)
{
    char tag[3];
    char type;
    void* value;

    //////////////////////////////////////////
    // Validate Record 1
    // Check the alignment end
    assert(samRecord.get0BasedAlignmentEnd() == expected0BasedAlignmentEnd);
    assert(samRecord.get1BasedAlignmentEnd() == expected1BasedAlignmentEnd);
    assert(samRecord.getAlignmentLength() == expectedAlignmentLength);
    assert(samRecord.get0BasedUnclippedStart() == expected0BasedUnclippedStart);
    assert(samRecord.get1BasedUnclippedStart() == expected1BasedUnclippedStart);
    assert(samRecord.get0BasedUnclippedEnd() == expected0BasedUnclippedEnd);
    assert(samRecord.get1BasedUnclippedEnd() == expected1BasedUnclippedEnd);

    // Check the accessors.
    assert(samRecord.getBlockSize() == expectedRecord.myBlockSize);
    assert(samRecord.getReferenceID() == expectedRecord.myReferenceID);
    assert(strcmp(samRecord.getReferenceName(), expectedReferenceName) == 0);
    assert(samRecord.get1BasedPosition() == expectedRecord.myPosition + 1);
    assert(samRecord.get0BasedPosition() == expectedRecord.myPosition);
    assert(samRecord.getReadNameLength() == 
           expectedRecord.myReadNameLength);
    assert(samRecord.getMapQuality() == expectedRecord.myMapQuality);
    assert(samRecord.getBin() == expectedRecord.myBin);
    assert(samRecord.getCigarLength() == expectedRecord.myCigarLength);
    assert(samRecord.getFlag() == expectedRecord.myFlag);
    assert(samRecord.getReadLength() == expectedRecord.myReadLength);
    assert(samRecord.getMateReferenceID() == 
           expectedRecord.myMateReferenceID);
    assert(strcmp(samRecord.getMateReferenceName(),
                  expectedMateReferenceName) == 0);
    assert(strcmp(samRecord.getMateReferenceNameOrEqual(),
                  expectedMateReferenceNameOrEqual) == 0);
    assert(samRecord.get1BasedMatePosition() ==
           expectedRecord.myMatePosition + 1);
    assert(samRecord.get0BasedMatePosition() == 
           expectedRecord.myMatePosition);
    assert(samRecord.getInsertSize() == expectedRecord.myInsertSize);
    assert(strcmp(samRecord.getReadName(), READ_NAMES[readIndex]) == 0);
    assert(strcmp(samRecord.getCigar(), expectedCigar) == 0);
    samRecord.setSequenceTranslation(SamRecord::BASES);
    assert(strcmp(samRecord.getSequence(), READ_SEQS_BASES[readIndex]) == 0);
    assert(strcmp(samRecord.getQuality(), expectedQuality) == 0);

    assert(samRecord.getSequence(0) == READ_SEQS_BASES[readIndex][0]);
    assert(samRecord.getQuality(0) == expectedQuality[0]);
    assert(samRecord.getSequence(1)== READ_SEQS_BASES[readIndex][1]);
    assert(samRecord.getQuality(1) == expectedQuality[1]);
    assert(samRecord.getSequence(2) == READ_SEQS_BASES[readIndex][2]);
    assert(samRecord.getQuality(2) == expectedQuality[2]);
    assert(samRecord.getSequence(3) == READ_SEQS_BASES[readIndex][3]);
    assert(samRecord.getQuality(3) == expectedQuality[3]);

    assert(strcmp(samRecord.getSequence(SamRecord::EQUAL), 
                  READ_SEQS_EQUALS[readIndex]) == 0);
    assert(samRecord.getSequence(0, SamRecord::EQUAL) == 
           READ_SEQS_EQUALS[readIndex][0]);
    assert(samRecord.getQuality(0) == expectedQuality[0]);
    assert(samRecord.getSequence(1, SamRecord::EQUAL) == 
           READ_SEQS_EQUALS[readIndex][1]);
    assert(samRecord.getQuality(1) == expectedQuality[1]);
    assert(samRecord.getSequence(2, SamRecord::EQUAL) == 
           READ_SEQS_EQUALS[readIndex][2]);
    assert(samRecord.getQuality(2) == expectedQuality[2]);
    assert(samRecord.getSequence(3, SamRecord::EQUAL) == 
           READ_SEQS_EQUALS[readIndex][3]);
    assert(samRecord.getQuality(3) == expectedQuality[3]);

    assert(strcmp(samRecord.getSequence(SamRecord::NONE), 
                  actualExpectedSequence) == 0);
    assert(samRecord.getSequence(0, SamRecord::NONE) == 
           actualExpectedSequence[0]);
    assert(samRecord.getQuality(0) == expectedQuality[0]);
    assert(samRecord.getSequence(1, SamRecord::NONE) == 
           actualExpectedSequence[1]);
    assert(samRecord.getQuality(1) == expectedQuality[1]);
    assert(samRecord.getSequence(2, SamRecord::NONE) == 
           actualExpectedSequence[2]);
    assert(samRecord.getQuality(2) == expectedQuality[2]);
    assert(samRecord.getSequence(3, SamRecord::NONE) == 
           actualExpectedSequence[3]);
    assert(samRecord.getQuality(3) == expectedQuality[3]);

    samRecord.setSequenceTranslation(SamRecord::NONE);
    assert(strcmp(samRecord.getSequence(), 
                  actualExpectedSequence) == 0);
    assert(samRecord.getSequence(0) == 
           actualExpectedSequence[0]);
    assert(samRecord.getQuality(0) == expectedQuality[0]);
    assert(samRecord.getSequence(1) == 
           actualExpectedSequence[1]);
    assert(samRecord.getQuality(1) == expectedQuality[1]);
    assert(samRecord.getSequence(2) == 
           actualExpectedSequence[2]);
    assert(samRecord.getQuality(2) == expectedQuality[2]);
    assert(samRecord.getSequence(3) == 
           actualExpectedSequence[3]);
    assert(samRecord.getQuality(3) == expectedQuality[3]);

    // No tags, should return false.
    
    assert(samRecord.getNextSamTag(tag, type, &value) == false);

    // Get the record ptr.   
    samRecord.setSequenceTranslation(SamRecord::BASES);
    validateEqReadBuffer(samRecord, READ_SEQS_BASES[readIndex]);
    samRecord.setSequenceTranslation(SamRecord::NONE);
    validateEqReadBuffer(samRecord, actualExpectedSequence);
    samRecord.setSequenceTranslation(SamRecord::EQUAL);
    validateEqReadBuffer(samRecord, READ_SEQS_EQUALS[readIndex]);
}


void EqualsTest::validateEqReadBuffer(SamRecord& samRecord, 
                                      const char* expectedSequence)
{
    const bamRecordStruct* bufferPtr;
    unsigned char* varPtr;

    bufferPtr = (const bamRecordStruct*)samRecord.getRecordBuffer();
    // Validate the buffers match.
    assert(bufferPtr->myBlockSize == expectedRecord.myBlockSize);
    assert(bufferPtr->myReferenceID == expectedRecord.myReferenceID);
    assert(bufferPtr->myPosition == expectedRecord.myPosition);
    assert(bufferPtr->myReadNameLength == expectedRecord.myReadNameLength);
    assert(bufferPtr->myMapQuality == expectedRecord.myMapQuality);
    assert(bufferPtr->myBin == expectedRecord.myBin);
    assert(bufferPtr->myCigarLength == expectedRecord.myCigarLength);
    assert(bufferPtr->myFlag == expectedRecord.myFlag);
    assert(bufferPtr->myReadLength == expectedRecord.myReadLength);
    assert(bufferPtr->myMateReferenceID == 
           expectedRecord.myMateReferenceID);
    assert(bufferPtr->myMatePosition == expectedRecord.myMatePosition);
    assert(bufferPtr->myInsertSize == expectedRecord.myInsertSize);

    // Validate the variable length fields in the buffer.
    // Set the pointer to the start of the variable fields.
    varPtr = (unsigned char*)(&(bufferPtr->myData[0]));

    // Validate the readname.
    for(int i = 0; i < expectedRecord.myReadNameLength; i++)
    {
        assert(*varPtr == samRecord.getReadName()[i]);
        varPtr++;
    }

    // Validate the cigar.
    for(int i = 0; i < expectedRecord.myCigarLength; i++)
    {
        assert(*(unsigned int*)varPtr == expectedCigarHex[i]);
        // Increment the varptr the size of an int.
        varPtr += 4;
    }

    // Validate the sequence.
    int expectedSeqHex = 0;
    for(int i = 0; i < expectedRecord.myReadLength; i++)
    {
        int hexChar = 0x0;
        switch(expectedSequence[i])
        {
            case '=':
                hexChar = 0x0;
                break;
            case 'A':
            case 'a':
                hexChar = 0x1;
                break;
            case 'C':
            case 'c':
                hexChar = 0x2;
                break;
            case 'G':
            case 'g':
                hexChar = 0x4;
                break;
            case 'T':
            case 't':
                hexChar = 0x8;
                break;
            case 'N':
            case 'n':
                hexChar = 0xF;
                break;
        }
        if(i%2 == 0)
        {
            expectedSeqHex = hexChar << 4;
        }
        else
        {
            expectedSeqHex |= hexChar;
            assert(*varPtr == expectedSeqHex);
            varPtr++;         
        }
    }
    if((expectedRecord.myReadLength%2) != 0)
    {
       // Odd number of sequences, so need to check the last one.
        assert(*varPtr == expectedSeqHex);
        varPtr++;
    }

    // Validate the Quality
    for(int i = 0; i < expectedRecord.myReadLength; i++)
    {
        assert(*varPtr == samRecord.getQuality()[i] - 33);
        varPtr++;
    }

}


