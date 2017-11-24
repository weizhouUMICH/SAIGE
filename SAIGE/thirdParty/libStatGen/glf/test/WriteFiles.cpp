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

#include "WriteFiles.h"
#include "Validate.h"

#include <assert.h>

void testWrite()
{
    TestWrite writeTest;
    writeTest.testWrite();
}

const std::string TestWrite::HEADER_TEXT1 = "This is my 1st test header.";
const std::string TestWrite::SEC1_REFNAME = "This is my 1st RefName";
const std::string TestWrite::SEC1REC2_INDELSEQ1 = "AC";
const std::string TestWrite::SEC1REC2_INDELSEQ2 = "TCA";
const std::string TestWrite::SEC2_REFNAME = "This is my 2nd RefName";
const std::string TestWrite::HEADER_TEXT2 = "This is my 2nd test header.";
const std::string TestWrite::HEADER_TEXT3 = "This is my 3rd test header.";

void TestWrite::testWrite()
{
    GlfFile glfOut;
    
    std::string testFile = "results/MyTestOut1.glf";

    assert(glfOut.openForWrite(testFile.c_str(), false));

    // Create a glf header.
    GlfHeader glfHeader;
    GlfRefSection glfSection;
    GlfRecord record;
  
    // Test writing refsection with no header - exception
    bool caughtException = false;
    try
    {
        assert(glfOut.writeRefSection(glfSection) == false);
    }
    catch (std::exception& e) 
    {
        caughtException = true;
    }
    assert(caughtException);

    // Test writing record with no header - exception.
    caughtException = false;
    try
    {
        assert(glfOut.writeRecord(record) == false);
    }
    catch (std::exception& e) 
    {
        caughtException = true;
    }
    assert(caughtException);

    // Write the header.
    writeHeader(glfOut, 1);

    // Test writing record with no refsection - exception.
    caughtException = false;
    try
    {
        assert(glfOut.writeRecord(record) == false);
    }
    catch (std::exception& e) 
    {
        caughtException = true;
    }
    assert(caughtException);


    //////////////////////////////////////////////
    writeRefSection1(glfOut);

    // Test writing header after refSection - exception
    caughtException = false;
    try
    {
        assert(glfOut.writeHeader(glfHeader) == false);
    }
    catch (std::exception& e) 
    {
        caughtException = true;
    }
    assert(caughtException);

    writeSec1Record1(glfOut);
    // Test writing header after record - exception
    caughtException = false;
    try
    {
        assert(glfOut.writeHeader(glfHeader) == false);
    }
    catch (std::exception& e) 
    {
        caughtException = true;
    }
    assert(caughtException);

    writeSec1Record2(glfOut);
    writeEndMarker(glfOut);

    writeRefSection2(glfOut);
    writeSec2Record1(glfOut);
    writeEndMarker(glfOut);

    ////////////////////
    // Close the file.
    glfOut.close();

    //////////////////////////////////////////////
    // Validate the just written file.
    GlfFile glfIn;
    assert(glfIn.openForRead(testFile.c_str()));

    readHeader(glfIn, 1);
    readRefSection1(glfIn);
    readSec1Record1(glfIn);
    readSec1Record2(glfIn);
    readEndMarker(glfIn);
    readRefSection2(glfIn);
    readSec2Record1(glfIn);
    readEndMarker(glfIn);
    checkEOF(glfIn);

    ////////////////////////////////
    // NEW FILE
    testFile = "results/MyTestOut2.glf";
    assert(glfOut.openForWrite(testFile.c_str()));
  
    writeHeader(glfOut, 2);
    writeRefSection1(glfOut);
    writeSec1Record1(glfOut);
    writeSec1Record2(glfOut);
    // Test writing new section without end of section marker - auto-added.
    writeRefSection2(glfOut);
    writeSec2Record1(glfOut);
    // Test closing file with no end of section marker - auto-added.
    glfOut.close();

    //////////////////////////////////////////////
    // Validate the just written file.
    assert(glfIn.openForRead(testFile.c_str()));

    readHeader(glfIn, 2);
    readRefSection1(glfIn);
    readSec1Record1(glfIn);
    readSec1Record2(glfIn);
    readEndMarker(glfIn);
    readRefSection2(glfIn);
    readSec2Record1(glfIn);
    readEndMarker(glfIn);
    checkEOF(glfIn);

 
    ////////////////////////////////
    // NEW FILE
    testFile = "results/MyTestOut3.glf";
    {
        GlfFile glfOutScoped;
        assert(glfOutScoped.openForWrite(testFile.c_str()));
        
        writeHeader(glfOutScoped, 3);
        writeRefSection1(glfOutScoped);
        writeSec1Record1(glfOutScoped);
        writeSec1Record2(glfOutScoped);
        // Test writing new section without end of section marker - auto-added.
        writeRefSection2(glfOutScoped);
        writeSec2Record1(glfOutScoped);
        // Test just letting the file go out of scope with no end
        // of section marker - auto added.
    }
    //////////////////////////////////////////////
    // Validate the just written file.
    assert(glfIn.openForRead(testFile.c_str()));
    
    // Test reading refsection with no header - exception.
    caughtException = false;
    try
    {
        assert(glfIn.getNextRefSection(glfSection) == false);
    }
    catch (std::exception& e) 
    {
        caughtException = true;
    }
    assert(caughtException);

    // Test reading record with no header - exception.
    caughtException = false;
    try
    {
        assert(glfIn.getNextRecord(record) == false);
    }
    catch (std::exception& e) 
    {
        caughtException = true;
    }
    assert(caughtException);

    readHeader(glfIn, 3);

    // Test reading record with no reference section - exception.
    caughtException = false;
    try
    {
        assert(glfIn.getNextRecord(record) == false);
    }
    catch (std::exception& e) 
    {
        caughtException = true;
    }
    assert(caughtException);

    // Test reading header after already read - exception
    caughtException = false;
    try
    {
        assert(glfIn.readHeader(glfHeader) == false);
    }
    catch (std::exception& e) 
    {
        caughtException = true;
    }
    assert(caughtException);

    readRefSection1(glfIn);
    readSec1Record1(glfIn);
    readSec1Record2(glfIn);
    readEndMarker(glfIn);
    readRefSection2(glfIn);
    readSec2Record1(glfIn);
    readEndMarker(glfIn);
    checkEOF(glfIn);
   

    // Read again, but text reading next refsection before 
    //end of current section - consumes the rest of the records.
    assert(glfIn.openForRead(testFile.c_str()));

    readHeader(glfIn, 3);
    readRefSection1(glfIn);
    readRefSection2(glfIn);
    readSec2Record1(glfIn);
    readEndMarker(glfIn);
    checkEOF(glfIn);
 

}


void TestWrite::writeHeader(GlfFile& glfOut, int headerNum)
{
    GlfHeader glfHeader;
    std::string headerString = "t";
    std::string expectedHeader = "";
    if(headerNum == 1)
    {
        expectedHeader = HEADER_TEXT1;
    }
    else if(headerNum == 2)
    {
        expectedHeader = HEADER_TEXT2;
    }
    else if(headerNum == 3)
    {
        expectedHeader = HEADER_TEXT3;
    }

    assert(glfHeader.getHeaderTextString(headerString));
    assert(headerString == "");
    assert(glfHeader.setHeaderTextString(expectedHeader));
    assert(glfHeader.getHeaderTextString(headerString));
    assert(headerString == expectedHeader);
    assert(glfOut.writeHeader(glfHeader));
}


void TestWrite::writeRefSection1(GlfFile& glfOut)
{
    GlfRefSection glfSection;

    ////////////////////////////////
    // Write the reference section.
    std::string refNameString = "";
    // Check the default settings (no data has been set yet).
    assert(glfSection.getName(refNameString));
    assert(refNameString == "");
    assert(glfSection.getRefLen() == 0);

    // Set the reference name.
    assert(glfSection.setName(SEC1_REFNAME));
    // Check properly set.
    assert(glfSection.getName(refNameString));
    assert(refNameString == SEC1_REFNAME);
    assert(glfSection.getRefLen() == 0);

    // Set the reference sequence length.
    assert(glfSection.setRefLen(SEC1_REFLEN));
    // Check properly set.
    assert(glfSection.getRefLen() == SEC1_REFLEN);
    assert(glfSection.getName(refNameString));
    assert(refNameString == SEC1_REFNAME);

    // Write the reference section
    assert(glfOut.writeRefSection(glfSection));
}


void TestWrite::writeSec1Record1(GlfFile& glfOut)
{
    GlfRecord record;
    assert(record.setRecordType(SEC1REC1_RECTYPE));
    assert(record.setRefBaseInt(SEC1REC1_REFBASE));
    assert(record.setOffset(SEC1REC1_OFFSET));
    assert(record.setMinLk(SEC1REC1_MINLK));
    assert(record.setReadDepth(SEC1REC1_READDEPTH));
    assert(record.setRmsMapQ(SEC1REC1_RMSMAPQ));
    assert(glfOut.writeRecord(record));
    
    // Verify the settings of record 1.
    assert(record.getRecordType() == SEC1REC1_RECTYPE);
    assert(record.getRefBase() == SEC1REC1_REFBASE);
    assert(record.getOffset() == SEC1REC1_OFFSET);
    assert(record.getMinLk() == SEC1REC1_MINLK);
    assert(record.getReadDepth() == SEC1REC1_READDEPTH);
    assert(record.getRmsMapQ() == SEC1REC1_RMSMAPQ);
}


void TestWrite::writeSec1Record2(GlfFile& glfOut)
{
    //////////////////////////////////////////////
    // Write a record of type 2.
    GlfRecord record;

    assert(record.setRecordType(SEC1REC2_RECTYPE));
    assert(record.setRefBaseInt(SEC1REC2_REFBASE));
    assert(record.setOffset(SEC1REC2_OFFSET));
    assert(record.setMinLk(SEC1REC2_MINLK));
    assert(record.setReadDepth(SEC1REC2_READDEPTH));
    assert(record.setRmsMapQ(SEC1REC2_RMSMAPQ));
    assert(record.setLkHom1(SEC1REC2_LKHOM1));
    assert(record.setLkHom2(SEC1REC2_LKHOM2));
    assert(record.setLkHet(SEC1REC2_LKHET));
    assert(record.setInsertionIndel1(SEC1REC2_INDELSEQ1));
    assert(record.setDeletionIndel2(SEC1REC2_INDELSEQ2));
    assert(glfOut.writeRecord(record));

    // Verify the settings of record 2.
    std::string indelSeq = "";
    assert(record.getRecordType() == SEC1REC2_RECTYPE);
    assert(record.getRefBase() == SEC1REC2_REFBASE);
    assert(record.getOffset() == SEC1REC2_OFFSET);
    assert(record.getMinLk() == SEC1REC2_MINLK);
    assert(record.getReadDepth() == SEC1REC2_READDEPTH);
    assert(record.getRmsMapQ() == SEC1REC2_RMSMAPQ);
    assert(record.getLkHom1() == SEC1REC2_LKHOM1);
    assert(record.getLkHom2() == SEC1REC2_LKHOM2);
    assert(record.getLkHet() == SEC1REC2_LKHET);
    assert(record.getIndel1(indelSeq) == SEC1REC2_INDELLEN1);
    assert(indelSeq == SEC1REC2_INDELSEQ1);
    assert(record.getIndel2(indelSeq) == SEC1REC2_INDELLEN2);
    assert(indelSeq == SEC1REC2_INDELSEQ2);
}


void TestWrite::writeEndMarker(GlfFile& glfOut)
{
    //////////////////////////////////////////////
    // Write a record of type 0.
    GlfRecord record;
    assert(glfOut.writeRecord(record));

    // Verify the settings of the types.
    assert(record.getRecordType() == 0);
    assert(record.getRefBase() == 0);
}


void TestWrite::writeRefSection2(GlfFile& glfOut)
{
    GlfRefSection glfSection;

    ////////////////////////////////
    // Write the reference section.
    std::string refNameString = "";
    // Check the default settings (no data has been set yet).
    assert(glfSection.getName(refNameString));
    assert(refNameString == "");
    assert(glfSection.getRefLen() == 0);

    // Set the reference name.
    assert(glfSection.setName(SEC2_REFNAME));
    // Check properly set.
    assert(glfSection.getName(refNameString));
    assert(refNameString == SEC2_REFNAME);
    assert(glfSection.getRefLen() == 0);

    // Set the reference sequence length.
    assert(glfSection.setRefLen(SEC2_REFLEN));
    // Check properly set.
    assert(glfSection.getRefLen() == SEC2_REFLEN);
    assert(glfSection.getName(refNameString));
    assert(refNameString == SEC2_REFNAME);

    // Write the reference section
    assert(glfOut.writeRefSection(glfSection));
}


void TestWrite::writeSec2Record1(GlfFile& glfOut)
{
    GlfRecord record;
    assert(record.setRecordType(SEC2REC1_RECTYPE));
    assert(record.setRefBaseInt(SEC2REC1_REFBASE));
    assert(record.setOffset(SEC2REC1_OFFSET));
    assert(record.setMinLk(SEC2REC1_MINLK));
    assert(record.setReadDepth(SEC2REC1_READDEPTH));
    assert(record.setRmsMapQ(SEC2REC1_RMSMAPQ));
    assert(glfOut.writeRecord(record));
    
    // Verify the settings of record 1.
    assert(record.getRecordType() == SEC2REC1_RECTYPE);
    assert(record.getRefBase() == SEC2REC1_REFBASE);
    assert(record.getOffset() == SEC2REC1_OFFSET);
    assert(record.getMinLk() == SEC2REC1_MINLK);
    assert(record.getReadDepth() == SEC2REC1_READDEPTH);
    assert(record.getRmsMapQ() == SEC2REC1_RMSMAPQ);
}


void TestWrite::readHeader(GlfFile& glfIn, int headerNum)
{
    GlfHeader glfHeader;
    std::string expectedHeader = "";
    std::string headerString;
    if(headerNum == 1)
    {
        expectedHeader = HEADER_TEXT1;
    }
    else if(headerNum == 2)
    {
        expectedHeader = HEADER_TEXT2;
    }
    else if(headerNum == 3)
    {
        expectedHeader = HEADER_TEXT3;
    }
    // Check the header string.
    assert(glfIn.readHeader(glfHeader));
    assert(glfHeader.getHeaderTextString(headerString));
    assert(headerString == expectedHeader);
}

void TestWrite::readRefSection1(GlfFile& glfIn)
{
    GlfRefSection glfSection;
    std::string refNameString;
    // Check the reference section.
    assert(glfIn.getNextRefSection(glfSection));
    assert(glfSection.getName(refNameString));
    assert(refNameString == SEC1_REFNAME);
    assert(glfSection.getRefLen() == SEC1_REFLEN);
}

void TestWrite::readSec1Record1(GlfFile& glfIn)
{
    GlfRecord record;
    // Check the record of type 1.
    assert(glfIn.getNextRecord(record));
    assert(record.getRecordType() == SEC1REC1_RECTYPE);
    assert(record.getRefBase() == SEC1REC1_REFBASE);
    assert(record.getOffset() == SEC1REC1_OFFSET);
    assert(record.getMinLk() == SEC1REC1_MINLK);
    assert(record.getReadDepth() == SEC1REC1_READDEPTH);
    assert(record.getRmsMapQ() == SEC1REC1_RMSMAPQ);
}

void TestWrite::readSec1Record2(GlfFile& glfIn)
{
    GlfRecord record;
    std::string indelSeq;
    //Check the record of type 2.
    assert(glfIn.getNextRecord(record));
    assert(record.getRecordType() == SEC1REC2_RECTYPE);
    assert(record.getRefBase() == SEC1REC2_REFBASE);
    assert(record.getOffset() == SEC1REC2_OFFSET);
    assert(record.getMinLk() == SEC1REC2_MINLK);
    assert(record.getReadDepth() == SEC1REC2_READDEPTH);
    assert(record.getRmsMapQ() == SEC1REC2_RMSMAPQ);
    assert(record.getLkHom1() == SEC1REC2_LKHOM1);
    assert(record.getLkHom2() == SEC1REC2_LKHOM2);
    assert(record.getLkHet() == SEC1REC2_LKHET);
    assert(record.getIndel1(indelSeq) == SEC1REC2_INDELLEN1);
    assert(indelSeq == SEC1REC2_INDELSEQ1);
    assert(record.getIndel2(indelSeq) == SEC1REC2_INDELLEN2);
    assert(indelSeq == SEC1REC2_INDELSEQ2);
}

void TestWrite::readEndMarker(GlfFile& glfIn)
{
    GlfRecord record;
    // Check the record of type 0.
    // False, since there are no more records in this section.
    assert(glfIn.getNextRecord(record) == false);
    assert(record.getRecordType() == 0);
    assert(record.getRefBase() == 0);
}

void TestWrite::readRefSection2(GlfFile& glfIn)
{
    GlfRefSection glfSection;
    std::string refNameString;
    // Check the reference section.
    assert(glfIn.getNextRefSection(glfSection));
    assert(glfSection.getName(refNameString));
    assert(refNameString == SEC2_REFNAME);
    assert(glfSection.getRefLen() == SEC2_REFLEN);
}


void TestWrite::readSec2Record1(GlfFile& glfIn)
{
    GlfRecord record;
    // Check the record of type 1.
    assert(glfIn.getNextRecord(record));
    assert(record.getRecordType() == SEC2REC1_RECTYPE);
    assert(record.getRefBase() == SEC2REC1_REFBASE);
    assert(record.getOffset() == SEC2REC1_OFFSET);
    assert(record.getMinLk() == SEC2REC1_MINLK);
    assert(record.getReadDepth() == SEC2REC1_READDEPTH);
    assert(record.getRmsMapQ() == SEC2REC1_RMSMAPQ);
}

void TestWrite::checkEOF(GlfFile& glfIn)
{
    GlfHeader glfHeader;
    GlfRefSection glfSection;
    GlfRecord record;
    // Check end of file - no more refsections
    assert(glfIn.getNextRefSection(glfSection) == false);
    assert(glfIn.isEOF());
}

