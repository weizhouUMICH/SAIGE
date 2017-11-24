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

#include <stdlib.h>
#include "GlfRecord.h"
#include "GlfException.h"
#include "StringBasics.h"

std::string GlfRecord::REF_BASE_CHAR = "XACMGRSVTQYHKSVN";

GlfRecord::GlfRecord()
{
    reset();
}


GlfRecord::~GlfRecord()
{
    reset();
}


// Reset the record for a new entry, clearing out previous values.
void GlfRecord::reset()
{
    myRecTypeRefBase = 0;
    myRec1Base.offset = 0;
    myRec1Base.min_depth = 0;
    myRec1Base.rmsMapQ = 0;
    for(int i = 0; i < 10; i++)
    {
        myRec1Base.lk[i] = 0;
    }

    myRec2Base.offset = 0;
    myRec2Base.min_depth = 0;
    myRec2Base.rmsMapQ = 0;
    myRec2Base.lkHom1 = 0;
    myRec2Base.lkHom2 = 0;
    myRec2Base.lkHet = 0;
    myRec2Base.indelLen1 = 0;
    myRec2Base.indelLen2 = 0;

    myIndelSeq1.reset();
    myIndelSeq2.reset();
}


// Read the record from the specified file.  Assumes the file is in
// the correct position for reading the record.
bool GlfRecord::read(IFILE filePtr)
{
    // Read the record type and reference base.
    int numRead = 0;
    int byteLen = sizeof(uint8_t);
    numRead = ifread(filePtr, &myRecTypeRefBase, byteLen);
    if(numRead != byteLen)
    {
        String errorMsg = "Failed to read the record type & reference base (";
        errorMsg += byteLen;
        errorMsg += " bytes).  Only read ";
        errorMsg += numRead;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
        return(false);
    }

    // TODO, split up by types of records...
    switch(getRecordType())
    {
        case 0:
            //  Last record.
            // Nothing more to read.
            break;
        case 1:
            // Read type 1.
            readType1(filePtr);
            break;
        case 2:
            // Read type 2.
            readType2(filePtr);
            break;
        default:
            String errorMsg = "Failed to read the record: unknown type: ";
            errorMsg += getRecordType();
            std::string errorString = errorMsg.c_str();
            throw(GlfException(GlfStatus::INVALID, errorString));
            return(false);
            break;
    };

    // Successfully read, return success.
    return(true);
}


// Write the record to the specified file.
bool GlfRecord::write(IFILE filePtr) const
{
    // TODO, split up by types of records...
    switch(getRecordType())
    {
        case 0:
            writeRtypeRef(filePtr);
            break;
        case 1:
            // write type 1.
            writeType1(filePtr);
            break;
        case 2:
            // write type 2.
            writeType2(filePtr);
            break;
        default:
            // unknown type, return error.
            String errorMsg = "Failed to write the record: unknown type: ";
            errorMsg += getRecordType();
            std::string errorString = errorMsg.c_str();
            throw(GlfException(GlfStatus::INVALID, errorString));
            return(false);
            break;
    };
    
    return(true);
}


void GlfRecord::print() const
{
    std::cout << "record_type: " << getRecordType()
              << "; ref_base: " << getRefBase()
              << "; ref_base_char: " << getRefBaseChar()
              << "\n";

    // TODO, split up by types of records...
    switch(getRecordType())
    {
        case 0:
            break;
        case 1:
            // print type 1.
            std::cout << "\toffset: " << myRec1Base.offset
                      << "; min_lk: " << (myRec1Base.min_depth >> 24)
                      << "; read_depth: " << (myRec1Base.min_depth & 0xFFFFFF)
                      << "; rmsMapQ: " << (int)myRec1Base.rmsMapQ;
            for(int i = 0; i < 10; ++i)
            {
                std::cout << "; lk[" << i << "]: " << (int)myRec1Base.lk[i];
            }

            std::cout << "\n";
            break;
        case 2:
            // print type 2.
            std::cout << "\toffset: " << myRec2Base.offset
                      << "; min_lk: " << (myRec2Base.min_depth >> 24)
                      << "; read_depth: " << (myRec2Base.min_depth & 0xFFFFF)
                      << "; rmsMapQ: " << (int)myRec2Base.rmsMapQ
                      << "; lkHom1: " << (int)myRec2Base.lkHom1
                      << "; lkHom2: " << (int)myRec2Base.lkHom2
                      << "; lkHet: " << (int)myRec2Base.lkHet
                      << "; indelLen1: " << myRec2Base.indelLen1
                      << "; indelLen2: " << myRec2Base.indelLen2
                      << "; myIndelSeq1: " << myIndelSeq1.c_str()
                      << "; myIndelSeq2: " << myIndelSeq2.c_str()
                      << "\n";
            break;
        default:
            break;
    };
}

bool GlfRecord::setRtypeRef(uint8_t rtypeRef)
{
    myRecTypeRefBase = rtypeRef;
    return(true);
}

bool GlfRecord::setRecordType(uint8_t recType)
{
    myRecTypeRefBase = 
        (myRecTypeRefBase & REF_BASE_MASK) | (recType << REC_TYPE_SHIFT);
    return(true);
}

bool GlfRecord::setRefBaseInt(uint8_t refBase)
{
    myRecTypeRefBase = 
        (myRecTypeRefBase & REC_TYPE_MASK) | (refBase & REF_BASE_MASK);
    return(true);
}

// bool GlfRecord::setRefBaseChar(char refBase)
// {

//     uint8_t refBaseInt = REF_BASE_CHAR_TO_INT[refBase];
//     return(setRefBaseInt(refBaseInt));
// }

bool GlfRecord::setOffset(uint32_t offset)
{
    myRec1Base.offset = offset;
    myRec2Base.offset = offset;
    return(true);
}

bool GlfRecord::setMinDepth(uint32_t minDepth)
{
    myRec1Base.min_depth = minDepth;
    myRec2Base.min_depth = minDepth;
    return(true);
}

bool GlfRecord::setMinLk(uint8_t minLk)
{
    setMinDepth((myRec1Base.min_depth & READ_DEPTH_MASK) |
                (minLk << MIN_LK_SHIFT));
    return(true);
}

bool GlfRecord::setReadDepth(uint32_t readDepth)
{
    setMinDepth((myRec1Base.min_depth & MIN_LK_MASK) |
                (readDepth & READ_DEPTH_MASK));
    return(true);
}

bool GlfRecord::setRmsMapQ(uint8_t rmsMapQ)
{
    myRec1Base.rmsMapQ = rmsMapQ;
    myRec2Base.rmsMapQ = rmsMapQ;
    return(true);
}

// Accessors to get the gneric values.
char GlfRecord::getRefBaseChar() const
{
    int index = myRecTypeRefBase & REF_BASE_MASK;
    if((index > REF_BASE_MAX) || (index < 0))
    {
        // TODO throw exception.
        return('N');
    }
    return(REF_BASE_CHAR[index]);
}


uint32_t GlfRecord::getOffset() const
{
    if(getRecordType() == 1)
    {
        return(myRec1Base.offset);
    }
    else if(getRecordType() == 2)
    {
        return(myRec2Base.offset);
    }
    throw(GlfException(GlfStatus::UNKNOWN, 
                       "Tried to call getOffset for Record not of type 1 or 2."));
    return(0);
}

uint32_t GlfRecord::getMinDepth() const
{
    if(getRecordType() == 1)
    {
        return(myRec1Base.min_depth);
    }
    else if(getRecordType() == 2)
    {
        return(myRec2Base.min_depth);
    }
    throw(GlfException(GlfStatus::UNKNOWN, 
                       "Tried to call getMinDepth for Record not of type 1 or 2."));
    return(0);
}

uint8_t GlfRecord::getMinLk() const
{
    if(getRecordType() == 1)
    {
        return(myRec1Base.min_depth >> MIN_LK_SHIFT);
    }
    else if(getRecordType() == 2)
    {
        return(myRec2Base.min_depth >> MIN_LK_SHIFT);
    }
    throw(GlfException(GlfStatus::UNKNOWN, 
                       "Tried to call getMinLk for Record not of type 1 or 2."));
    return(0);
}

uint32_t GlfRecord::getReadDepth() const
{
    if(getRecordType() == 1)
    {
        return(myRec1Base.min_depth & READ_DEPTH_MASK);
    }
    else if(getRecordType() == 2)
    {
        return(myRec2Base.min_depth & READ_DEPTH_MASK);
    }
    throw(GlfException(GlfStatus::UNKNOWN, 
                       "Tried to call getReadDepth for Record not of type 1 or 2."));
    return(0);
}

uint8_t GlfRecord::getRmsMapQ() const
{
    if(getRecordType() == 1)
    {
        return(myRec1Base.rmsMapQ);
    }
    else if(getRecordType() == 2)
    {
        return(myRec2Base.rmsMapQ);
    }
    throw(GlfException(GlfStatus::UNKNOWN, 
                       "Tried to call getRmsMapQ for Record not of type 1 or 2."));
    return(0);
}

    
// Accessors for getting record type 1
bool GlfRecord::setLk(int index, uint8_t value)
{
    if((index < 0) || (index >= NUM_REC1_LIKELIHOOD))
    {
        //  Out of range.
        throw(GlfException(GlfStatus::UNKNOWN, 
                           "Trying to set Record Type 1 likelihood position< 0 or >= 10."));
        return(false);
    }

    // In range.
    myRec1Base.lk[index] = value;
    return(true);
}

uint8_t GlfRecord::getLk(int index)
{
    if(getRecordType() != 1)
    {
        throw(GlfException(GlfStatus::UNKNOWN, 
                           "Tried to call getLk for Record not of type 1."));
        return(0);
    }
    if((index < 0) || (index >= NUM_REC1_LIKELIHOOD))
    {
        throw(GlfException(GlfStatus::UNKNOWN, 
                           "Tried to call getLk for index < 0 or >= 10."));
        return(0);
    }
    return(myRec1Base.lk[index]);
}

    
// Accessors for getting record type 2
bool GlfRecord::setLkHom1(uint8_t lk)
{
    myRec2Base.lkHom1 = lk;
    return(true);
}

bool GlfRecord::setLkHom2(uint8_t lk)
{
    myRec2Base.lkHom2 = lk;
    return(true);
}

bool GlfRecord::setLkHet(uint8_t lk)
{
    myRec2Base.lkHet = lk;
    return(true);
}

bool GlfRecord::setInsertionIndel1(const std::string& indelSeq)
{
    myRec2Base.indelLen1 = indelSeq.length();
    myIndelSeq1 = indelSeq;
    return(true);
}

bool GlfRecord::setDeletionIndel1(const std::string& indelSeq)
{
    myRec2Base.indelLen1 = 0 - (indelSeq.length());
    myIndelSeq1 = indelSeq;
    return(true);
}

bool GlfRecord::setInsertionIndel2(const std::string& indelSeq)
{
    myRec2Base.indelLen2 = indelSeq.length();
    myIndelSeq2 = indelSeq;
    return(true);
}

bool GlfRecord::setDeletionIndel2(const std::string& indelSeq)
{
    myRec2Base.indelLen2 = 0 - (indelSeq.length());
    myIndelSeq2 = indelSeq;
    return(true);
}

uint8_t GlfRecord::getLkHom1()
{
    if(getRecordType() != 2)
    {
        throw(GlfException(GlfStatus::UNKNOWN, 
                           "Tried to call getLkHom1 for Record not of type 2."));
        return(0);
    }
    return(myRec2Base.lkHom1);
}

uint8_t GlfRecord::getLkHom2()
{
    if(getRecordType() != 2)
    {
        throw(GlfException(GlfStatus::UNKNOWN, 
                           "Tried to call getLkHom2 for Record not of type 2."));
        return(0);
    }
    return(myRec2Base.lkHom2);
}

uint8_t GlfRecord::getLkHet()
{
    if(getRecordType() != 2)
    {
        throw(GlfException(GlfStatus::UNKNOWN, 
                           "Tried to call getLkHet for Record not of type 2."));
        return(0);
    }
    return(myRec2Base.lkHet);
}

int16_t GlfRecord::getIndel1(std::string& indelSeq)
{
    if(getRecordType() != 2)
    {
        throw(GlfException(GlfStatus::UNKNOWN, 
                           "Tried to call getIndel1 for Record not of type 2."));
        return(0);
    }
    indelSeq = myIndelSeq1.c_str();
    return(myRec2Base.indelLen1);
}

int16_t GlfRecord::getIndel2(std::string& indelSeq)
{
    if(getRecordType() != 2)
    {
        throw(GlfException(GlfStatus::UNKNOWN, 
                           "Tried to call getIndel2 for Record not of type 2."));
        return(0);
    }
    indelSeq = myIndelSeq2.c_str();
    return(myRec2Base.indelLen2);
}


void GlfRecord::readType1(IFILE filePtr)
{
    // Read record type 1 information.
    int numRead = 0;
    numRead = ifread(filePtr, &myRec1Base, REC1_BASE_SIZE);
    if(numRead != REC1_BASE_SIZE)
    {
        String errorMsg = "Failed to read record of type 1 (";
        errorMsg += REC1_BASE_SIZE;
        errorMsg += " bytes).  Only read ";
        errorMsg += numRead;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
    }

    // Record type 1 is fixed size and has no additional variable length
    // fields, so done reading.
}


void GlfRecord::readType2(IFILE filePtr)
{
    // Read record type 2 information.
    int numRead = 0;
    numRead = ifread(filePtr, &myRec2Base, REC2_BASE_SIZE);
    if(numRead != REC2_BASE_SIZE)
    {
        String errorMsg = "Failed to read record of type 2 base info (";
        errorMsg += REC2_BASE_SIZE;
        errorMsg += " bytes).  Only read ";
        errorMsg += numRead;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
    }

    // Record type 2 has 2 additional variable length fields.  Read those
    // fields.
    int16_t len = abs(myRec2Base.indelLen1);
    numRead = myIndelSeq1.readFromFile(filePtr, len);
    if(numRead != len)
    {
        String errorMsg = "Failed to read record of type 2, 1st indel sequence (";
        errorMsg += len;
        errorMsg += " bytes).  Only read ";
        errorMsg += numRead;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
    }
    len = abs(myRec2Base.indelLen2);
    numRead = myIndelSeq2.readFromFile(filePtr, len);
    if(numRead != len)
    {
        String errorMsg = "Failed to read record of type 2, 2nd indel sequence (";
        errorMsg += len;
        errorMsg += " bytes).  Only read ";
        errorMsg += numRead;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
    }
}


void GlfRecord::writeRtypeRef(IFILE filePtr) const
{
    int byteLen = sizeof(myRecTypeRefBase);
    int numWrite = 
        ifwrite(filePtr, &myRecTypeRefBase, byteLen);
    if(numWrite != byteLen)
    {
        String errorMsg = 
            "Failed to write the length of the record type and reference base (";
        errorMsg += byteLen;
        errorMsg += " bytes).  Only wrote ";
        errorMsg += numWrite;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
    }
}


void GlfRecord::writeType1(IFILE filePtr) const
{
    // Write the generic record field that all records have.
    writeRtypeRef(filePtr);
    
    // Record type 1 is fixed size and has no additional variable length
    // fields, so just write the base info.
    int numWrite = ifwrite(filePtr, &myRec1Base, REC1_BASE_SIZE);
    if(numWrite != REC1_BASE_SIZE)
    {
        // failed to write.
        String errorMsg = "Failed to write record of type 1 (";
        errorMsg += REC1_BASE_SIZE;
        errorMsg += " bytes).  Only wrote ";
        errorMsg += numWrite;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
    }
    

    // Done writing the record.
}


void GlfRecord::writeType2(IFILE filePtr) const
{
    // Write the generic record field that all records have.
    writeRtypeRef(filePtr);

    // Write the record type 2 base info.
    int numWrite = ifwrite(filePtr, &myRec2Base, REC2_BASE_SIZE);
    if(numWrite != REC2_BASE_SIZE)
    {
        // failed to write.
        String errorMsg = "Failed to write record of type 2 base info (";
        errorMsg += REC2_BASE_SIZE;
        errorMsg += " bytes).  Only wrote ";
        errorMsg += numWrite;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
    }
    
    // Record type 2 has 2 additional variable length fields.  Write those
    // fields.
    int len = myIndelSeq1.length();
    numWrite = ifwrite(filePtr, myIndelSeq1.c_str(), len);
    if(numWrite != len)
    {
        // failed to write.
        String errorMsg = "Failed to write record of type 2, 1st indel sequence (";
        errorMsg += len;
        errorMsg += " bytes).  Only wrote ";
        errorMsg += numWrite;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
    }
    len = myIndelSeq2.length();
    numWrite = ifwrite(filePtr, myIndelSeq2.c_str(), len);
    if(numWrite != len)
    {
        // failed to write.
        String errorMsg = "Failed to write record of type 2, 2nd indel sequence (";
        errorMsg += len;
        errorMsg += " bytes).  Only wrote ";
        errorMsg += numWrite;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
    }

    // Done writing the record.
}
