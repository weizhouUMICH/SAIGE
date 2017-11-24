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

#include "GlfRefSection.h"
#include "GlfException.h"
#include "StringBasics.h"

GlfRefSection::GlfRefSection()
    : myRefName()
{
    resetRefSection();
}


GlfRefSection::~GlfRefSection()
{
    resetRefSection();
}


// Copy Constructor   
GlfRefSection::GlfRefSection(const GlfRefSection& refSection)
    : myRefName()
{
    copy(refSection);
}


// Overload operator = to copy the passed in refSection into this refSection.
GlfRefSection & GlfRefSection::operator = (const GlfRefSection& refSection)
{
    copy(refSection);
    return(*this);
}


bool GlfRefSection::copy(const GlfRefSection& refSection)
{
    // Check to see if the passed in value is the same as this.
    if(this == &refSection)
    {
        return(true);
    }

    resetRefSection();

    // Copy the refSection.
    myRefName = refSection.myRefName;
    myRefLen = refSection.myRefLen;

    return(true);
}


// Reset the refSection for a new entry, clearing out previous values.
void GlfRefSection::resetRefSection()
{
    myRefName.reset();
    myRefLen = 0;
}


// Read the refSection from the specified file.  Assumes the file is in
// the correct position for reading the refSection.
bool GlfRefSection::read(IFILE filePtr)
{
    // Read the reference sequence name length
    int numRead = 0;
    int32_t refNameLen = 0;
    int byteLen = sizeof(int32_t);
    numRead = ifread(filePtr, &refNameLen, byteLen);
    if(numRead != byteLen)
    {
        // If no bytes were read and it is the end of the file, then return 
        // false, but do not throw an exception.  This is not an error, just
        // the end of the file.
        if((numRead == 0) && ifeof(filePtr))
        {
            return(false);
        }

        String errorMsg = 
            "Failed to read the length of the reference sequence name (";
        errorMsg += byteLen;
        errorMsg += " bytes).  Only read ";
        errorMsg += numRead;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
        return(false);
    }
       
    // Read the refSection from the file.
    numRead = myRefName.readFromFile(filePtr, refNameLen);
    if(numRead != refNameLen)
    {
        String errorMsg = "Failed to read the reference sequence name (";
        errorMsg += refNameLen;
        errorMsg += " bytes).  Only read ";
        errorMsg += numRead;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
        return(false);
    }

    // Read the ref length.
    byteLen = sizeof(uint32_t);
    numRead = ifread(filePtr, &myRefLen, byteLen);
    if(numRead != byteLen)
    {
        String errorMsg = "Failed to read the reference sequence length (";
        errorMsg += byteLen;
        errorMsg += " bytes).  Only read ";
        errorMsg += numRead;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
        return(false);
    }

    // Successfully read, return success.
    return(true);
}


// Write the refSection to the specified file.
bool GlfRefSection::write(IFILE filePtr) const
{
    int refNameLen = myRefName.length();
    int byteLen = sizeof(int32_t);
    int numWrite = ifwrite(filePtr, &refNameLen, byteLen);
    if(numWrite != byteLen)
    {
        String errorMsg = 
            "Failed to write the length of the reference sequence name (";
        errorMsg += byteLen;
        errorMsg += " bytes).  Only wrote ";
        errorMsg += numWrite;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
        return(false);
    }

    numWrite = ifwrite(filePtr, myRefName.c_str(), refNameLen);
    if(numWrite != refNameLen)
    {
        String errorMsg = "Failed to write the reference sequence name (";
        errorMsg += refNameLen;
        errorMsg += " bytes).  Only wrote ";
        errorMsg += numWrite;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
        return(false);
    }

    // Write the length of the reference sequence
    byteLen = sizeof(uint32_t);
    numWrite = ifwrite(filePtr, &myRefLen, byteLen);
    if(numWrite != byteLen)
    {
        String errorMsg = "Failed to write the reference sequence length (";
        errorMsg += byteLen;
        errorMsg += " bytes).  Only wrote ";
        errorMsg += numWrite;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
        return(false);
    }

    // Successfully wrote, return success.
    return(true);
}


bool GlfRefSection::getName(std::string& name) const
{
    name = myRefName.c_str();
    return(true);
}


uint32_t GlfRefSection::getRefLen() const
{
    return(myRefLen);
}


bool GlfRefSection::setName(const std::string& name)
{
    myRefName = name;
    return(true);
}


bool GlfRefSection::setRefLen(uint32_t refLen)
{
    myRefLen = refLen;
    return(true);
}


void GlfRefSection::print() const
{
    std::cout << "l_name: " << myRefName.length() 
              << "; name: " << myRefName.c_str()
              << "; ref_len: " << myRefLen
              << "\n";
}

