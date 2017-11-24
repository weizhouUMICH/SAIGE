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

#include "GlfHeader.h"
#include "GlfStatus.h"
#include "GlfException.h"
#include "StringBasics.h"

const std::string GlfHeader::GLF_MAGIC = "GLF\3";

GlfHeader::GlfHeader()
    : myText()
{
    resetHeader();
}


GlfHeader::~GlfHeader()
{
    resetHeader();
}


// Copy Constructor   
GlfHeader::GlfHeader(const GlfHeader& header)
    : myText()
{
    copy(header);
}


// Overload operator = to copy the passed in header into this header.
GlfHeader & GlfHeader::operator = (const GlfHeader& header)
{
    copy(header);
    return(*this);
}


bool GlfHeader::copy(const GlfHeader& header)
{
    // Check to see if the passed in value is the same as this.
    if(this == &header)
    {
        return(true);
    }

    resetHeader();

    // Copy the header.
    myText = header.myText;

    return(true);
}


// Reset the header for a new entry, clearing out previous values.
void GlfHeader::resetHeader()
{
   myText.reset();
}


// Read the header from the specified file.  Assumes the file is in
// the correct position for reading the header.
bool GlfHeader::read(IFILE filePtr)
{
    if((filePtr == NULL) || (filePtr->isOpen() == false))
    {
        // File is not open, return failure.
        std::string errorString = 
            "Failed to read the header since the file is not open.";
        throw(GlfException(GlfStatus::FAIL_ORDER, errorString));
        return(false);
    }

    // Read the magic
    int numRead = 0;
    char magic[GLF_MAGIC_LEN];
    numRead = ifread(filePtr, &magic, GLF_MAGIC_LEN);
    if(numRead != GLF_MAGIC_LEN)
    {
        String errorMsg = "Failed to read the magic number (";
        errorMsg += GLF_MAGIC_LEN;
        errorMsg += " bytes).  Only read ";
        errorMsg += numRead;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
        return(false);
    }
    // Read the header length.
    int32_t headerLen = 0;
    int byteLen = sizeof(int32_t);
    numRead = ifread(filePtr, &headerLen, byteLen);
    if(numRead != byteLen)
    {
        String errorMsg = "Failed to read the length of the header text (";
        errorMsg += byteLen;
        errorMsg += " bytes).  Only read ";
        errorMsg += numRead;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
        return(false);
    }
       
    // Read the header from the file.
    numRead = myText.readFromFile(filePtr, headerLen);
    if(numRead != headerLen)
    {
        String errorMsg = "Failed to read the header text (";
        errorMsg += headerLen;
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


// Write the header to the specified file.
bool GlfHeader::write(IFILE filePtr) const
{
    if((filePtr == NULL) || (filePtr->isOpen() == false))
    {
        // File is not open, return failure.
        std::string errorString = 
            "Failed to write the header since the file is not open.";
        throw(GlfException(GlfStatus::FAIL_ORDER, errorString));
        return(false);
    }

    int numWrite = 0;
    // Write the magic
    numWrite = ifwrite(filePtr, GLF_MAGIC.c_str(), GLF_MAGIC_LEN);
    if(numWrite != GLF_MAGIC_LEN)
    {
        String errorMsg = "Failed to write the magic number (";
        errorMsg += GLF_MAGIC_LEN;
        errorMsg += " bytes).  Only wrote ";
        errorMsg += numWrite;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
        return(false);
    }

    // Write the header length.
    int32_t headerLen = myText.length();
    int byteLen = sizeof(int32_t);
    numWrite = ifwrite(filePtr, &headerLen, byteLen);
    if(numWrite != byteLen)
    {
        String errorMsg = "Failed to write the length of the header text (";
        errorMsg += byteLen;
        errorMsg += " bytes).  Only wrote ";
        errorMsg += numWrite;
        errorMsg += " bytes.";
        std::string errorString = errorMsg.c_str();
        throw(GlfException(GlfStatus::FAIL_IO, errorString));
        return(false);
    }
       
    // Write the header to the file.
    numWrite = ifwrite(filePtr, myText.c_str(), headerLen);
    if(numWrite != headerLen)
    {
        String errorMsg = "Failed to write the header text (";
        errorMsg += headerLen;
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


// Set the passed in string to the text string stored in this header.
bool GlfHeader::getHeaderTextString(std::string& text)
{
    text = myText.c_str();
    return(true);
}


// Set the header to the passed in string.
bool GlfHeader::setHeaderTextString(const std::string& text)
{
    myText = text;
    return(true);
}
