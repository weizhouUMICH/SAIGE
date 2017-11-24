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
#include "CharBuffer.h"

CharBuffer::CharBuffer()
    : myBuffer(NULL)
{
    myBuffer = (char *) malloc(DEFAULT_BUFFER_SIZE);
    myBufferAllocatedLen = DEFAULT_BUFFER_SIZE;
    reset();
}


CharBuffer::CharBuffer(int32_t initialSize)
    : myBuffer(NULL)
{
    myBuffer = (char *) malloc(initialSize);
    myBufferAllocatedLen = DEFAULT_BUFFER_SIZE;

    reset();
}


CharBuffer::~CharBuffer()
{
    reset();
    if(myBuffer != NULL)
    {
        free(myBuffer);
        myBuffer = NULL;
    }
}


// Copy Constructor   
CharBuffer::CharBuffer(const CharBuffer& buffer)
    : myBuffer(NULL)
{
    myBuffer = 
        (char *) malloc(DEFAULT_BUFFER_SIZE);
    myBufferAllocatedLen = DEFAULT_BUFFER_SIZE;

    reset();

    copy(buffer);
}


// Overload operator = to copy the passed in buffer into this buffer.
CharBuffer& CharBuffer::operator = (const CharBuffer& buffer)
{
    copy(buffer);
    return(*this);
}


// Overload operator = to copy the passed in buffer into this buffer.
CharBuffer& CharBuffer::operator = (const std::string& stringBuffer)
{
    // First check lengh
    if(prepareNewLength(stringBuffer.length()))
    {
        memcpy(myBuffer, stringBuffer.c_str(), stringBuffer.length());
    }
    // TODO: on failure of prepareNewLength, should it throw an exception?
    
    return(*this);
}


bool CharBuffer::copy(const CharBuffer& buffer)
{
    // Check to see if the passed in value is the same as this.
    if(this == &buffer)
    {
        return(true);
    }

    // Copy the buffer.
    // First check lengh
    prepareNewLength(buffer.myBufferLen);

    memcpy(myBuffer, buffer.myBuffer, buffer.myBufferLen);
    myBufferLen = buffer.myBufferLen;

    return(true);
}


// Reset the buffer for a new entry, clearing out previous values.
void CharBuffer::reset()
{
    myBufferLen = 0;
    if(myBuffer != NULL)
    {
        myBuffer[0] = 0;
    }
}


// Read from a file into the buffer.  length is the amount of data to read.
// Returns the number of bytes read.
int CharBuffer::readFromFile(IFILE filePtr, int32_t length)
{
    if(filePtr == NULL)
    {
        return(0);
    }

    if(prepareNewLength(length))
    {
        return(ifread(filePtr, myBuffer, length));
    }
    // failed to setup the buffer, return false.
    return(false);
}


// newLen is the new length that this buffer needs to be.
bool CharBuffer::prepareNewLength(int32_t newLen)
{
    if(newLen < 0)
    {
        // Invalid length.
        return(false);
    }
    
    // myBufferAllocatedLen must be bigger than new length, because the
    // newLen position is set to 0.
    if(myBufferAllocatedLen <= newLen)
    {
        // Not enough space is allocated, so allocate more space.
        char* tmpBufferPtr = (char *)realloc(myBuffer, newLen);
        if(tmpBufferPtr == NULL)
        {
            // FAILED to allocate memory
            fprintf(stderr, "FAILED TO ALLOCATE MEMORY!!!");
            // myStatus.addError(GlfStatus::FAIL_MEM, "Failed Memory Allocation.");
            return(false);
        }
        // Successfully allocated memory, so set myRecordPtr.
        myBuffer = tmpBufferPtr;
        myBufferAllocatedLen = newLen;
    }
    myBufferLen = newLen;
    myBuffer[newLen] = 0;
    return(true);
}

