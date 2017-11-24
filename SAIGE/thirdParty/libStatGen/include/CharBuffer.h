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

#ifndef __BUFFER_H__
#define __BUFFER_H__

#include <stdint.h>
#include "InputFile.h"

class CharBuffer
{
public:
    CharBuffer();
    CharBuffer(int32_t initialSize);
    ~CharBuffer();

    // Copy Constructor
    CharBuffer(const CharBuffer& buffer);

    // Overload operator = to copy the passed in buffer into this buffer.
    CharBuffer& operator = (const CharBuffer& buffer);

    // Overload operator = to copy the passed in buffer into this buffer.
    CharBuffer& operator = (const std::string& stringBuffer);

    // Overload operator = to copy the passed in buffer into this buffer.
    bool copy(const CharBuffer& buffer);

    void reset();

    // Read from a file into the buffer.  length is the amount of data to read.
    // Returns the number of bytes read.
    int readFromFile(IFILE filePtr, int32_t length);

    inline const char* c_str() const
    {
        return(myBuffer);
    }

    inline int32_t length() const
    {
        return(myBufferLen);
    }

private:
    // newLen is the new length for the buffer.
    bool prepareNewLength(int32_t newLen);

    int32_t myBufferLen;
    char* myBuffer;

    int32_t myBufferAllocatedLen;

    static const int32_t DEFAULT_BUFFER_SIZE = 100;
};

#endif

