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

#include "GzipHeader.h"
#include <iostream>

#include <cstring>

// Constructor to initialize member data to 0.
GzipHeader::GzipHeader()
{
    // clear the union via memset:
    memset(headerBuffer, 0, sizeof(headerBuffer));
}


// Desctructor - nothing to do.
GzipHeader::~GzipHeader()
{
}


// Method to read the gzip header from a file.
// Returns true if the file is a gzip file, false, otherwise.
bool GzipHeader::readHeader(FILE* filePtr)
{
    bool isGzip = false;

    // If the file is not already open, return false.
    if (filePtr == NULL)
    {
        // File is not open, so return false - not a gzip file.
        return(false);
    }

    // Try to read a header from the file.
    //   if(144 == fread(buffer, 1, 144, filePtr))
    if (GZIP_HEADER_SIZE == fread(buffer, 1, GZIP_HEADER_SIZE, filePtr))
    {
        memcpy(headerBuffer, buffer, GZIP_HEADER_SIZE);

        // Successfully read enough bytes, so check to see if it is a GzipFile.
        if (isGzipFile())
        {
            // It is a gzip file.
            isGzip = true;
        }
    }

    return isGzip;
}


// Method to read the gzip header from a file.
// Returns true if the file is a gzip file, false, otherwise.
bool GzipHeader::readHeader(UncompressedFileType& file)
{
    bool isGzip = false;

    // If the file is not already open, return false.
    if (!file.isOpen())
    {
        // File is not open, so return false - not a gzip file.
        return(false);
    }

    // Try to read a header from the file.
    //   if(144 == file.read(buffer, 1, 144, filePtr))
    if ((int)GZIP_HEADER_SIZE == file.read(buffer, GZIP_HEADER_SIZE))
    {
        memcpy(headerBuffer, buffer, GZIP_HEADER_SIZE);

        // Successfully read enough bytes, so check to see if it is a GzipFile.
        if (isGzipFile())
        {
            // It is a gzip file.
            isGzip = true;
        }
    }

    return isGzip;
}


// Determine if the file is a gzip file.
bool GzipHeader::isGzipFile()
{
    if ((id1 == 31) && (id2 == 139))
    {
        return true;
    }
    return false;
}


// Determine if the file is a BGZF compressed file.
bool GzipHeader::isBgzfFile()
{
    if (isGzipFile() && (si1 == 66) && (si2 == 67))
    {
        return true;
    }
    return false;
}

