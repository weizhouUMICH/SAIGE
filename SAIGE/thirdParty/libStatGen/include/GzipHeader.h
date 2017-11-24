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

#ifndef __GZIPHEADER_H__
#define __GZIPHEADER_H__

#include <stdint.h>
#include <stdio.h>
#include "UncompressedFileType.h"

class GzipHeader
{
public:
    GzipHeader();
    ~GzipHeader();

    // Method to read the gzip header from a file.
    // Returns true if the file is a gzip file, false, otherwise.
    bool readHeader(FILE* filePtr);

    // Method to read the gzip header from a file of UncompresedFileType.
    // Returns true if the file is a gzip file, false, otherwise.
    bool readHeader(UncompressedFileType& file);

    // Determine if the file is a gzip file.
    bool isGzipFile();

    // Determine if the file is a BGZF compressed file.
    bool isBgzfFile();

private:

    static const unsigned int GZIP_HEADER_SIZE = 18;

    union
    {
        struct
        {
            uint8_t  id1;
            uint8_t  id2;
            uint8_t  cm;
            uint8_t  flg;
            uint32_t mtime;
            uint8_t  xfl;
            uint8_t  os;
            uint16_t xlen;
            uint8_t si1;
            uint8_t si2;
            uint16_t slen;
            uint16_t bsize;
        };
        char headerBuffer[GzipHeader::GZIP_HEADER_SIZE];
    };
    char buffer[GZIP_HEADER_SIZE];

};


#endif
