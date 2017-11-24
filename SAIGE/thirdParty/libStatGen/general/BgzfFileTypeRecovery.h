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

#ifndef __BGZFFILETYPERECOVERY_H__
#define __BGZFFILETYPERECOVERY_H__

#ifdef __ZLIB_AVAILABLE__

#include "FileType.h"
#include <stdio.h>  // for NULL

class BGZFReader;

class BgzfFileTypeRecovery : public FileType
{
public:
    BgzfFileTypeRecovery()
    {
        bgzfReader = NULL;
    }

    ~BgzfFileTypeRecovery()
    {
        close();
    }

    BgzfFileTypeRecovery(const char * filename, const char * mode);

    // these methods should not be used.  They are
    // misleading because the rhs could be anything,
    // (specifically not a BgzfFileTypeRecover object).
    bool operator == (void * rhs);

    bool operator != (void * rhs);

    // Close the file.
    int close();

    // Reset to the beginning of the file.
    inline void rewind()
    {
        // Just call rewind to move to the beginning of the file.
        seek(0LL, SEEK_SET);
    }

    // Check to see if we have reached the EOF.
    int eof();

    // Check to see if the file is open.
    bool isOpen()
    {
        return (bgzfReader != NULL);
    }

    // Write to the file
    unsigned int write(const void * buffer, unsigned int size);

    // Read into a buffer from the file.  Since the buffer is passed in and
    // this would bypass the fileBuffer used by this class, this method must
    // be protected.
    int read(void * buffer, unsigned int size);

    // Get current position in the file.
    // -1 return value indicates an error.
    int64_t tell();

    // Seek to the specified offset from the origin.
    // origin can be any of the following:
    // Note: not all are valid for all filetypes.
    //   SEEK_SET - Beginning of file
    //   SEEK_CUR - Current position of the file pointer
    //   SEEK_END - End of file
    // Returns true on successful seek and false on a failed seek.
    bool seek(int64_t offset, int origin);

    bool attemptRecoverySync(bool (*checkSignature)(void *data) , int length);

protected:
    // Read via BGZFReader
    BGZFReader* bgzfReader;

};

#endif
#endif
