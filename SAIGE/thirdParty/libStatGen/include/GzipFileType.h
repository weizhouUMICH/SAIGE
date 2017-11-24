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

#ifndef __GZFILETYPE_H__
#define __GZFILETYPE_H__

#ifdef  __ZLIB_AVAILABLE__

#if defined(_WIN32)
#include <stdio.h>  // for NULL!
#endif

#include <stdlib.h>
#include <zlib.h>
#include "FileType.h"

//#include <iostream>

class GzipFileType : public FileType
{
public:
    GzipFileType()
    {
        gzHandle = NULL;
    }

    virtual ~GzipFileType()
    {
    }

    GzipFileType(const char * filename, const char * mode);

    bool operator == (void * rhs)
    {
        // No two file pointers are the same, so if rhs is not NULL, then
        // the two pointers are different (false).
        if (rhs != NULL)
            return false;
        return (gzHandle == rhs);
    }

    bool operator != (void * rhs)
    {
        // No two file pointers are the same, so if rhs is not NULL, then
        // the two pointers are different (true).
        if (rhs != NULL)
            return true;
        return (gzHandle != rhs);
    }

    // Close the file.
    inline int close()
    {
        int result = gzclose(gzHandle);
        gzHandle = NULL;
        return result;
    }


    // Reset to the beginning of the file.
    inline void rewind()
    {
        // Just call rewind to move to the beginning of the file.
        gzrewind(gzHandle);
    }

    // Check to see if we have reached the EOF.
    inline int eof()
    {
        //  check the file for eof.
        return gzeof(gzHandle);
    }

    // Check to see if the file is open.
    virtual inline bool isOpen()
    {
        if (gzHandle != NULL)
        {
            // gzHandle is not null, so the file is open.
            return(true);
        }
        return(false);
    }

    // Write to the file
    inline unsigned int write(const void * buffer, unsigned int size)
    {
        return gzwrite(gzHandle, buffer, size);
    }

    // Read into a buffer from the file.  Since the buffer is passed in and
    // this would bypass the fileBuffer used by this class, this method must
    // be protected.
    inline int read(void * buffer, unsigned int size)
    {
        unsigned int numBytes = gzread(gzHandle, buffer, size);
//       if(numBytes != size)
//       {
//          std::cerr << "Error reading.  Read "<< numBytes << " instead of "<< size<<std::endl;
//          int error_code = 0;
//          const char* errorMsg = gzerror(gzHandle, &error_code);
//          std::cerr << "ERROR Code: " << error_code << ";  Error Msg: " << errorMsg << std::endl;
//       }
        return numBytes;
    }

    // Get current position in the file.
    // -1 return value indicates an error.
    virtual inline int64_t tell()
    {
        return gztell(gzHandle);
    }


    // Seek to the specified offset from the origin.
    // origin can be any of the following:
    // Note: not all are valid for all filetypes.
    //   SEEK_SET - Beginning of file
    //   SEEK_CUR - Current position of the file pointer
    //   SEEK_END - End of file
    // Returns true on successful seek and false on a failed seek.
    virtual inline bool seek(int64_t offset, int origin)
    {
        int64_t returnVal = gzseek(gzHandle, offset, origin);
        // Check for failure.
        if (returnVal == -1)
        {
            return false;
        }
        // Successful.
        return true;
    }


protected:
    // A gzFile is used.
    gzFile gzHandle;
};

#endif

#endif

