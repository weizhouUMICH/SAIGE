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

#ifndef __UNCOMPRESSEDFILETYPE_H__
#define __UNCOMPRESSEDFILETYPE_H__

#include <iostream>
#include <stdio.h>
#include "FileType.h"
#include "knetfile.h"

class UncompressedFileType : public FileType
{
public:
    UncompressedFileType()
    {
        filePtr = NULL;
        kfilePtr = NULL;
        keof = false;
    }

    virtual ~UncompressedFileType()
    {
        if((filePtr != NULL) || (kfilePtr != NULL))
        {
            close();
        }
    }

    UncompressedFileType(const char * filename, const char * mode);

    bool operator == (void * rhs)
    {
        // No two file pointers are the same, so if rhs is not NULL, then
        // the two pointers are different (false).
        if (rhs != NULL)
            return false;
        // rhs is NULL.  They are the same if both filePtr & kfilePtr are NULL.
        return((filePtr == rhs) && (kfilePtr == rhs));
    }

    bool operator != (void * rhs)
    {
        // No two file pointers are the same, so if rhs is not NULL, then
        // the two pointers are different (true).
        if (rhs != NULL)
            return true;
        // rhs is NULL.  They are the different if either filePtr or kfilePtr
        // are not NULL.
        return((filePtr != rhs) || (kfilePtr != rhs));
    }

    // Close the file.
    inline int close()
    {
        if(filePtr != NULL)
        {
            if((filePtr != stdout) && (filePtr != stdin))
            {
                int result = fclose(filePtr);
                filePtr = NULL;
                return result;
            }
            filePtr = NULL;
        }
        else if(kfilePtr != NULL)
        {
            int result = knet_close(kfilePtr);
            kfilePtr = NULL;
            return result;
        }
        return 0;
    }


    // Reset to the beginning of the file.
    inline void rewind()
    {
        // Just call rewind to move to the beginning of the file.
        if(filePtr != NULL)
        {
            ::rewind(filePtr);
        }
        else if (kfilePtr != NULL)
        {
            knet_seek(kfilePtr, 0, SEEK_SET);
        }
    }

    // Check to see if we have reached the EOF.
    inline int eof()
    {
        //  check the file for eof.
        if(kfilePtr != NULL)
        {
            return(keof);
        }
        else
        {
            return feof(filePtr);
        }
    }

    // Check to see if the file is open.
    virtual inline bool isOpen()
    {
        if((filePtr != NULL) || (kfilePtr != NULL))
        {
            // filePtr is not null, so the file is open.
            return(true);
        }
        return(false);
    }

    // Write to the file
    inline unsigned int write(const void * buffer, unsigned int size)
    {
        // knetfile is never used for writing.
        return fwrite(buffer, 1, size, filePtr);
    }

    // Read into a buffer from the file.  Since the buffer is passed in and
    // this would bypass the fileBuffer used by this class, this method must
    // be protected.
    inline int read(void * buffer, unsigned int size)
    {
        if(kfilePtr != NULL)
        {
            int bytesRead = knet_read(kfilePtr, buffer, size);
            if((bytesRead == 0) && (size != 0))
            {
                keof = true;
            }
            else if((bytesRead != (int)size) && (bytesRead >= 0))
            {
                // Less then the requested size was read and an error
                // was not returned (bgzf_read returns -1 on error).
                keof = true;
            }
            else
            {
                keof = false;
            }
            return(bytesRead);
        }
        return fread(buffer, 1, size, filePtr);
    }


    // Get current position in the file.
    // -1 return value indicates an error.
    virtual inline int64_t tell()
    {
        if(kfilePtr != NULL)
        {
            return knet_tell(kfilePtr);
        }
        return ftell(filePtr);
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
        int returnVal = 0;
        if(kfilePtr != NULL)
        {
            returnVal = knet_seek(kfilePtr, offset, origin);
            keof = false;
        }
        else
        {
            returnVal = fseek(filePtr, offset, origin);
        }
        // Check for success - 0 return value.
        if (returnVal == 0)
        {
            return true;
        }
        // Successful.
        return false;
    }


protected:
    // A FILE Pointer is used.
    FILE* filePtr;
    knetFile *kfilePtr;
    bool keof;
};

#endif


