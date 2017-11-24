/*
 *  Copyright (C) 2012-2013  Regents of the University of Michigan
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

#ifndef __TABIX_H__
#define __TABIX_H__

#include <stdint.h>
#include <vector>
#include <map>
#include <stdlib.h>

#include "IndexBase.h"

#include "InputFile.h"
#include "StatGenStatus.h"

class Tabix : public IndexBase
{
public:

    enum Format
        { 
            FORMAT_GENERIC = 0,
            FORMAT_SAM = 1,
            FORMAT_VCF = 2
        };

    Tabix();
    virtual ~Tabix();

    /// Reset the member data for a new index file.
    void resetIndex();

    // Read & parse the specified index file.
    /// \param filename the bam index file to be read.
    /// \return the status of the read.
    StatGenStatus::Status readIndex(const char* filename);

    /// Get the starting file offset to look for the specified start position.
    /// For an entire reference ID, set start to -1.
    /// To start at the beginning of the region, set start to 0/-1.
    bool getStartPos(const char* refName, int32_t start,
                     uint64_t& fileStartPos) const;

    /// Return the reference name at the specified index or
    /// throws an exception if out of range.
    const char* getRefName(unsigned int indexNum) const;

    // Get the format of this tabix file.
    inline int32_t getFormat() const { return myFormat.format; }

private:
    struct TabixFormat
    {
        int32_t format;
        int32_t col_seq;
        int32_t col_beg;
        int32_t col_end;
        int32_t meta; // character that starts header lines
        int32_t skip; // Number of lines to skip from putting into the index.
    };

    TabixFormat myFormat;

    char* myChromNamesBuffer;

    // vector pointing to the chromosome names.
    std::vector<const char*> myChromNamesVector;
};


#endif
