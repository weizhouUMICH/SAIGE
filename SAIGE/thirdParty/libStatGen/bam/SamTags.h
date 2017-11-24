/*
 *  Copyright (C) 2010-2011  Regents of the University of Michigan
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

#ifndef __SAM_TAGS_H__
#define __SAM_TAGS_H__

#include <stdint.h>
#include <stdexcept>
#include "SamRecord.h"

/// Class for parsing/creating/operating on SAM/BAM record tags.
class SamTags
{
public:
    ///////////////////////
    /// @name  Constants for parsing tags.
    //@{
    static const char* BQ_TAG;
    static const char BQ_TAG_TYPE;
    static const char* MD_TAG;
    static const char MD_TAG_TYPE;
    static const char* ORIG_POS_TAG;
    static const char ORIG_POS_TAG_TYPE;
    static const char* ORIG_CIGAR_TAG;
    static const char ORIG_CIGAR_TAG_TYPE;
    static const char* ORIG_QUAL_TAG;
    static const char ORIG_QUAL_TAG_TYPE;
    //@}

    /// Create the MD tag for the specified input record and the genome.
    /// \return returns true if an MD tag was created, false if one could not
    /// be created.
    static bool createMDTag(String& outputMDtag, SamRecord& inputRec, GenomeSequence& genome);
    /// Check to see if the MD tag in the record is accurate.
    static bool isMDTagCorrect(SamRecord& inputRec, GenomeSequence& genome);
    // Update/Add the MD tag in the inputRec.
    static bool updateMDTag(SamRecord& inputRec, GenomeSequence& genome);

private:
    SamTags();
};


#endif
