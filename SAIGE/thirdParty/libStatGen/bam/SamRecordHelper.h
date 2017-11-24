/*
 *  Copyright (C) 2012  Regents of the University of Michigan
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

#ifndef __SAM_RECORD_HELPER_H__
#define __SAM_RECORD_HELPER_H__

#include "SamRecord.h"

/// Class for extracting information from a SAM Flag.
class SamRecordHelper
{
public:

    /// Helper method that checks if the record's read sequence starting
    /// at the specified 0-based reference position matches the passed in
    /// sequence.
    /// \return returns -1 if it does not match,
    /// returns the cycle (read position) of pos0Based if it does match.
    static int checkSequence(SamRecord& record, int32_t pos0Based, 
                              const char* sequence);
    
    /// Helper to append the SAM string representation of all the tags to 
    /// the specified string.  Does NOT add a preceding delimiter before the
    /// first tag.
    /// \param record record whose tags to append.
    /// \param returnString string to append the tags to.
    /// \param delim delimiter to use to separate different tags.
    /// \return true on success, false on failure/partial generation.
    static bool genSamTagsString(SamRecord& record, String& returnString,
                                 char delim = '\t');


    /// Helper to append the SAM string representation of the specified tag to 
    /// the specified string.
    /// \param tag the tag name.
    /// \param vtype the vtype.
    /// \param value pointer to the value of the tag (will be cast
    /// to int, double, char, or string based on vtype).
    /// \param returnString string to append the tag to.
    /// \return true on success, false on failure/partial generation.
    static bool genSamTagString(const char* tag, char vtype, 
                                void* value, String& returnString);

private:
    SamRecordHelper();
};


#endif
