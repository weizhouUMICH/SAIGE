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

#ifndef __SAM_REFERENCE_INFO_H__
#define __SAM_REFERENCE_INFO_H__

#include "StringArray.h"
#include "StringHash.h"
#include "IntArray.h"

/// Class for tracking the reference information mapping between the
/// reference ids and the reference names.
class SamReferenceInfo
{
public:
    /// Constructor.
    SamReferenceInfo();
    /// Destructor.
    ~SamReferenceInfo();
    /// Add reference sequence name and reference sequence length.
    void add(const char* referenceSequenceName, 
             int32_t referenceSequenceLength);

    /// Get the reference ID for the specified name, if addID is set to true,
    /// a reference id will be created for the referenceName if one does not
    /// already exist, while if addID is set to false (default), it will return
    /// NO_REF_ID if the reference name does not exist.
    int getReferenceID(const String & referenceName, bool addID = false);
    /// Get the reference ID for the specified name, if addID is set to true,
    /// a reference id will be created for the referenceName if one does not
    /// already exist, while if addID is set to false (default), it will return
    /// NO_REF_ID if the reference name does not exist.
    int getReferenceID(const char* referenceName, bool addID = false);
    /// Get the reference name for the specified id, if the id is not found,
    /// return "*".
    const String & getReferenceLabel(int id) const;

    /// Get the number of entries contained here.
    int32_t getNumEntries() const;

    /// Return the reference name at the specified index, returning "" if the
    /// index is out of bounds.
    const char* getReferenceName(int index) const;
    
    /// Return the reference length at the specified index, returning 0 if the
    /// index is out of bounds.
    int32_t getReferenceLength(int index) const;

    /// Reset this reference info.
    void clear();

    /// Copy the reference information.
    SamReferenceInfo & operator = (const SamReferenceInfo & rhs);

    bool operator== (const SamReferenceInfo& rhs) const;
    bool operator!= (const SamReferenceInfo& rhs) const
    {
        return(!operator==(rhs));
    }

    /// Constant for the value returned if a reference id does not exist
    /// for a queried reference name.
    static const int NO_REF_ID = -3;

private:
    // Reference Name information
    StringArray    myReferenceContigs;
    StringIntHash  myReferenceHash;
    IntArray       myReferenceLengths;
};

#endif

