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

#ifndef __SAMHEADER_RECORD_H__
#define __SAMHEADER_RECORD_H__

#include "StringArray.h"
#include "StringHash.h"
#include "SamHeaderTag.h"

/// This class encapsulates the tag value pairs contained with a SAM Header
/// line with accessors for getting and setting the tags within this header.
class SamHeaderRecord
{
public:
    /// Specifies the Type for the sam header record (line).
    enum SamHeaderRecordType {
        HD, ///< Header
        SQ, ///< Sequence Dictionary
        RG, ///< Read Group
        PG  ///< Program
    };

    /// Constructor
    SamHeaderRecord();
   
    /// Destructor
    virtual ~SamHeaderRecord();

    /// Return a pointer to a newly created header record of the appropriate type
    /// that is a copy of this record. The newly created record will not be
    /// deleted by this class and it is the responsibility of the calling method
    /// to handle the deletion.
    /// Returns NULL on failure to copy.
    virtual SamHeaderRecord* createCopy() const = 0;
    
    /// Set the fields from the passed in line.
    /// Return true if successfully set.
    bool setFields(const StringArray& tokens);

    /// Check to see if the record is valid.
    bool isValid();

    /// Return the value associated with the specified tag.  Returns "" if it
    /// is not set.
    const char* getTagValue(const char* tag) const;

    /// Set the value of the specified tag to the specified value, deletes
    /// the tag when value is NULL.
    /// Returns whether or not it was successful, fails if tag is the key tag
    /// and the key tag already exists.
    bool setTag(const char* tag, const char* value);

    /// Reset this header record to an empty state with no tags.
    void reset();

    /// Appends the string representation of this header record
    /// to the passed in string.
    bool appendString(std::string& header);

    /// Add the key tag with the specified value (not for HD headers).
    bool addKey(const char* value);

    /// Get the value associated with the key tag.  Returns "" if it is not set.
    const char* getKeyValue() const;

    /// This record is active (true) if there is at least one tag set.
    bool isActiveHeaderRecord();

    /// Return the type of this header record (HD, SQ, RG, or PG) as a string.
    const char* getTypeString();

    /// Return the type of this header record (HD, SQ, RG, or PG) as an enum.
    SamHeaderRecordType getType();

protected:
    void addRequiredTag(const char* requiredTag);

    // Copy this record into the specified new one.
    virtual void internalCopy(SamHeaderRecord& newRec) const;

    // The type for this header record.
    std::string myTypeString;

    // The type for this header record.
    SamHeaderRecordType myType;

    // The TAG name that is the key for this record
    // Only applicable if more than one of this type
    // of record is allowed.
    std::string myKeyTag;

private:
    SamHeaderRecord(const SamHeaderRecord& samHeaderRecord);
    SamHeaderRecord& operator=(const SamHeaderRecord& samHeaderRecord);

    // hash from tag name to index into the tag values vector.
    StringIntHash myTagHash;
    std::vector<SamHeaderTag*> myTags;

    // The tags that are required for this record.
    std::vector<String> myRequiredTags;

    int myNumActiveTags;
};

#endif
