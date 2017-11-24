/*
 *  Copyright (C) 2011-2012  Regents of the University of Michigan,
 *                           Hyun Min Kang, Matthew Flickenger, Matthew Snyder,
 *                           and Goncalo Abecasis
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


#ifndef __VCF_RECORD_INFO_H__
#define __VCF_RECORD_INFO_H__

#include <list>
#include <utility>

#include "VcfRecordField.h"
#include "ReusableVector.h"

/// This header file provides interface to read/write VCF files.
class VcfRecordInfo : public VcfRecordField
{
public:
    /// Default Constructor, initializes the variables.
    VcfRecordInfo();

    /// Destructor
    virtual ~VcfRecordInfo();
    
    /// Read this info field from the file up until the next \t,\n, or EOF.
    /// \param filePtr IFILE to read from.
    /// \return true if a tab ended the field, false if it was \n or EOF.
    bool read(IFILE filePtr);

    /// Write the info field to the file, without printing the
    // starting/trailing '\t'.
    /// \param filePtr IFILE to write to.
    /// \return true if the field was successfully written to the specified
    ///  filePtr, false if not.
    bool write(IFILE filePtr);

    /// reset the field for a new entry.
    void reset();
    /// reset the field for a new entry.
    void clear() {reset();}

    int getNumInfoFields() const { return(myInfo.size()); }

    /// Set the string value associated with the specified key.  
    /// \param key to set the value for.
    /// \param const pointer to the string value for this key, NULL if
    /// the key was not found, a pointer to an empty string if the key
    /// was found, but does not have a value.
    void setString(const char* key, const char* stringVal);

    /// Get a pointer to the string containing the value associated with the
    /// specified key (the pointer will be invalid if the field is
    /// changed/reset).  
    /// \param key to find the value for.
    /// \return const pointer to the string value for this key, NULL if
    /// the key was not found, a pointer to an empty string if the key
    /// was found, but does not have a value.
    const std::string* getString(const char* key);

    /// Get a pointer to the string containing the value associated with the
    /// specified info index (the pointer will be invalid if the field is
    /// changed/reset).  
    /// \param index to get the value for.
    /// \return const pointer to the string value for this index, NULL if
    /// the index is out of range, a pointer to an empty string if the index
    /// is in range, but does not have a value.
    const std::string* getString(int index);

    /// Get a reference to the InfoElement containing the value associated with the
    /// specified info index (the pointer will be invalid if the field is
    /// changed/reset).
    /// \param index to get the value for.
    /// \return const references to the InfoElement for this index, index
    /// must be in range.
    std::pair<std::string, std::string> getInfoPair(int index) const;



protected:

private:
    VcfRecordInfo(const VcfRecordInfo& vcfRecordInfo);
    VcfRecordInfo& operator=(const VcfRecordInfo& vcfRecordInfo);

    static const char EMPTY_INFO = '.';

    class InfoElement
    {
    public:
        InfoElement() {key.clear(); value.clear();}

        std::string key;
        std::string value;

        void clear() {key.clear(); value.clear();}
    };

    ReusableVector<InfoElement> myInfo;
};


#endif
