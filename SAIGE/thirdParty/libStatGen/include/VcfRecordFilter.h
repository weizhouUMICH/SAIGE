/*
 *  Copyright (C) 2011  Regents of the University of Michigan,
 *                      Hyun Min Kang, Matthew Flickenger, Matthew Snyder,
 *                      and Goncalo Abecasis
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


#ifndef __VCF_RECORD_FILTER_H__
#define __VCF_RECORD_FILTER_H__

#include "InputFile.h"
#include "ReusableVector.h"

/// This header file provides interface to read/write VCF files.
class  VcfRecordFilter
{
public:
    /// Default Constructor, initializes the variables.
    VcfRecordFilter() {}
    /// Destructor
    virtual ~VcfRecordFilter() {}
    
    /// Whether or not to initially parse the filter as it is read
    /// from the file.
    static void parseRecordFilter(bool parse) {ourParseFilter = parse;}

    /// Read this field from the file up until the next \t,\n, or EOF.
    /// Reads the \t, \n, or EOF.
    /// \param filePtr IFILE to read from.
    /// \return true if the field was successfully read from the specified
    /// filePtr, false if not.
    bool read(IFILE filePtr);

    void reset();

    /// Return true if all filters were passed (contents is PASS).
    bool passedAllFilters();

    /// Get the entire filter string returned as a const reference.
    const std::string& getString();

    /// Get the number of filters in this string.
    int getNumFilters();

    /// Get the filter at the specified index (starting at 0).
    /// If the index is out of range, an empty string is returned.
    const std::string& getString(int index);

    void clear() {reset();}
    void setFilter(const char* filter);

    void addFilter(const char* filter);


protected: 

private:
    VcfRecordFilter(const VcfRecordFilter& vcfRecordFilter);
    VcfRecordFilter& operator=(const VcfRecordFilter& vcfRecordFilter);

    static bool ourParseFilter;
    static const char FILTER_DELIM = ';';

    std::string myFilterString;
    ReusableVector<std::string> myFilterVector;
};

#endif
