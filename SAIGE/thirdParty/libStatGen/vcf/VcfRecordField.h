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


#ifndef __VCF_RECORD_FIELD_H__
#define __VCF_RECORD_FIELD_H__

#include "InputFile.h"

/// This header file provides interface to read/write VCF files.
class  VcfRecordField
{
public:
    /// Default Constructor, initializes the variables.
    VcfRecordField() {}
    /// Destructor
    virtual ~VcfRecordField() {}
    
    /// Read this field from the file up until the next \t,\n, or EOF.
    /// Reads the \t, \n, or EOF.
    /// \param filePtr IFILE to read from.
    /// \return true if the field was successfully read from the specified
    /// filePtr, false if not.
    virtual bool read(IFILE filePtr) = 0;

    /// Write this field to the file, without printing the
    // starting/trailing '\t'.
    /// \return true if the field was successfully written to the specified
    /// filePtr, false if not.
    virtual bool write(IFILE filePtr) = 0;

protected: 

private:
    VcfRecordField(const VcfRecordField& field);
    VcfRecordField& operator=(const VcfRecordField& field);
};

#endif
