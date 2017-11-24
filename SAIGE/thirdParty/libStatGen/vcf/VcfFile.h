/*
 *  Copyright (C) 2010-2011  Regents of the University of Michigan,
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

#ifndef __VCF_FILE_H__
#define __VCF_FILE_H__

#include "InputFile.h"
#include "VcfHeader.h"

/// This header file provides interface to read/write VCF files.
class VcfFile {
public:
    /// Default Constructor, initializes the variables, but does not open
    /// any files.
    VcfFile();
    /// Destructor
    virtual ~VcfFile();

    /// Open the vcf file with the specified filename,
    /// overwritten by child classes for read/write.
    /// \param  filename the vcf file to open.
    /// \param header to be read/written from/to the file
    /// \return true = success; false = failure.
    virtual bool open(const char* filename, VcfHeader& header) = 0;
    
    /// Close the file if it is open.
    void close();

    /// When set to true, read only the first 8 columns, skipping the format
    /// and genotype fields, so when reading do not store them, and when
    /// writing do not write them.  Defaults to read/write all columns.
    /// This setting is maintained even when the file is reset/closed.
    /// \param siteOnly process only the first 8 columns
    void setSiteOnly(bool siteOnly) {mySiteOnly = siteOnly;}

    /// Get the number of VCF records that have been processed (read/written)
    /// so far including any filtered records.
    int getNumRecords() {return(myNumRecords);}

    // Get the Status of the last call that sets status.
    //    inline StatGenStatus::Status getStatus()
    //    {
    //        return(myStatus.getStatus());
    //    }

    /// Get the filename that is currently opened.
    /// \return filename associated with this class
    const char* getFileName() const
    {
        return(myFilePtr->getFileName());
    }


protected: 
    // Open the vcf file with the specified filename
    // with the specified mode.
    // \param  filename the vcf file to open.
    // \param  mode how to open (r/w).
    // \return true = success; false = failure.
    bool open(const char* filename, const char* mode,
              InputFile::ifileCompression compressionMode = InputFile::DEFAULT);
    
    void reset();
    virtual void resetFile() = 0;

    IFILE  myFilePtr;

    StatGenStatus myStatus;

    bool mySiteOnly;

    // Number of records read/written so far.  Child classes need to set this.
    int myNumRecords;

private:
    VcfFile(const VcfFile& vcfFile);
    VcfFile& operator=(const VcfFile& vcfFile);
};

#endif
