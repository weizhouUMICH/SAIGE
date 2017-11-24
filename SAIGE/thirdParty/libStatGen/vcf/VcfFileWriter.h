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

#ifndef __VCF_FILE_WRITER_H__
#define __VCF_FILE_WRITER_H__

#include "VcfFile.h"
#include "VcfRecord.h"

/// This header file provides interface to read/write VCF files.
class VcfFileWriter : public VcfFile
{
public:
    /// Default Constructor, initializes the variables, but does not open
    /// any files.
    VcfFileWriter();
    /// Destructor
    virtual ~VcfFileWriter();

    /// Open the vcf file with the specified filename for writing.
    /// \param filename the vcf file to open for writing.
    /// \param header to be written the file
    /// \param compressionMode type of compression to use for writing
    /// \return true = success; false = failure.
    bool open(const char* filename, VcfHeader& header,
              InputFile::ifileCompression compressionMode);

    /// Open the vcf file with the specified filename for writing using the
    /// default compression (BGZF).
    /// \param filename the vcf file to open for writing.
    /// \param header to be written the file
    /// \return true = success; false = failure.
    virtual bool open(const char* filename, VcfHeader& header);
    
    /// Write the VCF data line to the file.
    /// \param record record to write to the file.
    /// \return true if successfully wrote, false if not.
    bool writeRecord(VcfRecord& record);

protected: 
    virtual void resetFile() {}

private:
    VcfFileWriter(const VcfFileWriter& vcfFileWriter);
    VcfFileWriter& operator=(const VcfFileWriter& vcfFileWriter);
};

#endif
