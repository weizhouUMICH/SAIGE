/*
 *  Copyright (C) 2010-2012  Regents of the University of Michigan,
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

#include "VcfFileWriter.h"

VcfFileWriter::VcfFileWriter()
    : VcfFile()
{
}


VcfFileWriter::~VcfFileWriter() 
{
}


bool VcfFileWriter::open(const char* filename, VcfHeader& header,
                         InputFile::ifileCompression compressionMode)
{
    myStatus = StatGenStatus::SUCCESS;
    if(VcfFile::open(filename, "w", compressionMode))
    {
        // Successfully opened, so write the header.
        if(!header.write(myFilePtr))
        {
            // Failed, so copy the status.
            myStatus = header.getStatus();
            return(false);
        }
    }
    else
    {
        // Failed, status set by VcfFile::open.
        return(false);
    }

    // Successfully opened and read the header.
    return(true);
}


bool VcfFileWriter::open(const char* filename, VcfHeader& header)
{
    return(open(filename, header, InputFile::BGZF));
}


bool VcfFileWriter::writeRecord(VcfRecord& record)
{
    if(!record.write(myFilePtr, mySiteOnly))
    {
        myStatus = record.getStatus();
        return(false);
    }
    ++myNumRecords;
    return(true);
}
