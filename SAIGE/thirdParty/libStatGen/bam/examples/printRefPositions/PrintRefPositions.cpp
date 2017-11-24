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

//////////////////////////////////////////////////////////////////////////
#include "SamFile.h"

void printRefPositions(std::string inFile, std::string indexFile,
                       std::string rname, int startPosition, 
                       int endPosition)
{
    SamFileHeader header;
    // Open the bam file for reading and read the header.
    SamFile samIn(inFile.c_str(), SamFile::READ, &header);

    // Open the bam index file for reading.
    samIn.ReadBamIndex(indexFile.c_str());

    // Set the section to be read.
    samIn.SetReadSection(rname.c_str(), startPosition, endPosition);

    SamRecord record;
    // Keep reading BAM records until they aren't anymore.
    while(samIn.ReadRecord(header, record))
    {
        // Print the reference positions associated with this read.
        std::cout << "Read " << samIn.GetCurrentRecordCount() << ":";
        Cigar* cigar = record.getCigarInfo();
        for(int i = 0; i < record.getReadLength(); i++)
        {
            int refPos = 
                cigar->getRefPosition(i, record.get1BasedPosition());
            if(refPos != Cigar::INDEX_NA)
            {
                std::cout << "  " << refPos;
            }
        }
        std::cout << "\n";
    }
}
