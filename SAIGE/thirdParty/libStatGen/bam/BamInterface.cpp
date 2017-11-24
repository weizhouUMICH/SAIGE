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

#include "BamInterface.h"
#include "CharBuffer.h"

BamInterface::BamInterface()
{
}


BamInterface::~BamInterface()
{
}


// Read a BAM file's header.
bool BamInterface::readHeader(IFILE filePtr, SamFileHeader& header,
                              SamStatus& status)
{
    if(filePtr == NULL)
    {
        // File is not open, return false.
        status.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot read header since the file pointer is null");
        return(false);
    }
    if(filePtr->isOpen() == false)
    {
        status.setStatus(SamStatus::FAIL_ORDER, 
                         "Cannot read header since the file is not open");
        return(false);
    }

    // Clear the passed in header.
    header.resetHeader();

    int32_t headerLength;
    int readSize = ifread(filePtr, &headerLength, sizeof(headerLength));
    
    if(readSize != sizeof(headerLength))
    {
        String errMsg = "Failed to read the BAM header length, read ";
        errMsg += readSize;
        errMsg += " bytes instead of ";
        errMsg += (unsigned int)sizeof(headerLength);
        status.setStatus(SamStatus::FAIL_IO, errMsg.c_str());
        return(false);
    }

    String headerStr;
    if(headerLength > 0)
    {
        // Read the header.
        readSize = 
            ifread(filePtr, headerStr.LockBuffer(headerLength + 1), headerLength);
        headerStr[headerLength] = 0;
        headerStr.UnlockBuffer();
        if(readSize != headerLength)
        {
            // Failed to read the header.
            status.setStatus(SamStatus::FAIL_IO, "Failed to read the BAM header.");
            return(false);
        }
    }
    
    // Parse the header that was read.
    if(!header.addHeader(headerStr))
    {
        // Status is set in the method on failure.
        status.setStatus(SamStatus::FAIL_PARSE, header.getErrorMessage());
        return(false);
    }

    int referenceCount;
    // Read the number of references sequences.
    ifread(filePtr, &referenceCount, sizeof(int));

    // Get and clear the reference info so it can be set
    // from the bam reference table.
    SamReferenceInfo& refInfo = 
        header.getReferenceInfoForBamInterface();
    refInfo.clear();

    CharBuffer refName;
    
    // Read each reference sequence
    for (int i = 0; i < referenceCount; i++)
    {
        int nameLength;
        int rc;
        // Read the length of the reference name.
        rc = ifread(filePtr, &nameLength, sizeof(int));
        if(rc != sizeof(int))
        {
            status.setStatus(SamStatus::FAIL_IO, 
                             "Failed to read the BAM reference dictionary.");
            return(false);
        }
      
        // Read the name.
        refName.readFromFile(filePtr, nameLength);

        // Read the length of the reference sequence.
        int32_t refLen;
        rc = ifread(filePtr, &refLen, sizeof(int));

        if(rc != sizeof(int)) {
            status.setStatus(SamStatus::FAIL_IO, 
                             "Failed to read the BAM reference dictionary.");
            return(false);
        }

        refInfo.add(refName.c_str(), refLen);
    }

    // Successfully read the file.
    return(true);
}


bool BamInterface::writeHeader(IFILE filePtr, SamFileHeader& header,
                               SamStatus& status)
{
    if((filePtr == NULL) || (filePtr->isOpen() == false))
    {
        // File is not open, return false.
        status.setStatus(SamStatus::FAIL_ORDER, 
                         "Cannot write header since the file pointer is null");
        return(false);
    }

    char magic[4];
    magic[0] = 'B';
    magic[1] = 'A';
    magic[2] = 'M';
    magic[3] = 1;

    // Write magic to the file.
    ifwrite(filePtr, magic, 4);

    ////////////////////////////////
    // Write the header to the file.
    ////////////////////////////////
    // Construct a string containing the entire header.
    std::string headerString = "";
    header.getHeaderString(headerString);

    int32_t headerLen = headerString.length();
    int numWrite = 0;
    
    // Write the header length.
    numWrite = ifwrite(filePtr, &headerLen, sizeof(int32_t));
    if(numWrite != sizeof(int32_t))
    {
        status.setStatus(SamStatus::FAIL_IO, 
                         "Failed to write the BAM header length.");
        return(false);
    }
   
    // Write the header to the file.
    numWrite = ifwrite(filePtr, headerString.c_str(), headerLen);
    if(numWrite != headerLen)
    {
        status.setStatus(SamStatus::FAIL_IO, 
                         "Failed to write the BAM header.");
        return(false);
    }
    
    ////////////////////////////////////////////////////////
    // Write the Reference Information.
    const SamReferenceInfo& refInfo = header.getReferenceInfo();

    // Get the number of sequences.    
    int32_t numSeq = refInfo.getNumEntries();
    ifwrite(filePtr, &numSeq, sizeof(int32_t));

    // Write each reference sequence
    for (int i = 0; i < numSeq; i++)
    {
        const char* refName = refInfo.getReferenceName(i);
        // Add one for the null value.
        int32_t nameLength = strlen(refName) + 1;
        // Write the length of the reference name.
        ifwrite(filePtr, &nameLength, sizeof(int32_t));
      
        // Write the name.
        ifwrite(filePtr, refName, nameLength);
        // Write the length of the reference sequence.
        int32_t refLen = refInfo.getReferenceLength(i);
        ifwrite(filePtr, &refLen, sizeof(int32_t));
    }

    return(true);
}


void BamInterface::readRecord(IFILE filePtr, SamFileHeader& header,
                              SamRecord& record, 
                              SamStatus& samStatus)
{
    // TODO - need to validate there are @SQ lines in both sam/bam - MAYBE!

    // SetBufferFromFile will reset the record prior to reading a new one.
    if(record.setBufferFromFile(filePtr, header) != SamStatus::SUCCESS)
    {
        // Failed, so add the error message.
        samStatus.addError(record.getStatus());
    }
}

SamStatus::Status BamInterface::writeRecord(IFILE filePtr, 
                                            SamFileHeader& header,
                                            SamRecord& record,
                                            SamRecord::SequenceTranslation translation)
{
    // Write the file, returning the status.
    return(record.writeRecordBuffer(filePtr, translation));
}


