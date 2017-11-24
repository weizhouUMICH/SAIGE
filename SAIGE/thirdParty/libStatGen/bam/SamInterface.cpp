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

#include "SamInterface.h"
#include "SamRecordHelper.h"

#include <limits>
#include <stdint.h>

SamInterface::SamInterface()
    : myFirstRecord("")
{
}


SamInterface::~SamInterface()
{
}


// Read a SAM file's header.
bool SamInterface::readHeader(IFILE filePtr, SamFileHeader& header,
                              SamStatus& status)
{
    if(filePtr == NULL)
    {
        // File is not open.
        status.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot read header since the file pointer is null");
        return(false);
    }

    // Clear the passed in header.
    header.resetHeader();

    int numValid = 0;
    int numInvalid = 0;
    std::string errorMessages = "";

    do {
        StringIntHash tags;
        StringArray   values;
        buffer.ReadLine(filePtr);
      
        // Stop reading header lines if at the end of the file or
        // if the line is not blank and does not start with an @.
        if ( ifeof(filePtr) || 
             ((buffer.Length() != 0) && (buffer[0] != '@')) )
        {
            break;
        }
      
        // This is a header line, so add it to header.
        if(header.addHeaderLine(buffer.c_str()))
        {
            if(buffer.Length() != 0)
            {
                ++numValid;
            }
        }
        else
        {
            ++numInvalid;
            // Failed reading the header.
            errorMessages += header.getErrorMessage();
            // Skip further processing on this line since it was an error.
            continue;
        }
    } while (1);
   
    // Store the first record since it was read.
    myFirstRecord = buffer;

    if(numInvalid > 0)
    {
        if(numValid == 0)
        {
            std::cerr << "Failed to parse " << numInvalid << " header lines";
            std::cerr << ".  No valid header lines.\n";
            status.setStatus(SamStatus::FAIL_PARSE, errorMessages.c_str());
            return(false);
        }
    }
 
    // Successfully read.
    return(true);
}

bool SamInterface::writeHeader(IFILE filePtr, SamFileHeader& header,
                               SamStatus& status)
{
    if((filePtr == NULL) || (filePtr->isOpen() == false))
    {
        // File is not open, return failure.
        status.setStatus(SamStatus::FAIL_ORDER, 
                         "Cannot write header since the file pointer is null");
        return(false);
    }

    ////////////////////////////////
    // Write the header to the file.
    ////////////////////////////////
    // Construct a string containing the entire header.
    std::string headerString = "";
    header.getHeaderString(headerString);
    
    int32_t headerLen = headerString.length();
    int numWrite = 0;
    
    // Write the header to the file.
    numWrite = ifwrite(filePtr, headerString.c_str(), headerLen);
    if(numWrite != headerLen)
    {
        status.setStatus(SamStatus::FAIL_IO, 
                         "Failed to write the SAM header.");
        return(false);
    }
    return(true);
}


void SamInterface::readRecord(IFILE filePtr, SamFileHeader& header,
                              SamRecord& record, 
                              SamStatus& samStatus)
{
    // Initialize the status to success - will be set to false on failure.
    samStatus = SamStatus::SUCCESS;

    if((filePtr == NULL) || (filePtr->isOpen() == false))
    {
        // File is not open.
        samStatus.addError(SamStatus::FAIL_ORDER, 
                           "filePtr does not point to an open file.");
        return;
    }
    
    // If the first record has been set, use that and clear it,
    // otherwise read the record from the file.
    if(myFirstRecord.Length() != 0)
    {
        buffer = myFirstRecord;
        myFirstRecord.Clear();
    }
    else
    {
        // Read the next record.
        buffer.Clear();
        buffer.ReadLine(filePtr);
        // If the end of the file and nothing was read, return false.
        if ((ifeof(filePtr)) && (buffer.Length() == 0))
        {
            // end of the file and nothing to process.
            samStatus.addError(SamStatus::NO_MORE_RECS, 
                               "No more records in the file.");
            return;
        }
    }
    
    tokens.ReplaceColumns(buffer, '\t');
    
    
    // Error string for reporting a parsing failure.
    String errorString = "";
    
    if (tokens.Length() < 11)
    {
        errorString = "Too few columns (";
        errorString += tokens.Length();
        errorString += ") in the Record, expected at least 11.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
        return;
    }
        
    // Reset the record before setting any fields.
    record.resetRecord();

    if(!record.setReadName(tokens[0]))
    {
        samStatus.addError(record.getStatus());
    }
    
    long flagInt = 0;
    if(!tokens[1].AsInteger(flagInt))
    {
        errorString = "flag, ";
        errorString += tokens[1].c_str();
        errorString += ", is not an integer.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
    }
    else if((flagInt < 0) || (flagInt > UINT16_MAX))
    {
        errorString = "flag, ";
        errorString += tokens[1].c_str();
        errorString += ", is not between 0 and (2^16)-1 = 65535.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
    }
    else if(!record.setFlag(flagInt))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    if(!record.setReferenceName(header, tokens[2]))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    long posInt = 0;
    if(!tokens[3].AsInteger(posInt))
    {
        errorString = "position, ";
        errorString += tokens[3].c_str();
        errorString += ", is not an integer.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
    }
    else if((posInt < INT32_MIN) || (posInt > INT32_MAX))
    {
        // If it is not in this range, it cannot fit into a 32 bit int.
        errorString = "position, ";
        errorString += tokens[3].c_str();
        errorString += ", does not fit in a 32 bit signed int.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
    }
    else if(!record.set1BasedPosition(posInt))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    long mapInt = 0;
    if(!tokens[4].AsInteger(mapInt))
    {
        errorString = "map quality, ";
        errorString += tokens[4].c_str();
        errorString += ", is not an integer.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
    }
    else if((mapInt < 0) || (mapInt > UINT8_MAX))
    {
        errorString = "map quality, ";
        errorString += tokens[4].c_str();
        errorString += ", is not between 0 and (2^8)-1 = 255.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
    }
    else if(!record.setMapQuality(mapInt))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    if(!record.setCigar(tokens[5]))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    if(!record.setMateReferenceName(header, tokens[6]))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    long matePosInt = 0;
    if(!tokens[7].AsInteger(matePosInt))
    {
        errorString = "mate position, ";
        errorString += tokens[7].c_str();
        errorString += ", is not an integer.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
    }
    else if(!record.set1BasedMatePosition(matePosInt))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    long insertInt = 0;
    if(!tokens[8].AsInteger(insertInt))
    {
        errorString = "insert size, ";
        errorString += tokens[8].c_str();
        errorString += ", is not an integer.";
        samStatus.addError(SamStatus::FAIL_PARSE,
                           errorString.c_str());
    }
    else if(!record.setInsertSize(insertInt))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    if(!record.setSequence(tokens[9]))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }

    if(!record.setQuality(tokens[10]))
    {
        samStatus.addError(record.getStatus().getStatus(),
                           record.getStatus().getStatusMessage());
    }
    
    // Clear the tag fields.
    record.clearTags();
    
    // Add the tags to the record.
    for (int i = 11; i < tokens.Length(); i++)
    {
        String & nugget = tokens[i];
        
        if (nugget.Length() < 6 || nugget[2] != ':' || nugget[4] != ':')
        {
            // invalid tag format.
            errorString = "Invalid Tag Format: ";
            errorString += nugget.c_str();
            errorString += ", should be cc:c:x*.";
            samStatus.addError(SamStatus::FAIL_PARSE,
                               errorString.c_str());
            continue;
        }
        
        // Valid tag format.
        // Add the tag.
        if(!record.addTag((const char *)nugget, nugget[3],
                          (const char *)nugget + 5))
        {
            samStatus.addError(record.getStatus().getStatus(),
                               record.getStatus().getStatusMessage());
        }
    }

    return;
}


SamStatus::Status SamInterface::writeRecord(IFILE filePtr,
                                            SamFileHeader& header, 
                                            SamRecord& record,
                                            SamRecord::SequenceTranslation translation)
{
    // Store all the fields into a string, then write the string.
    String recordString = record.getReadName();
    recordString += "\t";
    recordString += record.getFlag();
    recordString += "\t";
    recordString += record.getReferenceName();
    recordString += "\t";
    recordString += record.get1BasedPosition();
    recordString += "\t";
    recordString += record.getMapQuality();
    recordString += "\t";
    recordString += record.getCigar();
    recordString += "\t";
    recordString += record.getMateReferenceNameOrEqual();
    recordString += "\t";
    recordString += record.get1BasedMatePosition();
    recordString += "\t";
    recordString += record.getInsertSize();
    recordString += "\t";
    recordString += record.getSequence(translation);
    recordString += "\t";
    recordString += record.getQuality();

    // If there are any tags, add a preceding tab.
    if(record.getTagLength() != 0)
    {
        recordString += "\t";
        SamRecordHelper::genSamTagsString(record, recordString);
    }

    recordString += "\n";
   
   
    // Write the record.
    ifwrite(filePtr, recordString.c_str(), recordString.Length());
    return(SamStatus::SUCCESS);
}


void SamInterface::ParseHeaderLine(StringIntHash & tags, StringArray & values)
{
    tags.Clear();
    values.Clear();

    tokens.AddColumns(buffer, '\t');

    for (int i = 1; i < tokens.Length(); i++)
    {
        tags.Add(tokens[i].Left(2), i - 1);
        values.Push(tokens[i].SubStr(3));
    }
}


bool SamInterface::isEOF(IFILE filePtr)
{

    if(myFirstRecord.Length() != 0)
    {
        // First record is set, so return false, not EOF, since we 
        // know the record still needs to be processed.
        return(false);
    }
    return(GenericSamInterface::isEOF(filePtr));
}
