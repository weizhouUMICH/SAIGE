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
#include <stdexcept>
#include "SamFile.h"
#include "SamFileHeader.h"
#include "SamRecord.h"
#include "BamInterface.h"
#include "SamInterface.h"
#include "BgzfFileType.h"

// Constructor, init variables.
SamFile::SamFile()
    : myStatus()
{
    init();
    resetFile();
}


// Constructor, init variables.
SamFile::SamFile(ErrorHandler::HandlingType errorHandlingType)
    : myStatus(errorHandlingType)
{
    init();
    resetFile();
}


// Constructor, init variables and open the specified file based on the
// specified mode (READ/WRITE).
SamFile::SamFile(const char* filename, OpenType mode)
    : myStatus()
{
    init(filename, mode, NULL);
}


// Constructor, init variables and open the specified file based on the
// specified mode (READ/WRITE).  Default is READ..
SamFile::SamFile(const char* filename, OpenType mode,
                 ErrorHandler::HandlingType errorHandlingType)
    : myStatus(errorHandlingType)
{
    init(filename, mode, NULL);
}


// Constructor, init variables and open the specified file based on the
// specified mode (READ/WRITE).
SamFile::SamFile(const char* filename, OpenType mode, SamFileHeader* header)
    : myStatus()
{
    init(filename, mode, header);
}


// Constructor, init variables and open the specified file based on the
// specified mode (READ/WRITE).  Default is READ..
SamFile::SamFile(const char* filename, OpenType mode,
                 ErrorHandler::HandlingType errorHandlingType, 
                 SamFileHeader* header)
    : myStatus(errorHandlingType)
{
    init(filename, mode, header);
}


SamFile::~SamFile()
{
    resetFile();
    if(myStatistics != NULL)
    {
        delete myStatistics;
    }
}


// Open a sam/bam file for reading with the specified filename.
bool SamFile::OpenForRead(const char * filename, SamFileHeader* header)
{
    // Reset for any previously operated on files.
    resetFile();

    int lastchar = 0;

    while (filename[lastchar] != 0) lastchar++;

    // If at least one character, check for '-'.
    if((lastchar >= 1) && (filename[0] == '-'))
    {
        // Read from stdin - determine type of file to read.
        // Determine if compressed bam.
        if(strcmp(filename, "-.bam") == 0)
        {
            // Compressed bam - open as bgzf.
            // -.bam is the filename, read compressed bam from stdin
            filename = "-";

            myFilePtr = new InputFile;
            // support recover mode - this switches in a reader
            // capable of recovering from bad BGZF compression blocks.
            myFilePtr->setAttemptRecovery(myAttemptRecovery);
            myFilePtr->openFile(filename, "rb", InputFile::BGZF);

            myInterfacePtr = new BamInterface;

            // Read the magic string.
            char magic[4];
            ifread(myFilePtr, magic, 4);
        }
        else if(strcmp(filename, "-.ubam") == 0)
        {
            // uncompressed BAM File.
            // -.ubam is the filename, read uncompressed bam from stdin.
            // uncompressed BAM is still compressed with BGZF, but using
            // compression level 0, so still open as BGZF since it has a
            // BGZF header.
            filename = "-";

            // Uncompressed, so do not require the eof block.
#ifdef __ZLIB_AVAILABLE__
            BgzfFileType::setRequireEofBlock(false);
#endif
            myFilePtr = ifopen(filename, "rb", InputFile::BGZF);
        
            myInterfacePtr = new BamInterface;

            // Read the magic string.
            char magic[4];
            ifread(myFilePtr, magic, 4);
        }
        else if((strcmp(filename, "-") == 0) || (strcmp(filename, "-.sam") == 0))
        {
            // SAM File.
            // read sam from stdin
            filename = "-";
            myFilePtr = ifopen(filename, "rb", InputFile::UNCOMPRESSED);
            myInterfacePtr = new SamInterface;
        }
        else
        {
            std::string errorMessage = "Invalid SAM/BAM filename, ";
            errorMessage += filename;
            errorMessage += ".  From stdin, can only be '-', '-.sam', '-.bam', or '-.ubam'";
            myStatus.setStatus(SamStatus::FAIL_IO, errorMessage.c_str());
            delete myFilePtr;
            myFilePtr = NULL;
            return(false);          
        }
    }
    else
    {
        // Not from stdin.  Read the file to determine the type.

        myFilePtr = new InputFile;

        // support recovery mode - this conditionally enables a reader
        // capable of recovering from bad BGZF compression blocks.
        myFilePtr->setAttemptRecovery(myAttemptRecovery);
        bool rc = myFilePtr->openFile(filename, "rb", InputFile::DEFAULT);

        if (rc == false)
        {
            std::string errorMessage = "Failed to Open ";
            errorMessage += filename;
            errorMessage += " for reading";
            myStatus.setStatus(SamStatus::FAIL_IO, errorMessage.c_str());
            delete myFilePtr;
            myFilePtr = NULL;
            return(false);
        }
        
        char magic[4];
        ifread(myFilePtr, magic, 4);
        
        if (magic[0] == 'B' && magic[1] == 'A' && magic[2] == 'M' &&
            magic[3] == 1)
        {
            myInterfacePtr = new BamInterface;
            // Set that it is a bam file open for reading.  This is needed to
            // determine if an index file can be used.
            myIsBamOpenForRead = true;
        }
        else
        {
            // Not a bam, so rewind to the beginning of the file so it
            // can be read.
            ifrewind(myFilePtr);
            myInterfacePtr = new SamInterface;
        }
    }

    // File is open for reading.
    myIsOpenForRead = true;

    // Read the header if one was passed in.
    if(header != NULL)
    {
        return(ReadHeader(*header));
    }

    // Successfully opened the file.
    myStatus = SamStatus::SUCCESS;
    return(true);
}


// Open a sam/bam file for writing with the specified filename.
bool SamFile::OpenForWrite(const char * filename, SamFileHeader* header)
{
    // Reset for any previously operated on files.
    resetFile();
    
    int lastchar = 0;
    while (filename[lastchar] != 0) lastchar++;   
    if (lastchar >= 4 && 
        filename[lastchar - 4] == 'u' &&
        filename[lastchar - 3] == 'b' &&
        filename[lastchar - 2] == 'a' &&
        filename[lastchar - 1] == 'm')
    {
        // BAM File.
        // if -.ubam is the filename, write uncompressed bam to stdout
        if((lastchar == 6) && (filename[0] == '-') && (filename[1] == '.'))
        {
            filename = "-";
        }

        myFilePtr = ifopen(filename, "wb0", InputFile::BGZF);

        myInterfacePtr = new BamInterface;
    }
    else if (lastchar >= 3 && 
             filename[lastchar - 3] == 'b' &&
             filename[lastchar - 2] == 'a' &&
             filename[lastchar - 1] == 'm')
    {
        // BAM File.
        // if -.bam is the filename, write compressed bam to stdout
        if((lastchar == 5) && (filename[0] == '-') && (filename[1] == '.'))
        {
            filename = "-";
        }
        myFilePtr = ifopen(filename, "wb", InputFile::BGZF);
        
        myInterfacePtr = new BamInterface;
    }
    else
    {
        // SAM File
        // if - (followed by anything is the filename,
        // write uncompressed sam to stdout
        if((lastchar >= 1) && (filename[0] == '-'))
        {
            filename = "-";
        }
        myFilePtr = ifopen(filename, "wb", InputFile::UNCOMPRESSED);
   
        myInterfacePtr = new SamInterface;
    }

    if (myFilePtr == NULL)
    {
        std::string errorMessage = "Failed to Open ";
        errorMessage += filename;
        errorMessage += " for writing";
        myStatus.setStatus(SamStatus::FAIL_IO, errorMessage.c_str());
        return(false);
    }
   
    myIsOpenForWrite = true;

    // Write the header if one was passed in.
    if(header != NULL)
    {
        return(WriteHeader(*header));
    }

    // Successfully opened the file.
    myStatus = SamStatus::SUCCESS;
    return(true);
}


// Read BAM Index file.
bool SamFile::ReadBamIndex(const char* bamIndexFilename)
{
    // Cleanup a previously setup index.
    if(myBamIndex != NULL)
    {
        delete myBamIndex;
        myBamIndex = NULL;
    }

    // Create a new bam index.
    myBamIndex = new BamIndex();
    SamStatus::Status indexStat = myBamIndex->readIndex(bamIndexFilename);

    if(indexStat != SamStatus::SUCCESS)
    {
        std::string errorMessage = "Failed to read the bam Index file: ";
        errorMessage += bamIndexFilename;
        myStatus.setStatus(indexStat, errorMessage.c_str());
        delete myBamIndex;
        myBamIndex = NULL;
        return(false);
    }
    myStatus = SamStatus::SUCCESS;
    return(true);
}


// Read BAM Index file.
bool SamFile::ReadBamIndex()
{
    if(myFilePtr == NULL)
    {
        // Can't read the bam index file because the BAM file has not yet been
        // opened, so we don't know the base filename for the index file.
        std::string errorMessage = "Failed to read the bam Index file -"
            " the BAM file needs to be read first in order to determine"
            " the index filename.";
        myStatus.setStatus(SamStatus::FAIL_ORDER, errorMessage.c_str());
        return(false);
    }

    const char* bamBaseName = myFilePtr->getFileName();
    
    std::string indexName = bamBaseName;
    indexName += ".bai";

    bool foundFile = true;
    try
    {
        if(ReadBamIndex(indexName.c_str()) == false)
        {
            foundFile = false;
        }
    }
    catch (std::exception&)
    {
        foundFile = false;
    }

    // Check to see if the index file was found.
    if(!foundFile)
    {
        // Not found - try without the bam extension.
        // Locate the start of the bam extension
        size_t startExt = indexName.find(".bam");
        if(startExt == std::string::npos)
        {
            // Could not find the .bam extension, so just return false since the
            // call to ReadBamIndex set the status.
            return(false);
        }
        // Remove ".bam" and try reading the index again.
        indexName.erase(startExt,  4);
        return(ReadBamIndex(indexName.c_str()));
    }
    return(true);
}


// Sets the reference to the specified genome sequence object.
void SamFile::SetReference(GenomeSequence* reference)
{
    myRefPtr = reference;
}


// Set the type of sequence translation to use when reading the sequence.
void SamFile::SetReadSequenceTranslation(SamRecord::SequenceTranslation translation)
{
    myReadTranslation = translation;
}


// Set the type of sequence translation to use when writing the sequence.
void SamFile::SetWriteSequenceTranslation(SamRecord::SequenceTranslation translation)
{
    myWriteTranslation = translation;
}

// Close the file if there is one open.
void SamFile::Close()
{
    // Resetting the file will close it if it is open, and
    // will reset all other variables.
    resetFile();
}


// Returns whether or not the file has been opened.
// return: int - true = open; false = not open.
bool SamFile::IsOpen()
{
    if (myFilePtr != NULL)
    {
        // File Pointer is set, so return if it is open.
        return(myFilePtr->isOpen());
    }
    // File pointer is not set, so return false, not open.
    return false;
}


// Returns whether or not the end of the file has been reached.
// return: int - true = EOF; false = not eof.
bool SamFile::IsEOF()
{
    if(myIsOpenForRead == false)
    {
        // Not open for read, return true.
        return(true);
    }
    return(myInterfacePtr->isEOF(myFilePtr));
}


// Returns whether or not the file is a stream.
// return: bool - true = stream; false = not stream/not open.
bool SamFile::IsStream()
{
    if (myFilePtr != NULL)
    {
        // File Pointer is set, so return if it is a stream.
        return((myFilePtr->getFileName())[0] == '-');
    }
    // File pointer is not set, so return false, not a stream.
    return false;
}


// Read the header from the currently opened file.
bool SamFile::ReadHeader(SamFileHeader& header)
{
    myStatus = SamStatus::SUCCESS;
    if(myIsOpenForRead == false)
    {
        // File is not open for read
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot read header since the file is not open for reading");
        return(false);
    }

    if(myHasHeader == true)
    {
        // The header has already been read.
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot read header since it has already been read.");
        return(false);
    }

    if(myInterfacePtr->readHeader(myFilePtr, header, myStatus))
    {
        // The header has now been successfully read.
        myHasHeader = true;
        return(true);
    }
    return(false);
}


// Write the header to the currently opened file.
bool SamFile::WriteHeader(SamFileHeader& header)
{
    myStatus = SamStatus::SUCCESS;
    if(myIsOpenForWrite == false)
    {
        // File is not open for write
        // -OR-
        // The header has already been written.
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot write header since the file is not open for writing");
        return(false);
    }

    if(myHasHeader == true)
    {
        // The header has already been written.
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot write header since it has already been written");
        return(false);
    }

    if(myInterfacePtr->writeHeader(myFilePtr, header, myStatus))
    {
        // The header has now been successfully written.
        myHasHeader = true;
        return(true);
    }

    // return the status.
    return(false);
}


// Read a record from the currently opened file.
bool SamFile::ReadRecord(SamFileHeader& header, 
                         SamRecord& record)
{
    myStatus = SamStatus::SUCCESS;

    if(myIsOpenForRead == false)
    {
        // File is not open for read
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot read record since the file is not open for reading");
        throw(std::runtime_error("SOFTWARE BUG: trying to read a SAM/BAM record prior to opening the file."));
        return(false);
    }

    if(myHasHeader == false)
    {
        // The header has not yet been read.
        // TODO - maybe just read the header.
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot read record since the header has not been read.");
        throw(std::runtime_error("SOFTWARE BUG: trying to read a SAM/BAM record prior to reading the header."));
        return(false);
    }

    // Check to see if a new region has been set.  If so, determine the
    // chunks for that region.
    if(myNewSection)
    {
        if(!processNewSection(header))
        {
            // Failed processing a new section.  Could be an 
            // order issue like the file not being open or the
            // indexed file not having been read.
            // processNewSection sets myStatus with the failure reason.
            return(false);
        }
    }

    // Read until a record is not successfully read or there are no more
    // requested records.
    while(myStatus == SamStatus::SUCCESS)
    {
        record.setReference(myRefPtr);
        record.setSequenceTranslation(myReadTranslation);

        // If reading by index, this method will setup to ensure it is in
        // the correct position for the next record (if not already there).
        // Sets myStatus if it could not move to a good section.
        // Just returns true if it is not setup to read by index.
        if(!ensureIndexedReadPosition())
        {
            // Either there are no more records in the section
            // or it failed to move to the right section, so there
            // is nothing more to read, stop looping.
            break;
        }
        
        // File is open for reading and the header has been read, so read the
        // next record.
        myInterfacePtr->readRecord(myFilePtr, header, record, myStatus);
        if(myStatus != SamStatus::SUCCESS)
        {
            // Failed to read the record, so break out of the loop.
            break;
        }

        // Successfully read a record, so check if we should filter it.
        // First check if it is out of the section.  Returns true
        // if not reading by sections, returns false if the record
        // is outside of the section.  Sets status to NO_MORE_RECS if
        // there is nothing left ot read in the section.
        if(!checkRecordInSection(record))
        {
            // The record is not in the section.
            // The while loop will detect if NO_MORE_RECS was set.
            continue;
        }

        // Check the flag for required/excluded flags.
        uint16_t flag = record.getFlag();
        if((flag & myRequiredFlags) != myRequiredFlags)
        {
            // The record does not conatain all required flags, so
            // continue looking.
            continue;
        }
        if((flag & myExcludedFlags) != 0)
        {
            // The record contains an excluded flag, so continue looking.
            continue;
        }

        //increment the record count.
        myRecordCount++;
        
        if(myStatistics != NULL)
        {
            // Statistics should be updated.
            myStatistics->updateStatistics(record);
        }
        
        // Successfully read the record, so check the sort order.
        if(!validateSortOrder(record, header))
        {
            // ValidateSortOrder sets the status on a failure.
            return(false);
        }
        return(true);

    } // End while loop that checks if a desired record is found or failure.

    // Return true if a record was found.
    return(myStatus == SamStatus::SUCCESS);
}



// Write a record to the currently opened file.
bool SamFile::WriteRecord(SamFileHeader& header, 
                          SamRecord& record)
{
    if(myIsOpenForWrite == false)
    {
        // File is not open for writing
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot write record since the file is not open for writing");
        return(false);
    }

    if(myHasHeader == false)
    {
        // The header has not yet been written.
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot write record since the header has not been written");
        return(false);
    }

    // Before trying to write the record, validate the sort order.
    if(!validateSortOrder(record, header))
    {
        // Not sorted like it is supposed to be, do not write the record
        myStatus.setStatus(SamStatus::INVALID_SORT, 
                           "Cannot write the record since the file is not properly sorted.");
        return(false);
    }

    if(myRefPtr != NULL)
    {
        record.setReference(myRefPtr);
    }

    // File is open for writing and the header has been written, so write the
    // record.
    myStatus = myInterfacePtr->writeRecord(myFilePtr, header, record,
                                           myWriteTranslation);

    if(myStatus == SamStatus::SUCCESS)
    {
        // A record was successfully written, so increment the record count.
        myRecordCount++;
        return(true);
    }
    return(false);
}


// Set the flag to validate that the file is sorted as it is read/written.
// Must be called after the file has been opened.
void SamFile::setSortedValidation(SortedType sortType)
{
    mySortedType = sortType;
}


// Return the number of records that have been read/written so far.
uint32_t SamFile::GetCurrentRecordCount()
{
    return(myRecordCount);
}


// Sets what part of the SamFile should be read.
bool SamFile::SetReadSection(int32_t refID)
{
    // No start/end specified, so set back to default -1.
    return(SetReadSection(refID, -1, -1));
}



// Sets what part of the SamFile should be read.
bool SamFile::SetReadSection(const char* refName)
{
    // No start/end specified, so set back to default -1.
    return(SetReadSection(refName, -1, -1));
}


// Sets what part of the BAM file should be read.
bool SamFile::SetReadSection(int32_t refID, int32_t start, int32_t end, 
                             bool overlap)
{
    // If there is not a BAM file open for reading, return failure.
    // Opening a new file clears the read section, so it must be
    // set after the file is opened.
    if(!myIsBamOpenForRead)
    {
        // There is not a BAM file open for reading.
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot set section since there is no bam file open");
        return(false);
    }

    myNewSection = true;
    myOverlapSection = overlap;
    myStartPos = start;
    myEndPos = end;
    myRefID = refID;
    myRefName.clear();
    myChunksToRead.clear();
    // Reset the end of the current chunk.  We are resetting our read, so
    // we no longer have a "current chunk" that we are reading.
    myCurrentChunkEnd = 0;
    myStatus = SamStatus::SUCCESS;

    // Reset the sort order criteria since we moved around in the file.    
    myPrevCoord = -1;
    myPrevRefID = 0;
    myPrevReadName.Clear();

    return(true);
}


// Sets what part of the BAM file should be read.
bool SamFile::SetReadSection(const char* refName, int32_t start, int32_t end,
                             bool overlap)
{
    // If there is not a BAM file open for reading, return failure.
    // Opening a new file clears the read section, so it must be
    // set after the file is opened.
    if(!myIsBamOpenForRead)
    {
        // There is not a BAM file open for reading.
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot set section since there is no bam file open");
        return(false);
    }

    myNewSection = true;
    myOverlapSection = overlap;
    myStartPos = start;
    myEndPos = end;
    if((strcmp(refName, "") == 0) || (strcmp(refName, "*") == 0))
    {
        // No Reference name specified, so read just the "-1" entries.
        myRefID = BamIndex::REF_ID_UNMAPPED;
    }
    else
    {
        // save the reference name and revert the reference ID to unknown
        // so it will be calculated later.
        myRefName = refName;
        myRefID = BamIndex::REF_ID_ALL;
    }
    myChunksToRead.clear();
    // Reset the end of the current chunk.  We are resetting our read, so
    // we no longer have a "current chunk" that we are reading.
    myCurrentChunkEnd = 0;
    myStatus = SamStatus::SUCCESS;
    
    // Reset the sort order criteria since we moved around in the file.    
    myPrevCoord = -1;
    myPrevRefID = 0;
    myPrevReadName.Clear();

    return(true);
}


void SamFile::SetReadFlags(uint16_t requiredFlags, uint16_t excludedFlags)
{
    myRequiredFlags = requiredFlags;
    myExcludedFlags = excludedFlags;
}


// Get the number of mapped reads in the specified reference id.  
// Returns -1 for out of range refIDs.
int32_t SamFile::getNumMappedReadsFromIndex(int32_t refID)
{
    // The bam index must have already been read.
    if(myBamIndex == NULL)
    {
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot get num mapped reads from the index until it has been read.");
        return(false);
    }
    return(myBamIndex->getNumMappedReads(refID));
}


// Get the number of unmapped reads in the specified reference id.  
// Returns -1 for out of range refIDs.
int32_t SamFile::getNumUnMappedReadsFromIndex(int32_t refID)
{
    // The bam index must have already been read.
    if(myBamIndex == NULL)
    {
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot get num unmapped reads from the index until it has been read.");
        return(false);
    }
    return(myBamIndex->getNumUnMappedReads(refID));
}


// Get the number of mapped reads in the specified reference id.  
// Returns -1 for out of range references.
int32_t SamFile::getNumMappedReadsFromIndex(const char* refName,
                                            SamFileHeader& header)
{
    // The bam index must have already been read.
    if(myBamIndex == NULL)
    {
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot get num mapped reads from the index until it has been read.");
        return(false);
    }
    int32_t refID = BamIndex::REF_ID_UNMAPPED;
    if((strcmp(refName, "") != 0) && (strcmp(refName, "*") != 0))
    {
        // Reference name specified, so read just the "-1" entries.
        refID =  header.getReferenceID(refName);
    }
    return(myBamIndex->getNumMappedReads(refID));
}


// Get the number of unmapped reads in the specified reference id.  
// Returns -1 for out of range refIDs.
int32_t SamFile::getNumUnMappedReadsFromIndex(const char* refName,
                                              SamFileHeader& header)
{
    // The bam index must have already been read.
    if(myBamIndex == NULL)
    {
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot get num unmapped reads from the index until it has been read.");
        return(false);
    }
    int32_t refID = BamIndex::REF_ID_UNMAPPED;
    if((strcmp(refName, "") != 0) && (strcmp(refName, "*") != 0))
    {
        // Reference name specified, so read just the "-1" entries.
        refID =  header.getReferenceID(refName);
    }
    return(myBamIndex->getNumUnMappedReads(refID));
}


// Returns the number of bases in the passed in read that overlap the
// region that is currently set.
uint32_t SamFile::GetNumOverlaps(SamRecord& samRecord)
{
    if(myRefPtr != NULL)
    {
        samRecord.setReference(myRefPtr);
    }
    samRecord.setSequenceTranslation(myReadTranslation);

    // Get the overlaps in the sam record for the region currently set
    // for this file.
    return(samRecord.getNumOverlaps(myStartPos, myEndPos));
}


void SamFile::GenerateStatistics(bool genStats)
{
    if(genStats)
    {
        if(myStatistics == NULL)
        {
            // Want to generate statistics, but do not yet have the
            // structure for them, so create one.
            myStatistics = new SamStatistics();
        }
    }
    else
    {
        // Do not generate statistics, so if myStatistics is not NULL, 
        // delete it.
        if(myStatistics != NULL)
        {
            delete myStatistics;
            myStatistics = NULL;
        }
    }

}


const BamIndex* SamFile::GetBamIndex()
{
    return(myBamIndex);
}


// initialize.
void SamFile::init()
{
    myFilePtr = NULL;
    myInterfacePtr = NULL;
    myStatistics = NULL;
    myBamIndex = NULL;
    myRefPtr = NULL;
    myReadTranslation = SamRecord::NONE;
    myWriteTranslation = SamRecord::NONE;
    myAttemptRecovery = false;
    myRequiredFlags = 0;
    myExcludedFlags = 0;
}


void SamFile::init(const char* filename, OpenType mode, SamFileHeader* header)
{
    init();
        
    resetFile();

    bool openStatus = true;
    if(mode == READ)
    {
        // open the file for read.
        openStatus = OpenForRead(filename, header);
    }
    else
    {
        // open the file for write.
        openStatus = OpenForWrite(filename, header);
    }
    if(!openStatus)
    {
        // Failed to open the file - print error and abort.
        fprintf(stderr, "%s\n", GetStatusMessage());
        std::cerr << "FAILURE - EXITING!!!" << std::endl;
        exit(-1);
    }
}


// Reset variables for each file.
void SamFile::resetFile()
{
    // Close the file.
    if (myFilePtr != NULL)
    {
        // If we already have an open file, close it.
        ifclose(myFilePtr);
        myFilePtr = NULL;
    }
    if(myInterfacePtr != NULL)
    {
        delete myInterfacePtr;
        myInterfacePtr = NULL;
    }

    myIsOpenForRead = false;
    myIsOpenForWrite = false;
    myHasHeader = false;
    mySortedType = UNSORTED;
    myPrevReadName.Clear();
    myPrevCoord = -1;
    myPrevRefID = 0;
    myRecordCount = 0;
    myStatus = SamStatus::SUCCESS;

    // Reset indexed bam values.
    myIsBamOpenForRead = false;
    myRefID = BamIndex::REF_ID_ALL;
    myStartPos = -1;
    myEndPos = -1;
    myNewSection = false;
    myOverlapSection = true;
    myCurrentChunkEnd = 0;
    myChunksToRead.clear();
    if(myBamIndex != NULL)
    {
        delete myBamIndex;
        myBamIndex = NULL;
    }

    // If statistics are being generated, reset them.
    if(myStatistics != NULL)
    {
        myStatistics->reset();
    }

    myRefName.clear();
}


// Validate that the record is sorted compared to the previously read record
// if there is one, according to the specified sort order.
// If the sort order is UNSORTED, true is returned.
bool SamFile::validateSortOrder(SamRecord& record, SamFileHeader& header)
{
    if(myRefPtr != NULL)
    {
        record.setReference(myRefPtr);
    }
    record.setSequenceTranslation(myReadTranslation);

    bool status = false;
    if(mySortedType == UNSORTED)
    {
        // Unsorted, so nothing to validate, just return true.
        status = true;
    }
    else 
    {
        // Check to see if mySortedType is based on the header.
        if(mySortedType == FLAG)
        {
            // Determine the sorted type from what was read out of the header.
            mySortedType = getSortOrderFromHeader(header);
        }

        if(mySortedType == QUERY_NAME)
        {
            // Validate that it is sorted by query name.
            // Get the query name from the record.
            const char* readName = record.getReadName();

            // Check if it is sorted either in samtools way or picard's way.
            if((myPrevReadName.Compare(readName) > 0) &&
               (strcmp(myPrevReadName.c_str(), readName) > 0))
            {
                // return false.
                String errorMessage = "ERROR: File is not sorted by read name at record ";
                errorMessage += myRecordCount;
                errorMessage += "\n\tPrevious record was ";
                errorMessage += myPrevReadName;
                errorMessage += ", but this record is ";
                errorMessage += readName;
                myStatus.setStatus(SamStatus::INVALID_SORT, 
                                   errorMessage.c_str());
                status = false;
            }
            else
            {
                myPrevReadName = readName;
                status = true;
            }
        }
        else 
        {
            // Validate that it is sorted by COORDINATES.
            // Get the leftmost coordinate and the reference index.
            int32_t refID = record.getReferenceID();
            int32_t coord = record.get0BasedPosition();
            // The unmapped reference id is at the end of a sorted file.
            if(refID == BamIndex::REF_ID_UNMAPPED)
            {
                // A new reference ID that is for the unmapped reads
                // is always valid.
                status = true;
                myPrevRefID = refID;
                myPrevCoord = coord;
            }
            else if(myPrevRefID == BamIndex::REF_ID_UNMAPPED)
            {
                // Previous reference ID was for unmapped reads, but the
                // current one is not, so this is not sorted.
                String errorMessage = "ERROR: File is not coordinate sorted at record ";
                errorMessage += myRecordCount;
                errorMessage += "\n\tPrevious record was unmapped, but this record is ";
                errorMessage += header.getReferenceLabel(refID) + ":" + coord;
                myStatus.setStatus(SamStatus::INVALID_SORT, 
                                   errorMessage.c_str());
                status = false;
            }
            else if(refID < myPrevRefID)
            {
                // Current reference id is less than the previous one, 
                //meaning that it is not sorted.
                String errorMessage = "ERROR: File is not coordinate sorted at record ";
                errorMessage += myRecordCount;
                errorMessage += "\n\tPrevious record was ";
                errorMessage += header.getReferenceLabel(myPrevRefID) + ":" + myPrevCoord;
                errorMessage += ", but this record is ";
                errorMessage += header.getReferenceLabel(refID) + ":" + coord;
                myStatus.setStatus(SamStatus::INVALID_SORT, 
                                   errorMessage.c_str());
                status = false;
            }
            else
            {
                // The reference IDs are in the correct order.
                if(refID > myPrevRefID)
                {
                    // New reference id, so set the previous coordinate to -1
                    myPrevCoord = -1;
                }
            
                // Check the coordinates.
                if(coord < myPrevCoord)
                {
                    // New Coord is less than the previous position.
                    String errorMessage = "ERROR: File is not coordinate sorted at record ";
                    errorMessage += myRecordCount;
                    errorMessage += "\n\tPreviousRecord was ";
                    errorMessage += header.getReferenceLabel(myPrevRefID) + ":" + myPrevCoord;
                    errorMessage += ", but this record is ";
                    errorMessage += header.getReferenceLabel(refID) + ":" + coord;
                    myStatus.setStatus(SamStatus::INVALID_SORT, 
                                       errorMessage.c_str());
                    status = false;
                }
                else
                {
                    myPrevRefID = refID;
                    myPrevCoord = coord;
                    status = true;
                }
            }
        }
    }

    return(status);
}


SamFile::SortedType SamFile::getSortOrderFromHeader(SamFileHeader& header)
{
    const char* tag = header.getSortOrder();
   
    // Default to unsorted since if it is not specified in the header
    // that is the value that should be used.
    SortedType headerSortOrder = UNSORTED;
    if(strcmp(tag, "queryname") == 0)
    {
        headerSortOrder = QUERY_NAME;
    }
    else if(strcmp(tag, "coordinate") == 0)
    {
        headerSortOrder = COORDINATE;
    }
    return(headerSortOrder);
}


 bool SamFile::ensureIndexedReadPosition()
 {
     // If no sections are specified, return true.
     if(myRefID == BamIndex::REF_ID_ALL)
     {
         return(true);
     }

     // Check to see if we have more to read out of the current chunk.
     // By checking the current position in relation to the current
     // end chunk.  If the current position is >= the end of the
     // current chunk, then we must see to the next chunk.
     uint64_t currentPos = iftell(myFilePtr);
     if(currentPos >= myCurrentChunkEnd)
     {
         // If there are no more chunks to process, return failure.
         if(myChunksToRead.empty())
         {
             myStatus = SamStatus::NO_MORE_RECS;
             return(false);
         }
         
         // There are more chunks left, so get the next chunk.
         Chunk newChunk = myChunksToRead.pop();
         
         // Move to the location of the new chunk if it is not adjacent
         // to the current chunk.
         if(newChunk.chunk_beg != currentPos)
         {
             // New chunk is not adjacent, so move to it.
             if(ifseek(myFilePtr, newChunk.chunk_beg, SEEK_SET) != true)
             {
                 // seek failed, cannot retrieve next record, return failure.
                 myStatus.setStatus(SamStatus::FAIL_IO, 
                                    "Failed to seek to the next record");
                 return(false);
             }
         }
         // Seek succeeded, set the end of the new chunk.
         myCurrentChunkEnd = newChunk.chunk_end;
     }
     return(true);
 }


bool SamFile::checkRecordInSection(SamRecord& record)
{
    bool recordFound = true;
    if(myRefID == BamIndex::REF_ID_ALL)
    {
        return(true);
    }
    // Check to see if it is in the correct reference/position.
    if(record.getReferenceID() != myRefID)
    {
        // Incorrect reference ID, return no more records.
        myStatus = SamStatus::NO_MORE_RECS;
        return(false);
    }
   
    // Found a record.
    recordFound = true;

    // If start/end position are set, verify that the alignment falls
    // within those.
    // If the alignment start is greater than the end of the region,
    // return NO_MORE_RECS.
    // Since myEndPos is Exclusive 0-based, anything >= myEndPos is outside
    // of the region.
    if((myEndPos != -1) && (record.get0BasedPosition() >= myEndPos))
    {
        myStatus = SamStatus::NO_MORE_RECS;
        return(false);
    }
        
    // We know the start is less than the end position, so the alignment
    // overlaps the region if the alignment end position is greater than the
    // start of the region.
    if((myStartPos != -1) && (record.get0BasedAlignmentEnd() < myStartPos))
    {
        // If it does not overlap the region, so go to the next
        // record...set recordFound back to false.
        recordFound = false;
    }

    if(!myOverlapSection)
    {
        // Needs to be fully contained.  Not fully contained if
        // 1) the record start position is < the region start position.
        // or
        // 2) the end position is specified and the record end position
        //    is greater than or equal to the region end position.
        //    (equal to since the region is exclusive.
        if((record.get0BasedPosition() < myStartPos) ||
           ((myEndPos != -1) && 
            (record.get0BasedAlignmentEnd() >= myEndPos)))
        {
            // This record is not fully contained, so move on to the next
            // record.
            recordFound = false;
        }
    }

    return(recordFound);
}
   

bool SamFile::processNewSection(SamFileHeader &header)
{
    myNewSection = false;

    // If there is no index file, return failure.
    if(myBamIndex == NULL)
    {
        // No bam index has been read.
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot read section since there is no index file open");
        throw(std::runtime_error("SOFTWARE BUG: trying to read a BAM record by section prior to opening the BAM Index file."));
        return(false);
    }

    // If there is not a BAM file open for reading, return failure.
    if(!myIsBamOpenForRead)
    {
        // There is not a BAM file open for reading.
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot read section since there is no bam file open");
        throw(std::runtime_error("SOFTWARE BUG: trying to read a BAM record by section without opening a BAM file."));
        return(false);
    }

    if(myHasHeader == false)
    {
        // The header has not yet been read.
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Cannot read record since the header has not been read.");
        throw(std::runtime_error("SOFTWARE BUG: trying to read a BAM record by section prior to opening the header."));
        return(false);
    }

    // Indexed Bam open for read, so disable read buffering because iftell
    // will be used.
    // Needs to be done here after we already know that the header has been
    // read.
    myFilePtr->disableBuffering();

    myChunksToRead.clear();
    // Reset the end of the current chunk.  We are resetting our read, so
    // we no longer have a "current chunk" that we are reading.
    myCurrentChunkEnd = 0;

    // Check to see if the read section was set based on the reference name
    // but not yet converted to reference id.
    if(!myRefName.empty())
    {
        myRefID = header.getReferenceID(myRefName.c_str());
        // Clear the myRefName length so this code is only executed once.
        myRefName.clear();

        // Check to see if a reference id was found.
        if(myRefID == SamReferenceInfo::NO_REF_ID)
        {
            myStatus = SamStatus::NO_MORE_RECS;
            return(false);
        }
    }

    // Get the chunks associated with this reference region.
    if(myBamIndex->getChunksForRegion(myRefID, myStartPos, myEndPos, 
                                      myChunksToRead) == true)
    {
        myStatus = SamStatus::SUCCESS;
    }
    else
    {
        String errorMsg = "Failed to get the specified region, refID = ";
        errorMsg += myRefID;
        errorMsg += "; startPos = ";
        errorMsg += myStartPos;
        errorMsg += "; endPos = ";
        errorMsg += myEndPos;
        myStatus.setStatus(SamStatus::FAIL_PARSE, 
                           errorMsg);
    }
    return(true);
}

//
// When the caller to SamFile::ReadRecord() catches an
// exception, it may choose to call this method to resync
// on the underlying binary stream.
//
// Arguments: a callback function that will requires length bytes
// of data to validate a record header.
//
// The expected use case is to re-sync on the next probably valid
// BAM record, so that we can resume reading even after detecting
// a corrupted BAM file.
//
bool SamFile::attemptRecoverySync(bool (*checkSignature)(void *data) , int length)
{
    if(myFilePtr==NULL) return false;
    // non-recovery aware objects will just return false:
    return myFilePtr->attemptRecoverySync(checkSignature, length);
}

// Default Constructor.
SamFileReader::SamFileReader()
    : SamFile()
{
}


// Constructor that opens the specified file for read.
SamFileReader::SamFileReader(const char* filename)
    : SamFile(filename, READ)
{
}




// Constructor that opens the specified file for read.
SamFileReader::SamFileReader(const char* filename,
                             ErrorHandler::HandlingType errorHandlingType)
    : SamFile(filename, READ, errorHandlingType)
{
}


// Constructor that opens the specified file for read.
SamFileReader::SamFileReader(const char* filename,
                             SamFileHeader* header)
    : SamFile(filename, READ, header)
{
}


// Constructor that opens the specified file for read.
SamFileReader::SamFileReader(const char* filename,
                             ErrorHandler::HandlingType errorHandlingType,
                             SamFileHeader* header)
    : SamFile(filename, READ, errorHandlingType, header)
{
}


SamFileReader::~SamFileReader()
{
}


// Default Constructor.
SamFileWriter::SamFileWriter()
    : SamFile()
{
}


// Constructor that opens the specified file for write.
SamFileWriter::SamFileWriter(const char* filename)
    : SamFile(filename, WRITE)
{
}


// Constructor that opens the specified file for write.
SamFileWriter::SamFileWriter(const char* filename,
                             ErrorHandler::HandlingType errorHandlingType)
    : SamFile(filename, WRITE, errorHandlingType)
{
}


// Constructor that opens the specified file for write.
SamFileWriter::SamFileWriter(const char* filename,
                             SamFileHeader* header)
    : SamFile(filename, WRITE, header)
{
}


// Constructor that opens the specified file for write.
SamFileWriter::SamFileWriter(const char* filename,
                             ErrorHandler::HandlingType errorHandlingType,
                             SamFileHeader* header)
    : SamFile(filename, WRITE, errorHandlingType, header)
{
}


SamFileWriter::~SamFileWriter()
{
}
