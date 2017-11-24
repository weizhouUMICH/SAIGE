/*
 *  Copyright (C) 2010-2012  Regents of the University of Michigan
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

#include <stdlib.h>
#include <limits>
#include <stdexcept>

#include "bam.h"

#include "SamRecord.h"
#include "SamValidation.h"

#include "BaseUtilities.h"
#include "SamQuerySeqWithRefHelper.h"

const char* SamRecord::DEFAULT_READ_NAME = "UNKNOWN";
const char* SamRecord::FIELD_ABSENT_STRING = "=";
int SamRecord::myNumWarns = 0;

SamRecord::SamRecord()
    : myStatus(),
      myRefPtr(NULL),
      mySequenceTranslation(NONE)
{
    int32_t defaultAllocSize = DEFAULT_BLOCK_SIZE + sizeof(int32_t);

    myRecordPtr = 
        (bamRecordStruct *) malloc(defaultAllocSize);

    myCigarTempBuffer = NULL;
    myCigarTempBufferAllocatedSize = 0;

    allocatedSize = defaultAllocSize;

    resetRecord();
}


SamRecord::SamRecord(ErrorHandler::HandlingType errorHandlingType)
    : myStatus(errorHandlingType),
      myRefPtr(NULL),
      mySequenceTranslation(NONE)
{
    int32_t defaultAllocSize = DEFAULT_BLOCK_SIZE + sizeof(int32_t);

    myRecordPtr = 
        (bamRecordStruct *) malloc(defaultAllocSize);

    myCigarTempBuffer = NULL;
    myCigarTempBufferAllocatedSize = 0;

    allocatedSize = defaultAllocSize;

    resetRecord();
}


SamRecord::~SamRecord()
{
    resetRecord();

    if(myRecordPtr != NULL)
    {
        free(myRecordPtr);
        myRecordPtr = NULL;
    }
    if(myCigarTempBuffer != NULL)
    {
        free(myCigarTempBuffer);
        myCigarTempBuffer = NULL;
        myCigarTempBufferAllocatedSize = 0;
    }
}


// Resets the fields of the record to a default value.
void SamRecord::resetRecord()
{
    myIsBufferSynced = true;

    myRecordPtr->myBlockSize = DEFAULT_BLOCK_SIZE;
    myRecordPtr->myReferenceID = -1;
    myRecordPtr->myPosition = -1;
    myRecordPtr->myReadNameLength = DEFAULT_READ_NAME_LENGTH;
    myRecordPtr->myMapQuality = 0;
    myRecordPtr->myBin = DEFAULT_BIN;
    myRecordPtr->myCigarLength = 0;
    myRecordPtr->myFlag = 0;
    myRecordPtr->myReadLength = 0;
    myRecordPtr->myMateReferenceID = -1;
    myRecordPtr->myMatePosition = -1;
    myRecordPtr->myInsertSize = 0;
   
    // Set the sam values for the variable length fields.
    // TODO - one way to speed this up might be to not set to "*" and just
    // clear them, and write out a '*' for SAM if it is empty.
    myReadName = DEFAULT_READ_NAME;
    myReferenceName = "*";
    myMateReferenceName = "*";
    myCigar = "*";
    mySequence = "*";
    mySeqWithEq.clear();
    mySeqWithoutEq.clear();
    myQuality = "*";
    myNeedToSetTagsFromBuffer = false;
    myNeedToSetTagsInBuffer = false;

    // Initialize the calculated alignment info to the uncalculated value.
    myAlignmentLength = -1;
    myUnclippedStartOffset = -1;
    myUnclippedEndOffset = -1;

    clearTags();

    // Set the bam values for the variable length fields.
    // Only the read name needs to be set, the others are a length of 0.
    // Set the read name.  The min size of myRecordPtr includes the size for
    // the default read name.
    memcpy(&(myRecordPtr->myData), myReadName.c_str(), 
           myRecordPtr->myReadNameLength);

    // Set that the variable length buffer fields are valid.
    myIsReadNameBufferValid = true;
    myIsCigarBufferValid = true;
    myPackedSequence = 
        (unsigned char *)myRecordPtr->myData + myRecordPtr->myReadNameLength +
        myRecordPtr->myCigarLength * sizeof(int);
    myIsSequenceBufferValid = true;
    myBufferSequenceTranslation = NONE;

    myPackedQuality = myPackedSequence;
    myIsQualityBufferValid = true;
    myIsTagsBufferValid = true;
    myIsBinValid = true;

    myCigarTempBufferLength = -1;

    myStatus = SamStatus::SUCCESS;

    NOT_FOUND_TAG_STRING = "";
    NOT_FOUND_TAG_INT = -1; // TODO - deprecate
}


// Returns whether or not the record is valid.
// Header is needed to perform some validation against it.
bool SamRecord::isValid(SamFileHeader& header)
{
    myStatus = SamStatus::SUCCESS;
    SamValidationErrors invalidSamErrors;
    if(!SamValidator::isValid(header, *this, invalidSamErrors))
    {
        // The record is not valid.
        std::string errorMessage = "";
        invalidSamErrors.getErrorString(errorMessage);
        myStatus.setStatus(SamStatus::INVALID, errorMessage.c_str());
        return(false);
    }
    // The record is valid.
    return(true);
}


void SamRecord::setReference(GenomeSequence* reference)
{
    myRefPtr = reference;
}


// Set the type of sequence translation to use when getting
// the sequence.  The default type (if this method is never called) is
// NONE (the sequence is left as-is).  This is used 
void SamRecord::setSequenceTranslation(SequenceTranslation translation)
{
    mySequenceTranslation = translation;
}


bool SamRecord::setReadName(const char* readName) 
{
    myReadName = readName;
    myIsBufferSynced = false;
    myIsReadNameBufferValid = false;
    myStatus = SamStatus::SUCCESS;

    // The read name must at least have some length, otherwise this is a parsing
    // error.
    if(myReadName.Length() == 0)
    {
        // Invalid - reset ReadName return false.
        myReadName = DEFAULT_READ_NAME;
        myRecordPtr->myReadNameLength = DEFAULT_READ_NAME_LENGTH;
        myStatus.setStatus(SamStatus::INVALID, "0 length Query Name.");
        return(false);
    }

    return true;
}


bool SamRecord::setFlag(uint16_t flag)
{
    myStatus = SamStatus::SUCCESS;
    myRecordPtr->myFlag = flag;
    return true;
}


bool SamRecord::setReferenceName(SamFileHeader& header,
                                 const char* referenceName) 
{
    myStatus = SamStatus::SUCCESS;

    myReferenceName = referenceName;
    // If the reference ID does not already exist, add it (pass true)
    myRecordPtr->myReferenceID = header.getReferenceID(referenceName, true);

    return true;
}


bool SamRecord::set1BasedPosition(int32_t position) 
{
    return(set0BasedPosition(position - 1));
}


bool SamRecord::set0BasedPosition(int32_t position)
{
    myStatus = SamStatus::SUCCESS;
    myRecordPtr->myPosition = position;
    myIsBinValid = false;
    return true;
}


bool SamRecord::setMapQuality(uint8_t mapQuality)
{
    myStatus = SamStatus::SUCCESS;
    myRecordPtr->myMapQuality = mapQuality;
    return true;
}


bool SamRecord::setCigar(const char* cigar)
{
    myStatus = SamStatus::SUCCESS;
    myCigar = cigar;
 
    myIsBufferSynced = false;
    myIsCigarBufferValid = false;
    myCigarTempBufferLength = -1;
    myIsBinValid = false;

    // Initialize the calculated alignment info to the uncalculated value.
    myAlignmentLength = -1;
    myUnclippedStartOffset = -1;
    myUnclippedEndOffset = -1;

    return true;
}


bool SamRecord::setCigar(const Cigar& cigar)
{
    myStatus = SamStatus::SUCCESS;
    cigar.getCigarString(myCigar);
 
    myIsBufferSynced = false;
    myIsCigarBufferValid = false;
    myCigarTempBufferLength = -1;
    myIsBinValid = false;

    // Initialize the calculated alignment info to the uncalculated value.
    myAlignmentLength = -1;
    myUnclippedStartOffset = -1;
    myUnclippedEndOffset = -1;

    return true;
}


bool SamRecord::setMateReferenceName(SamFileHeader& header,
                                     const char* mateReferenceName) 
{
    myStatus = SamStatus::SUCCESS;
    // Set the mate reference, if it is "=", set it to be equal
    // to myReferenceName.  This assumes that myReferenceName has already
    // been called.
    if(strcmp(mateReferenceName, FIELD_ABSENT_STRING) == 0)
    {
        myMateReferenceName = myReferenceName;
    }
    else
    {
        myMateReferenceName = mateReferenceName;
    }

    // Set the Mate Reference ID.
    // If the reference ID does not already exist, add it (pass true)
    myRecordPtr->myMateReferenceID = 
        header.getReferenceID(myMateReferenceName, true);

    return true;
}


bool SamRecord::set1BasedMatePosition(int32_t matePosition)
{
    return(set0BasedMatePosition(matePosition - 1));
}


bool SamRecord::set0BasedMatePosition(int32_t matePosition)
{
    myStatus = SamStatus::SUCCESS;
    myRecordPtr->myMatePosition = matePosition;
    return true;
}


bool SamRecord::setInsertSize(int32_t insertSize)
{
    myStatus = SamStatus::SUCCESS;
    myRecordPtr->myInsertSize = insertSize;
    return true;
}


bool SamRecord::setSequence(const char* seq) 
{
    myStatus = SamStatus::SUCCESS;
    mySequence = seq;
    mySeqWithEq.clear();
    mySeqWithoutEq.clear();
   
    myIsBufferSynced = false;
    myIsSequenceBufferValid = false;
    return true;
}


bool SamRecord::setQuality(const char* quality) 
{
    myStatus = SamStatus::SUCCESS;
    myQuality = quality;
    myIsBufferSynced = false;
    myIsQualityBufferValid = false;
    return true;
}


//Shift indels to the left
bool SamRecord::shiftIndelsLeft()
{
    // Check to see whether or not the Cigar has already been
    // set - this is determined by checking if alignment length
    // is set since alignment length and the cigar are set
    // at the same time.
    if(myAlignmentLength == -1)
    {
        // Not been set, so calculate it.
        parseCigar();
    }
    
    // Track whether or not there was a shift.
    bool shifted = false;

    // Cigar is set, so now myCigarRoller can be used.
    // Track where in the read we are.
    uint32_t currentPos = 0;

    // Since the loop starts at 1 because the first operation can't be shifted,
    // increment the currentPos past the first operation.
    if(Cigar::foundInQuery(myCigarRoller[0]))
    {
        // This op was found in the read, increment the current position.
        currentPos += myCigarRoller[0].count;
    }
   
    int numOps = myCigarRoller.size();
    
    // Loop through the cigar operations from the 2nd operation since
    // the first operation is already on the end and can't shift.
    for(int currentOp = 1; currentOp < numOps; currentOp++)
    {
        if(myCigarRoller[currentOp].operation == Cigar::insert)
        {
            // For now, only shift a max of 1 operation.
            int prevOpIndex = currentOp-1;
            // Track the next op for seeing if it is the same as the
            // previous for merging reasons.
            int nextOpIndex = currentOp+1;
            if(nextOpIndex == numOps)
            {
                // There is no next op, so set it equal to the current one.
                nextOpIndex = currentOp;
            }
            // The start of the previous operation, so we know when we hit it
            // so we don't shift past it.
            uint32_t prevOpStart = 
                currentPos - myCigarRoller[prevOpIndex].count;

            // We can only shift if the previous operation
            if(!Cigar::isMatchOrMismatch(myCigarRoller[prevOpIndex]))
            {
                // TODO - shift past pads
                // An insert is in the read, so increment the position.
                currentPos += myCigarRoller[currentOp].count;                 
                // Not a match/mismatch, so can't shift into it.
                continue;
            }
                    
            // It is a match or mismatch, so check to see if we can
            // shift into it.

            // The end of the insert is calculated by adding the size
            // of this insert minus 1 to the start of the insert.
            uint32_t insertEndPos = 
                currentPos + myCigarRoller[currentOp].count - 1;
                
            // The insert starts at the current position.
            uint32_t insertStartPos = currentPos;
                
            // Loop as long as the position before the insert start
            // matches the last character in the insert. If they match,
            // the insert can be shifted one index left because the
            // implied reference will not change.  If they do not match,
            // we can't shift because the implied reference would change.
            // Stop loop when insertStartPos = prevOpStart, because we 
            // don't want to move past that.
            while((insertStartPos > prevOpStart) && 
                  (getSequence(insertEndPos,BASES) == 
                   getSequence(insertStartPos - 1, BASES)))
            {
                // We can shift, so move the insert start & end one left.
                --insertEndPos;
                --insertStartPos;
            }

            // Determine if a shift has occurred.
            int shiftLen = currentPos - insertStartPos;
            if(shiftLen > 0)
            {
                // Shift occured, so adjust the cigar if the cigar will
                // not become more operations.
                // If the next operation is the same as the previous or
                // if the insert and the previous operation switch positions
                // then the cigar has the same number of operations.
                // If the next operation is different, and the shift splits
                // the previous operation in 2, then the cigar would
                // become longer, so we do not want to shift.
                if(myCigarRoller[nextOpIndex].operation == 
                   myCigarRoller[prevOpIndex].operation)
                {
                    // The operations are the same, so merge them by adding
                    // the length of the shift to the next operation.
                    myCigarRoller.IncrementCount(nextOpIndex, shiftLen);
                    myCigarRoller.IncrementCount(prevOpIndex, -shiftLen);

                    // If the previous op length is 0, just remove that
                    // operation.
                    if(myCigarRoller[prevOpIndex].count == 0)
                    {
                        myCigarRoller.Remove(prevOpIndex);
                    }
                    shifted = true;
                } 
                else
                {
                    // Can only shift if the insert shifts past the
                    // entire previous operation, otherwise an operation
                    // would need to be added.
                    if(insertStartPos == prevOpStart)
                    { 
                        // Swap the positions of the insert and the
                        // previous operation.
                        myCigarRoller.Update(currentOp,
                                             myCigarRoller[prevOpIndex].operation,
                                             myCigarRoller[prevOpIndex].count);
                        // Size of the previous op is the entire
                        // shift length.
                        myCigarRoller.Update(prevOpIndex, 
                                             Cigar::insert,
                                             shiftLen);
                        shifted = true;
                    }
                }
            }
            // An insert is in the read, so increment the position.
            currentPos += myCigarRoller[currentOp].count;                 
        }
        else if(Cigar::foundInQuery(myCigarRoller[currentOp]))
        {
            // This op was found in the read, increment the current position.
            currentPos += myCigarRoller[currentOp].count;
        }
    }
    if(shifted)
    {
        // TODO - setCigar is currently inefficient because later the cigar
        // roller will be recalculated, but for now it will work.
        setCigar(myCigarRoller);
    }
    return(shifted);
}


// Set the BAM record from the passeed in buffer of the specified size.
// Note: The size includes the block size.
SamStatus::Status SamRecord::setBuffer(const char* fromBuffer,
                                       uint32_t fromBufferSize,
                                       SamFileHeader& header)
{
    myStatus = SamStatus::SUCCESS;
    if((fromBuffer == NULL) || (fromBufferSize == 0))
    {
        // Buffer is empty.
        myStatus.setStatus(SamStatus::FAIL_PARSE,
                           "Cannot parse an empty file.");
        return(SamStatus::FAIL_PARSE);
    }

    // Clear the record.   
    resetRecord();

    // allocate space for the record size.
    if(!allocateRecordStructure(fromBufferSize))
    {
        // Failed to allocate space.
        return(SamStatus::FAIL_MEM);
    }
   
    memcpy(myRecordPtr, fromBuffer, fromBufferSize);

    setVariablesForNewBuffer(header);

    // Return the status of the record.
    return(SamStatus::SUCCESS);
}


// Read the BAM record from a file.
SamStatus::Status SamRecord::setBufferFromFile(IFILE filePtr, 
                                               SamFileHeader& header)
{
    myStatus = SamStatus::SUCCESS;
    if((filePtr == NULL) || (filePtr->isOpen() == false))
    {
        // File is not open, return failure.
        myStatus.setStatus(SamStatus::FAIL_ORDER, 
                           "Can't read from an unopened file.");
        return(SamStatus::FAIL_ORDER);
    }

    // Clear the record.
    resetRecord();

    // read the record size.
    int numBytes = 
        ifread(filePtr, &(myRecordPtr->myBlockSize), sizeof(int32_t));

    // Check to see if the end of the file was hit and no bytes were read.
    if(ifeof(filePtr) && (numBytes == 0))
    {
        // End of file, nothing was read, no more records.
            std::string statusMsg = "No more records left to read, ";
            statusMsg += filePtr->getFileName();
            statusMsg += ".";
        myStatus.setStatus(SamStatus::NO_MORE_RECS,
                           statusMsg.c_str());
        return(SamStatus::NO_MORE_RECS);
    }
    
    if(numBytes != sizeof(int32_t))
    {
        // Failed to read the entire block size.  Either the end of the file
        // was reached early or there was an error.
        if(ifeof(filePtr))
        {
            // Error: end of the file reached prior to reading the rest of the
            // record.
            std::string statusMsg = "EOF reached in the middle of a record, ";
            statusMsg += filePtr->getFileName();
            statusMsg += ".";
            myStatus.setStatus(SamStatus::FAIL_PARSE, 
                               statusMsg.c_str());
            return(SamStatus::FAIL_PARSE);
        }
        else
        {
            // Error reading.
            std::string statusMsg = "Failed to read the record size, ";
            statusMsg += filePtr->getFileName();
            statusMsg += ".";
            myStatus.setStatus(SamStatus::FAIL_IO, 
                               statusMsg.c_str());
            return(SamStatus::FAIL_IO);
        }
    }

    // allocate space for the record size.
    if(!allocateRecordStructure(myRecordPtr->myBlockSize + sizeof(int32_t)))
    {
        // Failed to allocate space.
        // Status is set by allocateRecordStructure.
        return(SamStatus::FAIL_MEM);
    }

    // Read the rest of the alignment block, starting at the reference id.
    if(ifread(filePtr, &(myRecordPtr->myReferenceID), myRecordPtr->myBlockSize)
       != (unsigned int)myRecordPtr->myBlockSize)
    {
        // Error reading the record.  Reset it and return failure.
        resetRecord();
        std::string statusMsg = "Failed to read the record, ";
        statusMsg += filePtr->getFileName();
        statusMsg += ".";
        myStatus.setStatus(SamStatus::FAIL_IO, 
                           statusMsg.c_str());
        return(SamStatus::FAIL_IO);
    }

    setVariablesForNewBuffer(header);

    // Return the status of the record.
    return(SamStatus::SUCCESS);
}


// Add the specified tag to the record.
// Returns true if the tag was successfully added, false otherwise.
bool SamRecord::addIntTag(const char* tag, int32_t value)
{
    myStatus = SamStatus::SUCCESS;
    int key = 0;
    int index = 0;
    char bamvtype;

    int tagBufferSize = 0;

    // First check to see if the tags need to be synced to the buffer.
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read tags from the buffer, so cannot add new ones.
            return(false);
        }
    }

    // Ints come in as int.  But it can be represented in fewer bits.
    // So determine a more specific type that is in line with the
    // types for BAM files.
    // First check to see if it is a negative.
    if(value < 0)
    {
        // The int is negative, so it will need to use a signed type.
        // See if it is greater than the min value for a char.
        if(value > ((std::numeric_limits<char>::min)()))
        {
            // It can be stored in a signed char.
            bamvtype = 'c';
            tagBufferSize += 4;
        }
        else if(value > ((std::numeric_limits<short>::min)()))
        {
            // It fits in a signed short.
            bamvtype = 's';
            tagBufferSize += 5;
        }
        else
        {
            // Just store it as a signed int.
            bamvtype = 'i';
            tagBufferSize += 7;
        }
    }
    else
    {
        // It is positive, so an unsigned type can be used.
        if(value < ((std::numeric_limits<unsigned char>::max)()))
        {
            // It is under the max of an unsigned char.
            bamvtype = 'C';
            tagBufferSize += 4;
        }
        else if(value < ((std::numeric_limits<unsigned short>::max)()))
        {
            // It is under the max of an unsigned short.
            bamvtype = 'S';
            tagBufferSize += 5;
        }
        else
        {
            // Just store it as an unsigned int.
            bamvtype = 'I';
            tagBufferSize += 7;
        }
    }

    // Check to see if the tag is already there.
    key = MAKEKEY(tag[0], tag[1], bamvtype);
    unsigned int hashIndex = extras.Find(key);
    if(hashIndex != LH_NOTFOUND)
    {
        // Tag was already found.
        index = extras[hashIndex];
        
        // Since the tagBufferSize was already updated with the new value,
        // subtract the size for the previous tag (even if they are the same).
        switch(intType[index])
        {
            case 'c':
            case 'C':
            case 'A':
                tagBufferSize -= 4;
                break;
            case 's':
            case 'S':
                tagBufferSize -= 5;
                break;
            case 'i':
            case 'I':
                tagBufferSize -= 7;
                break;
            default:
                myStatus.setStatus(SamStatus::INVALID, 
                                   "unknown tag inttype type found.\n");
                return(false);              
        }
            
        // Tag already existed, print message about overwriting.
        // WARN about dropping duplicate tags.
        if(myNumWarns++ < myMaxWarns)
        {
            String newVal;
            String origVal;
            appendIntArrayValue(index, origVal);
            appendIntArrayValue(bamvtype, value, newVal);
            fprintf(stderr, "WARNING: Duplicate Tags, overwritting %c%c:%c:%s with %c%c:%c:%s\n",
                    tag[0], tag[1], intType[index], origVal.c_str(), tag[0], tag[1], bamvtype, newVal.c_str());
            if(myNumWarns == myMaxWarns)
            {
                fprintf(stderr, "Suppressing rest of Duplicate Tag warnings.\n");
            }
        }

        // Update the integer value and type.
        integers[index] = value;
        intType[index] = bamvtype;
    }
    else
    {
        // Tag is not already there, so add it.
        index = integers.Length();
        
        integers.Push(value);
        intType.push_back(bamvtype);

        extras.Add(key, index);
    }

    // The buffer tags are now out of sync.
    myNeedToSetTagsInBuffer = true;
    myIsTagsBufferValid = false;
    myIsBufferSynced = false;
    myTagBufferSize += tagBufferSize;

    return(true);
}


// Add the specified tag to the record, replacing it if it is already there and
// is different from the previous value.
// Returns true if the tag was successfully added (or was already there), false otherwise.
bool SamRecord::addTag(const char* tag, char vtype, const char* valuePtr)
{
    if(vtype == 'i')
    {
        // integer type.  Call addIntTag to handle it.
        int intVal = atoi(valuePtr);
        return(addIntTag(tag, intVal));
    }

    // Non-int type.
    myStatus = SamStatus::SUCCESS;
    bool status = true; // default to successful.
    int key = 0;
    int index = 0;

    int tagBufferSize = 0;

    // First check to see if the tags need to be synced to the buffer.
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read tags from the buffer, so cannot add new ones.
            return(false);
        }
    }

    // First check to see if the tag is already there.
    key = MAKEKEY(tag[0], tag[1], vtype);
    unsigned int hashIndex = extras.Find(key);
    if(hashIndex != LH_NOTFOUND)
    {
        // The key was found in the hash, so get the lookup index.
        index = extras[hashIndex];

        String origTag;
        char origType = vtype;

        // Adjust the currently pointed to value to the new setting.
        switch (vtype)
        {
            case 'A' :
                // First check to see if the value changed.
                if((integers[index] == (const int)*(valuePtr)) &&
                   (intType[index] == vtype))
                {
                    // The value & type has not changed, so do nothing.
                    return(true);
                }
                else
                {
                    // Tag buffer size changes if type changes, so subtract & add.
                    origType = intType[index];
                    appendIntArrayValue(index, origTag);
                    tagBufferSize -= getNumericTagTypeSize(intType[index]);
                    tagBufferSize += getNumericTagTypeSize(vtype);
                    integers[index] = (const int)*(valuePtr);
                    intType[index] = vtype;
                }
                break;
            case 'Z' :
                // First check to see if the value changed.
                if(strings[index] == valuePtr)
                {
                    // The value has not changed, so do nothing.
                    return(true);
                }
                else
                {
                    // Adjust the tagBufferSize by removing the size of the old string.
                    origTag = strings[index];
                    tagBufferSize -= strings[index].Length();
                    strings[index] = valuePtr;
                    // Adjust the tagBufferSize by adding the size of the new string.
                    tagBufferSize += strings[index].Length();
                }
                break;
            case 'B' :
                // First check to see if the value changed.
                if(strings[index] == valuePtr)
                {
                    // The value has not changed, so do nothing.
                    return(true);
                }
                else
                {
                    // Adjust the tagBufferSize by removing the size of the old field.
                    origTag = strings[index];
                    tagBufferSize -= getBtagBufferSize(strings[index]);
                    strings[index] = valuePtr;
                    // Adjust the tagBufferSize by adding the size of the new field.
                    tagBufferSize += getBtagBufferSize(strings[index]);
                }
                break;
            case 'f' :
                // First check to see if the value changed.
                if(floats[index] == (float)atof(valuePtr))
                {
                    // The value has not changed, so do nothing.
                    return(true);
                }
                else
                {
                    // Tag buffer size doesn't change between different 'f' entries.
                    origTag.appendFullFloat(floats[index]);
                    floats[index] = (float)atof(valuePtr);
                }
                break;
            default :
                fprintf(stderr,
                        "samRecord::addTag() - Unknown custom field of type %c\n",
                        vtype);
                myStatus.setStatus(SamStatus::FAIL_PARSE, 
                                   "Unknown custom field in a tag");
                status = false;
                break;
        }

        // Duplicate tag in this record.
        // Tag already existed, print message about overwriting.
        // WARN about dropping duplicate tags.
        if(myNumWarns++ < myMaxWarns)
        {
            fprintf(stderr, "WARNING: Duplicate Tags, overwritting %c%c:%c:%s with %c%c:%c:%s\n",
                    tag[0], tag[1], origType, origTag.c_str(), tag[0], tag[1], vtype, valuePtr);
            if(myNumWarns == myMaxWarns)
            {
                fprintf(stderr, "Suppressing rest of Duplicate Tag warnings.\n");
            }
        }
    }
    else
    {
        // The key was not found in the hash, so add it.
        switch (vtype)
        {
            case 'A' :
                index = integers.Length();
                integers.Push((const int)*(valuePtr));
                intType.push_back(vtype);
                tagBufferSize += 4;
                break;
            case 'Z' :
                index = strings.Length();
                strings.Push(valuePtr);
                tagBufferSize += 4 + strings.Last().Length();
                break;
            case 'B' :
                index = strings.Length();
                strings.Push(valuePtr);
                tagBufferSize += 3 + getBtagBufferSize(strings[index]);
                break;
            case 'f' :
                index = floats.size();
                floats.push_back((float)atof(valuePtr));
                tagBufferSize += 7;
                break;
            default :
                fprintf(stderr,
                        "samRecord::addTag() - Unknown custom field of type %c\n",
                        vtype);
                myStatus.setStatus(SamStatus::FAIL_PARSE, 
                                   "Unknown custom field in a tag");
                status = false;
                break;
        }
        if(status)
        {
            // If successful, add the key to extras.
            extras.Add(key, index);
        }
    }

    // Only add the tag if it has so far been successfully processed.
    if(status)
    {
        // The buffer tags are now out of sync.
        myNeedToSetTagsInBuffer = true;
        myIsTagsBufferValid = false;
        myIsBufferSynced = false;
        myTagBufferSize += tagBufferSize;
    }
    return(status);
}


void SamRecord::clearTags()
{
    if(extras.Entries() != 0)
    {
        extras.Clear();
    }
    strings.Clear();
    integers.Clear();
    intType.clear();
    floats.clear();
    myTagBufferSize = 0;
    resetTagIter();
}


bool SamRecord::rmTag(const char* tag, char type)
{
    // Check the length of tag.
    if(strlen(tag) != 2)
    {
        // Tag is the wrong length.
        myStatus.setStatus(SamStatus::INVALID, 
                           "rmTag called with tag that is not 2 characters\n");
        return(false);
    }

    myStatus = SamStatus::SUCCESS;
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.
            return(false);
        }
    }

    // Construct the key.
    int key = MAKEKEY(tag[0], tag[1], type);
    // Look to see if the key exsists in the hash.
    int offset = extras.Find(key);

    if(offset < 0)
    {
        // Not found, so return true, successfully removed since
        // it is not in tag.
        return(true);
    }

    // Offset is set, so the key was found.
    // First if it is an integer, determine the actual type of the int.
    char vtype;
    getTypeFromKey(key, vtype);
    if(vtype == 'i')
    {
        vtype = getIntegerType(offset);
    }

    // Offset is set, so recalculate the buffer size without this entry.
    // Do NOT remove from strings, integers, or floats because then
    // extras would need to be updated for all entries with the new indexes
    // into those variables.
    int rmBuffSize = 0;
    switch(vtype)
    {
        case 'A':
        case 'c':
        case 'C':
            rmBuffSize = 4;
            break;
        case 's':
        case 'S':
            rmBuffSize = 5;
            break;
        case 'i':
        case 'I':
            rmBuffSize = 7;
            break;
        case 'f':
            rmBuffSize = 7;
            break;
        case 'Z':
            rmBuffSize = 4 + getString(offset).Length();
            break;
        case 'B':
            rmBuffSize = 3 + getBtagBufferSize(getString(offset));
            break;
        default:
            myStatus.setStatus(SamStatus::INVALID, 
                               "rmTag called with unknown type.\n");
            return(false);
            break;
    };

    // The buffer tags are now out of sync.
    myNeedToSetTagsInBuffer = true;
    myIsTagsBufferValid = false;
    myIsBufferSynced = false;
    myTagBufferSize -= rmBuffSize;

    // Remove from the hash.
    extras.Delete(offset);
    return(true);
}


bool SamRecord::rmTags(const char* tags)
{
    const char* currentTagPtr = tags;

    myStatus = SamStatus::SUCCESS;
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.
            return(false);
        }
    }
    
    bool returnStatus = true;

    int rmBuffSize = 0;
    while(*currentTagPtr != '\0')
    {

        // Tags are formatted as: XY:Z
        // Where X is [A-Za-z], Y is [A-Za-z], and
        // Z is A,i,f,Z,H (cCsSI are also excepted)
        if((currentTagPtr[0] == '\0') || (currentTagPtr[1] == '\0') ||
           (currentTagPtr[2] != ':') || (currentTagPtr[3] == '\0'))
        {
            myStatus.setStatus(SamStatus::INVALID, 
                               "rmTags called with improperly formatted tags.\n");
            returnStatus = false;
            break;
        }

        // Construct the key.
        int key = MAKEKEY(currentTagPtr[0], currentTagPtr[1], 
                          currentTagPtr[3]);
        // Look to see if the key exsists in the hash.
        int offset = extras.Find(key);

        if(offset >= 0)
        {
            // Offset is set, so the key was found.
            // First if it is an integer, determine the actual type of the int.
            char vtype;
            getTypeFromKey(key, vtype);
            if(vtype == 'i')
            {
                vtype = getIntegerType(offset);
            }
            
            // Offset is set, so recalculate the buffer size without this entry.
            // Do NOT remove from strings, integers, or floats because then
            // extras would need to be updated for all entries with the new indexes
            // into those variables.
            switch(vtype)
            {
                case 'A':
                case 'c':
                case 'C':
                    rmBuffSize += 4;
                    break;
                case 's':
                case 'S':
                    rmBuffSize += 5;
                    break;
                case 'i':
                case 'I':
                    rmBuffSize += 7;
                    break;
                case 'f':
                    rmBuffSize += 7;
                    break;
                case 'Z':
                    rmBuffSize += 4 + getString(offset).Length();
                    break;
                case 'B':
                    rmBuffSize += 3 + getBtagBufferSize(getString(offset));
                    break;
                default:
                    myStatus.setStatus(SamStatus::INVALID, 
                                       "rmTag called with unknown type.\n");
                    returnStatus = false;
                    break;
            };
            
            // Remove from the hash.
            extras.Delete(offset);
        }
        // Increment to the next tag.
        if((currentTagPtr[4] == ';') || (currentTagPtr[4] == ','))
        {
            // Increment once more.
            currentTagPtr += 5;
        }
        else if(currentTagPtr[4] != '\0')
        {
            // Invalid tag format. 
            myStatus.setStatus(SamStatus::INVALID, 
                               "rmTags called with improperly formatted tags.\n");
            returnStatus = false;
            break;
        }
        else
        {
            // Last Tag.
            currentTagPtr += 4;
        }
    }

    // The buffer tags are now out of sync.
    myNeedToSetTagsInBuffer = true;
    myIsTagsBufferValid = false;
    myIsBufferSynced = false;
    myTagBufferSize -= rmBuffSize;
    

    return(returnStatus);
}


// Get methods for record fields.
const void* SamRecord::getRecordBuffer()
{
    return(getRecordBuffer(mySequenceTranslation));
}


// Get methods for record fields.
const void* SamRecord::getRecordBuffer(SequenceTranslation translation) 
{
    myStatus = SamStatus::SUCCESS;
    bool status = true;
    // If the buffer is not synced or the sequence in the buffer is not
    // properly translated, fix the buffer.
    if((myIsBufferSynced == false) ||
       (myBufferSequenceTranslation != translation))
    {
        status &= fixBuffer(translation);
    }
    // If the buffer is synced, check to see if the tags need to be synced.
    if(myNeedToSetTagsInBuffer)
    {
        status &= setTagsInBuffer();
    }
    if(!status)
    {
        return(NULL);
    }
    return (const void *)myRecordPtr;
}


// Write the record as a buffer into the file using the class's 
// sequence translation setting.
SamStatus::Status SamRecord::writeRecordBuffer(IFILE filePtr)
{
    return(writeRecordBuffer(filePtr, mySequenceTranslation));
}


// Write the record as a buffer into the file using the specified translation.
SamStatus::Status SamRecord::writeRecordBuffer(IFILE filePtr, 
                                               SequenceTranslation translation)
{
    myStatus = SamStatus::SUCCESS;
    if((filePtr == NULL) || (filePtr->isOpen() == false))
    {
        // File is not open, return failure.
        myStatus.setStatus(SamStatus::FAIL_ORDER,
                           "Can't write to an unopened file.");
        return(SamStatus::FAIL_ORDER);
    }

    if((myIsBufferSynced == false) ||
       (myBufferSequenceTranslation != translation))
    {
        if(!fixBuffer(translation))
        {
            return(myStatus.getStatus());
        }
    }

    // Write the record.
    unsigned int numBytesToWrite = myRecordPtr->myBlockSize + sizeof(int32_t);
    unsigned int numBytesWritten = 
        ifwrite(filePtr, myRecordPtr, numBytesToWrite);

    // Return status based on if the correct number of bytes were written.
    if(numBytesToWrite == numBytesWritten)
    {
        return(SamStatus::SUCCESS);
    }
    // The correct number of bytes were not written.
    myStatus.setStatus(SamStatus::FAIL_IO, "Failed to write the entire record.");
    return(SamStatus::FAIL_IO);
}


int32_t SamRecord::getBlockSize() 
{
    myStatus = SamStatus::SUCCESS;
    // If the buffer isn't synced, sync the buffer to determine the
    // block size.
    if(myIsBufferSynced == false)
    {
        // Since this just returns the block size, the translation of
        // the sequence does not matter, so just use the currently set
        // value.
        fixBuffer(myBufferSequenceTranslation);
    }
    return myRecordPtr->myBlockSize;
}


// This method returns the reference name.
const char* SamRecord::getReferenceName()
{
    myStatus = SamStatus::SUCCESS;
    return myReferenceName.c_str();
}


int32_t SamRecord::getReferenceID()
{
    myStatus = SamStatus::SUCCESS;
    return myRecordPtr->myReferenceID;
}


int32_t SamRecord::get1BasedPosition()
{
    myStatus = SamStatus::SUCCESS;
    return (myRecordPtr->myPosition + 1);
}


int32_t SamRecord::get0BasedPosition()
{
    myStatus = SamStatus::SUCCESS;
    return myRecordPtr->myPosition;
}


uint8_t SamRecord::getReadNameLength() 
{
    myStatus = SamStatus::SUCCESS;
    // If the buffer is valid, return the size from there, otherwise get the 
    // size from the string length + 1 (ending null).
    if(myIsReadNameBufferValid)
    {
        return(myRecordPtr->myReadNameLength);
    }
   
    return(myReadName.Length() + 1);
}


uint8_t SamRecord::getMapQuality()
{
    myStatus = SamStatus::SUCCESS;
    return myRecordPtr->myMapQuality;
}


uint16_t SamRecord::getBin()
{
    myStatus = SamStatus::SUCCESS;
    if(!myIsBinValid)
    {
        // The bin that is set in the record is not valid, so
        // reset it.
        myRecordPtr->myBin = 
            bam_reg2bin(myRecordPtr->myPosition, get1BasedAlignmentEnd());      
        myIsBinValid = true;
    }
    return(myRecordPtr->myBin);
}


uint16_t SamRecord::getCigarLength()
{
    myStatus = SamStatus::SUCCESS;
    // If the cigar buffer is valid
    // then get the length from there.
    if(myIsCigarBufferValid)
    {
        return myRecordPtr->myCigarLength;      
    }

    if(myCigarTempBufferLength == -1)
    {
        // The cigar buffer is not valid and the cigar temp buffer is not set,
        // so parse the string.
        parseCigarString();
    }
   
    // The temp buffer is now set, so return the size.
    return(myCigarTempBufferLength);
}


uint16_t SamRecord::getFlag()
{
    myStatus = SamStatus::SUCCESS;
    return myRecordPtr->myFlag;
}


int32_t SamRecord::getReadLength() 
{
    myStatus = SamStatus::SUCCESS;
    if(myIsSequenceBufferValid == false)
    {
        // If the sequence is "*", then return 0.
        if((mySequence.Length() == 1) && (mySequence[0] == '*'))
        {
            return(0);
        }
        // Do not add 1 since it is not null terminated.
        return(mySequence.Length());
    }
    return(myRecordPtr->myReadLength);
}


// This method returns the mate reference name.  If it is equal to the
// reference name, it still returns the reference name.
const char* SamRecord::getMateReferenceName()
{
    myStatus = SamStatus::SUCCESS;
    return myMateReferenceName.c_str();
}


// This method returns the mate reference name.  If it is equal to the
// reference name, it returns "=", unless they are both "*" in which case
// "*" is returned.
const char* SamRecord::getMateReferenceNameOrEqual()
{
    myStatus = SamStatus::SUCCESS;
    if(myMateReferenceName == "*")
    {
        return(myMateReferenceName);
    }
    if(myMateReferenceName == getReferenceName())
    {
        return(FIELD_ABSENT_STRING);
    }
    else
    {
        return(myMateReferenceName);
    }
}


int32_t SamRecord::getMateReferenceID()
{
    myStatus = SamStatus::SUCCESS;
    return myRecordPtr->myMateReferenceID;
}


int32_t SamRecord::get1BasedMatePosition()
{
    myStatus = SamStatus::SUCCESS;
    return (myRecordPtr->myMatePosition + 1);
}


int32_t SamRecord::get0BasedMatePosition()
{
    myStatus = SamStatus::SUCCESS;
    return myRecordPtr->myMatePosition;
}


int32_t SamRecord::getInsertSize()
{
    myStatus = SamStatus::SUCCESS;
    return myRecordPtr->myInsertSize;
}


// Returns the inclusive rightmost position of the clipped sequence.
int32_t SamRecord::get0BasedAlignmentEnd()
{
    myStatus = SamStatus::SUCCESS;
    if(myAlignmentLength == -1)
    {
        // Alignment end has not been set, so calculate it.
        parseCigar();
    }
    // If alignment length > 0, subtract 1 from it to get the end.
    if(myAlignmentLength == 0)
    {
        // Length is 0, just return the start position.
        return(myRecordPtr->myPosition);
    }
    return(myRecordPtr->myPosition + myAlignmentLength - 1);
}


// Returns the inclusive rightmost position of the clipped sequence.
int32_t SamRecord::get1BasedAlignmentEnd()
{
    return(get0BasedAlignmentEnd() + 1);
}

   
// Return the length of the alignment.
int32_t SamRecord::getAlignmentLength()
{
    myStatus = SamStatus::SUCCESS;
    if(myAlignmentLength == -1)
    {
        // Alignment end has not been set, so calculate it.
        parseCigar();
    }
    // Return the alignment length.
    return(myAlignmentLength);
}

// Returns the inclusive left-most position adjust for clipped bases.
int32_t SamRecord::get0BasedUnclippedStart()
{
    myStatus = SamStatus::SUCCESS;
    if(myUnclippedStartOffset == -1)
    {
        // Unclipped has not yet been calculated, so parse the cigar to get it
        parseCigar();
    }
    return(myRecordPtr->myPosition - myUnclippedStartOffset);
}


// Returns the inclusive left-most position adjust for clipped bases.
int32_t SamRecord::get1BasedUnclippedStart()
{
    return(get0BasedUnclippedStart() + 1);
}


// Returns the inclusive right-most position adjust for clipped bases.
int32_t SamRecord::get0BasedUnclippedEnd()
{
    // myUnclippedEndOffset will be set by get0BasedAlignmentEnd if the 
    // cigar has not yet been parsed, so no need to check it here.
    return(get0BasedAlignmentEnd() + myUnclippedEndOffset);
}


// Returns the inclusive right-most position adjust for clipped bases.
int32_t SamRecord::get1BasedUnclippedEnd()
{
    return(get0BasedUnclippedEnd() + 1);
}


// Get the read name.
const char* SamRecord::getReadName()
{
    myStatus = SamStatus::SUCCESS;
    if(myReadName.Length() == 0)
    {
        // 0 Length, means that it is in the buffer, but has not yet
        // been synced to the string, so do the sync.
        myReadName = (char*)&(myRecordPtr->myData);
    }
    return myReadName.c_str();
}


const char* SamRecord::getCigar()
{
    myStatus = SamStatus::SUCCESS;
    if(myCigar.Length() == 0)
    {
        // 0 Length, means that it is in the buffer, but has not yet
        // been synced to the string, so do the sync.
        parseCigarBinary();
    }
    return myCigar.c_str();
}


const char* SamRecord::getSequence()
{
    return(getSequence(mySequenceTranslation));
}


const char* SamRecord::getSequence(SequenceTranslation translation)
{
    myStatus = SamStatus::SUCCESS;
    if(mySequence.Length() == 0)
    {
        // 0 Length, means that it is in the buffer, but has not yet
        // been synced to the string, so do the sync.
        setSequenceAndQualityFromBuffer();
    }

    // Determine if translation needs to be done.
    if((translation == NONE) || (myRefPtr == NULL))
    {
        return mySequence.c_str();
    }
    else if(translation == EQUAL)
    {
        if(mySeqWithEq.length() == 0)
        {
            // Check to see if the sequence is defined.
            if(mySequence == "*")
            {
                // Sequence is undefined, so no translation necessary.
                mySeqWithEq = '*';
            }
            else
            {
                // Sequence defined, so translate it.
                SamQuerySeqWithRef::seqWithEquals(mySequence.c_str(), 
                                                  myRecordPtr->myPosition,
                                                  *(getCigarInfo()),
                                                  getReferenceName(),
                                                  *myRefPtr,
                                                  mySeqWithEq);
            }
        }
        return(mySeqWithEq.c_str());
    }
    else
    {
        // translation == BASES
        if(mySeqWithoutEq.length() == 0)
        {
            if(mySequence == "*")
            {
                // Sequence is undefined, so no translation necessary.
                mySeqWithoutEq = '*';
            }
            else
            {
                // Sequence defined, so translate it.
                SamQuerySeqWithRef::seqWithoutEquals(mySequence.c_str(), 
                                                     myRecordPtr->myPosition,
                                                     *(getCigarInfo()),
                                                     getReferenceName(),
                                                     *myRefPtr,
                                                     mySeqWithoutEq);
            }
        }
        return(mySeqWithoutEq.c_str());
    }
}


const char* SamRecord::getQuality() 
{
    myStatus = SamStatus::SUCCESS;
    if(myQuality.Length() == 0)
    {
        // 0 Length, means that it is in the buffer, but has not yet
        // been synced to the string, so do the sync.
        setSequenceAndQualityFromBuffer();      
    }
    return myQuality.c_str();
}


char SamRecord::getSequence(int index)
{
    return(getSequence(index, mySequenceTranslation));
}


char SamRecord::getSequence(int index, SequenceTranslation translation)
{
    static const char * asciiBases = "=AC.G...T......N";

    // Determine the read length.
    int32_t readLen = getReadLength();

    // If the read length is 0, this method should not be called.
    if(readLen == 0)
    {
        String exceptionString = "SamRecord::getSequence(";
        exceptionString += index;
        exceptionString += ") is not allowed since sequence = '*'";
        throw std::runtime_error(exceptionString.c_str());
    }
    else if((index < 0) || (index >= readLen))
    {
        // Only get here if the index was out of range, so thow an exception.
        String exceptionString = "SamRecord::getSequence(";
        exceptionString += index;
        exceptionString += ") is out of range. Index must be between 0 and ";
        exceptionString += (readLen - 1);
        throw std::runtime_error(exceptionString.c_str());
    }

    // Determine if translation needs to be done.
    if((translation == NONE) || (myRefPtr == NULL))
    {
        // No translation needs to be done.
        if(mySequence.Length() == 0)
        {
            // Parse BAM sequence.
            if(myIsSequenceBufferValid)
            {
                return(index & 1 ?
                       asciiBases[myPackedSequence[index / 2] & 0xF] :
                       asciiBases[myPackedSequence[index / 2] >> 4]);
            }
            else
            {
                String exceptionString = "SamRecord::getSequence(";
                exceptionString += index;
                exceptionString += ") called with no sequence set";
                throw std::runtime_error(exceptionString.c_str());
            }
        }
        // Already have string.
        return(mySequence[index]);
    }
    else
    {
        // Need to translate the sequence either to have '=' or to not
        // have it.
        // First check to see if the sequence has been set.
        if(mySequence.Length() == 0)
        {
            // 0 Length, means that it is in the buffer, but has not yet
            // been synced to the string, so do the sync.
            setSequenceAndQualityFromBuffer();
        }

        // Check the type of translation.
        if(translation == EQUAL)
        {
            // Check whether or not the string has already been 
            // retrieved that has the '=' in it.
            if(mySeqWithEq.length() == 0)
            {
                // The string with '=' has not yet been determined,
                // so get the string.
                // Check to see if the sequence is defined.
                if(mySequence == "*")
                {
                    // Sequence is undefined, so no translation necessary.
                    mySeqWithEq = '*';
                }
                else
                {
                    // Sequence defined, so translate it.
                    SamQuerySeqWithRef::seqWithEquals(mySequence.c_str(), 
                                                      myRecordPtr->myPosition, 
                                                      *(getCigarInfo()),
                                                      getReferenceName(),
                                                      *myRefPtr,
                                                      mySeqWithEq);
                }
            }
            // Sequence is set, so return it.
            return(mySeqWithEq[index]);
        }
        else
        {
            // translation == BASES
            // Check whether or not the string has already been 
            // retrieved that does not have the '=' in it.
            if(mySeqWithoutEq.length() == 0)
            {
                // The string with '=' has not yet been determined,
                // so get the string.
                // Check to see if the sequence is defined.
                if(mySequence == "*")
                {
                    // Sequence is undefined, so no translation necessary.
                    mySeqWithoutEq = '*';
                }
                else
                {
                    // Sequence defined, so translate it.
                    // The string without '=' has not yet been determined,
                    // so get the string.
                    SamQuerySeqWithRef::seqWithoutEquals(mySequence.c_str(), 
                                                         myRecordPtr->myPosition, 
                                                         *(getCigarInfo()),
                                                         getReferenceName(),
                                                         *myRefPtr,
                                                         mySeqWithoutEq);
                }
            }
            // Sequence is set, so return it.
            return(mySeqWithoutEq[index]);
        }
    }
}


char SamRecord::getQuality(int index)
{
    // Determine the read length.
    int32_t readLen = getReadLength();

    // If the read length is 0, return ' ' whose ascii code is below
    // the minimum ascii code for qualities.
    if(readLen == 0)
    {
        return(BaseUtilities::UNKNOWN_QUALITY_CHAR);
    }
    else if((index < 0) || (index >= readLen))
    {
        // Only get here if the index was out of range, so thow an exception.
        String exceptionString = "SamRecord::getQuality(";
        exceptionString += index;
        exceptionString += ") is out of range. Index must be between 0 and ";
        exceptionString += (readLen - 1);
        throw std::runtime_error(exceptionString.c_str());
    }

    if(myQuality.Length() == 0) 
    {
        // Parse BAM Quality.
        // Know that myPackedQuality is correct since readLen != 0.
        return(myPackedQuality[index] + 33);
    }
    else
    {
        // Already have string.
        if((myQuality.Length() == 1) && (myQuality[0] == '*'))
        {
            // Return the unknown quality character.
            return(BaseUtilities::UNKNOWN_QUALITY_CHAR);
        }
        else if(index >= myQuality.Length())
        {
            // Only get here if the index was out of range, so thow an exception.
            // Technically the myQuality string is not guaranteed to be the same length
            // as the sequence, so this catches that error.
            String exceptionString = "SamRecord::getQuality(";
            exceptionString += index;
            exceptionString += ") is out of range. Index must be between 0 and ";
            exceptionString += (myQuality.Length() - 1);
            throw std::runtime_error(exceptionString.c_str());
        }
        else
        {
            return(myQuality[index]);
        }
    }
}

   
Cigar* SamRecord::getCigarInfo()
{
    // Check to see whether or not the Cigar has already been
    // set - this is determined by checking if alignment length
    // is set since alignment length and the cigar are set
    // at the same time.
    if(myAlignmentLength == -1)
    {
        // Not been set, so calculate it.
        parseCigar();
    }
    return(&myCigarRoller);
}


// Return the number of bases in this read that overlap the passed in
// region.  (start & end are 0-based)
uint32_t SamRecord::getNumOverlaps(int32_t start, int32_t end)
{
    // Determine whether or not the cigar has been parsed, which sets up
    // the cigar roller.  This is determined by checking the alignment length.
    if(myAlignmentLength == -1)
    {
        parseCigar();
    }
    return(myCigarRoller.getNumOverlaps(start, end, get0BasedPosition()));
}


// Returns the values of all fields except the tags.
bool SamRecord::getFields(bamRecordStruct& recStruct, String& readName, 
                          String& cigar, String& sequence, String& quality)
{
    return(getFields(recStruct, readName, cigar, sequence, quality,
                     mySequenceTranslation));
}


// Returns the values of all fields except the tags.
bool SamRecord::getFields(bamRecordStruct& recStruct, String& readName, 
                          String& cigar, String& sequence, String& quality,
                          SequenceTranslation translation)
{
    myStatus = SamStatus::SUCCESS;
    if(myIsBufferSynced == false)
    {
        if(!fixBuffer(translation))
        {
            // failed to set the buffer, return false.
            return(false);
        }
    }
    memcpy(&recStruct, myRecordPtr, sizeof(bamRecordStruct));

    readName = getReadName();
    // Check the status.
    if(myStatus != SamStatus::SUCCESS)
    {
        // Failed to set the fields, return false.
        return(false);
    }
    cigar = getCigar();
    // Check the status.
    if(myStatus != SamStatus::SUCCESS)
    {
        // Failed to set the fields, return false.
        return(false);
    }
    sequence = getSequence(translation);
    // Check the status.
    if(myStatus != SamStatus::SUCCESS)
    {
        // Failed to set the fields, return false.
        return(false);
    }
    quality = getQuality();
    // Check the status.
    if(myStatus != SamStatus::SUCCESS)
    {
        // Failed to set the fields, return false.
        return(false);
    }
    return(true);
}


// Returns the reference pointer.
GenomeSequence* SamRecord::getReference()
{
    return(myRefPtr);
}


uint32_t SamRecord::getTagLength()
{
    myStatus = SamStatus::SUCCESS;
    if(myNeedToSetTagsFromBuffer)
    {
        // Tags are only set in the buffer, so the size of the tags is 
        // the length of the record minus the starting location of the tags.
        unsigned char * tagStart = 
            (unsigned char *)myRecordPtr->myData 
            + myRecordPtr->myReadNameLength 
            + myRecordPtr->myCigarLength * sizeof(int)
            + (myRecordPtr->myReadLength + 1) / 2 + myRecordPtr->myReadLength;
      
        // The non-tags take up from the start of the record to the tag start.
        // Do not include the block size part of the record since it is not
        // included in the size.
        uint32_t nonTagSize = 
            tagStart - (unsigned char*)&(myRecordPtr->myReferenceID);
        // Tags take up the size of the block minus the non-tag section.
        uint32_t tagSize = myRecordPtr->myBlockSize - nonTagSize;
        return(tagSize);
    }

    // Tags are stored outside the buffer, so myTagBufferSize is set.
    return(myTagBufferSize);
}


// Returns true if there is another tag and sets tag and vtype to the
// appropriate values, and returns a pointer to the value.
// Sets the Status to SUCCESS when a tag is successfully returned or
// when there are no more tags.  Otherwise the status is set to describe
// why it failed (parsing, etc).
bool SamRecord::getNextSamTag(char* tag, char& vtype, void** value)
{
    myStatus = SamStatus::SUCCESS;
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.
            return(false);
        }
    }

    // Increment the tag index to start looking at the next tag.
    // At the beginning, it is set to -1.
    myLastTagIndex++;
    int maxTagIndex = extras.Capacity();
    if(myLastTagIndex >= maxTagIndex)
    {
        // Hit the end of the tags, return false, no more tags.
        // Status is still success since this is not an error, 
        // it is just the end of the list.
        return(false);
    }

    bool tagFound = false;
    // Loop until a tag is found or the end of extras is hit.
    while((tagFound == false) && (myLastTagIndex < maxTagIndex))
    {
        if(extras.SlotInUse(myLastTagIndex))
        {
            // Found a slot to use.
            int key = extras.GetKey(myLastTagIndex);
            getTag(key, tag);
            getTypeFromKey(key, vtype);
            tagFound = true;
            // Get the value associated with the key based on the vtype.
            switch (vtype)
            {
                case 'f' :
                    *value = getFloatPtr(myLastTagIndex);
                    break;
                case 'i' :
                    *value = getIntegerPtr(myLastTagIndex, vtype);
                    if(vtype != 'A')
                    {
                        // Convert all int types to 'i'
                        vtype = 'i';
                    }
                    break;
                case 'Z' :
                case 'B' :
                    *value = getStringPtr(myLastTagIndex);
                    break;
                default:
                    myStatus.setStatus(SamStatus::FAIL_PARSE,
                                       "Unknown tag type");
                    tagFound = false;
                    break;
            }
        }
        if(!tagFound)
        {
            // Increment the index since a tag was not found.
            myLastTagIndex++;
        }
    }
    return(tagFound);
}


// Reset the tag iterator to the beginning of the tags.
void SamRecord::resetTagIter()
{
    myLastTagIndex = -1;
}


bool SamRecord::isIntegerType(char vtype)
{
    if((vtype == 'c') || (vtype == 'C') ||
       (vtype == 's') || (vtype == 'S') ||
       (vtype == 'i') || (vtype == 'I'))
    {
        return(true);
    }
    return(false);
}


bool SamRecord::isFloatType(char vtype)
{
    if(vtype == 'f')
    {
        return(true);
    }
    return(false);
}


bool SamRecord::isCharType(char vtype)
{
    if(vtype == 'A')
    {
        return(true);
    }
    return(false);
}


bool SamRecord::isStringType(char vtype)
{
    if((vtype == 'Z') || (vtype == 'B'))
    {
        return(true);
    }
    return(false);
}


bool SamRecord::getTagsString(const char* tags, String& returnString, char delim)
{
    const char* currentTagPtr = tags;

    returnString.Clear();
    myStatus = SamStatus::SUCCESS;
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.
            return(false);
        }
    }
    
    bool returnStatus = true;

    while(*currentTagPtr != '\0')
    {
        // Tags are formatted as: XY:Z
        // Where X is [A-Za-z], Y is [A-Za-z], and
        // Z is A,i,f,Z,H (cCsSI are also excepted)
        if((currentTagPtr[0] == '\0') || (currentTagPtr[1] == '\0') ||
           (currentTagPtr[2] != ':') || (currentTagPtr[3] == '\0'))
        {
            myStatus.setStatus(SamStatus::INVALID, 
                               "getTagsString called with improperly formatted tags.\n");
            returnStatus = false;
            break;
        }

        // Construct the key.
        int key = MAKEKEY(currentTagPtr[0], currentTagPtr[1], 
                          currentTagPtr[3]);
        // Look to see if the key exsists in the hash.
        int offset = extras.Find(key);

        if(offset >= 0)
        {
            // Offset is set, so the key was found.
            if(!returnString.IsEmpty())
            {
                returnString += delim;
            }
            returnString += currentTagPtr[0];
            returnString += currentTagPtr[1];
            returnString += ':';
            returnString += currentTagPtr[3];
            returnString += ':';

            // First if it is an integer, determine the actual type of the int.
            char vtype;
            getTypeFromKey(key, vtype);

            switch(vtype)
            {
                case 'i':
                    returnString += *(int*)getIntegerPtr(offset, vtype);
                    break;
                case 'f':
                    returnString += *(float*)getFloatPtr(offset);
                    break;
                case 'Z':
                case 'B':
                    returnString += *(String*)getStringPtr(offset);
                    break;
                default:
                    myStatus.setStatus(SamStatus::INVALID, 
                                       "rmTag called with unknown type.\n");
                    returnStatus = false;
                    break;
            };
        }
        // Increment to the next tag.
        if((currentTagPtr[4] == ';') || (currentTagPtr[4] == ','))
        {
            // Increment once more.
            currentTagPtr += 5;
        }
        else if(currentTagPtr[4] != '\0')
        {
            // Invalid tag format. 
            myStatus.setStatus(SamStatus::INVALID, 
                               "rmTags called with improperly formatted tags.\n");
            returnStatus = false;
            break;
        }
        else
        {
            // Last Tag.
            currentTagPtr += 4;
        }
    }
    return(returnStatus);
}


const String* SamRecord::getStringTag(const char * tag)
{
    // Parse the buffer if necessary.
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.  setTagsFromBuffer set the errors,
            // so just return null.
            return(NULL);
        }
    }
    
    int key = MAKEKEY(tag[0], tag[1], 'Z');
    int offset = extras.Find(key);

    int value;
    if (offset < 0)
    {
        // Check for 'B' tag.
        key = MAKEKEY(tag[0], tag[1], 'B');
        offset = extras.Find(key);
        if(offset < 0)
        {
            // Tag not found.
            return(NULL);
        }
    }

    // Offset is valid, so return the tag.
    value = extras[offset];
    return(&(strings[value]));
}


int* SamRecord::getIntegerTag(const char * tag)
{
    // Init to success.
    myStatus = SamStatus::SUCCESS;
    // Parse the buffer if necessary.
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.  setTagsFromBuffer set the errors,
            // so just return NULL.
            return(NULL);
        }
    }
    
    int key = MAKEKEY(tag[0], tag[1], 'i');
    int offset = extras.Find(key);

    int value;
    if (offset < 0)
    {
        // Failed to find the tag.
        return(NULL);
    }
    else
        value = extras[offset];

    return(&(integers[value]));
}


bool SamRecord::getIntegerTag(const char * tag, int& tagVal)
{
    // Init to success.
    myStatus = SamStatus::SUCCESS;
    // Parse the buffer if necessary.
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.  setTagsFromBuffer set the errors,
            // so just return false.
            return(false);
        }
    }
    
    int key = MAKEKEY(tag[0], tag[1], 'i');
    int offset = extras.Find(key);

    int value;
    if (offset < 0)
    {
        // Failed to find the tag.
        return(false);
    }
    else
        value = extras[offset];

    tagVal = integers[value];
    return(true);
}


bool SamRecord::getFloatTag(const char * tag, float& tagVal)
{
    // Init to success.
    myStatus = SamStatus::SUCCESS;
    // Parse the buffer if necessary.
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.  setTagsFromBuffer set the errors,
            // so just return false.
            return(false);
        }
    }
    
    int key = MAKEKEY(tag[0], tag[1], 'f');
    int offset = extras.Find(key);

    int value;
    if (offset < 0)
    {
        // Failed to find the tag.
        return(false);
    }
    else
        value = extras[offset];

    tagVal = floats[value];
    return(true);
}


const String & SamRecord::getString(const char * tag)
{
    // Init to success.
    myStatus = SamStatus::SUCCESS;
    // Parse the buffer if necessary.
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.
            // TODO - what do we want to do on failure?            
        }
    }
    
    int key = MAKEKEY(tag[0], tag[1], 'Z');
    int offset = extras.Find(key);

    int value;
    if (offset < 0)
    {
    
        key = MAKEKEY(tag[0], tag[1], 'B');
        offset = extras.Find(key);
        if (offset < 0)
        {
            // TODO - what do we want to do on failure?
            return(NOT_FOUND_TAG_STRING);
        }
    }
    value = extras[offset];

    return strings[value];
}


int & SamRecord::getInteger(const char * tag)
{
    // Init to success.
    myStatus = SamStatus::SUCCESS;
    // Parse the buffer if necessary.
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.  setTagsFromBuffer set the error.
            // TODO - what do we want to do on failure?
        }
    }
    
    int key = MAKEKEY(tag[0], tag[1], 'i');
    int offset = extras.Find(key);
    
    int value;
    if (offset < 0)
    {
        // TODO - what do we want to do on failure?
        return NOT_FOUND_TAG_INT;
    }
    else
        value = extras[offset];
    
    return integers[value];
}


bool SamRecord::checkTag(const char * tag, char type)
{
    // Init to success.
    myStatus = SamStatus::SUCCESS;
    // Parse the buffer if necessary.
    if(myNeedToSetTagsFromBuffer)
    {
        if(!setTagsFromBuffer())
        {
            // Failed to read the tags from the buffer, so cannot
            // get tags.  setTagsFromBuffer set the error.
            return("");
        }
    }
    
    int key = MAKEKEY(tag[0], tag[1], type);

    return (extras.Find(key) != LH_NOTFOUND);
}


// Return the error after a failed SamRecord call.
const SamStatus& SamRecord::getStatus()
{
    return(myStatus);
}


// Allocate space for the record - does a realloc.  
// The passed in size is the size of the entire record including the
// block size field.
bool SamRecord::allocateRecordStructure(int size)
{
    if (allocatedSize < size)
    {
        bamRecordStruct* tmpRecordPtr = 
            (bamRecordStruct *)realloc(myRecordPtr, size);
        if(tmpRecordPtr == NULL)
        {
            // FAILED to allocate memory
            fprintf(stderr, "FAILED TO ALLOCATE MEMORY!!!");
            myStatus.addError(SamStatus::FAIL_MEM, "Failed Memory Allocation.");
            return(false);
        }
        // Successfully allocated memory, so set myRecordPtr.
        myRecordPtr = tmpRecordPtr;

        // Reset the pointers into the record.
        if(myIsSequenceBufferValid)
        {
            myPackedSequence = (unsigned char *)myRecordPtr->myData + 
                myRecordPtr->myReadNameLength +
                myRecordPtr->myCigarLength * sizeof(int);
        }
        if(myIsQualityBufferValid)
        {
            myPackedQuality = (unsigned char *)myRecordPtr->myData + 
                myRecordPtr->myReadNameLength +
                myRecordPtr->myCigarLength * sizeof(int) + 
                (myRecordPtr->myReadLength + 1) / 2;
        }

        allocatedSize = size;
    }
    return(true);
}


// Index is the index into the strings array.
void* SamRecord::getStringPtr(int index)
{
    int value = extras[index];

    return &(strings[value]);
}

void* SamRecord::getIntegerPtr(int offset, char& type)
{
    int value = extras[offset];

    type = intType[value];

    return &(integers[value]);
}

void* SamRecord::getFloatPtr(int offset)
{
    int value = extras[offset];

    return &(floats[value]);
}


// Fixes the buffer to match the variable length fields if they are set.
bool SamRecord::fixBuffer(SequenceTranslation translation)
{
    // Check to see if the buffer is already synced.
    if(myIsBufferSynced &&
       (myBufferSequenceTranslation == translation))
    {
        // Already synced, nothing to do.
        return(true);
    }
   
    // Set the bin if necessary.
    if(!myIsBinValid)
    {
        // The bin that is set in the record is not valid, so
        // reset it.
        myRecordPtr->myBin = 
            bam_reg2bin(myRecordPtr->myPosition, get1BasedAlignmentEnd());      
        myIsBinValid = true;
    }

    // Not synced.
    bool status = true;

    // First determine the size the buffer needs to be.
    uint8_t newReadNameLen = getReadNameLength();
    uint16_t newCigarLen = getCigarLength();
    int32_t newReadLen = getReadLength();
    uint32_t newTagLen = getTagLength();
    uint32_t bamSequenceLen = (newReadLen+1)/2;

    // The buffer size extends from the start of the record to data
    // plus the length of the variable fields,
    // Multiply the cigar length by 4 since it is the number of uint32_t fields.
    int newBufferSize = 
        ((unsigned char*)(&(myRecordPtr->myData)) - 
         (unsigned char*)myRecordPtr) +
        newReadNameLen + ((newCigarLen)*4) +
        newReadLen + bamSequenceLen + newTagLen;
   
    if(!allocateRecordStructure(newBufferSize))
    {
        // Failed to allocate space.
        return(false);
    }

    // Now that space has been added to the buffer, check to see what if
    // any fields need to be extracted from the buffer prior to starting to
    // overwrite it.  Fields need to be extracted from the buffer if the 
    // buffer is valid for the field and a previous variable length field has
    // changed length.
    bool readNameLenChange = (newReadNameLen != myRecordPtr->myReadNameLength);
    bool cigarLenChange = (newCigarLen != myRecordPtr->myCigarLength);
    bool readLenChange = (newReadLen != myRecordPtr->myReadLength);

    // If the tags are still stored in the buffer and any other fields changed
    // lengths, they need to be extracted.
    if(myIsTagsBufferValid &&
       (readNameLenChange | cigarLenChange | readLenChange))
    {
        status &= setTagsFromBuffer();
        // The tag buffer will not be valid once the other fields
        // are written, so set it to not valid.
        myIsTagsBufferValid = false;
    }

    // If the sequence or quality strings are still stored in the buffer, and
    // any of the previous fields have changed length, extract it from the
    // current buffer.
    if((myIsQualityBufferValid | myIsSequenceBufferValid) && 
       (readNameLenChange | cigarLenChange | readLenChange))
    {
        setSequenceAndQualityFromBuffer();
        // The quality and sequence buffers will not be valid once the other
        // fields are written, so set them to not valid.
        myIsQualityBufferValid = false;
        myIsSequenceBufferValid = false;
    }

    // If the cigar is still stored in the buffer, and any of the
    // previous fields have changed length, extract it from the current buffer.
    if((myIsCigarBufferValid) && 
       (readNameLenChange))
    {
        status &= parseCigarBinary();
        myIsCigarBufferValid = false;
    }

    // Set each value in the buffer if it is not already valid.
    if(!myIsReadNameBufferValid)
    {
        memcpy(&(myRecordPtr->myData), myReadName.c_str(), 
               newReadNameLen);
   
        // Set the new ReadNameLength.
        myRecordPtr->myReadNameLength = newReadNameLen;
        myIsReadNameBufferValid = true;
    }

    unsigned char * readNameEnds = (unsigned char*)(&(myRecordPtr->myData)) +
        myRecordPtr->myReadNameLength;
   
    // Set the Cigar.  Need to reformat from the string to 
    unsigned int * packedCigar = (unsigned int *) (void *) readNameEnds;
      
    if(!myIsCigarBufferValid)
    {
        // The cigar was already parsed when it was set, so just copy
        // data from the temporary buffer.
        myRecordPtr->myCigarLength = newCigarLen;
        memcpy(packedCigar, myCigarTempBuffer, 
               myRecordPtr->myCigarLength * sizeof(uint32_t));
      
        myIsCigarBufferValid = true;
    }

    unsigned char * packedSequence = readNameEnds + 
        myRecordPtr->myCigarLength * sizeof(int);
    unsigned char * packedQuality = packedSequence + bamSequenceLen;
   
    if(!myIsSequenceBufferValid || !myIsQualityBufferValid || 
       (myBufferSequenceTranslation != translation))
    {
        myRecordPtr->myReadLength = newReadLen;
        // Determine if the quality needs to be set and is just a * and needs to
        // be set to 0xFF.
        bool noQuality = false;
        if((myQuality.Length() == 1) && (myQuality[0] == '*'))
        {
            noQuality = true;
        }
      
        const char* translatedSeq = NULL;
        // If the sequence is not valid in the buffer or it is not
        // properly translated, get the properly translated sequence
        // that needs to be put into the buffer.
        if((!myIsSequenceBufferValid) ||
           (translation != myBufferSequenceTranslation))
        {
            translatedSeq = getSequence(translation);
        }

        for (int i = 0; i < myRecordPtr->myReadLength; i++) 
        {
            if((!myIsSequenceBufferValid) ||
               (translation != myBufferSequenceTranslation))
            {
                // Sequence buffer is not valid, so set the sequence.
                int seqVal = 0;
                switch(translatedSeq[i])
                {
                    case '=':
                        seqVal = 0;
                        break;
                    case 'A':
                    case 'a':
                        seqVal = 1;
                        break;
                    case 'C':
                    case 'c':
                        seqVal = 2;
                        break;
                    case 'G':
                    case 'g':
                        seqVal = 4;
                        break;
                    case 'T':
                    case 't':
                        seqVal = 8;
                        break;
                    case 'N':
                    case 'n':
                    case '.':
                        seqVal = 15;
                        break;
                    default:
                        myStatus.addError(SamStatus::FAIL_PARSE,
                                          "Unknown Sequence character found.");
                        status = false;
                        break;
                };
            
                if(i & 1)
                {
                    // Odd number i's go in the lower 4 bits, so OR in the
                    // lower bits
                    packedSequence[i/2] |= seqVal;
                }
                else
                {
                    // Even i's go in the upper 4 bits and are always set first.
                    packedSequence[i/2] = seqVal << 4;
                }
            }

            if(!myIsQualityBufferValid)
            {
                // Set the quality.
                if((noQuality) || (myQuality.Length() <= i))
                {
                    // No quality or the quality is smaller than the sequence,
                    // so set it to 0xFF
                    packedQuality[i] = 0xFF;
                }
                else
                {
                    // Copy the quality string.
                    packedQuality[i] = myQuality[i] - 33;
                }
            }
        }
        myPackedSequence = (unsigned char *)myRecordPtr->myData + 
            myRecordPtr->myReadNameLength + 
            myRecordPtr->myCigarLength * sizeof(int);
        myPackedQuality = myPackedSequence + 
            (myRecordPtr->myReadLength + 1) / 2;
        myIsSequenceBufferValid = true;
        myIsQualityBufferValid = true;
        myBufferSequenceTranslation = translation;
    }

    if(!myIsTagsBufferValid)
    {
        status &= setTagsInBuffer();
    }

    // Set the lengths in the buffer.
    myRecordPtr->myReadNameLength = newReadNameLen;
    myRecordPtr->myCigarLength = newCigarLen;
    myRecordPtr->myReadLength = newReadLen;

    // Set the buffer block size to the size of the buffer minus the
    // first field.
    myRecordPtr->myBlockSize = newBufferSize - sizeof(int32_t);

    if(status)
    {
        myIsBufferSynced = true;
    }

    return(status);
}


// Sets the Sequence and Quality strings from the buffer.
// They are done together in one method because they require the same
// loop, so might as well be done at the same time.
void SamRecord::setSequenceAndQualityFromBuffer()
{
    // NOTE: If the sequence buffer is not valid, do not set the sequence
    // string from the buffer.
    // NOTE: If the quality buffer is not valid, do not set the quality string
    // from the buffer.

    // Extract the sequence if the buffer is valid and the string's length is 0.
    bool extractSequence = false;
    if(myIsSequenceBufferValid && (mySequence.Length() == 0))
    {
        extractSequence = true;
    }

    // Extract the quality if the buffer is valid and the string's length is 0.
    bool extractQuality = false;
    if(myIsQualityBufferValid && (myQuality.Length() == 0))
    {
        extractQuality = true;
    }

    // If neither the quality nor the sequence need to be extracted,
    // just return.
    if(!extractSequence && !extractQuality)
    {
        return;
    }

    // Set the sequence and quality strings..
    if(extractSequence)
    {
        mySequence.SetLength(myRecordPtr->myReadLength);
    }
    if(extractQuality)
    {
        myQuality.SetLength(myRecordPtr->myReadLength);
    }
   
    const char * asciiBases = "=AC.G...T......N";

    // Flag to see if the quality is specified - the quality contains a value
    // other than 0xFF.  If all values are 0xFF, then there is no quality.
    bool qualitySpecified = false;

    for (int i = 0; i < myRecordPtr->myReadLength; i++)
    {
        if(extractSequence)
        {
            mySequence[i] = i & 1 ?
                asciiBases[myPackedSequence[i / 2] & 0xF] :
                asciiBases[myPackedSequence[i / 2] >> 4];
        }

        if(extractQuality)
        {
            if(myPackedQuality[i] != 0xFF)
            {
                // Quality is specified, so mark the flag.
                qualitySpecified = true;
            }

            myQuality[i] = myPackedQuality[i] + 33;
        }
    }

    // If the read length is 0, then set the sequence and quality to '*'
    if(myRecordPtr->myReadLength == 0)
    {
        if(extractSequence)
        {
            mySequence = "*";
        }
        if(extractQuality)
        {
            myQuality = "*";
        }
    }
    else if(extractQuality && !qualitySpecified)
    {
        // No quality was specified, so set it to "*"
        myQuality = "*";
    }
}


// Parse the cigar to calculate the alignment/unclipped end.
bool SamRecord::parseCigar()
{
    // Determine if the cigar string or cigar binary needs to be parsed.
    if(myCigar.Length() == 0)
    {
        // The cigar string is not yet set, so parse the binary.
        return(parseCigarBinary());
    }
    return(parseCigarString());
}

// Parse the cigar to calculate the alignment/unclipped end.
bool SamRecord::parseCigarBinary()
{
    // Only need to parse if the string is not already set.
    // The length of the cigar string is set to zero when the 
    // record is read from a file into the buffer.
    if(myCigar.Length() != 0)
    {
        // Already parsed.
        return(true);
    }

    unsigned char * readNameEnds = 
        (unsigned char *)myRecordPtr->myData + myRecordPtr->myReadNameLength;
   
    unsigned int * packedCigar = (unsigned int *) (void *) readNameEnds;

    myCigarRoller.Set(packedCigar, myRecordPtr->myCigarLength);
    
    myCigarRoller.getCigarString(myCigar);

    myAlignmentLength = myCigarRoller.getExpectedReferenceBaseCount();

    myUnclippedStartOffset = myCigarRoller.getNumBeginClips();
    myUnclippedEndOffset = myCigarRoller.getNumEndClips();

    // if the cigar length is 0, then set the cigar string to "*"
    if(myRecordPtr->myCigarLength == 0)
    {
        myCigar = "*";
        return(true);
    }

    // Copy the cigar into a temporary buffer.
    int newBufferSize = myRecordPtr->myCigarLength * sizeof(uint32_t);
    if(newBufferSize > myCigarTempBufferAllocatedSize)
    {
        uint32_t* tempBufferPtr = 
            (uint32_t*)realloc(myCigarTempBuffer, newBufferSize);
        if(tempBufferPtr == NULL)
        {
            // Failed to allocate memory.
            // Do not parse, just return.
            fprintf(stderr, "FAILED TO ALLOCATE MEMORY!!!");
            myStatus.addError(SamStatus::FAIL_MEM,
                              "Failed to Allocate Memory.");
            return(false);
        }
        myCigarTempBuffer = tempBufferPtr;
        myCigarTempBufferAllocatedSize = newBufferSize;
    }

    memcpy(myCigarTempBuffer, packedCigar, 
           myRecordPtr->myCigarLength * sizeof(uint32_t));

    // Set the length of the temp buffer.
    myCigarTempBufferLength = myRecordPtr->myCigarLength;

    return(true);
}

// Parse the cigar string to calculate the cigar length and alignment end.
bool SamRecord::parseCigarString()
{
    myCigarTempBufferLength = 0;
    if(myCigar == "*")
    {
        // Cigar is empty, so initialize the variables.
        myAlignmentLength = 0;
        myUnclippedStartOffset = 0;
        myUnclippedEndOffset = 0;
        myCigarRoller.clear();
        return(true);
    }

    myCigarRoller.Set(myCigar);
    
    myAlignmentLength = myCigarRoller.getExpectedReferenceBaseCount();
    
    myUnclippedStartOffset = myCigarRoller.getNumBeginClips();
    myUnclippedEndOffset = myCigarRoller.getNumEndClips();

    // Check to see if the Temporary Cigar Buffer is large enough to contain
    // this cigar.  If we make it the size of the length of the cigar string,
    // it will be more than large enough.
    int newBufferSize = myCigar.Length() * sizeof(uint32_t);
    if(newBufferSize > myCigarTempBufferAllocatedSize)
    {
        uint32_t* tempBufferPtr = 
            (uint32_t*)realloc(myCigarTempBuffer, newBufferSize);
        if(tempBufferPtr == NULL)
        {
            // Failed to allocate memory.
            // Do not parse, just return.
            fprintf(stderr, "FAILED TO ALLOCATE MEMORY!!!");
            myStatus.addError(SamStatus::FAIL_MEM,
                              "Failed to Allocate Memory.");
            return(false);
        }
        myCigarTempBuffer = tempBufferPtr;
        myCigarTempBufferAllocatedSize = newBufferSize;
    }

    // Track if there were any errors.
    bool status = true;

    // Track the index into the cigar string that is being parsed.
    char *cigarOp;
    const char* cigarEntryStart = myCigar.c_str();
    int opLen = 0;
    int op = 0;

    unsigned int * packedCigar = myCigarTempBuffer;
    // TODO - maybe one day make a cigar list... or maybe make a 
    // reference cigar string for ease of lookup....
    const char* endCigarString = cigarEntryStart + myCigar.Length();
    while(cigarEntryStart < endCigarString)
    {
        bool validCigarEntry = true;
        // Get the opLen from the string.  cigarOp will then point to 
        // the operation.
        opLen = strtol(cigarEntryStart, &cigarOp, 10);
        // Switch on the type of operation.
        switch(*cigarOp)
        {
            case('M'):
                op = 0;
                break;
            case('I'):
                // Insert into the reference position, so do not increment the
                // reference end position.
                op = 1;
                break;
            case('D'):
                op = 2;
                break;
            case('N'):
                op = 3;
                break;
            case('S'):
                op = 4;
                break;
            case('H'):
                op = 5;
                break;
            case('P'):
                op = 6;
                break;
            default:
                fprintf(stderr, "ERROR parsing cigar\n");
                validCigarEntry = false;
                status = false;
                myStatus.addError(SamStatus::FAIL_PARSE,
                                  "Unknown operation found when parsing the Cigar.");
                break;
        };
        if(validCigarEntry)
        {
            // Increment the cigar length.
            ++myCigarTempBufferLength;
            *packedCigar = (opLen << 4) | op;
            packedCigar++;
        }
        // The next Entry starts at one past the cigar op, so set the start.
        cigarEntryStart = ++cigarOp;
    }

    // Check clipLength to adjust the end position.
    return(status);
}


bool SamRecord::setTagsFromBuffer()
{
    // If the tags do not need to be set from the buffer, return true.
    if(myNeedToSetTagsFromBuffer == false)
    {
        // Already been set from the buffer.
        return(true);
    }

    // Mark false, as they are being set now.
    myNeedToSetTagsFromBuffer = false;

    unsigned char * extraPtr = myPackedQuality + myRecordPtr->myReadLength;

    // Default to success, will be changed to false on failure.
    bool status = true;

    // Clear any previously set tags.
    clearTags();
    while (myRecordPtr->myBlockSize + 4 - 
           (extraPtr - (unsigned char *)myRecordPtr) > 0)
    {
        int key = 0;
        int value = 0;
        void * content = extraPtr + 3;
        int tagBufferSize = 0;

        key = MAKEKEY(extraPtr[0], extraPtr[1], extraPtr[2]);

        // First check if the tag already exists.
        unsigned int location = extras.Find(key);
        int origIndex = 0;
        String* duplicate = NULL;
        String* origTag = NULL;
        if(location != LH_NOTFOUND)
        {
            duplicate = new String;
            origTag = new String;
            origIndex = extras[location];

            *duplicate = (char)(extraPtr[0]);
            *duplicate += (char)(extraPtr[1]);
            *duplicate += ':';

            *origTag = *duplicate;
            *duplicate += (char)(extraPtr[2]);
            *duplicate += ':';
        }

        switch (extraPtr[2])
        {
            case 'A' :
                if(duplicate != NULL)
                {
                    *duplicate += (* (char *) content);
                    *origTag += intType[origIndex];
                    *origTag += ':';
                    appendIntArrayValue(origIndex, *origTag);
                    tagBufferSize -= getNumericTagTypeSize(intType[origIndex]);
                    integers[origIndex] = *(char *)content;
                    intType[origIndex] = extraPtr[2];
                    tagBufferSize += getNumericTagTypeSize(intType[origIndex]);
                }
                else
                {
                    value = integers.Length();
                    integers.Push(* (char *) content);
                    intType.push_back(extraPtr[2]);
                    tagBufferSize += 4;
                }
                extraPtr += 4;
                break;
            case 'c' :
                if(duplicate != NULL)
                {
                    *duplicate += (* (char *) content);
                    *origTag += intType[origIndex];
                    *origTag += ':';
                    appendIntArrayValue(origIndex, *origTag);
                    tagBufferSize -= getNumericTagTypeSize(intType[origIndex]);
                    integers[origIndex] = *(char *)content;
                    intType[origIndex] = extraPtr[2];
                    tagBufferSize += getNumericTagTypeSize(intType[origIndex]);
                }
                else
                {
                    value = integers.Length();
                    integers.Push(* (char *) content);
                    intType.push_back(extraPtr[2]);
                    tagBufferSize += 4;
                }
                extraPtr += 4;
                break;
            case 'C' :
                if(duplicate != NULL)
                {
                    *duplicate += (* (unsigned char *) content);
                    *origTag += intType[origIndex];
                    *origTag += ':';
                    appendIntArrayValue(origIndex, *origTag);
                    tagBufferSize -= getNumericTagTypeSize(intType[origIndex]);
                    integers[origIndex] = *(unsigned char *)content;
                    intType[origIndex] = extraPtr[2];
                    tagBufferSize += getNumericTagTypeSize(intType[origIndex]);
                }
                else
                {
                    value = integers.Length();
                    integers.Push(* (unsigned char *) content);
                    intType.push_back(extraPtr[2]);
                    tagBufferSize += 4;
                }
                extraPtr += 4;
                break;
            case 's' :
               if(duplicate != NULL)
                {
                    *duplicate += (* (short *) content);
                    *origTag += intType[origIndex];
                    *origTag += ':';
                    appendIntArrayValue(origIndex, *origTag);
                    tagBufferSize -= getNumericTagTypeSize(intType[origIndex]);
                    integers[origIndex] = *(short *)content;
                    intType[origIndex] = extraPtr[2];
                    tagBufferSize += getNumericTagTypeSize(intType[origIndex]);
                }
                else
                {
                    value = integers.Length();
                    integers.Push(* (short *) content);
                    intType.push_back(extraPtr[2]);
                    tagBufferSize += 5;
                }
                extraPtr += 5;
                break;
            case 'S' :
                if(duplicate != NULL)
                {
                    *duplicate += (* (unsigned short *) content);
                    *origTag += intType[origIndex];
                    *origTag += ':';
                    appendIntArrayValue(origIndex, *origTag);
                    tagBufferSize -= getNumericTagTypeSize(intType[origIndex]);
                    integers[origIndex] = *(unsigned short *)content;
                    intType[origIndex] = extraPtr[2];
                    tagBufferSize += getNumericTagTypeSize(intType[origIndex]);
                }
                else
                {
                    value = integers.Length();
                    integers.Push(* (unsigned short *) content);
                    intType.push_back(extraPtr[2]);
                    tagBufferSize += 5;
                }
                extraPtr += 5;
                break;
            case 'i' :
                if(duplicate != NULL)
                {
                    *duplicate += (* (int *) content);
                    *origTag += intType[origIndex];
                    *origTag += ':';
                    appendIntArrayValue(origIndex, *origTag);
                    tagBufferSize -= getNumericTagTypeSize(intType[origIndex]);
                    integers[origIndex] = *(int *)content;
                    intType[origIndex] = extraPtr[2];
                    tagBufferSize += getNumericTagTypeSize(intType[origIndex]);
                }
                else
                {
                    value = integers.Length();
                    integers.Push(* (int *) content);
                    intType.push_back(extraPtr[2]);
                    tagBufferSize += 7;
                }
                extraPtr += 7;
                break;
            case 'I' :
               if(duplicate != NULL)
                {
                    *duplicate += (* (unsigned int *) content);
                    *origTag += intType[origIndex];
                    *origTag += ':';
                    appendIntArrayValue(origIndex, *origTag);
                    tagBufferSize -= getNumericTagTypeSize(intType[origIndex]);
                    integers[origIndex] = *(unsigned int *)content;
                    intType[origIndex] = extraPtr[2];
                    tagBufferSize += getNumericTagTypeSize(intType[origIndex]);
                }
                else
                {
                    value = integers.Length();
                    integers.Push((int) * (unsigned int *) content);
                    intType.push_back(extraPtr[2]);
                    tagBufferSize += 7;
                }
                extraPtr += 7;
                break;
            case 'Z' :
                if(duplicate != NULL)
                {
                    *duplicate += ((const char *) content);
                    *origTag += 'Z';
                    *origTag += ':';
                    *origTag += (char*)(strings[origIndex]);
                    tagBufferSize -= strings[origIndex].Length();
                    strings[origIndex] = (const char *) content;
                    extraPtr += 4 + strings[origIndex].Length();
                    tagBufferSize += strings[origIndex].Length();
                }
                else
                {
                    value = strings.Length();
                    strings.Push((const char *) content);
                    tagBufferSize += 4 + strings.Last().Length();
                    extraPtr += 4 + strings.Last().Length();
                }
                break;
            case 'B' :
                if(duplicate != NULL)
                {
                    *origTag += 'B';
                    *origTag += ':';
                    *origTag += (char*)(strings[origIndex]);
                    tagBufferSize -= 
                        getBtagBufferSize(strings[origIndex]);
                    int bufferSize = 
                        getStringFromBtagBuffer((unsigned char*)content,
                                                strings[origIndex]);
                    *duplicate += (char *)(strings[origIndex]);
                    tagBufferSize += bufferSize;
                    extraPtr += 3 + bufferSize;
                }
                else
                {
                    value = strings.Length();
                    String tempBTag;
                    int bufferSize = 
                        getStringFromBtagBuffer((unsigned char*)content,
                                                tempBTag);
                    strings.Push(tempBTag);
                    tagBufferSize += 3 + bufferSize;
                    extraPtr += 3 + bufferSize;
                }
                break;
            case 'f' :
                if(duplicate != NULL)
                {
                    duplicate->appendFullFloat(* (float *) content);
                    *origTag += 'f';
                    *origTag += ':';
                    origTag->appendFullFloat(floats[origIndex]);
                    floats[origIndex] = *(float *)content;
                }
                else
                {
                    value = floats.size();
                    floats.push_back(* (float *) content);
                    tagBufferSize += 7;
                }
                extraPtr += 7;
                break;
            default :
                fprintf(stderr, 
                        "parsing BAM - Unknown custom field of type %c%c:%c\n",
                        extraPtr[0], extraPtr[1], extraPtr[2]);
                fprintf(stderr, "BAM Tags: \n");

                unsigned char* tagInfo = myPackedQuality + myRecordPtr->myReadLength;

                fprintf(stderr, "\n\n");
                tagInfo = myPackedQuality + myRecordPtr->myReadLength;
                while(myRecordPtr->myBlockSize + 4 - 
                      (tagInfo - (unsigned char *)myRecordPtr) > 0)
                {
                        fprintf(stderr, "%02x",tagInfo[0]);
                        ++tagInfo;
                }
                fprintf(stderr, "\n");

                // Failed on read.
                // Increment extraPtr just by the size of the 3 known fields
                extraPtr += 3;
                myStatus.addError(SamStatus::FAIL_PARSE,
                                  "Unknown tag type.");
                status = false;
        }

        if(duplicate != NULL)
        {
            // Duplicate tag in this record.
            // Tag already existed, print message about overwriting.
            // WARN about dropping duplicate tags.
            if(myNumWarns++ < myMaxWarns)
            {
                fprintf(stderr, "WARNING: Duplicate Tags, overwritting %s with %s\n",
                        origTag->c_str(), duplicate->c_str());
                if(myNumWarns == myMaxWarns)
                {
                    fprintf(stderr, "Suppressing rest of Duplicate Tag warnings.\n");
                }
            }

            continue;
        }

        // Only add the tag if it has so far been successfully processed.
        if(status)
        {
            // Add the tag.
            extras.Add(key, value);
            myTagBufferSize += tagBufferSize;
        }
    }
    return(status);
}


bool SamRecord::setTagsInBuffer()
{
    // The buffer size extends from the start of the record to data
    // plus the length of the variable fields,
    // Multiply the cigar length by 4 since it is the number of uint32_t fields.
    int bamSequenceLength = (myRecordPtr->myReadLength+1)/2;
    int newBufferSize = ((unsigned char*)(&(myRecordPtr->myData)) - 
                         (unsigned char*)myRecordPtr) +  
        myRecordPtr->myReadNameLength + ((myRecordPtr->myCigarLength)*4) +
        myRecordPtr->myReadLength + bamSequenceLength + myTagBufferSize;

    // Make sure the buffer is big enough.
    if(!allocateRecordStructure(newBufferSize))
    {
        // Failed to allocate space.
        return(false);
    }

    char * extraPtr = (char*)myPackedQuality + myRecordPtr->myReadLength;

    bool status = true;

    // Set the tags in the buffer.
    if (extras.Entries())
    {
        for (int i = 0; i < extras.Capacity(); i++)
        {
            if (extras.SlotInUse(i))
            {
                int key = extras.GetKey(i);
                getTag(key, extraPtr);
                extraPtr += 2;
                char vtype;
                getTypeFromKey(key, vtype);
 
                if(vtype == 'i')
                {
                    vtype = getIntegerType(i);
                }

                extraPtr[0] = vtype;

                // increment the pointer to where the value is.
                extraPtr += 1;

                switch (vtype)
                {
                    case 'A' :
                        *(char*)extraPtr = (char)getInteger(i);
                        // sprintf(extraPtr, "%d", getInteger(i));
                        extraPtr += 1;
                        break;
                    case 'c' :
                        *(int8_t*)extraPtr = (int8_t)getInteger(i);
                        // sprintf(extraPtr, "%.4d", getInteger(i));
                        extraPtr += 1;
                        break;
                    case 'C' :
                        *(uint8_t*)extraPtr = (uint8_t)getInteger(i);
                        // sprintf(extraPtr, "%.4d", getInteger(i));
                        extraPtr += 1;
                        break;
                    case 's' :
                        *(int16_t*)extraPtr = (int16_t)getInteger(i);
                        // sprintf(extraPtr, "%.4d", getInteger(i));
                        extraPtr += 2;
                        break;
                    case 'S' :
                        *(uint16_t*)extraPtr = (uint16_t)getInteger(i);
                        // sprintf(extraPtr, "%.4d", getInteger(i));
                        extraPtr += 2;
                        break;
                    case 'i' :
                        *(int32_t*)extraPtr = (int32_t)getInteger(i);
                        // sprintf(extraPtr, "%.4d", getInteger(i));
                        extraPtr += 4;
                        break;
                    case 'I' :
                        *(uint32_t*)extraPtr = (uint32_t)getInteger(i);
                        // sprintf(extraPtr, "%.4d", getInteger(i));
                        extraPtr += 4;
                        break;
                    case 'Z' :
                        sprintf(extraPtr, "%s", getString(i).c_str());
                        extraPtr += getString(i).Length() + 1;
                        break;
                    case 'B' :
                        extraPtr += setBtagBuffer(getString(i), extraPtr);
                        //--TODO-- Set buffer with correct B tag
                        //sprintf(extraPtr, "%s", getString(i).c_str());
                        // extraPtr += getBtagBufferSize(getString(i));
                        break;
                    case 'f' :
                        *(float*)extraPtr = getFloat(i);
                        extraPtr += 4;
                        break;
                    default :
                        myStatus.addError(SamStatus::FAIL_PARSE,
                                          "Unknown tag type.");
                        status = false;
                        break;
                }
            }
        }
    }

    // Validate that the extra pointer is at the end of the allocated buffer.
    // If not then there was a problem.
    if(extraPtr != (char*)myRecordPtr + newBufferSize)
    {
        fprintf(stderr, "ERROR updating the buffer.  Incorrect size.");
        myStatus.addError(SamStatus::FAIL_PARSE,
                          "ERROR updating the buffer.  Incorrect size.");
        status = false;
    }


    // The buffer tags are now in sync.
    myNeedToSetTagsInBuffer = false;
    myIsTagsBufferValid = true;
    return(status);
}


// Reset the variables for a newly set buffer.  The buffer must be set first
// since this looks up the reference ids in the buffer to set the reference
// names.
void SamRecord::setVariablesForNewBuffer(SamFileHeader& header)
{
    // Lookup the reference name & mate reference name associated with this
    // record.
    myReferenceName = 
        header.getReferenceLabel(myRecordPtr->myReferenceID);
    myMateReferenceName = 
        header.getReferenceLabel(myRecordPtr->myMateReferenceID);      

    // Clear the SAM Strings that are now not in-sync with the buffer.
    myReadName.SetLength(0);
    myCigar.SetLength(0);
    mySequence.SetLength(0);
    mySeqWithEq.clear();
    mySeqWithoutEq.clear();
    myQuality.SetLength(0);
    myNeedToSetTagsFromBuffer = true;
    myNeedToSetTagsInBuffer = false;

    //Set that the buffer is valid.
    myIsBufferSynced = true;
    // Set that the variable length buffer fields are valid.
    myIsReadNameBufferValid = true;
    myIsCigarBufferValid = true;
    myPackedSequence = (unsigned char *)myRecordPtr->myData + 
        myRecordPtr->myReadNameLength + myRecordPtr->myCigarLength * sizeof(int);
    myIsSequenceBufferValid = true;
    myBufferSequenceTranslation = NONE;
    myPackedQuality = myPackedSequence + 
        (myRecordPtr->myReadLength + 1) / 2;
    myIsQualityBufferValid = true;
    myIsTagsBufferValid = true;
    myIsBinValid = true;
}


// Extract the vtype from the key.
void SamRecord::getTypeFromKey(int key, char& type) const
{
    // Extract the type from the key.
    type = (key >> 16) & 0xFF;
}


// Extract the tag from the key.
void SamRecord::getTag(int key, char* tag) const
{
    // Extract the tag from the key.
    tag[0] = key & 0xFF;
    tag[1] = (key >> 8) & 0xFF;
    tag[2] = 0;
}


// Index is the index into the strings array.
String & SamRecord::getString(int index)
{
    int value = extras[index];

    return strings[value];
}

int & SamRecord::getInteger(int offset)
{
    int value = extras[offset];

    return integers[value];
}

const char & SamRecord::getIntegerType(int offset) const
{
    int value = extras[offset];

    return intType[value];
}

float & SamRecord::getFloat(int offset)
{
    int value = extras[offset];

    return floats[value];
 }
 

void SamRecord::appendIntArrayValue(char type, int value, String& strVal) const
{
    switch(type)
    {
        case 'A':
            strVal += (char)value;
            break;
        case 'c':
        case 's':
        case 'i':
        case 'C':
        case 'S':
        case 'I':
            strVal += value;
            break;
        default:
            // Do nothing.
            ;
    }
}


int SamRecord::getBtagBufferSize(String& tagStr)
{
    if(tagStr.Length() < 1)
    {
        // ERROR, needs at least the type.
        myStatus.addError(SamStatus::FAIL_PARSE, 
                          "SamRecord::getBtagBufferSize no tag subtype specified");
        return(0);
    }
    char type = tagStr[0];
    int elementSize = getNumericTagTypeSize(type);
    if(elementSize <= 0)
    {
        // ERROR, 'B' tag subtype must be numeric, so should be non-zero
        String errorMsg = "SamRecord::getBtagBufferSize invalid tag subtype, ";
        errorMsg += type;
        myStatus.addError(SamStatus::FAIL_PARSE, errorMsg.c_str());
        return(0);
    }

    // Separated by ',', so count the number of commas.
    int numElements = 0;
    int index = tagStr.FastFindChar(',', 0);
    while(index > 0)
    {
        ++numElements;
        index = tagStr.FastFindChar(',', index+1);
    }
    // Add 5 bytes: type & numElements
    return(numElements * elementSize + 5);
}


int SamRecord::setBtagBuffer(String& tagStr, char* extraPtr)
{
    if(tagStr.Length() < 1)
    {
        // ERROR, needs at least the type.
        myStatus.addError(SamStatus::FAIL_PARSE, 
                          "SamRecord::getBtagBufferSize no tag subtype specified");
        return(0);
    }
    char type = tagStr[0];
    int elementSize = getNumericTagTypeSize(type);
    if(elementSize <= 0)
    {
        // ERROR, 'B' tag subtype must be numeric, so should be non-zero
        String errorMsg = "SamRecord::getBtagBufferSize invalid tag subtype, ";
        errorMsg += type;
        myStatus.addError(SamStatus::FAIL_PARSE, errorMsg.c_str());
        return(0);
    }

    int totalInc = 0;
    // Write the type.
    *(char*)extraPtr = type;
    ++extraPtr;
    ++totalInc;

    // Get the number of elements by counting ','s
    uint32_t numElements = 0;
    int index = tagStr.FastFindChar(',', 0);
    while(index > 0)
    {
        ++numElements;
        index = tagStr.FastFindChar(',', index+1);
    }
    *(uint32_t*)extraPtr = numElements;
    extraPtr += 4;
    totalInc += 4;

    const char* stringPtr = tagStr.c_str();
    const char* endPtr = stringPtr + tagStr.Length();
    // increment past the subtype and ','.
    stringPtr += 2;

    char* newPtr = NULL;
    while(stringPtr < endPtr)
    {
        switch(type)
        {
            case 'f':
                *(float*)extraPtr = (float)(strtod(stringPtr, &newPtr));
                break;
            case 'c':
                *(int8_t*)extraPtr = (int8_t)strtol(stringPtr, &newPtr, 0);
                break;
            case 's':
                *(int16_t*)extraPtr = (int16_t)strtol(stringPtr, &newPtr, 0);
                break;
            case 'i':
                *(int32_t*)extraPtr = (int32_t)strtol(stringPtr, &newPtr, 0);
                break;
            case 'C':
                *(uint8_t*)extraPtr = (uint8_t)strtoul(stringPtr, &newPtr, 0);
                break;
            case 'S':
                *(uint16_t*)extraPtr = (uint16_t)strtoul(stringPtr, &newPtr, 0);
                break;
            case 'I':
                *(uint32_t*)extraPtr = (uint32_t)strtoul(stringPtr, &newPtr, 0);
                break;
            default :
                myStatus.addError(SamStatus::FAIL_PARSE,
                                  "Unknown 'B' tag subtype.");
                break;
        }
        extraPtr += elementSize;
        totalInc += elementSize;
        stringPtr = newPtr + 1; // skip the ','
    }
    return(totalInc);
}


int SamRecord::getStringFromBtagBuffer(unsigned char* buffer, 
                                       String& tagStr)
{
    tagStr.Clear();

    int bufferSize = 0;

    // 1st byte is the type.
    char type = *buffer;
    ++buffer;
    ++bufferSize;
    tagStr = type;

    // 2nd-5th bytes are the size
    unsigned int numEntries = *(unsigned int *)buffer;
    buffer += 4;
    bufferSize += 4;
    // Num Entries is not included in the string.

    int subtypeSize = getNumericTagTypeSize(type);

    for(unsigned int i = 0; i < numEntries; i++)
    {
        tagStr += ',';
        switch(type)
        {
            case 'f':
                tagStr.appendFullFloat(*(float *)buffer);
                break;
            case 'c':
                tagStr += *(int8_t *)buffer;
                break;
            case 's':
                tagStr += *(int16_t *)buffer;
                break;
            case 'i':
                 tagStr += *(int32_t *)buffer;
                break;
            case 'C':
                tagStr += *(uint8_t *)buffer;
                break;
            case 'S':
                tagStr += *(uint16_t *)buffer;
                break;
            case 'I':
                tagStr += *(uint32_t *)buffer;
                break;
            default :
                myStatus.addError(SamStatus::FAIL_PARSE,
                                  "Unknown 'B' tag subtype.");
                break;
        }
        buffer += subtypeSize;
        bufferSize += subtypeSize;
    }
    return(bufferSize);
}
