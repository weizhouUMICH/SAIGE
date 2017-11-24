/*
 *  Copyright (C) 2010-2015  Regents of the University of Michigan
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

#include <iostream>

#include "SamValidation.h"
#include "CigarRoller.h"
#include "SamTags.h"

const char* SamValidationError::enumSeverityString[] = {
    "WARNING", "ERROR"};

const char* SamValidationError::enumTypeString[] = {
    "INVALID_QNAME",
    "INVALID_REF_ID",
    "INVALID_RNAME",
    "INVALID_POS",
    "INVALID_MAPQ",
    "INVALID_CIGAR",
    "INVALID_MRNM",
    "INVALID_QUAL",
    "INVALID_TAG"
};

const char* SamValidationError::getTypeString(Type type)
{
    return(enumTypeString[type]);
}


SamValidationError::SamValidationError(Type type, Severity severity, 
                                       std::string message)
{
    myType = type;
    mySeverity = severity;
    myMessage = message;
}


SamValidationError::Type SamValidationError::getType() const
{
    return(myType);
}


SamValidationError::Severity SamValidationError::getSeverity() const
{
    return(mySeverity);
}


const char* SamValidationError::getMessage() const
{
    return(myMessage.c_str());
}


const char* SamValidationError::getTypeString() const
{
    return(enumTypeString[myType]);
}


const char* SamValidationError::getSeverityString() const
{
    return(enumSeverityString[mySeverity]);
}


void SamValidationError::getErrorString(std::string& errorString) const
{
    errorString = getTypeString();
    errorString += " (";
    errorString += getSeverityString();
    errorString += ") : ";
    errorString += getMessage();
    errorString += "\n";
}


void SamValidationError::printError() const
{
    std::cerr << this;
}



// Constructor.
SamValidationErrors::SamValidationErrors()
    : myValidationErrors()
{
    myErrorIter = myValidationErrors.begin();
}


// Destructor
SamValidationErrors::~SamValidationErrors()
{
    clear();
}


void SamValidationErrors::clear()
{
    // Clear the errors.
    std::list<const SamValidationError*>::iterator errorIter;
    for(errorIter = myValidationErrors.begin(); 
        errorIter != myValidationErrors.end(); ++errorIter)
    {
        delete *errorIter;
        *errorIter = NULL;
    }
    myValidationErrors.clear();
    myErrorIter = myValidationErrors.end();
}


void SamValidationErrors::addError(SamValidationError::Type newType, 
                                   SamValidationError::Severity newSeverity, 
                                   const char* newMessage)
{
    myValidationErrors.push_back(new SamValidationError(newType, 
                                                        newSeverity, 
                                                        newMessage));

    // If this is the first element in the list, set the iterator.
    if(myValidationErrors.size() == 1)
    {
        // set the iterator to the first element.
        myErrorIter = myValidationErrors.begin();
    }
}



// Return the number of validation errors that are contained in this object.
unsigned int SamValidationErrors::numErrors()
{
    return(myValidationErrors.size());
}


// Return a pointer to the next error.  It does not remove it from the list.
// Returns null once all errors have been retrieved until resetErrorIter
// is called.
const SamValidationError* SamValidationErrors::getNextError()
{
    if(myErrorIter == myValidationErrors.end())
    {
        // at the end of the list, return null.
        return(NULL);
    }
    // Not at the end of the list, return the last element and increment.
    return(*myErrorIter++);
}

   
// Resets the iterator to the begining of the errors.
void SamValidationErrors::resetErrorIter()
{
    myErrorIter = myValidationErrors.begin();
}


// Appends the error messages to the passed in string.
void SamValidationErrors::getErrorString(std::string& errorString) const
{
    for(std::list<const SamValidationError*>::
            const_iterator validationErrorIter = 
            myValidationErrors.begin(); 
        validationErrorIter != myValidationErrors.end(); 
        validationErrorIter++)
    {
        std::string error = "";
        (*validationErrorIter)->getErrorString(error);
        errorString += error;
    }
}


bool SamValidator::isValid(SamFileHeader& samHeader, SamRecord& samRecord, 
                           SamValidationErrors& validationErrors)
{
    bool status = true;
    status &= isValidQname(samRecord.getReadName(), 
                           samRecord.getReadNameLength(), 
                           validationErrors);

    status &= isValidFlag(samRecord.getFlag(), 
                          validationErrors);

    // Validate the RName including validating it against the header.
    status &= isValidRname(samHeader,
                           samRecord.getReferenceName(), 
                           validationErrors);

    status &= isValidRefID(samRecord.getReferenceID(), 
                           samHeader.getReferenceInfo(), 
                           validationErrors);

    status &= isValid1BasedPos(samRecord.get1BasedPosition(),
                               validationErrors);
    
    status &= isValidMapQuality(samRecord.getMapQuality(), validationErrors);

    status &= isValidSequence(samRecord, validationErrors);

    status &= isValidCigar(samRecord, validationErrors);
    
    status &= isValidQuality(samRecord, validationErrors);

    status &= isValidTags(samRecord, validationErrors);

    return(status);
}


// qname is the query (read) name - result of SamRecord::getReadName().
// readNameLen is the length of the read name including the null (the result
// of SamRecord::getReadNameLength()).
// For some invalid records, the getReadNameLength may be different than the 
// length of qname.
// NOTE: Query Name and Read Name both refer to the same field.
bool SamValidator::isValidQname(const char* qname, uint8_t readNameLen, 
                                SamValidationErrors& validationErrors)
{
    // Validation for QNAME is:
    //   a) length of the qname string is the same as the read name length
    //   b) length is between 1 and 254.
    //   c) [ \t\n\r] are not allowed in the name.

    bool status = true;

    // Get the length of the qname string.
    int32_t qnameLenNull = strlen(qname) + 1;

    ////////////////////////////////////
    //   a) length of the qname string is the same as the read name length
    if(qnameLenNull != readNameLen)
    {
        // This results from a poorly formatted bam file, where the null
        // terminated read_name field is not the same length as specified by 
        // read_name_len.
        String message = "Invalid Query Name - the string length (";
        message += qnameLenNull;
        message += ") does not match the specified query name length (";
        message += readNameLen;
        message += ").";

        validationErrors.addError(SamValidationError::INVALID_QNAME,
                                  SamValidationError::ERROR, 
                                  message.c_str());
        status = false;
    }

    ////////////////////////////////////
    //   b) length is between 1 and 254
    // The length with the terminating null must be between 2 & 255,
    if((qnameLenNull < 2) || (qnameLenNull > 255))
    {
        String message = "Invalid Query Name (QNAME) length: ";
        message += qnameLenNull;
        message += ".  Length with the terminating null must be between 2 & 255.";
      
        validationErrors.addError(SamValidationError::INVALID_QNAME,
                                  SamValidationError::WARNING, 
                                  message.c_str());
        status = false;
    }

    ////////////////////////////////////
    // Loop through and validate they all characters are valid.
    //   c) [ \t\n\r] are not allowed in the name.
    String message;
    for(int i = 0; i < qnameLenNull; ++i)
    {
        switch(qname[i])
        {
            case ' ':
                // Invalid character.
                message = "Invalid character in the Query Name (QNAME): ' ' at position ";
                message += i;
                message += ".";
                validationErrors.addError(SamValidationError::INVALID_QNAME,
                                          SamValidationError::WARNING, 
                                          message.c_str());
                status = false;
                break;
            case '\t':
                // Invalid character.
                message = "Invalid character in the Query Name (QNAME): '\t' at position ";
                message += i;
                message += ".";
                validationErrors.addError(SamValidationError::INVALID_QNAME,
                                          SamValidationError::WARNING, 
                                          message.c_str());
                status = false;
                break;
            case '\n':
                // Invalid character.
                message = "Invalid character in the Query Name (QNAME): '\n' at position ";
                message += i;
                message += ".";
                validationErrors.addError(SamValidationError::INVALID_QNAME,
                                          SamValidationError::WARNING, 
                                          message.c_str());
                status = false;
                break;
            case '\r':
                // Invalid character.
                message = "Invalid character in the Query Name (QNAME): '\r' at position ";
                message += i;
                message += ".";
                validationErrors.addError(SamValidationError::INVALID_QNAME,
                                          SamValidationError::WARNING, 
                                          message.c_str());
                status = false;
                break;
        }
    }

    return(status);
}


bool SamValidator::isValidFlag(uint16_t flag,
                               SamValidationErrors& validationErrors)
{
    // All values in a uint16_t are valid, so return true.
    return(true);
}


bool SamValidator::isValidRname(SamFileHeader& samHeader,
                                const char* rname,
                                SamValidationErrors& validationErrors)
{
    bool status = true;

    // Cross validate the rname and the header.
    // If the rname is not '*'
    // AND there are any SQ records in the header,
    // Then the rname must be in one of them.
    if((strcmp(rname, "*") != 0) &&
       (samHeader.getNumSQs() != 0) && 
       (samHeader.getSQ(rname) == NULL))
    {
        // There are SQ fields, but the ref name is not in it.
        status = false;
        std::string message = "RNAME, ";
        message += rname;
        message += ", was not found in a SAM Header SQ record";
        validationErrors.addError(SamValidationError::INVALID_RNAME,
                                  SamValidationError::WARNING,
                                  message.c_str());
    }
    status &= isValidRname(rname, validationErrors);
    return(status);
}


bool SamValidator::isValidRname(const char* rname,
                                SamValidationErrors& validationErrors)
{
    // Validation for RNAME is:
    //   a) cannot be 0 length.
    //   b) [ \t\n\r@=] are not allowed in the name.

    bool status = true;

    // Get the length of the rname string.
    int32_t rnameLen = strlen(rname);

    String message;

    if(rnameLen == 0)
    {
        validationErrors.addError(SamValidationError::INVALID_RNAME,
                                  SamValidationError::WARNING, 
                                  "Reference Sequence Name (RNAME) cannot have 0 length.");
        status = false;
    }

    ////////////////////////////////////
    ////////////////////////////////////
    // Loop through and validate they all characters are valid.
    //   b) [ \t\n\r] are not allowed in the name.
    for(int i = 0; i < rnameLen; ++i)
    {
        switch(rname[i])
        {
            case ' ':
                // Invalid character.
                message = "Invalid character in the Reference Sequence Name (RNAME): ' ' at position ";
                message += i;
                message += ".";
                validationErrors.addError(SamValidationError::INVALID_RNAME,
                                          SamValidationError::WARNING, 
                                          message.c_str());
                status = false;
                break;
            case '\t':
                // Invalid character.
                message = "Invalid character in the Reference Sequence Name (RNAME): '\t' at position ";
                message += i;
                message += ".";
                validationErrors.addError(SamValidationError::INVALID_RNAME,
                                          SamValidationError::WARNING, 
                                          message.c_str());
                status = false;
                break;
            case '\n':
                // Invalid character.
                message = "Invalid character in the Reference Sequence Name (RNAME): '\n' at position ";
                message += i;
                message += ".";
                validationErrors.addError(SamValidationError::INVALID_RNAME,
                                          SamValidationError::WARNING, 
                                          message.c_str());
                status = false;
                break;
            case '\r':
                // Invalid character.
                message = "Invalid character in the Reference Sequence Name (RNAME): '\r' at position ";
                message += i;
                message += ".";
                validationErrors.addError(SamValidationError::INVALID_RNAME,
                                          SamValidationError::WARNING, 
                                          message.c_str());
                status = false;
                break;
            case '@':
                // Invalid character.
                message = "Invalid character in the Reference Sequence Name (RNAME): '@' at position ";
                message += i;
                message += ".";
                validationErrors.addError(SamValidationError::INVALID_RNAME,
                                          SamValidationError::WARNING, 
                                          message.c_str());
                status = false;
                break;
            case '=':
                // Invalid character.
                message = "Invalid character in the Reference Sequence Name (RNAME): '=' at position ";
                message += i;
                message += ".";
                validationErrors.addError(SamValidationError::INVALID_RNAME,
                                          SamValidationError::WARNING, 
                                          message.c_str());
                status = false;
                break;
            default:
                // Allowed character.
                break;
        }
    }

    return(status);
}


bool SamValidator::isValidRefID(int32_t refID, 
                                const SamReferenceInfo& refInfo,
                                SamValidationErrors& validationErrors)
{
    // Validation for rID is:
    //   a) must be between -1 and the number of refInfo.
    //      -1 is allowed, and otherwise it must properly index into the array.

    bool status = true;
    if((refID < -1) || (refID >= refInfo.getNumEntries()))
    {
        // Reference ID is too large or too small.
        String message = "Invalid Reference ID, out of range (";
        message += refID;
        message += ") must be between -1 and ";
        message += refInfo.getNumEntries() - 1;
        message += ".";

        validationErrors.addError(SamValidationError::INVALID_REF_ID,
                                  SamValidationError::WARNING, 
                                  message.c_str());
        status = false;
    }

    return(status);
}


bool SamValidator::isValid1BasedPos(int32_t pos, 
                                    SamValidationErrors& validationErrors)
{
    // Validation for pos is:
    //   a) must be between 0 and (2^29)-1.

    bool status = true;

    if((pos < 0) || (pos > 536870911))
    {
        String message = "POS out of range (";
        message += pos;
        message += ") must be between 0 and (2^29)-1.";

        validationErrors.addError(SamValidationError::INVALID_POS,
                                  SamValidationError::WARNING, 
                                  message.c_str());
        status = false;
    }

    return(status);
}


bool SamValidator::isValidMapQuality(uint8_t mapQuality,
                                     SamValidationErrors& validationErrors)
{
    // All values in a uint8_t are valid, so return true.
    return(true);
}


bool SamValidator::isValidSequence(SamRecord& samRecord,
                                   SamValidationErrors& validationErrors)
{
    return(true);
}


bool SamValidator::isValidCigar(SamRecord& samRecord,
                                SamValidationErrors& validationErrors)
{
    return(isValidCigar(samRecord.getCigar(), 
                        samRecord.getReadLength(),
                        validationErrors));
}

bool SamValidator::isValidCigar(const char* cigar,
                                const char* sequence,
                                SamValidationErrors& validationErrors)
{
    if(strcmp(sequence, "*") != 0)
    {
        return(isValidCigar(cigar, strlen(sequence), validationErrors));
    }
    // Sequence is '*', so the length is 0.
    return(isValidCigar(cigar, 0, validationErrors));
}

bool SamValidator::isValidCigar(const char* cigar,
                                int seqLen,
                                SamValidationErrors& validationErrors)
{
    // Validation for CIGAR is:
    //   a) cannot be 0 length.
    // if not "*", validate the following:
    //   b) must have an integer length for each operator (if not "*"). TODO
    //   c) all operators must be valid (if not "*"). TODO
    //   d) evaluates to the same read length as the sequence string if not '*'.
    bool status = true;
    String message;

    int32_t cigarLen = strlen(cigar);

    //   a) cannot be 0 length.
    if(cigarLen == 0)
    {
        validationErrors.addError(SamValidationError::INVALID_CIGAR,
                                  SamValidationError::WARNING,
                                  "Cigar must not be blank.");
        status = false;
    }
    
    if(strcmp(cigar, "*") != 0)
    {
        // The cigar is not "*", so validate it.
        CigarRoller cigarRoller(cigar);
        
        //   b) must have an integer length for each operator.
        // TODO
        //   c) all operators must be valid.
        // TODO

        //   d) is the same length as the sequence string.
        int cigarSeqLen = cigarRoller.getExpectedQueryBaseCount();
        if((cigarSeqLen != seqLen) && (seqLen != 0))
        {
            message = "CIGAR does not evaluate to the same length as SEQ, (";
            message += cigarSeqLen;
            message += " != ";
            message += seqLen;
            message += ").";
            validationErrors.addError(SamValidationError::INVALID_CIGAR,
                                      SamValidationError::WARNING, 
                                      message.c_str());
            status = false;
        }
    }
    return(status);
}


bool SamValidator::isValidQuality(SamRecord& samRecord,
                                  SamValidationErrors& validationErrors)
{
    return(isValidQuality(samRecord.getQuality(), 
                          samRecord.getReadLength(),
                          validationErrors));
}


bool SamValidator::isValidQuality(const char* quality,
                                  const char* sequence,
                                  SamValidationErrors& validationErrors)
{
    // Determine the length of the sequence.
    int seqLen = strlen(sequence);

    // Check if the sequence is '*' since then the seqLength is 0.
    if(strcmp(sequence, "*") == 0)
    {
        seqLen = 0;
    }
    return(isValidQuality(quality, seqLen, validationErrors));
}


bool SamValidator::isValidQuality(const char* quality,
                                  int seqLength,
                                  SamValidationErrors& validationErrors)
{
    bool status = true;

    // If the quality or the sequence are non-"*", validate that the quality
    // and sequence have the same length.
    if((seqLength != 0) && (strcmp(quality, "*") != 0))
    {
        int qualLen = strlen(quality);
        // Both the sequence and the quality are not "*", so validate
        // that they are the same length.
        if(seqLength != qualLen)
        {
            // Both fields are specified but are different lengths.
            
            String message = "QUAL is not the same length as SEQ, (";
            message += qualLen;
            message += " != ";
            message += seqLength;
            message += ").";
            
            validationErrors.addError(SamValidationError::INVALID_QUAL,
                                      SamValidationError::WARNING, 
                                      message.c_str());
        status = false;
        }
    }
    return(status);
}


bool SamValidator::isValidTags(SamRecord& samRecord,
                               SamValidationErrors& validationErrors)
{
    bool status = true;

    GenomeSequence* reference = samRecord.getReference();
    // If the reference is not null, check the MD tag.
    if(reference != NULL)
    {
        const String* recordMD = samRecord.getStringTag(SamTags::MD_TAG);
        if(recordMD != NULL)
        {
            // The record has an MD tag so check to see if it is correct.
            if(!SamTags::isMDTagCorrect(samRecord, *reference))
            {
                // Invalid MD tags.
                String correctMD;
                if(!SamTags::createMDTag(correctMD, samRecord, *reference))
                {
                    // Failed to get the MD tag, so indicate that it is unknown.
                    correctMD = "UNKNOWN";
                }
                String message = "Incorrect MD Tag, ";
                message += *recordMD;
                message += ", should be ";
                message += correctMD;
                message += ".";
                
                validationErrors.addError(SamValidationError::INVALID_TAG,
                                          SamValidationError::WARNING, 
                                          message.c_str());
                
                status = false;
            }
        }
    }

    return(status);
}
