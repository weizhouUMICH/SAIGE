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

#ifndef __SAM_VALIDATION_H__
#define __SAM_VALIDATION_H__

#include "SamFile.h"
#include <list>

// On windows, ERROR and WARNING are pre-defined macros, so undefine them.
#ifdef WARNING
#undef WARNING
#endif
#ifdef ERROR
#undef ERROR
#endif

/// The SamValidationError class describes a validation error that occured,
/// containing the error type, severity, and textual error message.
class SamValidationError
{
public:
    /// Severity of the error.
    enum Severity 
        {
            WARNING, ///< Warning is used if it is just an invalid value.
            ERROR    ///< Error is used if parsing could not succeed.
        };

    /// Type of the error.
    /// TODO: NOT ALL INVALID TYPES HAVE BEEN ADDED SINCE NOT ALL VALIDATION
    /// IS COMPLETE YET
    enum Type
        {
            INVALID_QNAME, ///< Invalid read/query name
            INVALID_REF_ID, ///< Invalid reference id
            INVALID_RNAME, ///< Invalid reference name
            INVALID_POS, ///< Invalid position
            INVALID_MAPQ, ///< Invalid mapping quality
            INVALID_CIGAR, ///< Invalid CIGAR
            INVALID_MRNM, ///< Invalid mate/next fragment reference name
            INVALID_QUAL, ///< Invalid base quality
            INVALID_TAG ///< Invalid tag
        };

    /// Get the string representing the specified type of validation error.
    static const char* getTypeString(Type type);

    /// Constructor that sets the type, severity, and message for the
    /// validation error.
    SamValidationError(Type type, Severity severity, std::string Message);
   
    /// Return the type enum of this validation error object.
    Type getType() const;

    /// Return the severity enum of this validation error object.
    Severity getSeverity() const;

    /// Return the error message of this validation error object.
    const char* getMessage() const;

    /// Return the string representing this object's type of validation error.
    const char* getTypeString() const;

    /// Return the string representing this object's severity of validation
    /// error.
    const char* getSeverityString() const;

    /// Get the error string representing this object's error.
    void getErrorString(std::string& errorString) const;

    /// Print a formatted output of the error to cerr.
    void printError() const;

private:
    SamValidationError();

    static const char* enumTypeString[];
    static const char* enumSeverityString[];

    Type myType;
    Severity mySeverity;
    std::string myMessage;

};


/// stream output for validation failure information
inline std::ostream &operator << (std::ostream &stream, 
                                  const SamValidationError &error)
{
    std::string errorMessage;
    error.getErrorString(errorMessage);
    stream << errorMessage;
    return stream;
}


/// The SamValidationErrors class is a container class that holds
/// SamValidationError Objects, allowing a validation method to return all
/// of the invalid errors rather than just one.
class SamValidationErrors
{
public:
    /// Constructor.
    SamValidationErrors();
    /// Destructor
    ~SamValidationErrors();

    /// Remove all the errors from the container.
    void clear();

    /// Add the specified error to this container.
    void addError(SamValidationError::Type newType, 
                  SamValidationError::Severity newSeverity,
                  const char* newMessage);

    /// Return the number of validation errors contained in this object.
    unsigned int numErrors();

    /// Return a pointer to the next error without removing it from the 
    /// container, and returning null once all errors have been retrieved
    /// until resetErrorIter is called.
    const SamValidationError* getNextError();
   
    /// Reset the iterator to the begining of the errors.
    void resetErrorIter();

    /// Append the error messages contained in this container to the passed
    /// in string.
    void getErrorString(std::string& errorString) const;

private:
    std::list<const SamValidationError*> myValidationErrors;
    std::list<const SamValidationError*>::const_iterator myErrorIter;
};


/// stream output for all validation failures information
inline std::ostream& operator << (std::ostream& stream,
                                  const SamValidationErrors& errors)
{
    std::string errorString = "";
    errors.getErrorString(errorString);
    stream << errorString;
    return stream;
}


/// The SamValidator class contains static methods for validating the SAM/BAM
/// Record and each of its fields. The generic isValid method performs all of
/// the other validations. The SamValidator methods return whether or not what
/// is being validated is valid. True means it is valid, false means it is not.
/// The specifics of the invalid value(s) are contained in the
/// SamValidationErrors object that is passed in (by reference) to the method.
/// The specific errors can be pulled out of that object.
/// TODO: VALIDATION METHODS STILL NEED TO BE ADDED, and isValid does not yet
/// validate all fields!!!
class SamValidator
{
public:

    /// Validates whether or not the specified SamRecord is valid, calling
    /// all of the other validations.
    /// TODO: more validation needs to be added.
    /// \param samHeader header associated with the record to be validated.
    /// \param samRecord record to be validated.
    /// \param validationErrors status  to append any errors too.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValid(SamFileHeader& samHeader, SamRecord& samRecord, 
                        SamValidationErrors& validationErrors);

    /// Determines whether or not the specified qname is valid.
    /// Validation for QNAME is:
    ///   a) length of the qname string is the same as the read name length
    ///   b) length is between 1 and 254.
    ///   c) [ \t\n\r] are not allowed in the name.
    /// \param qname the read/query name.
    /// \param qnameLen length of the read including the null (result of 
    /// SamRecord::getReadNameLength().
    /// \param validationErrors status  to append any errors too.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidQname(const char* qname, uint8_t qnameLen, 
                             SamValidationErrors& validationErrors);

    /// Determines whether or not the flag is valid.
    /// TODO: currently no validation is done on the flag.
    /// \param flag flag to be validated.
    /// \param validationErrors status  to append any errors too.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidFlag(uint16_t flag,
                            SamValidationErrors& validationErrors);

    /// Validate the reference name including validating against the header.
    /// 1) Cross validate the rname and the header.
    /// 2) perform the validation in the method that doesn't take the header.
    /// \param samHeader header associated with the rname to be validated.
    /// \param rname reference name to be validated.
    /// \param validationErrors status  to append any errors too.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidRname(SamFileHeader& samHeader, 
                             const char* rname,
                             SamValidationErrors& validationErrors);
    /// Validate the rname without validating against the header.
    /// Validation for RNAME is:
    ///   a) cannot be 0 length.
    ///   b) [ \t\n\r@=] are not allowed in the name.
    /// \param rname reference name to be validated.
    /// \param validationErrors status  to append any errors too.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidRname(const char* rname,
                             SamValidationErrors& validationErrors);

    /// Validate whether or not the specified reference id is valid.
    /// Validation for rID is:
    ///  a) must be between -1 and the number of refInfo.
    ///     -1 is allowed, and otherwise it must properly index into the array.
    /// \param refID reference id to be validated.
    /// \param refInfo sam reference information containing the mapping
    /// from reference id to reference name for this refID.
    /// \param validationErrors status  to append any errors too.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidRefID(int32_t refID, const SamReferenceInfo& refInfo, 
                             SamValidationErrors& validationErrors);

    /// Validate the refeference position.
    /// Validation for pos is:
    ///   a) must be between 0 and (2^29)-1.
    /// \param pos position to be validated.
    /// \param validationErrors status  to append any errors too.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValid1BasedPos(int32_t pos, 
                                 SamValidationErrors& validationErrors);

    /// Validate the mapping quality.
    /// TODO: currently no validation is done on the mapping quality.
    /// \param mapQuality mapping quality to be validated.
    /// \param validationErrors status  to append any errors too.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidMapQuality(uint8_t mapQuality,
                                  SamValidationErrors& validationErrors);

    /// Validate the sequence, but not against the cigar or quality string.
    /// Validation against cigar is done in isValidCigar.
    /// Validation against the quality string is done in isValidQuality.
    /// TODO: currently no validation is done in this method.
    /// \param samRecord record whose sequence should be validated.
    /// \param validationErrors status  to append any errors too.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidSequence(SamRecord& samRecord,
                                SamValidationErrors& validationErrors);

    /// Validate the cigar.  Cigar validation depends on sequence.
    /// Validation for CIGAR is:
    ///   a) cannot be 0 length.
    /// if not "*", validate the following:
    ///   b) must have an integer length for each operator (if not "*"). TODO
    ///   c) all operators must be valid (if not "*"). TODO
    ///   d) evaluates to the same read length as the sequence string.
    /// \param samRecord record whose cigar should be validated.
    /// \param validationErrors status  to append any errors too.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidCigar(SamRecord& samRecord,
                             SamValidationErrors& validationErrors);

    /// Validate the cigar.  Cigar validation depends on sequence.
    /// Validation for CIGAR is:
    ///   a) cannot be 0 length.
    /// if not "*", validate the following:
    ///   b) must have an integer length for each operator (if not "*"). TODO
    ///   c) all operators must be valid (if not "*"). TODO
    ///   d) evaluates to the same read length as the sequence string.
    /// \param cigar cigar string to be validated.
    /// \param sequence sequence to check the cigar against.
    /// \param validationErrors status  to append any errors too.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidCigar(const char* cigar, const char* sequence,
                             SamValidationErrors& validationErrors);

    /// Validate the cigar.  Cigar validation depends on sequence.
    /// Validation for CIGAR is:
    ///   a) cannot be 0 length.
    /// if not "*", validate the following:
    ///   b) TODO: must have an integer length for each operator (if not "*").
    ///   c) TODO: all operators must be valid (if not "*").
    ///   d) evaluates to the same read length as the sequence string.
    /// \param cigar cigar string to be validated.
    /// \param seqLen sequence length to check the cigar against.
    /// \param validationErrors status  to append any errors too.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidCigar(const char* cigar,
                             int seqLen,
                             SamValidationErrors& validationErrors);

    /// TODO: validate the mate/next fragment's reference name.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidMrnm();

    /// TODO: validate the mate/next fragment's position.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidMpos();

    /// TODO: validate the insertion size/observed template length.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidIsize();

    /// TODO, validate the sequence.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidSeq();

    /// Validate the base quality.
    /// Quality validation depends on sequence.
    /// Validation for quality is:
    ///   a) quality & sequence are the same length if both are specified.
    /// TODO: more validation.
    /// \param samRecord record whose quality should be validated.
    /// \param validationErrors status  to append any errors too.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidQuality(SamRecord& samRecord,
                               SamValidationErrors& validationErrors);

    /// Validate the base quality.
    /// Quality validation depends on sequence.
    /// Validation for quality is:
    ///   a) quality & sequence are the same length if both are specified.
    /// TODO: more validation.
    /// \param quality quality string to be validated.
    /// \param seqLen sequence length to check the quality against.
    /// \param validationErrors status  to append any errors too.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidQuality(const char* quality, const char* sequence,
                               SamValidationErrors& validationErrors);

    /// Validate the base quality.
    /// Quality validation depends on sequence.
    /// Validation for quality is:
    ///   a) quality & sequence are the same length if both are specified.
    /// TODO: more validation.
    /// \param quality quality string to be validated.
    /// \param seqLen sequence length to check the quality against.
    /// \param validationErrors status  to append any errors too.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    bool static isValidQuality(const char* quality,
                               int seqLength,
                               SamValidationErrors& validationErrors);

    /// Validate the tags.
    /// Validation for tags is:
    ///   a) check that the "MD" tag is correct if it is present.
    /// TODO: more validation.
    /// \param samRecord record whose tags should be validated.
    /// \param validationErrors status  to append any errors too.
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidTags(SamRecord& samRecord,
                            SamValidationErrors& validationErrors);

    /// TODO validate the tag vtype
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidVtype();

    /// TODO validate the tag vtype
    /// \return true if it is valid, false and appends to SamValidationErrors
    /// if it is not
    static bool isValidValue();
};


#endif
