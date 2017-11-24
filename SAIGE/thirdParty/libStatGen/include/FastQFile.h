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

#ifndef __FASTQ_VALIDATOR_H__
#define __FASTQ_VALIDATOR_H__

#include <iostream>
#include <map>
#include "StringBasics.h"
#include "InputFile.h"
#include "BaseComposition.h"
#include "FastQStatus.h"

/// Class for reading/validating a fastq file.
class FastQFile
{
 public:
    /// Constructor.
    /// /param minReadLength The minimum length that a base sequence must be for
    ///                      it to be valid.
    /// \param numPrintableErrors The maximum number of errors that should be reported
    ///                           in detail before suppressing the errors.
    FastQFile(int minReadLength = 10, int numPrintableErrors = 20);
    
    /// Disable messages - do not write to cout.
    void disableMessages();
    
    /// Enable messages - write to cout.
    void enableMessages();
    
    /// Disable Unique Sequence ID checking
    /// (Unique Sequence ID checking is enabled by default).
    void disableSeqIDCheck();
    
    /// Enable Unique Sequence ID checking.
    /// (Unique Sequence ID checking is enabled by default).
    void enableSeqIDCheck();
    
    /// Interleaved.
    void interleaved();
    
    /// Set the number of errors after which to quit reading/validating a file,
    /// defaults to -1.
    /// \param maxErrors # of errors before quitting, 
    /// -1 indicates to not quit until the entire file has been read/validated (default), 
    /// 0 indicates to quit without reading/validating anything.
    void setMaxErrors(int maxErrors);

    /// Open a FastQFile.
    /// Use the specified SPACE_TYPE to determine BASE, COLOR, or UNKNOWN.
    FastQStatus::Status openFile(const char* fileName,
                                 BaseAsciiMap::SPACE_TYPE spaceType = BaseAsciiMap::UNKNOWN);

    /// Close a FastQFile.
    FastQStatus::Status closeFile();

    /// Check to see if the file is open.
    bool isOpen();

    /// Check to see if the file is at the end of the file.
    bool isEof();

    /// Returns whether or not to keep reading the file,
    /// it stops reading (false) if eof or there is a problem reading the file.
    bool keepReadingFile();

    /// Validate the specified fastq file
    /// \param filename fastq file to be validated.
    /// \param printBaseComp whether or not to print the base composition for the file.
    ///                      true means print it, false means do not.
    /// \param spaceType the spaceType to use for validation - BASE_SPACE, COLOR_SPACE,
    ///                  or UNKNOWN (UNKNOWN means to determine the spaceType to
    ///                  validate against from the first character of the first
    ///                  sequence).
    /// \param printQualAvg  whether or not to print the quality averages for the file.
    ///                      true means to print it, false (default) means do not.
    /// \return the fastq validation status,  SUCCESS on a successfully
    /// validated fastq file.
    FastQStatus::Status validateFastQFile(const String &filename,  
                                          bool printBaseComp,
                                          BaseAsciiMap::SPACE_TYPE spaceType,
                                          bool printQualAvg = false);

    /// Read 1 FastQSequence, validating it.
    FastQStatus::Status readFastQSequence();

    ///////////////////////
    /// @name  Public Sequence Line variables.
    /// Keep public variables for a sequence's line so they can be accessed
    /// without having to do string copies.
    //@{
    String myRawSequence;
    String mySequenceIdLine;
    String mySequenceIdentifier;
    String myPlusLine;
    String myQualityString;
    //@}

    /// Get the space type used for this file.
    inline BaseAsciiMap::SPACE_TYPE getSpaceType()
    {
        return(myBaseComposition.getSpaceType());
    }
    
private:
    // Validates a single fastq sequence from myFile.
    bool validateFastQSequence();

    // Reads and validates the sequence identifier line of a fastq sequence.
    bool validateSequenceIdentifierLine();

    // Reads and validates the raw sequence line(s) and the plus line.  Both are
    // included in one method since it is unknown when the raw sequence line
    // ends until you find the plus line that divides it from the quality
    // string.  Since this method will read the plus line to know when the
    // raw sequence ends, it also validates that line.
    bool validateRawSequenceAndPlusLines();

    // Reads and validates the quality string line(s).
    bool validateQualityStringLines();

    // Method to validate a line that contains part of the raw sequence.
    // offset specifies where in the sequence to start validating.
    bool validateRawSequence(int offset);

    // Method to validate the "+" line that seperates the raw sequence and the
    // quality string.
    bool validateSequencePlus();

    // Method to validate the quality string.
    // offset specifies where in the quality string to start validating.
    bool validateQualityString(int offset);

    // Helper method to read a line from the input file into a string.
    // It also tracks the line number.
    void readLine();

    // Helper method for printing the contents of myErrorString.  It will
    // only print the errors until the maximum number of reportable errors is
    // reached.
    void reportErrorOnLine();

    // Reset the member data for each fastq file.
    void reset();

    // Reset the member data for each sequence.
    void resetForEachSequence();

    // Log the specified message if enabled.
    void logMessage(const char* message);

    // Determine if it is time to quit by checking if we are to quit after a
    // certain number of errors and that many errors have been encountered.
    bool isTimeToQuit();

    void printAvgQual();

    //////////////////////////////////////////////////////////////////////
    // Following member data elements are reset for each validated sequence.
    //

    // Buffer for storing the contents of the line read.
    // Stored as member data so memory allocation is only done once.
    String myLineBuffer;

    // Buffer for storing the error string.  This prevents the reallocation of
    // the string buffer for each error.
    String myErrorString;

    String myTempPartialQuality;

    //////////////////////////////////////////////////////////////////////
    // Following member data elements are reset for each validated file.
    //
    IFILE myFile; // Input file to be read.
    String myFileName; // Name of file being processed.
    int myNumErrors;   // Tracks the number of errors.
    unsigned int myLineNum;    // Track the line number - used for reporting errors.
    BaseComposition myBaseComposition;  // Tracks the base composition.
    std::vector<int> myQualPerCycle;  // Tracks the quality by cycle.
    std::vector<int> myCountPerCycle;  // Tracks the number of entries by cycle.

    // Whether or not to check the sequence identifier for uniqueness.
    // Checking may use up a lot of memory.
    bool myCheckSeqID;

    // Whether or not to check that the file is interleaved.
    // Disabled by myCheckSeqID
    bool myInterleaved;

    // Previous sequence id for checking interleaved.
    std::string myPrevSeqID;

    // Map to track which identifiers have appeared in the file.
    std::map<std::string, unsigned int> myIdentifierMap;
 
    //////////////////////////////////////////////////////////////////////
    // Following member data do not change for each call to the validator.
    //
    int myMinReadLength; // Min Length for a read.
    int myNumPrintableErrors;  // Max number of errors to print the details of.

    // Number of errors after which to quit reading/validating a file.
    // Defaults to -1.
    //   -1 indicates to not quit until the entire file has been read/validated.
    //    0 indicates to quit without reading/validating anything.
    int myMaxErrors;

    // Whether or not messages should be printed.  
    // Defaulted to false (they should be printed).
    bool myDisableMessages;

    // Track if there is a problem reading the file.  If there are read
    // problems, stop reading the file.
    bool myFileProblem;
};

#endif
