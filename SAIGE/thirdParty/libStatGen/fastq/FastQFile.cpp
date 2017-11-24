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

#include <iostream>

#include "InputFile.h"
#include "FastQFile.h"
#include "BaseUtilities.h"

// Constructor.
// minReadLength - The minimum length that a base sequence must be for
//                 it to be valid.
// numPrintableErrors - The maximum number of errors that should be reported
//                      in detail before suppressing the errors.
// 
FastQFile::FastQFile(int minReadLength, int numPrintableErrors)
   : myFile(NULL),
     myBaseComposition(),
     myQualPerCycle(),
     myCountPerCycle(),
     myCheckSeqID(true),
     myInterleaved(false),
     myPrevSeqID(""),
     myMinReadLength(minReadLength),
     myNumPrintableErrors(numPrintableErrors),
     myMaxErrors(-1),
     myDisableMessages(false),
     myFileProblem(false)
{
   // Reset the member data.
   reset();
}


void FastQFile::disableMessages()
{
   myDisableMessages = true;
}


void FastQFile::enableMessages()
{
   myDisableMessages = false;
}


// Disable Unique Sequence ID checking.  
// Unique Sequence ID checking is enabled by default.
void FastQFile::disableSeqIDCheck()
{
    myCheckSeqID = false;
}


// Enable Unique Sequence ID checking.
// Unique Sequence ID checking is enabled by default.
void FastQFile::enableSeqIDCheck()
{
    myCheckSeqID = true;
}


/// Interleaved.
void FastQFile::interleaved()
{
    myInterleaved = true;
}   


// Set the number of errors after which to quit reading/validating a file.
void FastQFile::setMaxErrors(int maxErrors)
{
   myMaxErrors = maxErrors;
}


// Open a FastQFile.
FastQStatus::Status FastQFile::openFile(const char* fileName,
                                        BaseAsciiMap::SPACE_TYPE spaceType)
{
   // reset the member data.
   reset();

   myBaseComposition.resetBaseMapType();
   myBaseComposition.setBaseMapType(spaceType);
   myQualPerCycle.clear();
   myCountPerCycle.clear();

   FastQStatus::Status status = FastQStatus::FASTQ_SUCCESS;

   // Close the file if there is already one open - checked by close.
   status = closeFile();
   if(status == FastQStatus::FASTQ_SUCCESS)
   {
      // Successfully closed a previously opened file if there was one.
      
      // Open the file
      myFile = ifopen(fileName, "rt");
      myFileName = fileName;
      
      if(myFile == NULL)
      {
         // Failed to open the file.
         status = FastQStatus::FASTQ_OPEN_ERROR;
      }
   }

   if(status != FastQStatus::FASTQ_SUCCESS)
   {
      // Failed to open the file.
      std::string errorMessage = "ERROR: Failed to open file: ";
      errorMessage += fileName;
      logMessage(errorMessage.c_str());
   }
   return(status);
}


// Close a FastQFile.
FastQStatus::Status FastQFile::closeFile()
{
   int closeStatus = 0; // Success.

   // If a file has been opened, close it.
   if(myFile != NULL)
   {
      // Close the file.
      closeStatus = ifclose(myFile);
      myFile = NULL;
   }
   if(closeStatus == 0)
   {
      // Success - either there wasn't a file to close or it was closed
      // successfully.
      return(FastQStatus::FASTQ_SUCCESS);
   }
   else
   {
      std::string errorMessage = "Failed to close file: ";
      errorMessage += myFileName.c_str();
      logMessage(errorMessage.c_str());
      return(FastQStatus::FASTQ_CLOSE_ERROR);
   }
}


// Check to see if the file is open.
bool FastQFile::isOpen()
{
   // Check to see if the file is open.
   if((myFile != NULL) && (myFile->isOpen()))
   {
      // File pointer exists and the file is open.
      return true;
   }

   // File is not open.
   return false;
}


// Check to see if the file is at the end of the file.
bool FastQFile::isEof()
{
   // Check to see if the file is open.
   if((myFile != NULL) && (ifeof(myFile)))
   {
      // At EOF.
      return true;
   }

   // Not at EOF.
   return false;
}


// Returns whether or not to keep reading the file.
// Stop reading (false) if eof or there is a problem reading the file.
bool FastQFile::keepReadingFile()
{
   if(isEof() || myFileProblem)
   {
      return(false);
   }
   return(true);
}


// Validate the specified fastq file
FastQStatus::Status FastQFile::validateFastQFile(const String& filename,
                                                 bool printBaseComp,
                                                 BaseAsciiMap::SPACE_TYPE spaceType,
                                                 bool printQualAvg)
{
   // Open the fastqfile.
   if(openFile(filename, spaceType) != FastQStatus::FASTQ_SUCCESS)
   {
      // Failed to open the specified file.
      return(FastQStatus::FASTQ_OPEN_ERROR);
   }

   // Track the total number of sequences that were validated.
   int numSequences = 0;

   // Keep reading the file until there are no more fastq sequences to process
   // and not configured to quit after a certain number of errors or there
   // has not yet been that many errors.
   // Or exit if there is a problem reading the file.
   FastQStatus::Status status = FastQStatus::FASTQ_SUCCESS;
   while (keepReadingFile() &&
          ((myMaxErrors == -1) || (myMaxErrors > myNumErrors)))
   {
      // Validate one sequence.  This call will read all the lines for 
      // one sequence.
      status = readFastQSequence();
      if((status == FastQStatus::FASTQ_SUCCESS) || (status == FastQStatus::FASTQ_INVALID))
      {
         // Read a sequence and it is either valid or invalid, but
         // either way, a sequence was read, so increment the sequence count.
         ++numSequences;
      }
      else
      {
         // Other error, so break out of processing.
         break;
      }
   }
   
   // Report Base Composition Statistics.
   if(printBaseComp)
   {
      myBaseComposition.print();
   }

   if(printQualAvg)
   {
      printAvgQual();
   }

   std::string finishMessage = "Finished processing ";
   finishMessage += myFileName.c_str();
   char buffer[100];
   if(sprintf(buffer, 
              " with %u lines containing %d sequences.", 
              myLineNum, numSequences) > 0)
   {
      finishMessage += buffer;
      logMessage(finishMessage.c_str());
   }
   if(sprintf(buffer, 
              "There were a total of %d errors.", 
              myNumErrors) > 0)
   {
      logMessage(buffer);
   }

   // Close the input file.
   FastQStatus::Status closeStatus = closeFile();

   if((status != FastQStatus::FASTQ_SUCCESS) && (status != FastQStatus::FASTQ_INVALID) &&
      (status != FastQStatus::FASTQ_NO_SEQUENCE_ERROR))
   {
      // Stopped validating due to some error other than invalid, so
      // return that error.
      return(status);
   }
   else if(myNumErrors == 0)
   {
      // No errors, check to see if there were any sequences.
      // Finished processing all of the sequences in the file.
      // If there are no sequences, report an error.
      if(numSequences == 0)
      {
         // Empty file, return error.
         logMessage("ERROR: No FastQSequences in the file.");
         return(FastQStatus::FASTQ_NO_SEQUENCE_ERROR);
      }
      return(FastQStatus::FASTQ_SUCCESS);
   }
   else
   {
      // The file is invalid.  But check the close status.  If the close
      // failed, it means there is a problem with the file itself not just
      // with validation, so the close failure should be returned.
      if(closeStatus != FastQStatus::FASTQ_SUCCESS)
      {
         return(closeStatus);
      }
      return(FastQStatus::FASTQ_INVALID);
   }
}


// Reads and validates a single fastq sequence from myFile.
FastQStatus::Status FastQFile::readFastQSequence()
{
   // First verify that a file is open, if not, return failure.
   if(!isOpen())
   {
      std::string message = 
         "ERROR: Trying to read a fastq file but no file is open.";
      logMessage(message.c_str());
      return(FastQStatus::FASTQ_ORDER_ERROR);
   }

   // Reset variables for each sequence.
   resetForEachSequence();
   
   bool valid = true;

   // No sequence was read.
   if(isTimeToQuit())
   {
      return(FastQStatus::FASTQ_NO_SEQUENCE_ERROR);
   }

   // The first line is the sequence identifier, so validate that.
   valid = validateSequenceIdentifierLine();
   
   if(myFileProblem)
   {
      return(FastQStatus::FASTQ_READ_ERROR);
   }
    
   // If we are at the end of the file, check to see if it is a partial
   // sequence or just an empty line at the end.
   if(ifeof(myFile))
   {
      // If the sequence identifier line was empty and we are at the
      // end of the file, there is nothing more to validate.
      if(mySequenceIdLine.Length() != 0)
      { 
         // There was a sequence identifier line, so this is an incomplete 
         // sequence.
         myErrorString = "Incomplete Sequence.\n";
         reportErrorOnLine();

         valid = false;
      }
      if(valid)
      {
         // Return failure - no sequences were left to read.  At the end
         // of the file.  It wasn't invalid and it wasn't really an error.
         return(FastQStatus::FASTQ_NO_SEQUENCE_ERROR);
      }
      else
      {
         return(FastQStatus::FASTQ_INVALID);
      }
   }

   // If enough errors, quit before reading any more.
   if(isTimeToQuit())
   {
      // Means there was an error, so mark it as invalid.
      return(FastQStatus::FASTQ_INVALID);
   }

   // Validate the Raw Sequence Line(s) and the "+" line.
   valid &= validateRawSequenceAndPlusLines();

   if(myFileProblem)
   {
      return(FastQStatus::FASTQ_READ_ERROR);
   }
    
   // If enough errors, quit before reading any more.
   if(isTimeToQuit())
   {
      return(FastQStatus::FASTQ_INVALID);
   }

   // If it is the end of a file, it is missing the quality string.
   if(ifeof(myFile))
   {
      // There was a sequence identifier line, so this is an incomplete 
      // sequence.
      myErrorString = "Incomplete Sequence, missing Quality String.";
      reportErrorOnLine();
      valid = false;
      return(FastQStatus::FASTQ_INVALID);
   }
    
   // All that is left is to validate the quality string line(s).
   valid &= validateQualityStringLines();

   if(myFileProblem)
   {
      return(FastQStatus::FASTQ_READ_ERROR);
   }
    
   if(valid)
   {
      return(FastQStatus::FASTQ_SUCCESS);
   }
   return(FastQStatus::FASTQ_INVALID);
}


// Reads and validates the sequence identifier line of a fastq sequence.
bool FastQFile::validateSequenceIdentifierLine()
{
   // Read the first line of the sequence.
   int readStatus = mySequenceIdLine.ReadLine(myFile);

   // Check to see if the read was successful.
   if(readStatus <= 0)
   {
      // If EOF, not an error.
      if(ifeof(myFile))
      {
         return true;
      }
      myFileProblem = true;
      myErrorString = "Failure trying to read sequence identifier line";
      reportErrorOnLine();
      return false;
   }

   // If the line is 0 length and it is the end of the file, just
   // return since this is the eof - no error.
   if((mySequenceIdLine.Length() == 0) && (ifeof(myFile)))
   {
      // Not an error, just a new line at the end of the file.
      return true;
   }

   // Increment the line number.
   myLineNum++;
   
   // Verify that the line has at least 2 characters: '@' and at least
   // one character for the sequence identifier.
   if(mySequenceIdLine.Length() < 2)
   {
      // Error. Sequence Identifier line not long enough.
      myErrorString = "The sequence identifier line was too short.";
      reportErrorOnLine();
      return false;
   }
   
   // The sequence identifier line must start wtih a '@'
   if(mySequenceIdLine[0] != '@')
   {
      // Error - sequence identifier line does not begin with an '@'.
      myErrorString = "First line of a sequence does not begin with @";
      reportErrorOnLine();
      return false;
   }

   // Valid Sequence Identifier Line.

   // The sequence identifier ends at the first space or at the end of the
   // line if there is no space.
   // Use fast find since this is a case insensitive search.
   // Start at 1 since we know that 0 is '@'
   int endSequenceIdentifier = mySequenceIdLine.FastFindChar(' ', 1);
   
   // Check if a " " was found.
   if(endSequenceIdentifier == -1)
   {
      // Did not find a ' ', so the identifier is the rest of the line.
      // It starts at 1 since @ is at offset 0.
      mySequenceIdentifier = (mySequenceIdLine.SubStr(1)).c_str();
   }
   else
   {
      // Found a ' ', so the identifier ends just before that.
      // The sequence identifier must be at least 1 character long, 
      // therefore the endSequenceIdentifier must be greater than 1.
      if(endSequenceIdentifier <= 1)
      {
         myErrorString = 
            "No Sequence Identifier specified before the comment.";
         reportErrorOnLine();
         return false;
      }

      mySequenceIdentifier = 
         (mySequenceIdLine.SubStr(1, endSequenceIdentifier - 1)).c_str();
   }

   // If myInterleaved, validate matches the previous seqID.
   if(myInterleaved && (myPrevSeqID != ""))
   {
       // Valid if the sequence identifiers are identical or if
       // the only difference is a trailing 1 or 2.
       if(myPrevSeqID.compare(mySequenceIdentifier) != 0)
       {
           // Compare all but the last characters, then check the last characters for 1 or 2.
           if((myPrevSeqID.compare(0, myPrevSeqID.length()-1, mySequenceIdentifier.c_str(), mySequenceIdentifier.Length()-1) != 0) || 
              (((myPrevSeqID[myPrevSeqID.length()-1] != '1') || (mySequenceIdentifier[mySequenceIdentifier.Length()-1] != '2')) && 
               (myPrevSeqID[myPrevSeqID.length()-1] != mySequenceIdentifier[mySequenceIdentifier.Length()-1])))
           {
               myErrorString = "Interleaved: consecutive reads do not have matching sequence identifiers: ";
               myErrorString += mySequenceIdentifier.c_str();
               myErrorString += " and ";
               myErrorString += myPrevSeqID.c_str();
               reportErrorOnLine();
               myPrevSeqID.clear();
               return(false);
           }
       }
       myPrevSeqID.clear();
   }
   else
   {
       if(myInterleaved)
       {
           myPrevSeqID = mySequenceIdentifier.c_str();
       }

       // Check if sequence identifier should be validated for uniqueness if it is 
       // not the 2nd in an interleaved pair.
       if(myCheckSeqID)
       {
           // Check to see if the sequenceIdentifier is a repeat by adding
           // it to the set and seeing if it already existed.
           std::pair<std::map<std::string, unsigned int>::iterator,bool> insertResult;
           insertResult = 
               myIdentifierMap.insert(std::make_pair(mySequenceIdentifier.c_str(), 
                                                     myLineNum));
           
           if(insertResult.second == false)
           {
               // Sequence Identifier is a repeat.
               myErrorString = "Repeated Sequence Identifier: ";
               myErrorString += mySequenceIdentifier.c_str();
               myErrorString += " at Lines ";
               myErrorString += insertResult.first->second;
               myErrorString += " and ";
               myErrorString += myLineNum;
               reportErrorOnLine();
               return(false);
           }
       }
   }

   // Valid, return true.
   return(true);
}


// Reads and validates the raw sequence line(s) and the plus line.  Both are
// included in one method since it is unknown when the raw sequence line
// ends until you find the plus line that divides it from the quality
// string.  Since this method will read the plus line to know when the
// raw sequence ends, it also validates that line.
bool FastQFile::validateRawSequenceAndPlusLines()
{
   // Read the raw sequence.
   int readStatus = myRawSequence.ReadLine(myFile);

   myLineNum++;

   if(readStatus <= 0)
   {
      myFileProblem = true;
      myErrorString = "Failure trying to read sequence line";
      reportErrorOnLine();
      return false;
   }

   // Offset into the raw sequence to be validated.
   int offset = 0;
   
   // Validate the raw sequence.
   bool valid = validateRawSequence(offset);

   // Increment the offset for what was just read.
   offset = myRawSequence.Length();

   // The next line is either a continuation of the raw sequence or it starts
   // with a '+'
   // Keep validating til the '+' line or the end of file is found.
   bool stillRawLine = true;
   while(stillRawLine && 
         !ifeof(myFile))
   {
      // If enough errors, quit before reading any more.
      if(isTimeToQuit())
      {
         return(false);
      }

      // Read the next line.
      // Read into the plus line, but if it isn't a plus line, then
      // it will be copied into the raw sequence line.
      readStatus = myPlusLine.ReadLine(myFile);
      myLineNum++;

      if(readStatus <= 0)
      {
         myFileProblem = true;
         myErrorString = "Failure trying to read sequence/plus line";
         reportErrorOnLine();
         return false;
      }

      // Check if the next line is blank
      if(myPlusLine.Length() == 0)
      {
         // The next line is blank.  Assume it is part of the raw sequence and
         // report an error since there are no valid characters on the line.
         myErrorString = 
            "Looking for continuation of Raw Sequence or '+' instead found a blank line, assuming it was part of Raw Sequence.";
         reportErrorOnLine();
      }
      // Check for the plus line.
      else if(myPlusLine[0] == '+')
      {
         // This is the + line.
         valid &= validateSequencePlus();
         stillRawLine = false;
      }
      else
      {
         // Not a plus line, so assume this is a continuation of the Raw
         // Sequence.
         // Copy from the plus line to the raw sequence line.
         myRawSequence += myPlusLine;
         myPlusLine.SetLength(0);
         valid &= validateRawSequence(offset);
         
         // Increment the offset.
         offset = myRawSequence.Length();
      }
   }
   
   // If enough errors, quit before reading any more.
   if(isTimeToQuit())
   {
      return(false);
   }
   
   // Now that the entire raw sequence has been obtained, check its length
   // against the minimum allowed length.
   if(myRawSequence.Length() < myMinReadLength)
   {
      // Raw sequence is not long enough - error.
      myErrorString = "Raw Sequence is shorter than the min read length: ";
      myErrorString += myRawSequence.Length();
      myErrorString += " < ";
      myErrorString += myMinReadLength;
      reportErrorOnLine();
      valid = false;
   }

   // If enough errors, quit before reading any more.
   if(isTimeToQuit())
   {
      return(false);
   }

   // if the flag still indicates it is processing the raw sequence that means
   // we reached the end of the file without a '+' line.  So report that error.
   if(stillRawLine)
   {
      myErrorString = "Reached the end of the file without a '+' line.";
      reportErrorOnLine();
      valid = false;
   }

   return(valid);
}


// Reads and validates the quality string line(s).
bool FastQFile::validateQualityStringLines()
{
   // Read the quality string.
   int readStatus = myQualityString.ReadLine(myFile);
   myLineNum++;

   if(readStatus <= 0)
   {
      myFileProblem = true;
      myErrorString = "Failure trying to read quality line";
      reportErrorOnLine();
      return false;
   }

   // track the offset into the quality string to validate.
   int offset = 0;

   // Validate this line of the quality string.
   bool valid = validateQualityString(offset);

   offset = myQualityString.Length();

   // Keep reading quality string lines until the length of the 
   // raw sequence has been hit or the end of the file is reached.
   while((myQualityString.Length() < myRawSequence.Length()) && 
         (!ifeof(myFile)))
   {
      // If enough errors, quit before reading any more.
      if(isTimeToQuit())
      {
         return(false);
      }

      // Read another line of the quality string.
      readStatus = myTempPartialQuality.ReadLine(myFile);
      myLineNum++;

      if(readStatus <= 0)
      {
         myFileProblem = true;
         myErrorString = "Failure trying to read quality line";
         reportErrorOnLine();
         return false;
      }

      myQualityString += myTempPartialQuality;
      myTempPartialQuality.Clear();

      // Validate this line of the quality string.
      valid &= validateQualityString(offset);
      offset = myQualityString.Length();
   }

   // If enough errors, quit before reading any more.
   if(isTimeToQuit())
   {
      return(false);
   }

   // Validate that the quality string length is the same as the
   // raw sequence length.
   if(myQualityString.Length() != myRawSequence.Length()) 
   {
      myErrorString = "Quality string length (";
      myErrorString += myQualityString.Length();
      myErrorString += ") does not equal raw sequence length (";
      myErrorString += myRawSequence.Length();
      myErrorString += ")";
      reportErrorOnLine();
      valid = false;
   }
   return(valid);
}


// Method to validate a line that contains part of the raw sequence.
bool FastQFile::validateRawSequence(int offset)
{
   bool validBase = false; 
   bool valid = true;

   // Loop through validating each character is valid for the raw sequence.
   for(int sequenceIndex = offset; sequenceIndex < myRawSequence.Length(); 
       sequenceIndex++)
   {
      // Update the composition for this position.  Returns false if the
      // character was not a valid base.
      validBase = 
         myBaseComposition.updateComposition(sequenceIndex, 
                                             myRawSequence[sequenceIndex]);
      // Check the return
      if(!validBase)
      {
         // Error, found a value that is not a valid base character.
         myErrorString = "Invalid character ('";
         myErrorString += myRawSequence[sequenceIndex];
         myErrorString += "') in base sequence.";
         reportErrorOnLine();
         valid = false;
         // If enough errors, quit before reading any more.
         if(isTimeToQuit())
         {
            return(false);
         }
      }
   }
   return(valid);
}


// Method to validate the "+" line that seperates the raw sequence and the
// quality string.
bool FastQFile::validateSequencePlus()
{
   // Validate that optional sequence identifier is the same
   // as the one on the @ line.

   // Check to see if there is more to the line than just the plus
   int lineLength = myPlusLine.Length();
   
   // If the line is only 1 character or the second character is a space,
   // then there is no sequence identifier on this line and there is nothing
   // further to validate.
   if((lineLength == 1) || (myPlusLine[1] == ' '))
   {
      // No sequence identifier, so just return valid.
      return true;
   }
   
   // There is a sequence identifier on this line, so validate that
   // it matches the one from the associated @ line.
   // The read in line must be at least 1 more character ('+') than the
   // sequence identifier read from the '@' line.
   // If it is not longer than the sequence identifier, then we know that it
   // cannot be the same.
   int sequenceIdentifierLength = mySequenceIdentifier.Length();
   if(lineLength <= sequenceIdentifierLength)
   {
      myErrorString = 
         "Sequence Identifier on '+' line does not equal the one on the '@' line.";
      reportErrorOnLine();
      return false;
   }

   bool same = true;
   int seqIndex = 0;
   int lineIndex = 1;  // Start at 1 since offset 0 has '+'

   // Loop through the sequence index and the line buffer verifying they
   // are the same until a difference is found or the end of the sequence
   // identifier is found.
   while((same == true) && (seqIndex < sequenceIdentifierLength))
   {
      if(myPlusLine[lineIndex] != mySequenceIdentifier[seqIndex])
      {
         myErrorString = 
            "Sequence Identifier on '+' line does not equal the one on the '@' line.";
         reportErrorOnLine();
         same = false;
      }
      lineIndex++;
      seqIndex++;
   }
   return(same);
}


// Method to validate the quality string.
bool FastQFile::validateQualityString(int offset)
{
   bool valid = true;
   if(myQualityString.Length() > (int)(myQualPerCycle.size()))
   {
       myQualPerCycle.resize(myQualityString.Length());
       myCountPerCycle.resize(myQualityString.Length());
   }
   // For each character in the line, verify that it is ascii > 32.
   for(int i=offset; i < myQualityString.Length(); i++)
   {
      if(myQualityString[i] <= 32)
      {
         myErrorString = "Invalid character ('";
         myErrorString += myQualityString[i];
         myErrorString += "') in quality string.";
         reportErrorOnLine();
         valid = false;
         // If enough errors, quit before reading any more.
         if(isTimeToQuit())
         {
            return(false);
         }
      }
      else
      {
          myQualPerCycle[i] += BaseUtilities::getPhredBaseQuality(myQualityString[i]);
          myCountPerCycle[i] += 1;
      }
   }
   return(valid);
}


// Helper method for printing the contents of myErrorString.  It will
// only print the errors until the maximum number of reportable errors is
// reached.
void FastQFile::reportErrorOnLine()
{
   // Increment the total number of errors.
   myNumErrors++;
   
   // Only display the first X number of errors.
   if(myNumErrors <= myNumPrintableErrors)
   {
      // Write the error with the line number.
      char buffer[100];
      sprintf(buffer, "ERROR on Line %u: ", myLineNum);
      std::string message = buffer;
      message += myErrorString.c_str();
      logMessage(message.c_str());
   }
}


// Reset member data that is unique for each fastQFile.
void FastQFile::reset()
{
   // Each fastq file processing needs to also reset the member data for each
   // sequence.
   resetForEachSequence();
   myNumErrors = 0;  // per fastqfile
   myLineNum = 0;    // per fastqfile
   myFileName.SetLength(0);  // reset the filename string.
   myIdentifierMap.clear(); // per fastqfile
   myBaseComposition.clear(); // clear the base composition.
   myQualPerCycle.clear();
   myCountPerCycle.clear();
   myFileProblem = false;
}


// Reset the member data that is unique for each sequence.
void FastQFile::resetForEachSequence()
{
   myLineBuffer.SetLength(0);
   myErrorString.SetLength(0);
   myRawSequence.SetLength(0);
   mySequenceIdLine.SetLength(0);
   mySequenceIdentifier.SetLength(0);
   myPlusLine.SetLength(0);
   myQualityString.SetLength(0);
   myTempPartialQuality.SetLength(0);
}


void FastQFile::logMessage(const char* logMessage)
{
   // Write the message if they are not disabled.
   if(!myDisableMessages)
   {
      std::cout << logMessage << std::endl;
   }
}


// Determine if it is time to quit by checking if we are to quit after a
// certain number of errors and that many errors have been encountered.
bool FastQFile::isTimeToQuit()
{
   // It is time to quit if we are to quit after a certain number of errors
   // and that many errors have been encountered.
   if((myMaxErrors != -1) && (myNumErrors >= myMaxErrors))
   {
      return(true);
   }
   return(false);
}


void FastQFile::printAvgQual()
{
   std::cout << std::endl << "Average Phred Quality by Read Index (starts at 0):" << std::endl;
   std::cout.precision(2);
   std::cout << std::fixed << "Read Index\tAverage Quality" 
             << std::endl;
   if(myQualPerCycle.size() != myCountPerCycle.size())
   {
       // This is a code error and should NEVER happen.
       std::cerr << "ERROR calculating the average Qualities per cycle\n";
   }

   double sumQual = 0;
   double count = 0;
   double avgQual = 0;
   for(unsigned int i = 0; i < myQualPerCycle.size(); i++)
   {
       avgQual = 0;
       if(myCountPerCycle[i] != 0)
       {
           avgQual = myQualPerCycle[i] / (double)(myCountPerCycle[i]);
       }
       std::cout << i << "\t" << avgQual << "\n";
       sumQual += myQualPerCycle[i];
       count += myCountPerCycle[i];
   }
   std::cout << std::endl;
   avgQual = 0;
   if(count != 0)
   {
       avgQual = sumQual / count;
   }
   std::cout << "Overall Average Phred Quality = " << avgQual << std::endl;
}
