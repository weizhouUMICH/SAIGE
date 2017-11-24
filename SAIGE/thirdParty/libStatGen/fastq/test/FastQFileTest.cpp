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

#include "FastQFile.h"
#include <assert.h>

const String FIRST_SEQID_LINE = "@Valid with comment";
const String FIRST_SEQID = "Valid";
const String FIRST_RAW_SEQ = "ACTGNactng.0123";
const String FIRST_PLUS_LINE = "+";
const String FIRST_QUALITY = "!#\"$%&'()*+,-./";
const String SECOND_SEQID_LINE = "@Valid1 with comment";
const String SECOND_SEQID = "Valid1";
const String SECOND_RAW_SEQ = "ACTGACTNactngaac";
const String SECOND_PLUS_LINE = "+";
const String SECOND_QUALITY = "0123456789:;<=>@";
const String THIRD_SEQID_LINE = "@Valid2";
const String THIRD_SEQID = "Valid2";
const String THIRD_RAW_SEQ = "A123.0321.011";
const String THIRD_PLUS_LINE = "+";
const String THIRD_QUALITY = "?@ABCDEFGHIJK";
const String FOURTH_SEQID_LINE = "@Valid3";
const String FOURTH_SEQID = "Valid3";
const String FOURTH_RAW_SEQ = "ACTGACTNactngACTGACTNactng";
const String FOURTH_PLUS_LINE = "+";
const String FOURTH_QUALITY = "LMNOPQRSTUVWXYZ[\\]^_'abcde";
const String FIFTH_SEQID_LINE = "@Valid4";
const String FIFTH_SEQID = "Valid4";
const String FIFTH_RAW_SEQ = "ACTGACTNactngACTGACTNactng";
const String FIFTH_PLUS_LINE = "+";
const String FIFTH_QUALITY = "fghijklmnopqrstuvwxyz{|}~~";
const String SIXTH_SEQID_LINE = "@";
const String SIXTH_SEQID = "";
const String SIXTH_RAW_SEQ = "ACTGACTNactng";
const String SIXTH_PLUS_LINE = "+";
const String SIXTH_QUALITY = "?@ABCDEFGHIJK";
const String SEVENTH_SEQID_LINE = "Line no start with @";
const String SEVENTH_SEQID = "";
const String SEVENTH_RAW_SEQ = "ACTGACTNactng";
const String SEVENTH_PLUS_LINE = "+";
const String SEVENTH_QUALITY = "LMNOPQRSTUVWX";
const String EIGHTH_SEQID_LINE = "@ a";
const String EIGHTH_SEQID = "";
const String EIGHTH_RAW_SEQ = "ACTGACTNactng";
const String EIGHTH_PLUS_LINE = "+";
const String EIGHTH_QUALITY = "YZ[\\]^_'abcde";
const String NINTH_SEQID_LINE = "@ ";
const String NINTH_SEQID = "";
const String NINTH_RAW_SEQ = "ACTGACTNactng";
const String NINTH_PLUS_LINE = "+";
const String NINTH_QUALITY = "fghijklmnopqr";
const String TENTH_SEQID_LINE = "@Valid";
const String TENTH_SEQID = "Valid";
const String TENTH_RAW_SEQ = "ACTGNactng";
const String TENTH_PLUS_LINE = "+";
const String TENTH_QUALITY = "!#\"$%&'()*";
const String ELEVENTH_SEQID_LINE = "@RawError1";
const String ELEVENTH_SEQID = "RawError1";
const String ELEVENTH_RAW_SEQ = "ACTNaHtng0aBZa";
const String ELEVENTH_PLUS_LINE = "+";
const String ELEVENTH_QUALITY = "ACTNactng0aBaZ";
const String TWELFTH_SEQID_LINE = "@RawError2";
const String TWELFTH_SEQID = "RawError2";
const String TWELFTH_RAW_SEQ = "aaa";
const String TWELFTH_PLUS_LINE = "+";
const String TWELFTH_QUALITY = "aaa";
const String THIRTEENTH_SEQID_LINE = "@RawError3";
const String THIRTEENTH_SEQID = "RawError3";
const String THIRTEENTH_RAW_SEQ = "ACTGACTNactng";
const String THIRTEENTH_PLUS_LINE = "+";
const String THIRTEENTH_QUALITY = "ACTGACTNactng";
const String FOURTEENTH_SEQID_LINE = "@QualityError1";
const String FOURTEENTH_SEQID = "QualityError1";
const String FOURTEENTH_RAW_SEQ = "ACTGCacgnc";
const String FOURTEENTH_PLUS_LINE = "+";
const String FOURTEENTH_QUALITY = "ac gcacg n";
const String FIFTEENTH_SEQID_LINE = "@QualityError2";
const String FIFTEENTH_SEQID = "QualityError2";
const String FIFTEENTH_RAW_SEQ = "ACTGCacgnc";
const String FIFTEENTH_PLUS_LINE = "+";
const String FIFTEENTH_QUALITY = "actgc@cgnc";
const String SIXTEENTH_SEQID_LINE = "@QualityError3";
const String SIXTEENTH_SEQID = "QualityError3";
const String SIXTEENTH_RAW_SEQ = "ACTGCacgnc";
const String SIXTEENTH_PLUS_LINE = "+";
const String SIXTEENTH_QUALITY = "actgc77acgnc";
const String SEVENTEENTH_SEQID_LINE = "@PlusValid1";
const String SEVENTEENTH_SEQID = "PlusValid1";
const String SEVENTEENTH_RAW_SEQ = "ACTGCacgnc";
const String SEVENTEENTH_PLUS_LINE = "+PlusValid1";
const String SEVENTEENTH_QUALITY = "actgcacgnc";
const String EIGHTEENTH_SEQID_LINE = "@PlusValid2";
const String EIGHTEENTH_SEQID = "PlusValid2";
const String EIGHTEENTH_RAW_SEQ = "ACTGCacgnc";
const String EIGHTEENTH_PLUS_LINE = "+PlusValid2 Added comment";
const String EIGHTEENTH_QUALITY = "actgcacgnc";
const String NINETEENTH_SEQID_LINE = "@PlusError1";
const String NINETEENTH_SEQID = "PlusError1";
const String NINETEENTH_RAW_SEQ = "ACTGCacgnc";
const String NINETEENTH_PLUS_LINE = "+PlusError2";
const String NINETEENTH_QUALITY = "actgcacgnc";

const String TWENTIETH_SEQID_LINE = "@InvalidColor";
const String TWENTIETH_SEQID = "InvalidColor";
const String TWENTIETH_RAW_SEQ = "0123.0321.011";
const String TWENTIETH_PLUS_LINE = "+";
const String TWENTIETH_QUALITY = "0123.0321.011";


const String TWENTY_FIRST_SEQID_LINE = "@PlusError2";
const String TWENTY_FIRST_SEQID = "PlusError2";
const String TWENTY_FIRST_RAW_SEQ = "ACTGCacgnc";
const String TWENTY_FIRST_PLUS_LINE = "";
const String TWENTY_FIRST_QUALITY = "";

void testReadUnOpenedFile()
{
   FastQFile fastqfile;

   assert(fastqfile.isOpen() == false);
   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_ORDER_ERROR);
   assert(fastqfile.isOpen() == false);
}

void testOpenFile()
{
   FastQFile fastqfile;

   // Test for non-existent file.
   assert(fastqfile.isOpen() == false);
   assert(fastqfile.openFile("noexist.txt", BaseAsciiMap::UNKNOWN) == FastQStatus::FASTQ_OPEN_ERROR);
   assert(fastqfile.isOpen() == false);


}


void testCloseFile()
{
   FastQFile fastqfile;

   // Test closing a file even though there isn't one open - counts as success.
   assert(fastqfile.isOpen() == false);
   assert(fastqfile.closeFile() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.isOpen() == false);
}


void testReadSequence()
{
   FastQFile fastqfile;
   
   assert(fastqfile.isOpen() == false);
   assert(fastqfile.openFile("testFile.txt") == FastQStatus::FASTQ_SUCCESS);

   assert(fastqfile.isOpen() == true);

   assert(fastqfile.getSpaceType() == BaseAsciiMap::UNKNOWN);

   // Read Sequence from test file.
   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == FIRST_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FIRST_SEQID);
   assert(fastqfile.myRawSequence == FIRST_RAW_SEQ);
   assert(fastqfile.myPlusLine == FIRST_PLUS_LINE);
   assert(fastqfile.myQualityString == FIRST_QUALITY);
   assert(fastqfile.getSpaceType() == BaseAsciiMap::BASE_SPACE);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == SECOND_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SECOND_SEQID);
   assert(fastqfile.myRawSequence == SECOND_RAW_SEQ);
   assert(fastqfile.myPlusLine == SECOND_PLUS_LINE);
   assert(fastqfile.myQualityString == SECOND_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == THIRD_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == THIRD_SEQID);
   assert(fastqfile.myRawSequence == THIRD_RAW_SEQ);
   assert(fastqfile.myPlusLine == THIRD_PLUS_LINE);
   assert(fastqfile.myQualityString == THIRD_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == FOURTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FOURTH_SEQID);
   assert(fastqfile.myRawSequence == FOURTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FOURTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FOURTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == FIFTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FIFTH_SEQID);
   assert(fastqfile.myRawSequence == FIFTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FIFTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FIFTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SIXTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SIXTH_SEQID);
   assert(fastqfile.myRawSequence == SIXTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SIXTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SIXTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SEVENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SEVENTH_SEQID);
   assert(fastqfile.myRawSequence == SEVENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SEVENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SEVENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == EIGHTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == EIGHTH_SEQID);
   assert(fastqfile.myRawSequence == EIGHTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == EIGHTH_PLUS_LINE);
   assert(fastqfile.myQualityString == EIGHTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == NINTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == NINTH_SEQID);
   assert(fastqfile.myRawSequence == NINTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == NINTH_PLUS_LINE);
   assert(fastqfile.myQualityString == NINTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TENTH_SEQID);
   assert(fastqfile.myRawSequence == TENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == TENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == TENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == ELEVENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == ELEVENTH_SEQID);
   assert(fastqfile.myRawSequence == ELEVENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == ELEVENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == ELEVENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TWELFTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TWELFTH_SEQID);
   assert(fastqfile.myRawSequence == TWELFTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == TWELFTH_PLUS_LINE);
   assert(fastqfile.myQualityString == TWELFTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == THIRTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == THIRTEENTH_SEQID);
   assert(fastqfile.myRawSequence == THIRTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == THIRTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == THIRTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == FOURTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FOURTEENTH_SEQID);
   assert(fastqfile.myRawSequence == FOURTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FOURTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FOURTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == FIFTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FIFTEENTH_SEQID);
   assert(fastqfile.myRawSequence == FIFTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FIFTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FIFTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SIXTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SIXTEENTH_SEQID);
   assert(fastqfile.myRawSequence == SIXTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SIXTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SIXTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == SEVENTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SEVENTEENTH_SEQID);
   assert(fastqfile.myRawSequence == SEVENTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SEVENTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SEVENTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == EIGHTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == EIGHTEENTH_SEQID);
   assert(fastqfile.myRawSequence == EIGHTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == EIGHTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == EIGHTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == NINETEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == NINETEENTH_SEQID);
   assert(fastqfile.myRawSequence == NINETEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == NINETEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == NINETEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TWENTIETH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TWENTIETH_SEQID);
   assert(fastqfile.myRawSequence == TWENTIETH_RAW_SEQ);
   assert(fastqfile.myPlusLine == TWENTIETH_PLUS_LINE);
   assert(fastqfile.myQualityString == TWENTIETH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TWENTY_FIRST_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TWENTY_FIRST_SEQID);
   assert(fastqfile.myRawSequence == TWENTY_FIRST_RAW_SEQ);
   assert(fastqfile.myPlusLine == TWENTY_FIRST_PLUS_LINE);
   assert(fastqfile.myQualityString == TWENTY_FIRST_QUALITY);

   // Close the file, and verify isOpen = false;
   assert(fastqfile.closeFile() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.isOpen() == false);
   

   //////////////////////////////////
   // Repeat test specifying base space
   assert(fastqfile.isOpen() == false);
   assert(fastqfile.openFile("testFile.txt", BaseAsciiMap::BASE_SPACE) == FastQStatus::FASTQ_SUCCESS);

   assert(fastqfile.isOpen() == true);

   assert(fastqfile.getSpaceType() == BaseAsciiMap::BASE_SPACE);

   // Read Sequence from test file.
   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == FIRST_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FIRST_SEQID);
   assert(fastqfile.myRawSequence == FIRST_RAW_SEQ);
   assert(fastqfile.myPlusLine == FIRST_PLUS_LINE);
   assert(fastqfile.myQualityString == FIRST_QUALITY);
   assert(fastqfile.getSpaceType() == BaseAsciiMap::BASE_SPACE);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == SECOND_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SECOND_SEQID);
   assert(fastqfile.myRawSequence == SECOND_RAW_SEQ);
   assert(fastqfile.myPlusLine == SECOND_PLUS_LINE);
   assert(fastqfile.myQualityString == SECOND_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == THIRD_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == THIRD_SEQID);
   assert(fastqfile.myRawSequence == THIRD_RAW_SEQ);
   assert(fastqfile.myPlusLine == THIRD_PLUS_LINE);
   assert(fastqfile.myQualityString == THIRD_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == FOURTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FOURTH_SEQID);
   assert(fastqfile.myRawSequence == FOURTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FOURTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FOURTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == FIFTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FIFTH_SEQID);
   assert(fastqfile.myRawSequence == FIFTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FIFTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FIFTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SIXTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SIXTH_SEQID);
   assert(fastqfile.myRawSequence == SIXTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SIXTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SIXTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SEVENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SEVENTH_SEQID);
   assert(fastqfile.myRawSequence == SEVENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SEVENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SEVENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == EIGHTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == EIGHTH_SEQID);
   assert(fastqfile.myRawSequence == EIGHTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == EIGHTH_PLUS_LINE);
   assert(fastqfile.myQualityString == EIGHTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == NINTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == NINTH_SEQID);
   assert(fastqfile.myRawSequence == NINTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == NINTH_PLUS_LINE);
   assert(fastqfile.myQualityString == NINTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TENTH_SEQID);
   assert(fastqfile.myRawSequence == TENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == TENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == TENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == ELEVENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == ELEVENTH_SEQID);
   assert(fastqfile.myRawSequence == ELEVENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == ELEVENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == ELEVENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TWELFTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TWELFTH_SEQID);
   assert(fastqfile.myRawSequence == TWELFTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == TWELFTH_PLUS_LINE);
   assert(fastqfile.myQualityString == TWELFTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == THIRTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == THIRTEENTH_SEQID);
   assert(fastqfile.myRawSequence == THIRTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == THIRTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == THIRTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == FOURTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FOURTEENTH_SEQID);
   assert(fastqfile.myRawSequence == FOURTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FOURTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FOURTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == FIFTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FIFTEENTH_SEQID);
   assert(fastqfile.myRawSequence == FIFTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FIFTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FIFTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SIXTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SIXTEENTH_SEQID);
   assert(fastqfile.myRawSequence == SIXTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SIXTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SIXTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == SEVENTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SEVENTEENTH_SEQID);
   assert(fastqfile.myRawSequence == SEVENTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SEVENTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SEVENTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == EIGHTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == EIGHTEENTH_SEQID);
   assert(fastqfile.myRawSequence == EIGHTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == EIGHTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == EIGHTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == NINETEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == NINETEENTH_SEQID);
   assert(fastqfile.myRawSequence == NINETEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == NINETEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == NINETEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TWENTIETH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TWENTIETH_SEQID);
   assert(fastqfile.myRawSequence == TWENTIETH_RAW_SEQ);
   assert(fastqfile.myPlusLine == TWENTIETH_PLUS_LINE);
   assert(fastqfile.myQualityString == TWENTIETH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TWENTY_FIRST_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TWENTY_FIRST_SEQID);
   assert(fastqfile.myRawSequence == TWENTY_FIRST_RAW_SEQ);
   assert(fastqfile.myPlusLine == TWENTY_FIRST_PLUS_LINE);
   assert(fastqfile.myQualityString == TWENTY_FIRST_QUALITY);

   // Close the file, and verify isOpen = false;
   assert(fastqfile.closeFile() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.isOpen() == false);
   

   ////////////////////////////////
   // Repeat test specifying color space
   assert(fastqfile.isOpen() == false);
   assert(fastqfile.openFile("testFile.txt", BaseAsciiMap::COLOR_SPACE) == FastQStatus::FASTQ_SUCCESS);

   assert(fastqfile.isOpen() == true);

   assert(fastqfile.getSpaceType() == BaseAsciiMap::COLOR_SPACE);

   // Read Sequence from test file.
   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == FIRST_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FIRST_SEQID);
   assert(fastqfile.myRawSequence == FIRST_RAW_SEQ);
   assert(fastqfile.myPlusLine == FIRST_PLUS_LINE);
   assert(fastqfile.myQualityString == FIRST_QUALITY);
   assert(fastqfile.getSpaceType() == BaseAsciiMap::COLOR_SPACE);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SECOND_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SECOND_SEQID);
   assert(fastqfile.myRawSequence == SECOND_RAW_SEQ);
   assert(fastqfile.myPlusLine == SECOND_PLUS_LINE);
   assert(fastqfile.myQualityString == SECOND_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == THIRD_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == THIRD_SEQID);
   assert(fastqfile.myRawSequence == THIRD_RAW_SEQ);
   assert(fastqfile.myPlusLine == THIRD_PLUS_LINE);
   assert(fastqfile.myQualityString == THIRD_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == FOURTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FOURTH_SEQID);
   assert(fastqfile.myRawSequence == FOURTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FOURTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FOURTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == FIFTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FIFTH_SEQID);
   assert(fastqfile.myRawSequence == FIFTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FIFTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FIFTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SIXTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SIXTH_SEQID);
   assert(fastqfile.myRawSequence == SIXTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SIXTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SIXTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SEVENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SEVENTH_SEQID);
   assert(fastqfile.myRawSequence == SEVENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SEVENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SEVENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == EIGHTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == EIGHTH_SEQID);
   assert(fastqfile.myRawSequence == EIGHTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == EIGHTH_PLUS_LINE);
   assert(fastqfile.myQualityString == EIGHTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == NINTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == NINTH_SEQID);
   assert(fastqfile.myRawSequence == NINTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == NINTH_PLUS_LINE);
   assert(fastqfile.myQualityString == NINTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TENTH_SEQID);
   assert(fastqfile.myRawSequence == TENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == TENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == TENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == ELEVENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == ELEVENTH_SEQID);
   assert(fastqfile.myRawSequence == ELEVENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == ELEVENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == ELEVENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TWELFTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TWELFTH_SEQID);
   assert(fastqfile.myRawSequence == TWELFTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == TWELFTH_PLUS_LINE);
   assert(fastqfile.myQualityString == TWELFTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == THIRTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == THIRTEENTH_SEQID);
   assert(fastqfile.myRawSequence == THIRTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == THIRTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == THIRTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == FOURTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FOURTEENTH_SEQID);
   assert(fastqfile.myRawSequence == FOURTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FOURTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FOURTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == FIFTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FIFTEENTH_SEQID);
   assert(fastqfile.myRawSequence == FIFTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FIFTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FIFTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SIXTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SIXTEENTH_SEQID);
   assert(fastqfile.myRawSequence == SIXTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SIXTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SIXTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SEVENTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SEVENTEENTH_SEQID);
   assert(fastqfile.myRawSequence == SEVENTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SEVENTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SEVENTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == EIGHTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == EIGHTEENTH_SEQID);
   assert(fastqfile.myRawSequence == EIGHTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == EIGHTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == EIGHTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == NINETEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == NINETEENTH_SEQID);
   assert(fastqfile.myRawSequence == NINETEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == NINETEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == NINETEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TWENTIETH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TWENTIETH_SEQID);
   assert(fastqfile.myRawSequence == TWENTIETH_RAW_SEQ);
   assert(fastqfile.myPlusLine == TWENTIETH_PLUS_LINE);
   assert(fastqfile.myQualityString == TWENTIETH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TWENTY_FIRST_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TWENTY_FIRST_SEQID);
   assert(fastqfile.myRawSequence == TWENTY_FIRST_RAW_SEQ);
   assert(fastqfile.myPlusLine == TWENTY_FIRST_PLUS_LINE);
   assert(fastqfile.myQualityString == TWENTY_FIRST_QUALITY);

   // Close the file, and verify isOpen = false;
   assert(fastqfile.closeFile() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.isOpen() == false);   

   ////////////////////////////////
   // Repeat test specifying Unknown space
   assert(fastqfile.isOpen() == false);
   assert(fastqfile.openFile("testFile.txt", BaseAsciiMap::UNKNOWN) == FastQStatus::FASTQ_SUCCESS);

   assert(fastqfile.isOpen() == true);

   assert(fastqfile.getSpaceType() == BaseAsciiMap::UNKNOWN);

   // Read Sequence from test file.
   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == FIRST_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FIRST_SEQID);
   assert(fastqfile.myRawSequence == FIRST_RAW_SEQ);
   assert(fastqfile.myPlusLine == FIRST_PLUS_LINE);
   assert(fastqfile.myQualityString == FIRST_QUALITY);
   assert(fastqfile.getSpaceType() == BaseAsciiMap::BASE_SPACE);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == SECOND_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SECOND_SEQID);
   assert(fastqfile.myRawSequence == SECOND_RAW_SEQ);
   assert(fastqfile.myPlusLine == SECOND_PLUS_LINE);
   assert(fastqfile.myQualityString == SECOND_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == THIRD_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == THIRD_SEQID);
   assert(fastqfile.myRawSequence == THIRD_RAW_SEQ);
   assert(fastqfile.myPlusLine == THIRD_PLUS_LINE);
   assert(fastqfile.myQualityString == THIRD_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == FOURTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FOURTH_SEQID);
   assert(fastqfile.myRawSequence == FOURTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FOURTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FOURTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == FIFTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FIFTH_SEQID);
   assert(fastqfile.myRawSequence == FIFTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FIFTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FIFTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SIXTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SIXTH_SEQID);
   assert(fastqfile.myRawSequence == SIXTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SIXTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SIXTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SEVENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SEVENTH_SEQID);
   assert(fastqfile.myRawSequence == SEVENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SEVENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SEVENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == EIGHTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == EIGHTH_SEQID);
   assert(fastqfile.myRawSequence == EIGHTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == EIGHTH_PLUS_LINE);
   assert(fastqfile.myQualityString == EIGHTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == NINTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == NINTH_SEQID);
   assert(fastqfile.myRawSequence == NINTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == NINTH_PLUS_LINE);
   assert(fastqfile.myQualityString == NINTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TENTH_SEQID);
   assert(fastqfile.myRawSequence == TENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == TENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == TENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == ELEVENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == ELEVENTH_SEQID);
   assert(fastqfile.myRawSequence == ELEVENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == ELEVENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == ELEVENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TWELFTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TWELFTH_SEQID);
   assert(fastqfile.myRawSequence == TWELFTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == TWELFTH_PLUS_LINE);
   assert(fastqfile.myQualityString == TWELFTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == THIRTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == THIRTEENTH_SEQID);
   assert(fastqfile.myRawSequence == THIRTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == THIRTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == THIRTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == FOURTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FOURTEENTH_SEQID);
   assert(fastqfile.myRawSequence == FOURTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FOURTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FOURTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == FIFTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FIFTEENTH_SEQID);
   assert(fastqfile.myRawSequence == FIFTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FIFTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FIFTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SIXTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SIXTEENTH_SEQID);
   assert(fastqfile.myRawSequence == SIXTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SIXTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SIXTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == SEVENTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SEVENTEENTH_SEQID);
   assert(fastqfile.myRawSequence == SEVENTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SEVENTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SEVENTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == EIGHTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == EIGHTEENTH_SEQID);
   assert(fastqfile.myRawSequence == EIGHTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == EIGHTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == EIGHTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == NINETEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == NINETEENTH_SEQID);
   assert(fastqfile.myRawSequence == NINETEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == NINETEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == NINETEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TWENTIETH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TWENTIETH_SEQID);
   assert(fastqfile.myRawSequence == TWENTIETH_RAW_SEQ);
   assert(fastqfile.myPlusLine == TWENTIETH_PLUS_LINE);
   assert(fastqfile.myQualityString == TWENTIETH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TWENTY_FIRST_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TWENTY_FIRST_SEQID);
   assert(fastqfile.myRawSequence == TWENTY_FIRST_RAW_SEQ);
   assert(fastqfile.myPlusLine == TWENTY_FIRST_PLUS_LINE);
   assert(fastqfile.myQualityString == TWENTY_FIRST_QUALITY);

   // Close the file, and verify isOpen = false;
   assert(fastqfile.closeFile() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.isOpen() == false);
   

   ////////////////////////////////
   // Repeat test specifying to not check for unique sequence id.
   fastqfile.disableSeqIDCheck();
   assert(fastqfile.isOpen() == false);
   assert(fastqfile.openFile("testFile.txt", BaseAsciiMap::UNKNOWN) == FastQStatus::FASTQ_SUCCESS);

   assert(fastqfile.isOpen() == true);

   assert(fastqfile.getSpaceType() == BaseAsciiMap::UNKNOWN);

   // Read Sequence from test file.
   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == FIRST_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FIRST_SEQID);
   assert(fastqfile.myRawSequence == FIRST_RAW_SEQ);
   assert(fastqfile.myPlusLine == FIRST_PLUS_LINE);
   assert(fastqfile.myQualityString == FIRST_QUALITY);
   assert(fastqfile.getSpaceType() == BaseAsciiMap::BASE_SPACE);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == SECOND_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SECOND_SEQID);
   assert(fastqfile.myRawSequence == SECOND_RAW_SEQ);
   assert(fastqfile.myPlusLine == SECOND_PLUS_LINE);
   assert(fastqfile.myQualityString == SECOND_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == THIRD_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == THIRD_SEQID);
   assert(fastqfile.myRawSequence == THIRD_RAW_SEQ);
   assert(fastqfile.myPlusLine == THIRD_PLUS_LINE);
   assert(fastqfile.myQualityString == THIRD_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == FOURTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FOURTH_SEQID);
   assert(fastqfile.myRawSequence == FOURTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FOURTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FOURTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == FIFTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FIFTH_SEQID);
   assert(fastqfile.myRawSequence == FIFTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FIFTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FIFTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SIXTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SIXTH_SEQID);
   assert(fastqfile.myRawSequence == SIXTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SIXTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SIXTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SEVENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SEVENTH_SEQID);
   assert(fastqfile.myRawSequence == SEVENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SEVENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SEVENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == EIGHTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == EIGHTH_SEQID);
   assert(fastqfile.myRawSequence == EIGHTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == EIGHTH_PLUS_LINE);
   assert(fastqfile.myQualityString == EIGHTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == NINTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == NINTH_SEQID);
   assert(fastqfile.myRawSequence == NINTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == NINTH_PLUS_LINE);
   assert(fastqfile.myQualityString == NINTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == TENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TENTH_SEQID);
   assert(fastqfile.myRawSequence == TENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == TENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == TENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == ELEVENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == ELEVENTH_SEQID);
   assert(fastqfile.myRawSequence == ELEVENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == ELEVENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == ELEVENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TWELFTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TWELFTH_SEQID);
   assert(fastqfile.myRawSequence == TWELFTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == TWELFTH_PLUS_LINE);
   assert(fastqfile.myQualityString == TWELFTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == THIRTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == THIRTEENTH_SEQID);
   assert(fastqfile.myRawSequence == THIRTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == THIRTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == THIRTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == FOURTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FOURTEENTH_SEQID);
   assert(fastqfile.myRawSequence == FOURTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FOURTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FOURTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == FIFTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == FIFTEENTH_SEQID);
   assert(fastqfile.myRawSequence == FIFTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == FIFTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == FIFTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == SIXTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SIXTEENTH_SEQID);
   assert(fastqfile.myRawSequence == SIXTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SIXTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SIXTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == SEVENTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == SEVENTEENTH_SEQID);
   assert(fastqfile.myRawSequence == SEVENTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == SEVENTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == SEVENTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.mySequenceIdLine == EIGHTEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == EIGHTEENTH_SEQID);
   assert(fastqfile.myRawSequence == EIGHTEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == EIGHTEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == EIGHTEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == NINETEENTH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == NINETEENTH_SEQID);
   assert(fastqfile.myRawSequence == NINETEENTH_RAW_SEQ);
   assert(fastqfile.myPlusLine == NINETEENTH_PLUS_LINE);
   assert(fastqfile.myQualityString == NINETEENTH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TWENTIETH_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TWENTIETH_SEQID);
   assert(fastqfile.myRawSequence == TWENTIETH_RAW_SEQ);
   assert(fastqfile.myPlusLine == TWENTIETH_PLUS_LINE);
   assert(fastqfile.myQualityString == TWENTIETH_QUALITY);

   assert(fastqfile.readFastQSequence() == FastQStatus::FASTQ_INVALID);
   assert(fastqfile.mySequenceIdLine == TWENTY_FIRST_SEQID_LINE);
   assert(fastqfile.mySequenceIdentifier == TWENTY_FIRST_SEQID);
   assert(fastqfile.myRawSequence == TWENTY_FIRST_RAW_SEQ);
   assert(fastqfile.myPlusLine == TWENTY_FIRST_PLUS_LINE);
   assert(fastqfile.myQualityString == TWENTY_FIRST_QUALITY);

   // Close the file, and verify isOpen = false;
   assert(fastqfile.closeFile() == FastQStatus::FASTQ_SUCCESS);
   assert(fastqfile.isOpen() == false);
   

}

int main(int argc, char ** argv)
{   
   testReadUnOpenedFile();
   testOpenFile();
   testCloseFile();
   testReadSequence();
}

