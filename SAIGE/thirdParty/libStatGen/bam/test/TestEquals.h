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

#include "SamFile.h"
void testSeqEquals();

class EqualsTest
{
public:
    enum FileType{SAM, BAM};
    static void testEq(FileType inputType);
    
private:
    static void reset();

    static void validateEqRead(SamRecord& samRecord, 
                               int readIndex,
                               const char* actualExpectedSequence);
    static void validateEqReadBuffer(SamRecord& samRecord, 
                                     const char* expectedSequence);

    static const char* READ_NAMES[];
    static const char* READ_SEQS_BASES[];
    static const char* READ_SEQS_EQUALS[];
    static const char* READ_SEQS_MIXED[];

    static const char* expectedReferenceName;
    static const char* expectedMateReferenceName;
    static const char* expectedMateReferenceNameOrEqual;
    static const char* expectedCigar;
    static const char* expectedQuality;

    static std::vector<unsigned int> expectedCigarHex;

    static int expected0BasedAlignmentEnd;
    static int expected1BasedAlignmentEnd;
    static int expectedAlignmentLength;
    static int expected0BasedUnclippedStart;
    static int expected1BasedUnclippedStart;
    static int expected0BasedUnclippedEnd;
    static int expected1BasedUnclippedEnd;

    static bamRecordStruct expectedRecord;
};
