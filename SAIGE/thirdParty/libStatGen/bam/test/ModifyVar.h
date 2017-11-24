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

void modifyFirstBase();
void modifyFirstBaseLong();
void testModifyVar();

class modifyVar
{
public:
    void testModifyVar(const char* filename, bool valBufFirst);

private:
    void testModifyReadNameOnlySameLength();
    void testModifyCigarOnlySameLength();
    void testModifySequenceOnlySameLength();
    void testModifyQualityOnlySameLength();
    void testRemoveQuality();
    void testShortenQuality();
    void testLengthenQuality();
   
    void testShortenReadName();
    void testShortenCigar();
    void testShortenSequence();

    void testLengthenReadName();
    void testLengthenCigar();
    void testLengthenSequence();
   
    void testRemoveCigar();
    void testRemoveSequence();
   
    void testLengthenSequenceAndQuality();

    void validate();

    void validateReadName(const bamRecordStruct* recordBuffer);
    void validateCigar(const bamRecordStruct* recordBuffer);
    void validateSequence(const bamRecordStruct* recordBuffer);
    void validateQuality(const bamRecordStruct* recordBuffer);
    void validateTags(const bamRecordStruct* recordBuffer);
   
    void validateReadNameString();
    void validateCigarString();
    void validateSequenceString();
    void validateQualityString();
    void validateTagsString();

    // Open and read the first record.
    void openAndRead1Rec();
    void resetExpected();

    // Variables.
    const char* myFilename;
    bool myValBufFirst;   

    // Rather than passing around all these variables, just store them in the class.
    SamFile samIn;
    SamFileHeader samHeader;
    SamRecord samRecord;
    const bamRecordStruct* recordBuffer;

    // Expected values.
    int expectedCigarBufLen;
    unsigned int expectedCigarBuffer[100];
    unsigned char expectedSequenceBuffer[100];
    int expectedTagsLen;
    unsigned char expectedTagsBuffer[100];

    // Expected values for the strings.
    std::string expectedReadNameString;
    std::string expectedCigarString;
    std::string expectedSequenceString;
    std::string expectedQualityString;
};

