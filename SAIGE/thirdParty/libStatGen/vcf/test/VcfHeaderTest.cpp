/*
 *  Copyright (C) 2011  Regents of the University of Michigan
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

#include "VcfHeaderTest.h"
#include "VcfHeader.h"
#include <assert.h>

//extern const std::string HEADER_LINE="#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003";
//extern const std::string SAMPLES[NUM_SAMPLES] = {"NA00001","NA00002","NA00003"};
//extern const std::string META_LINES[NUM_META_LINES]

void testVcfHeader()
{
    VcfHeader header;

    // Test accessing the header without having read anything.
    assert(header.getNumMetaLines() == 0);
    assert(header.getMetaLine(0) == NULL);
    assert(header.getMetaLine(2) == NULL);
    assert(header.getHeaderLine() == NULL);
    assert(header.getNumSamples() == 0);
    assert(header.getSampleName(2) == NULL);
    assert(header.getSampleName(0) == NULL);
    assert(header.getSampleName(1) == NULL);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == -1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == -1);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == -1);


    IFILE filePtr = NULL; // Input File
    IFILE outputFile = NULL; // Output File.

    // Try reading without opening.
    bool caughtException = false;
    try
    {
        assert(header.read(filePtr) == false);
    }
    catch (std::exception& e) 
    {
        caughtException = true;
    }
    assert(caughtException);

    // Open the file, then read.
    filePtr = ifopen("testFiles/vcfFile.vcf", "r");
    assert(header.read(filePtr));
    assert(header.getNumMetaLines() == NUM_META_LINES);
    assert(header.getMetaLine(0) == META_LINES[0]);
    assert(header.getMetaLine(2) == META_LINES[2]);
    assert(header.getMetaLine(23) == NULL);
    assert(header.getHeaderLine() == HEADER_LINE);
    assert(header.getNumSamples() == NUM_SAMPLES);
    assert(header.getSampleName(2) == SAMPLES[2]);
    assert(header.getSampleName(0) == SAMPLES[0]);
    assert(header.getSampleName(1) == SAMPLES[1]);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == 1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == 2);


    // Reset and verify it is empty.
    header.reset();
    assert(header.getNumMetaLines() == 0);
    assert(header.getMetaLine(0) == NULL);
    assert(header.getMetaLine(2) == NULL);
    assert(header.getHeaderLine() == NULL);
    assert(header.getNumSamples() == 0);
    assert(header.getSampleName(2) == NULL);
    assert(header.getSampleName(0) == NULL);
    assert(header.getSampleName(1) == NULL);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == -1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == -1);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == -1);

    // Close the file and read again.
    ifclose(filePtr);
    filePtr = ifopen("testFiles/vcfFile.vcf", "r");
    assert(header.read(filePtr));
    assert(header.getNumMetaLines() == NUM_META_LINES);
    assert(header.getMetaLine(0) == META_LINES[0]);
    assert(header.getMetaLine(2) == META_LINES[2]);
    assert(header.getMetaLine(23) == NULL);
    assert(header.getHeaderLine() == HEADER_LINE);
    assert(header.getNumSamples() == NUM_SAMPLES);
    assert(header.getSampleName(2) == SAMPLES[2]);
    assert(header.getSampleName(0) == SAMPLES[0]);
    assert(header.getSampleName(1) == SAMPLES[1]);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == 1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == 2);

    // Try writing without opening.
    caughtException = false;
    try
    {
        assert(header.write(outputFile) == false);
    }
    catch (std::exception& e) 
    {
        caughtException = true;
    }
    assert(caughtException);
    caughtException = false;

    // write.
    outputFile = ifopen("results/vcfHeader.vcf", "w");
    assert(header.write(outputFile));

    ////////////////////////////////
    // Test creating a new header starting with the header line.
    VcfHeader newHeader;
    // Header starts empty.
    assert(newHeader.getNumMetaLines() == 0);
    assert(newHeader.getMetaLine(0) == NULL);
    assert(newHeader.getHeaderLine() == NULL);
    assert(newHeader.getNumSamples() == 0);
    assert(newHeader.getSampleName(0) == NULL);
    assert(newHeader.getSampleIndex(SAMPLES[0].c_str()) == -1);

    // Try adding a header line first.
    newHeader.addHeaderLine(HEADER_LINE.c_str());
    assert(newHeader.getNumMetaLines() == 0);
    assert(newHeader.getMetaLine(0) == HEADER_LINE);
    assert(newHeader.getMetaLine(1) == NULL);
    assert(newHeader.getHeaderLine() == HEADER_LINE);
    assert(newHeader.getNumSamples() == NUM_SAMPLES);
    assert(newHeader.getSampleName(2) == SAMPLES[2]);
    assert(newHeader.getSampleName(0) == SAMPLES[0]);
    assert(newHeader.getSampleName(1) == SAMPLES[1]);
    assert(newHeader.getSampleIndex(SAMPLES[1].c_str()) == 1);
    assert(newHeader.getSampleIndex(SAMPLES[0].c_str()) == 0);
    assert(newHeader.getSampleIndex(SAMPLES[2].c_str()) == 2);

    // Add an invalid meta line.
    assert(newHeader.appendMetaLine("# bad line") == false);
    assert(newHeader.getNumMetaLines() == 0);
    assert(newHeader.getMetaLine(0) == HEADER_LINE);
    assert(newHeader.getMetaLine(1) == NULL);


    // Add the meta lines.
    for(int i = 1; i <= NUM_META_LINES; i++)
    {
        assert(newHeader.appendMetaLine(META_LINES[i-1].c_str()));
        assert(newHeader.getNumMetaLines() == i);
        for(int j = 0; j < i; j++)
        {
            assert(newHeader.getMetaLine(j) == META_LINES[j]);
        }
        assert(newHeader.getMetaLine(i) == HEADER_LINE);
        for(int k = i+1; k <= NUM_META_LINES; k++)
        {
            assert(newHeader.getMetaLine(k) == NULL);
        }
    }
    // write.
    outputFile = ifopen("results/vcfHeaderAddedFirst.vcf", "w");
    assert(newHeader.write(outputFile));

    ////////////////////////////////
    // Test creating a new header ending with the header line.
    VcfHeader newHeader2;
    // Header starts empty.
    assert(newHeader2.getNumMetaLines() == 0);
    assert(newHeader2.getMetaLine(0) == NULL);
    assert(newHeader2.getHeaderLine() == NULL);
    assert(newHeader2.getNumSamples() == 0);
    assert(newHeader2.getSampleName(0) == NULL);
    assert(newHeader2.getSampleIndex(SAMPLES[0].c_str()) == -1);

    // Add an invalid meta line.
    assert(newHeader2.appendMetaLine("# bad line") == false);
    assert(newHeader2.getNumMetaLines() == 0);
    assert(newHeader2.getMetaLine(0) == NULL);
    assert(newHeader2.getMetaLine(1) == NULL);


    // Add the meta lines.
    for(int i = 1; i <= NUM_META_LINES; i++)
    {
        assert(newHeader2.appendMetaLine(META_LINES[i-1].c_str()));
        assert(newHeader2.getNumMetaLines() == i);
        for(int j = 0; j < i; j++)
        {
            assert(newHeader2.getMetaLine(j) == META_LINES[j]);
        }
        assert(newHeader2.getMetaLine(i) == NULL);
        for(int k = i+1; k <= NUM_META_LINES; k++)
        {
            assert(newHeader2.getMetaLine(k) == NULL);
        }
    }
 
    // Try adding a header line last.
    newHeader2.addHeaderLine(HEADER_LINE.c_str());
    assert(newHeader2.getNumMetaLines() == NUM_META_LINES);
    for(int i = 0; i < NUM_META_LINES; i++)
    {
        assert(newHeader2.getMetaLine(i) == META_LINES[i]);
    }
    assert(newHeader2.getMetaLine(NUM_META_LINES) == HEADER_LINE);
    assert(newHeader2.getHeaderLine() == HEADER_LINE);
    assert(newHeader2.getNumSamples() == NUM_SAMPLES);
    assert(newHeader2.getSampleName(2) == SAMPLES[2]);
    assert(newHeader2.getSampleName(0) == SAMPLES[0]);
    assert(newHeader2.getSampleName(1) == SAMPLES[1]);
    assert(newHeader2.getSampleIndex(SAMPLES[1].c_str()) == 1);
    assert(newHeader2.getSampleIndex(SAMPLES[0].c_str()) == 0);
    assert(newHeader2.getSampleIndex(SAMPLES[2].c_str()) == 2);

    // write.
    outputFile = ifopen("results/vcfHeaderAddedLast.vcf", "w");
    assert(newHeader2.write(outputFile));


    ////////////////////////////////
    // Test creating a new header adding the header line in the middle.
    VcfHeader newHeader3;
    // Header starts empty.
    assert(newHeader3.getNumMetaLines() == 0);
    assert(newHeader3.getMetaLine(0) == NULL);
    assert(newHeader3.getHeaderLine() == NULL);
    assert(newHeader3.getNumSamples() == 0);
    assert(newHeader3.getSampleName(0) == NULL);
    assert(newHeader3.getSampleIndex(SAMPLES[0].c_str()) == -1);

    // Add an invalid meta line.
    assert(newHeader3.appendMetaLine("# bad line") == false);
    assert(newHeader3.getNumMetaLines() == 0);
    assert(newHeader3.getMetaLine(0) == NULL);
    assert(newHeader3.getMetaLine(1) == NULL);


    // Add the meta lines.
    int subVal = 5;
    for(int i = 1; i <= NUM_META_LINES-subVal; i++)
    {
        assert(newHeader3.appendMetaLine(META_LINES[i-1].c_str()));
        assert(newHeader3.getNumMetaLines() == i);
        for(int j = 0; j < i; j++)
        {
            assert(newHeader3.getMetaLine(j) == META_LINES[j]);
        }
        assert(newHeader3.getMetaLine(i) == NULL);
        for(int k = i+1; k <= NUM_META_LINES; k++)
        {
            assert(newHeader3.getMetaLine(k) == NULL);
        }
    }
 
    // Try adding a header line.
    newHeader3.addHeaderLine(HEADER_LINE.c_str());
    assert(newHeader3.getNumMetaLines() == NUM_META_LINES - subVal);
    for(int i = 0; i < NUM_META_LINES-subVal; i++)
    {
        assert(newHeader3.getMetaLine(i) == META_LINES[i]);
    }
    assert(newHeader3.getMetaLine(NUM_META_LINES-subVal) == HEADER_LINE);
    assert(newHeader3.getHeaderLine() == HEADER_LINE);
    assert(newHeader3.getNumSamples() == NUM_SAMPLES);
    assert(newHeader3.getSampleName(2) == SAMPLES[2]);
    assert(newHeader3.getSampleName(0) == SAMPLES[0]);
    assert(newHeader3.getSampleName(1) == SAMPLES[1]);
    assert(newHeader3.getSampleIndex(SAMPLES[1].c_str()) == 1);
    assert(newHeader3.getSampleIndex(SAMPLES[0].c_str()) == 0);
    assert(newHeader3.getSampleIndex(SAMPLES[2].c_str()) == 2);

    // Add the rest of the meta lines.
    for(int i = NUM_META_LINES - subVal + 1; i <= NUM_META_LINES; i++)
    {
        assert(newHeader3.appendMetaLine(META_LINES[i-1].c_str()));
        assert(newHeader3.getNumMetaLines() == i);
        for(int j = 0; j < i; j++)
        {
            assert(newHeader3.getMetaLine(j) == META_LINES[j]);
        }
        assert(newHeader3.getMetaLine(i) == HEADER_LINE);
        for(int k = i+1; k <= NUM_META_LINES; k++)
        {
            assert(newHeader3.getMetaLine(k) == NULL);
        }
    }
 
    assert(newHeader3.getNumMetaLines() == NUM_META_LINES);
    for(int i = 0; i < NUM_META_LINES; i++)
    {
        assert(newHeader3.getMetaLine(i) == META_LINES[i]);
    }
    assert(newHeader3.getMetaLine(NUM_META_LINES) == HEADER_LINE);
    assert(newHeader3.getHeaderLine() == HEADER_LINE);
    assert(newHeader3.getNumSamples() == NUM_SAMPLES);
    assert(newHeader3.getSampleName(2) == SAMPLES[2]);
    assert(newHeader3.getSampleName(0) == SAMPLES[0]);
    assert(newHeader3.getSampleName(1) == SAMPLES[1]);
    assert(newHeader3.getSampleIndex(SAMPLES[1].c_str()) == 1);
    assert(newHeader3.getSampleIndex(SAMPLES[0].c_str()) == 0);
    assert(newHeader3.getSampleIndex(SAMPLES[2].c_str()) == 2);

    // write.
    outputFile = ifopen("results/vcfHeaderAddedMiddle.vcf", "w");
    assert(newHeader3.write(outputFile));

}
