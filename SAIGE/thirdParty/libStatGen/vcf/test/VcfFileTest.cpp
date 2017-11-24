/*
 *  Copyright (C) 2011-2013  Regents of the University of Michigan
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

#include "VcfFileTest.h"
#include "VcfFileReader.h"
#include "VcfFileWriter.h"
#include "VcfHeaderTest.h"
#include <assert.h>

const std::string HEADER_LINE_SUBSET1="#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002";
const std::string HEADER_LINE_SUBSET2="#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00002	NA00003";
const int NUM_SAMPLES_SUBSET1 = 2;
const int NUM_SAMPLES_SUBSET2 = 2;
const std::string HEADER_LINE_EXCLUDE_SUBSET1="#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00003";
const int NUM_SAMPLES_EXCLUDE2 = 2;
const std::string HEADER_LINE_EXCLUDE2="#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00003";

void testVcfFile()
{
    testVcfReadFile();
    testVcfWriteFile();

    testVcfReadSection();
    testVcfReadSectionNoIndex();
    testVcfReadSectionBadIndex();
}


void testVcfReadFile()
{
    //    VcfFileHeader header;

    // Test open for read via the constructor with return.
    VcfFileReader reader;
    VcfHeader header;
    VcfRecord record;
    // Try reading without opening.
    bool caughtException = false;
    try
    {
        assert(reader.readRecord(record) == false);
    }
    catch (std::exception& e) 
    {
        caughtException = true;
    }

    assert(caughtException);


    // Try opening a file that doesn't exist.
    caughtException = false;
    try
    {
        assert(reader.open("fileDoesNotExist.txt", header) == false);
    }
    catch (std::exception& e) 
    {
        caughtException = true;
    }
    assert(caughtException);
// "testFiles/testVcf.vcf");
//     assert(vcfInConstructorReadDefault.WriteHeader(header) == false);
//     assert(vcfInConstructorReadDefault.ReadHeader(header) == true);

//     // Test open for write via the constructor.
//     VcfFile vcfInConstructorWrite("results/newWrite.vcf", VcfFile::WRITE,
//                                   ErrorHandler::RETURN);
//     assert(vcfInConstructorWrite.ReadHeader(header) == false);
//     assert(vcfInConstructorWrite.WriteHeader(header) == true);

//     // Test open for read via the constructor
//     VcfFile vcfInConstructorRead("testFiles/testVcf.vcf", VcfFile::READ);
//     bool caughtException = false;
//     try
//     {
//         assert(vcfInConstructorRead.WriteHeader(header) == false);
//     }
//     catch (std::exception& e) 
//     {
//         caughtException = true;
//     }
//     assert(caughtException);
//     assert(vcfInConstructorRead.ReadHeader(header) == true);

//     // Test open for write via child class.
//     VcfFileWriter vcfWriteConstructor("results/newWrite1.vcf");
//     caughtException = false;
//     try
//     {
//         assert(vcfWriteConstructor.ReadHeader(header) == false);
//     }
//     catch (std::exception& e) 
//     {
//         caughtException = true;
//     }
//     assert(caughtException);
//     assert(vcfWriteConstructor.WriteHeader(header) == true);

//     // Test open for read via child class.
//     VcfFileReader vcfReadConstructor("testFiles/testVcf.vcf");
//     caughtException = false;
//     try
//     {
//         assert(vcfReadConstructor.WriteHeader(header) == false);
//     }
//     catch (std::exception& e) 
//     {
//         caughtException = true;
//     }
//     assert(caughtException);
//     assert(vcfReadConstructor.ReadHeader(header) == true);


    ////////////////////////////////
    // Test the subset logic.
    VcfRecordGenotype* sampleInfo = NULL;

    reader.open("testFiles/vcfFile.vcf", header, "testFiles/subset1.txt",
                NULL, NULL, ";");

    assert(header.getHeaderLine() == HEADER_LINE_SUBSET1);
    assert(header.getNumSamples() == NUM_SAMPLES_SUBSET1);
    assert(header.getSampleName(2) == NULL);
    assert(header.getSampleName(0) == SAMPLES[0]);
    assert(header.getSampleName(1) == SAMPLES[1]);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == 1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == -1);

    assert(strcmp(reader.getFileName(), "testFiles/vcfFile.vcf") == 0);
 
    // Read the records to make sure they were subset.
    assert(reader.readRecord(record));
    assert(record.getGT(0,0) == 0);
    assert(record.getGT(1,1) == 0);
    assert(record.getGT(1,0) == 1);
    assert(record.getGT(0,1) == 0);
    assert(record.getGT(1,2) == VcfGenotypeSample::INVALID_GT);
    assert(record.getGT(2,0) == VcfGenotypeSample::INVALID_GT);
    assert(strcmp(record.getAlleles(0), "G") == 0);
    assert(strcmp(record.getAlleles(1), "A") == 0);
    assert(record.getIntAllele(0) == 3);
    assert(record.getIntAllele(1) == 1);
    assert(record.getNumAlts() == 1);
    try
    {
        caughtException = false;
        assert(record.getIntAllele(2) == 0);
    }
    catch (std::exception& e)
    {
        caughtException = true;
    }
    try
    {
        caughtException = false;
        assert(record.getAlleles(2) == NULL);
    }
    catch (std::exception& e)
    {
        caughtException = true;
    }
    assert(caughtException);
    assert(record.allPhased() == true);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "1|0");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 1);

    assert(reader.readRecord(record));
    assert(record.getGT(0,0) == 0);
    assert(record.getGT(1,1) == 1);
    assert(record.getGT(1,0) == 0);
    assert(record.getGT(0,1) == 0);
    assert(record.getGT(1,2) == VcfGenotypeSample::INVALID_GT);
    assert(record.getGT(2,0) == VcfGenotypeSample::INVALID_GT);
    assert(record.getNumAlts() == 1);
    assert(strcmp(record.getAlleles(0), "T") == 0);
    assert(strcmp(record.getAlleles(1), "A") == 0);
    try
    {
        caughtException = false;
        assert(record.getAlleles(2) == NULL);
    }
    catch (std::exception& e)
    {
        caughtException = true;
    }
    assert(caughtException);
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0/0");
    assert(*(sampleInfo->getString("GT", 1)) == "0|1");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);
    assert(record.passedAllFilters() == false);
    assert(record.getNumAlts() == 1);

    assert(reader.readRecord(record));
    assert(record.getGT(0,0) == 1);
    assert(record.getGT(1,1) == 1);
    assert(record.getGT(1,0) == 2);
    assert(record.getGT(0,1) == 2);
    assert(record.getGT(1,2) == VcfGenotypeSample::INVALID_GT);
    assert(record.getGT(2,0) == VcfGenotypeSample::INVALID_GT);
    assert(strcmp(record.getAlleles(0), "A") == 0);
    assert(strcmp(record.getAlleles(1), "G") == 0);
    assert(strcmp(record.getAlleles(2), "T") == 0);
    assert(record.getIntAllele(2) == 4);
    assert(record.getNumAlts() == 2);
    try
    {
        caughtException = false;
        assert(record.getAlleles(3) == NULL);
    }
    catch (std::exception& e)
    {
        caughtException = true;
    }
    assert(caughtException);
    assert(record.allPhased() == true);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "1|2");
    assert(*(sampleInfo->getString("GT", 1)) == "2|1");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 2);

    assert(reader.readRecord(record));
    assert(strcmp(record.getAlleles(0), "T") == 0);
    try
    {
        caughtException = false;
        assert(record.getAlleles(1) == NULL);
    }
    catch (std::exception& e)
    {
        caughtException = true;
    }
    assert(caughtException);
    assert(record.getNumAlts() == 0);
    assert(record.allPhased() == true);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "0|0");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 0);

    assert(reader.readRecord(record));
    assert(strcmp(record.getAlleles(0), "GTC") == 0);
    assert(record.getIntAllele(0) == 3);
    assert(strcmp(record.getAlleles(1), "G") == 0);
    assert(record.getIntAllele(1) == 3);
    assert(record.getNumAlts() == 2);
    assert(strcmp(record.getAlleles(2), "GTCT") == 0);
    assert(record.getIntAllele(2) == 3);
    try
    {
        caughtException = false;
        assert(record.getAlleles(3) == NULL);
    }
    catch (std::exception& e)
    {
        caughtException = true;
    }
    assert(caughtException);
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0/1");
    assert(*(sampleInfo->getString("GT", 1)) == "0/2");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 2);

    assert(reader.readRecord(record));
    assert(strcmp(record.getAlleles(0), "GTC") == 0);
    assert(record.getNumAlts() == 2);
    assert(strcmp(record.getAlleles(1), "G") == 0);
    assert(strcmp(record.getAlleles(2), "GTCT") == 0);
    try
    {
        caughtException = false;
        assert(record.getAlleles(3) == NULL);
    }
    catch (std::exception& e)
    {
        caughtException = true;
    }
    assert(caughtException);
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == false);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(sampleInfo->getString("GT", 0) == NULL);
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 2);

    assert(reader.readRecord(record));
    assert(record.getGT(0,0) == 0);
    assert(record.getGT(1,1) == VcfGenotypeSample::MISSING_GT);
    assert(record.getGT(1,0) == 0);
    assert(record.getGT(0,1) == 1);
    assert(record.getGT(1,2) == VcfGenotypeSample::INVALID_GT);
    assert(record.getGT(2,0) == VcfGenotypeSample::INVALID_GT);
    assert(record.getNumAlts() == 1);
    assert(strcmp(record.getAlleles(0), "GTC") == 0);
    assert(strcmp(record.getAlleles(1), "G") == 0);
    try
    {
        caughtException = false;
        assert(record.getAlleles(2) == NULL);
    }
    catch (std::exception& e)
    {
        caughtException = true;
    }
    assert(caughtException);
    assert(record.allPhased() == true);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == false);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|1");
    assert(*(sampleInfo->getString("GT", 1)) == "0|.");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 1);

    assert(reader.readRecord(record) == false);

    assert(reader.getNumKeptRecords() == 7);
    assert(reader.getNumRecords() == 7);

    reader.close();

    //////////////////////////
    // Subset with a different file.
    reader.open("testFiles/vcfFile.vcf", header, "testFiles/subset2.txt", 
                NULL, NULL);
    
    assert(header.getHeaderLine() == HEADER_LINE_SUBSET2);
    assert(header.getNumSamples() == NUM_SAMPLES_SUBSET2);
    assert(header.getSampleName(2) == NULL);
    assert(header.getSampleName(0) == SAMPLES[1]);
    assert(header.getSampleName(1) == SAMPLES[2]);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == 0);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == -1);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == 1);

    // Read the records to make sure they were subset.
    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "1|0");
    assert(*(sampleInfo->getString("GT", 1)) == "1/1");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|1");
    assert(*(sampleInfo->getString("GT", 1)) == "0/0");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "2|1");
    assert(*(sampleInfo->getString("GT", 1)) == "2/2");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "0/0");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0/2");
    assert(*(sampleInfo->getString("GT", 1)) == "1/1");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == false);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(sampleInfo->getString("GT", 0) == NULL);
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == false);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == true);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == false);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|.");
    assert(*(sampleInfo->getString("GT", 1)) == "1|1");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == true);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == false);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record) == false);

    assert(reader.getNumKeptRecords() == 7);
    assert(reader.getNumRecords() == 7);


    //////////////////////////
    // Subset using an exclude file
    reader.open("testFiles/vcfFile.vcf", header, NULL, NULL,
                "testFiles/subset1.txt", ";");
    
    assert(header.getHeaderLine() == HEADER_LINE_EXCLUDE_SUBSET1);
    assert(header.getNumSamples() == NUM_SAMPLES - NUM_SAMPLES_SUBSET1);
    assert(header.getSampleName(2) == NULL);
    assert(header.getSampleName(0) == SAMPLES[2]);
    assert(header.getSampleName(1) == NULL);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == -1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == -1);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == 0);

    // Read the records to make sure they were subset.
    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "1/1");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "0/0");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "2/2");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "0/0");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "1/1");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == false);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(sampleInfo->getString("GT", 0) == NULL);
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == false);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == true);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "1|1");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == true);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record) == false);

    assert(reader.getNumKeptRecords() == 7);
    assert(reader.getNumRecords() == 7);


    //////////////////////////
    // Subset with a different exclude.
    reader.open("testFiles/vcfFile.vcf", header, NULL, NULL,
                "testFiles/exclude2.txt");
    
    assert(header.getHeaderLine() == HEADER_LINE_EXCLUDE2);
    assert(header.getNumSamples() == NUM_SAMPLES_EXCLUDE2);
    assert(header.getSampleName(2) == NULL);
    assert(header.getSampleName(0) == SAMPLES[0]);
    assert(header.getSampleName(1) == SAMPLES[2]);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == -1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == 1);

    // Read the records to make sure they were subset.
    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "1/1");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0/0");
    assert(*(sampleInfo->getString("GT", 1)) == "0/0");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "1|2");
    assert(*(sampleInfo->getString("GT", 1)) == "2/2");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "0/0");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0/1");
    assert(*(sampleInfo->getString("GT", 1)) == "1/1");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == false);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(sampleInfo->getString("GT", 0) == NULL);
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == false);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == true);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|1");
    assert(*(sampleInfo->getString("GT", 1)) == "1|1");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == true);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record) == false);

    assert(reader.getNumKeptRecords() == 7);
    assert(reader.getNumRecords() == 7);


    //////////////////////////
    // Subset with an exclude sample.
    reader.open("testFiles/vcfFile.vcf", header, NULL, "NA00002",
                NULL);
    
    assert(header.getHeaderLine() == HEADER_LINE_EXCLUDE2);
    assert(header.getNumSamples() == NUM_SAMPLES_EXCLUDE2);
    assert(header.getSampleName(2) == NULL);
    assert(header.getSampleName(0) == SAMPLES[0]);
    assert(header.getSampleName(1) == SAMPLES[2]);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == -1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == 1);

    // Read the records to make sure they were subset.
    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "1/1");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0/0");
    assert(*(sampleInfo->getString("GT", 1)) == "0/0");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "1|2");
    assert(*(sampleInfo->getString("GT", 1)) == "2/2");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "0/0");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0/1");
    assert(*(sampleInfo->getString("GT", 1)) == "1/1");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == false);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(sampleInfo->getString("GT", 0) == NULL);
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == false);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == true);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|1");
    assert(*(sampleInfo->getString("GT", 1)) == "1|1");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == true);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record) == false);

    assert(reader.getNumKeptRecords() == 7);
    assert(reader.getNumRecords() == 7);


    //////////////////////////
    // Subset using an exclude file and exclude sample.
    reader.open("testFiles/vcfFile.vcf", header, NULL, "NA00001",
                "testFiles/exclude2.txt");
    
    assert(header.getHeaderLine() == HEADER_LINE_EXCLUDE_SUBSET1);
    assert(header.getNumSamples() == NUM_SAMPLES - NUM_SAMPLES_SUBSET1);
    assert(header.getSampleName(2) == NULL);
    assert(header.getSampleName(0) == SAMPLES[2]);
    assert(header.getSampleName(1) == NULL);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == -1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == -1);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == 0);

    // Read the records to make sure they were subset.
    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "1/1");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "0/0");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "2/2");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "0/0");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "1/1");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == false);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(sampleInfo->getString("GT", 0) == NULL);
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == false);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(record.allPhased() == true);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "1|1");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == true);
    assert(sampleInfo->allUnphased() == false);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record) == false);

    assert(reader.getNumKeptRecords() == 7);
    assert(reader.getNumRecords() == 7);


    //////////////////////////
    // Subset using an exclude file and exclude sample.
    // Add variant discard
    reader.open("testFiles/vcfFile.vcf", header, NULL, "NA00001",
                "testFiles/exclude2.txt");
    reader.setExcludeIDs("testFiles/excludeIDs.txt");

    assert(header.getHeaderLine() == HEADER_LINE_EXCLUDE_SUBSET1);
    assert(header.getNumSamples() == NUM_SAMPLES - NUM_SAMPLES_SUBSET1);
    assert(header.getSampleName(2) == NULL);
    assert(header.getSampleName(0) == SAMPLES[2]);
    assert(header.getSampleName(1) == NULL);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == -1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == -1);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == 0);

    // Read the records to make sure they were subset.
    assert(reader.readRecord(record));
    assert(strcmp(record.getIDStr(), ".") == 0);
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "0/0");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(strcmp(record.getIDStr(), "rs6040355") == 0);
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "2/2");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(strcmp(record.getIDStr(), ".") == 0);
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "0/0");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record) == false);

    assert(reader.getNumKeptRecords() == 3);
    assert(reader.getNumRecords() == 7);

    //////////////////////////
    // Subset using an exclude file and exclude sample.
    // Add variant discard
    reader.open("testFiles/vcfFile.vcf", header, NULL, "NA00001",
                "testFiles/exclude2.txt");
    reader.setIncludeIDs("testFiles/includeIDs.txt");

    assert(header.getHeaderLine() == HEADER_LINE_EXCLUDE_SUBSET1);
    assert(header.getNumSamples() == NUM_SAMPLES - NUM_SAMPLES_SUBSET1);
    assert(header.getSampleName(2) == NULL);
    assert(header.getSampleName(0) == SAMPLES[2]);
    assert(header.getSampleName(1) == NULL);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == -1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == -1);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == 0);

    // Read the records to make sure they were subset.
    assert(reader.readRecord(record));
    assert(strcmp(record.getIDStr(), ".") == 0);
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "0/0");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(strcmp(record.getIDStr(), "rs6040355") == 0);
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "2/2");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    assert(strcmp(record.getIDStr(), ".") == 0);
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 1);
    assert(*(sampleInfo->getString("GT", 0)) == "0/0");
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->allPhased() == false);
    assert(sampleInfo->allUnphased() == true);
    assert(sampleInfo->hasAllGenotypeAlleles() == true);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record) == false);

    assert(reader.getNumKeptRecords() == 3);
    assert(reader.getNumRecords() == 7);

    //////////////////////////
    // Add in discarding non-phased.
    reader.setDiscardRules(VcfFileReader::DISCARD_NON_PHASED);
    reader.open("testFiles/vcfFile.vcf", header, "testFiles/subset1.txt", 
                NULL, NULL, ";");

    assert(header.getHeaderLine() == HEADER_LINE_SUBSET1);
    assert(header.getNumSamples() == NUM_SAMPLES_SUBSET1);
    assert(header.getSampleName(2) == NULL);
    assert(header.getSampleName(0) == SAMPLES[0]);
    assert(header.getSampleName(1) == SAMPLES[1]);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == 1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == -1);
 
    // Read the records to make sure they were subset.
    assert(reader.readRecord(record));
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "1|0");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "1|2");
    assert(*(sampleInfo->getString("GT", 1)) == "2|1");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "0|0");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|1");
    assert(*(sampleInfo->getString("GT", 1)) == "0|.");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record) == false);

    assert(reader.getNumKeptRecords() == 4);
    assert(reader.getNumRecords() == 7);

    reader.close();

    //////////////////////////
    // Discard missing GTs.
    reader.setDiscardRules(VcfFileReader::DISCARD_MISSING_GT);
    reader.open("testFiles/vcfFile.vcf", header, "testFiles/subset1.txt", 
                NULL, NULL, ";");

    assert(header.getHeaderLine() == HEADER_LINE_SUBSET1);
    assert(header.getNumSamples() == NUM_SAMPLES_SUBSET1);
    assert(header.getSampleName(2) == NULL);
    assert(header.getSampleName(0) == SAMPLES[0]);
    assert(header.getSampleName(1) == SAMPLES[1]);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == 1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == -1);
 
    // Read the records to make sure they were subset.
    assert(reader.readRecord(record));
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "1|0");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0/0");
    assert(*(sampleInfo->getString("GT", 1)) == "0|1");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "1|2");
    assert(*(sampleInfo->getString("GT", 1)) == "2|1");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "0|0");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0/1");
    assert(*(sampleInfo->getString("GT", 1)) == "0/2");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record) == false);

    assert(reader.getNumKeptRecords() == 5);
    assert(reader.getNumRecords() == 7);

    //////////////////////////
    // Discard missing GTs & non-Phased.
    reader.setDiscardRules(VcfFileReader::DISCARD_MISSING_GT | 
                           VcfFileReader::DISCARD_NON_PHASED);
    reader.open("testFiles/vcfFile.vcf", header, "testFiles/subset1.txt", 
                NULL, NULL, ";");

    assert(header.getHeaderLine() == HEADER_LINE_SUBSET1);
    assert(header.getNumSamples() == NUM_SAMPLES_SUBSET1);
    assert(header.getSampleName(2) == NULL);
    assert(header.getSampleName(0) == SAMPLES[0]);
    assert(header.getSampleName(1) == SAMPLES[1]);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == 1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == -1);
 
    // Read the records to make sure they were subset.
    assert(reader.readRecord(record));
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "1|0");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "1|2");
    assert(*(sampleInfo->getString("GT", 1)) == "2|1");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record));
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "0|0");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);

    assert(reader.readRecord(record) == false);

    assert(reader.getNumKeptRecords() == 3);
    assert(reader.getNumRecords() == 7);

    reader.close();


    ////////////////////////////////
    // Test Discarding filtered without subsetting.
    reader.open("testFiles/vcfFile.vcf", header);

    reader.setDiscardRules(VcfFileReader::DISCARD_FILTERED); 

    assert(header.getHeaderLine() == HEADER_LINE);
    assert(header.getNumSamples() == NUM_SAMPLES);
    assert(header.getSampleName(2) == SAMPLES[2]);
    assert(header.getSampleName(0) == SAMPLES[0]);
    assert(header.getSampleName(1) == SAMPLES[1]);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == 1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == 2);
 
    // Read the records.
    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 3);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "1|0");
    assert(*(sampleInfo->getString("GT", 2)) == "1/1");
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == true);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 1);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 3);
    assert(*(sampleInfo->getString("GT", 0)) == "1|2");
    assert(*(sampleInfo->getString("GT", 1)) == "2|1");
    assert(*(sampleInfo->getString("GT", 2)) == "2/2");
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == true);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 2);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 3);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "0|0");
    assert(*(sampleInfo->getString("GT", 2)) == "0/0");
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == true);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 0);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == true);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 3);
    assert(*(sampleInfo->getString("GT", 0)) == "0/1");
    assert(*(sampleInfo->getString("GT", 1)) == "0/2");
    assert(*(sampleInfo->getString("GT", 2)) == "1/1");
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == true);
    assert(sampleInfo->isUnphased(2) == true);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 2);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == false);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 3);
    assert(sampleInfo->getString("GT", 0) == NULL);
    assert(sampleInfo->getString("GT", 1) == NULL);
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == false);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 2);

    assert(reader.readRecord(record));
    assert(record.allPhased() == true);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == false);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 3);
    assert(*(sampleInfo->getString("GT", 0)) == "0|1");
    assert(*(sampleInfo->getString("GT", 1)) == "0|.");
    assert(*(sampleInfo->getString("GT", 2)) == "1|1");
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == true);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 1);

    assert(reader.readRecord(record) == false);

    assert(reader.getNumKeptRecords() == 6);
    assert(reader.getNumRecords() == 7);

    reader.close();

    ////////////////////////////////
    // Test Discarding multiple Alts without subsetting.
    reader.open("testFiles/vcfFile.vcf", header);

    reader.setDiscardRules(VcfFileReader::DISCARD_MULTIPLE_ALTS); 

    assert(header.getHeaderLine() == HEADER_LINE);
    assert(header.getNumSamples() == NUM_SAMPLES);
    assert(header.getSampleName(2) == SAMPLES[2]);
    assert(header.getSampleName(0) == SAMPLES[0]);
    assert(header.getSampleName(1) == SAMPLES[1]);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == 1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == 2);
 
    // Read the records.
    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 3);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "1|0");
    assert(*(sampleInfo->getString("GT", 2)) == "1/1");
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == true);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 1);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 3);
    assert(*(sampleInfo->getString("GT", 0)) == "0/0");
    assert(*(sampleInfo->getString("GT", 1)) == "0|1");
    assert(*(sampleInfo->getString("GT", 2)) == "0/0");
    assert(sampleInfo->isPhased(0) == false);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == true);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == true);
    assert(record.passedAllFilters() == false);
    assert(record.getNumAlts() == 1);

    assert(reader.readRecord(record));
    assert(record.allPhased() == false);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 3);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "0|0");
    assert(*(sampleInfo->getString("GT", 2)) == "0/0");
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == true);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 0);

    assert(reader.readRecord(record));
    assert(record.getGT(0,0) == 0);
    assert(record.getGT(1,1) == VcfGenotypeSample::MISSING_GT);
    assert(record.getGT(1,0) == 0);
    assert(record.getGT(0,1) == 1);
    assert(record.getGT(2,0) == 1);
    assert(record.getGT(2,1) == 1);
    assert(record.getGT(1,2) == VcfGenotypeSample::INVALID_GT);
    assert(record.getGT(3,0) == VcfGenotypeSample::INVALID_GT);
    assert(record.allPhased() == true);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == false);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 3);
    assert(*(sampleInfo->getString("GT", 0)) == "0|1");
    assert(*(sampleInfo->getString("GT", 1)) == "0|.");
    assert(*(sampleInfo->getString("GT", 2)) == "1|1");
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == true);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 1);

    assert(reader.readRecord(record) == false);

    assert(reader.getNumKeptRecords() == 4);
    assert(reader.getNumRecords() == 7);

    reader.close();

    ////////////////////////////////
    // Test subsetting and discarding multiple Alts, filter failures,
    // non-phased, and missing genotypes.
    reader.open("testFiles/vcfFile.vcf", header);

    reader.setDiscardRules(VcfFileReader::DISCARD_MULTIPLE_ALTS|
                           VcfFileReader::DISCARD_FILTERED |
                           VcfFileReader::DISCARD_MISSING_GT | 
                           VcfFileReader::DISCARD_NON_PHASED);
    reader.open("testFiles/vcfFile.vcf", header, "testFiles/subset1.txt", 
                NULL, NULL, ";");

    assert(header.getHeaderLine() == HEADER_LINE_SUBSET1);
    assert(header.getNumSamples() == NUM_SAMPLES_SUBSET1);
    assert(header.getSampleName(2) == NULL);
    assert(header.getSampleName(0) == SAMPLES[0]);
    assert(header.getSampleName(1) == SAMPLES[1]);
    assert(header.getSampleIndex(SAMPLES[1].c_str()) == 1);
    assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
    assert(header.getSampleIndex(SAMPLES[2].c_str()) == -1);
 
    // Read the records to make sure they were subset.
    assert(reader.readRecord(record));
    assert(record.allPhased() == true);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "1|0");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 1);

    assert(reader.readRecord(record));
    assert(record.allPhased() == true);
    assert(record.allUnphased() == false);
    assert(record.hasAllGenotypeAlleles() == true);
    sampleInfo = &(record.getGenotypeInfo());
    assert(sampleInfo->getNumSamples() == 2);
    assert(*(sampleInfo->getString("GT", 0)) == "0|0");
    assert(*(sampleInfo->getString("GT", 1)) == "0|0");
    assert(sampleInfo->getString("GT", 2) == NULL);
    assert(sampleInfo->isPhased(0) == true);
    assert(sampleInfo->isPhased(1) == true);
    assert(sampleInfo->isPhased(2) == false);
    assert(sampleInfo->isUnphased(0) == false);
    assert(sampleInfo->isUnphased(1) == false);
    assert(sampleInfo->isUnphased(2) == false);
    assert(record.passedAllFilters() == true);
    assert(record.getNumAlts() == 0);

    assert(reader.readRecord(record) == false);

    assert(reader.getNumKeptRecords() == 2);
    assert(reader.getNumRecords() == 7);

    reader.close();

   //////////////////////////
   // Discard missing GTs & non-Phased and filtering 
   // AND discard without at least 2 alternates with no additional subsetting.
   reader.setDiscardRules(VcfFileReader::DISCARD_MISSING_GT | 
                          VcfFileReader::DISCARD_NON_PHASED);
   reader.addDiscardMinAltAlleleCount(2, NULL);
   reader.open("testFiles/vcfFile.vcf", header, "testFiles/subset1.txt", 
               NULL, NULL, ";");
   
   assert(header.getHeaderLine() == HEADER_LINE_SUBSET1);
   assert(header.getNumSamples() == NUM_SAMPLES_SUBSET1);
   assert(header.getSampleName(2) == NULL);
   assert(header.getSampleName(0) == SAMPLES[0]);
   assert(header.getSampleName(1) == SAMPLES[1]);
   assert(header.getSampleIndex(SAMPLES[1].c_str()) == 1);
   assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
   assert(header.getSampleIndex(SAMPLES[2].c_str()) == -1);
   
   // Read the records to make sure they were subset.
   assert(reader.readRecord(record));
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 2);
   assert(*(sampleInfo->getString("GT", 0)) == "1|2");
   assert(*(sampleInfo->getString("GT", 1)) == "2|1");
   assert(sampleInfo->getString("GT", 2) == NULL);
   assert(sampleInfo->isPhased(0) == true);
   assert(sampleInfo->isPhased(1) == true);
   assert(sampleInfo->isPhased(2) == false);
   assert(sampleInfo->isUnphased(0) == false);
   assert(sampleInfo->isUnphased(1) == false);
   assert(sampleInfo->isUnphased(2) == false);

   assert(reader.readRecord(record) == false);
   
   assert(reader.getNumKeptRecords() == 1);
   assert(reader.getNumRecords() == 7);
   
   reader.close();

   ////////////////////////////////
   // Test Discarding multiple Alts without subsetting
   // and discard any without at least 3 alts.
   reader.open("testFiles/vcfFile.vcf", header);
   reader.addDiscardMinAltAlleleCount(3, NULL);
   reader.setDiscardRules(VcfFileReader::DISCARD_MULTIPLE_ALTS); 

   assert(header.getHeaderLine() == HEADER_LINE);
   assert(header.getNumSamples() == NUM_SAMPLES);
   assert(header.getSampleName(2) == SAMPLES[2]);
   assert(header.getSampleName(0) == SAMPLES[0]);
   assert(header.getSampleName(1) == SAMPLES[1]);
   assert(header.getSampleIndex(SAMPLES[1].c_str()) == 1);
   assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
   assert(header.getSampleIndex(SAMPLES[2].c_str()) == 2);
 
   // Read the records.
   assert(reader.readRecord(record));
   assert(record.allPhased() == false);
   assert(record.allUnphased() == false);
   assert(record.hasAllGenotypeAlleles() == true);
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 3);
   assert(*(sampleInfo->getString("GT", 0)) == "0|0");
   assert(*(sampleInfo->getString("GT", 1)) == "1|0");
   assert(*(sampleInfo->getString("GT", 2)) == "1/1");
   assert(sampleInfo->isPhased(0) == true);
   assert(sampleInfo->isPhased(1) == true);
   assert(sampleInfo->isPhased(2) == false);
   assert(sampleInfo->isUnphased(0) == false);
   assert(sampleInfo->isUnphased(1) == false);
   assert(sampleInfo->isUnphased(2) == true);
   assert(record.passedAllFilters() == true);
   assert(record.getNumAlts() == 1);

   assert(reader.readRecord(record));
   assert(record.getGT(0,0) == 0);
   assert(record.getGT(1,1) == VcfGenotypeSample::MISSING_GT);
   assert(record.getGT(1,0) == 0);
   assert(record.getGT(0,1) == 1);
   assert(record.getGT(2,0) == 1);
   assert(record.getGT(2,1) == 1);
   assert(record.getGT(1,2) == VcfGenotypeSample::INVALID_GT);
   assert(record.getGT(3,0) == VcfGenotypeSample::INVALID_GT);
   assert(record.allPhased() == true);
   assert(record.allUnphased() == false);
   assert(record.hasAllGenotypeAlleles() == false);
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 3);
   assert(*(sampleInfo->getString("GT", 0)) == "0|1");
   assert(*(sampleInfo->getString("GT", 1)) == "0|.");
   assert(*(sampleInfo->getString("GT", 2)) == "1|1");
   assert(sampleInfo->isPhased(0) == true);
   assert(sampleInfo->isPhased(1) == true);
   assert(sampleInfo->isPhased(2) == true);
   assert(sampleInfo->isUnphased(0) == false);
   assert(sampleInfo->isUnphased(1) == false);
   assert(sampleInfo->isUnphased(2) == false);
   assert(record.passedAllFilters() == true);
   assert(record.getNumAlts() == 1);

   assert(reader.readRecord(record) == false);

   assert(reader.getNumKeptRecords() == 2);
   assert(reader.getNumRecords() == 7);

   reader.close();

   ////////////////////////////////
   // Test Discarding multiple Alts without subsetting
   // and discard any without at least 3 alts and only samples 1 & 2.
   reader.open("testFiles/vcfFile.vcf", header);
   VcfSubsetSamples minAltAlleleSubset;
   minAltAlleleSubset.init(header, true);
   minAltAlleleSubset.addExcludeSample("NA00002");
   reader.addDiscardMinAltAlleleCount(3, &minAltAlleleSubset);
   reader.setDiscardRules(VcfFileReader::DISCARD_MULTIPLE_ALTS); 

   assert(header.getHeaderLine() == HEADER_LINE);
   assert(header.getNumSamples() == 3);
   assert(header.getSampleName(0) == SAMPLES[0]);
   assert(header.getSampleName(2) == SAMPLES[2]);
   assert(header.getSampleName(1) == SAMPLES[1]);
   assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
   assert(header.getSampleIndex(SAMPLES[1].c_str()) == 1);
   assert(header.getSampleIndex(SAMPLES[2].c_str()) == 2);
 
   // Read the records.
   assert(reader.readRecord(record));
   assert(record.getGT(0,0) == 0);
   assert(record.getGT(1,1) == VcfGenotypeSample::MISSING_GT);
   assert(record.getGT(1,0) == 0);
   assert(record.getGT(0,1) == 1);
   assert(record.getGT(2,0) == 1);
   assert(record.getGT(2,1) == 1);
   assert(record.getGT(1,2) == VcfGenotypeSample::INVALID_GT);
   assert(record.getGT(3,0) == VcfGenotypeSample::INVALID_GT);
   assert(record.allPhased() == true);
   assert(record.allUnphased() == false);
   assert(record.hasAllGenotypeAlleles() == false);
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 3);
   assert(*(sampleInfo->getString("GT", 0)) == "0|1");
   assert(*(sampleInfo->getString("GT", 1)) == "0|.");
   assert(*(sampleInfo->getString("GT", 2)) == "1|1");
   assert(sampleInfo->isPhased(0) == true);
   assert(sampleInfo->isPhased(1) == true);
   assert(sampleInfo->isPhased(2) == true);
   assert(sampleInfo->isUnphased(0) == false);
   assert(sampleInfo->isUnphased(1) == false);
   assert(sampleInfo->isUnphased(2) == false);
   assert(record.passedAllFilters() == true);
   assert(record.getNumAlts() == 1);

   assert(reader.readRecord(record) == false);

   assert(reader.getNumKeptRecords() == 1);
   assert(reader.getNumRecords() == 7);

   reader.close();

   
   ////////////////////////////////
   // Test Discarding minor allele count < 1 without subsetting.
   reader.rmDiscardMinAltAlleleCount();
   reader.setDiscardRules(0);
   reader.addDiscardMinMinorAlleleCount(1, NULL);

   reader.open("testFiles/vcfFile.vcf", header);
   
   assert(header.getHeaderLine() == HEADER_LINE);
   assert(header.getNumSamples() == NUM_SAMPLES);
   assert(header.getSampleName(2) == SAMPLES[2]);
   assert(header.getSampleName(0) == SAMPLES[0]);
   assert(header.getSampleName(1) == SAMPLES[1]);
   assert(header.getSampleIndex(SAMPLES[1].c_str()) == 1);
   assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
   assert(header.getSampleIndex(SAMPLES[2].c_str()) == 2);
 
   // Read the records.
   assert(reader.readRecord(record));
   assert(record.allPhased() == false);
   assert(record.allUnphased() == false);
   assert(record.hasAllGenotypeAlleles() == true);
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 3);
   assert(*(sampleInfo->getString("GT", 0)) == "0|0");
   assert(*(sampleInfo->getString("GT", 1)) == "1|0");
   assert(*(sampleInfo->getString("GT", 2)) == "1/1");
   assert(sampleInfo->isPhased(0) == true);
   assert(sampleInfo->isPhased(1) == true);
   assert(sampleInfo->isPhased(2) == false);
   assert(sampleInfo->isUnphased(0) == false);
   assert(sampleInfo->isUnphased(1) == false);
   assert(sampleInfo->isUnphased(2) == true);
   assert(record.passedAllFilters() == true);
   assert(record.getNumAlts() == 1);

   assert(reader.readRecord(record));
   assert(record.allPhased() == false);
   assert(record.allUnphased() == false);
   assert(record.hasAllGenotypeAlleles() == true);
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 3);
   assert(*(sampleInfo->getString("GT", 0)) == "0/0");
   assert(*(sampleInfo->getString("GT", 1)) == "0|1");
   assert(*(sampleInfo->getString("GT", 2)) == "0/0");
   assert(sampleInfo->isPhased(0) == false);
   assert(sampleInfo->isPhased(1) == true);
   assert(sampleInfo->isPhased(2) == false);
   assert(sampleInfo->isUnphased(0) == true);
   assert(sampleInfo->isUnphased(1) == false);
   assert(sampleInfo->isUnphased(2) == true);
   assert(record.passedAllFilters() == false);
   assert(record.getNumAlts() == 1);

   assert(reader.readRecord(record));
   assert(record.allPhased() == false);
   assert(record.allUnphased() == false);
   assert(record.hasAllGenotypeAlleles() == true);
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 3);
   assert(*(sampleInfo->getString("GT", 0)) == "0|0");
   assert(*(sampleInfo->getString("GT", 1)) == "0|0");
   assert(*(sampleInfo->getString("GT", 2)) == "0/0");
   assert(sampleInfo->isPhased(0) == true);
   assert(sampleInfo->isPhased(1) == true);
   assert(sampleInfo->isPhased(2) == false);
   assert(sampleInfo->isUnphased(0) == false);
   assert(sampleInfo->isUnphased(1) == false);
   assert(sampleInfo->isUnphased(2) == true);
   assert(record.passedAllFilters() == true);
   assert(record.getNumAlts() == 0);

   assert(reader.readRecord(record));
   assert(record.allPhased() == false);
   assert(record.allUnphased() == true);
   assert(record.hasAllGenotypeAlleles() == true);
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 3);
   assert(*(sampleInfo->getString("GT", 0)) == "0/1");
   assert(*(sampleInfo->getString("GT", 1)) == "0/2");
   assert(*(sampleInfo->getString("GT", 2)) == "1/1");
   assert(sampleInfo->isPhased(0) == false);
   assert(sampleInfo->isPhased(1) == false);
   assert(sampleInfo->isPhased(2) == false);
   assert(sampleInfo->isUnphased(0) == true);
   assert(sampleInfo->isUnphased(1) == true);
   assert(sampleInfo->isUnphased(2) == true);
   assert(record.passedAllFilters() == true);
   assert(record.getNumAlts() == 2);

   assert(reader.readRecord(record));
   assert(record.allPhased() == true);
   assert(record.allUnphased() == false);
   assert(record.hasAllGenotypeAlleles() == false);
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 3);
   assert(*(sampleInfo->getString("GT", 0)) == "0|1");
   assert(*(sampleInfo->getString("GT", 1)) == "0|.");
   assert(*(sampleInfo->getString("GT", 2)) == "1|1");
   assert(sampleInfo->isPhased(0) == true);
   assert(sampleInfo->isPhased(1) == true);
   assert(sampleInfo->isPhased(2) == true);
   assert(sampleInfo->isUnphased(0) == false);
   assert(sampleInfo->isUnphased(1) == false);
   assert(sampleInfo->isUnphased(2) == false);
   assert(record.passedAllFilters() == true);
   assert(record.getNumAlts() == 1);

   assert(reader.readRecord(record) == false);

   assert(reader.getNumKeptRecords() == 5);
   assert(reader.getNumRecords() == 7);

   reader.close();

   ////////////////////////////////
   // Test Discarding minor allele count < 2 without subsetting.
   reader.addDiscardMinMinorAlleleCount(2, NULL);

   reader.open("testFiles/vcfFile.vcf", header);
   
   assert(header.getHeaderLine() == HEADER_LINE);
   assert(header.getNumSamples() == NUM_SAMPLES);
   assert(header.getSampleName(2) == SAMPLES[2]);
   assert(header.getSampleName(0) == SAMPLES[0]);
   assert(header.getSampleName(1) == SAMPLES[1]);
   assert(header.getSampleIndex(SAMPLES[1].c_str()) == 1);
   assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
   assert(header.getSampleIndex(SAMPLES[2].c_str()) == 2);
 
   // Read the records.
   assert(reader.readRecord(record));
   assert(record.allPhased() == false);
   assert(record.allUnphased() == false);
   assert(record.hasAllGenotypeAlleles() == true);
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 3);
   assert(*(sampleInfo->getString("GT", 0)) == "0|0");
   assert(*(sampleInfo->getString("GT", 1)) == "1|0");
   assert(*(sampleInfo->getString("GT", 2)) == "1/1");
   assert(sampleInfo->isPhased(0) == true);
   assert(sampleInfo->isPhased(1) == true);
   assert(sampleInfo->isPhased(2) == false);
   assert(sampleInfo->isUnphased(0) == false);
   assert(sampleInfo->isUnphased(1) == false);
   assert(sampleInfo->isUnphased(2) == true);
   assert(record.passedAllFilters() == true);
   assert(record.getNumAlts() == 1);

   assert(reader.readRecord(record));
   assert(record.allPhased() == false);
   assert(record.allUnphased() == false);
   assert(record.hasAllGenotypeAlleles() == true);
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 3);
   assert(*(sampleInfo->getString("GT", 0)) == "0|0");
   assert(*(sampleInfo->getString("GT", 1)) == "0|0");
   assert(*(sampleInfo->getString("GT", 2)) == "0/0");
   assert(sampleInfo->isPhased(0) == true);
   assert(sampleInfo->isPhased(1) == true);
   assert(sampleInfo->isPhased(2) == false);
   assert(sampleInfo->isUnphased(0) == false);
   assert(sampleInfo->isUnphased(1) == false);
   assert(sampleInfo->isUnphased(2) == true);
   assert(record.passedAllFilters() == true);
   assert(record.getNumAlts() == 0);

   assert(reader.readRecord(record));
   assert(record.allPhased() == true);
   assert(record.allUnphased() == false);
   assert(record.hasAllGenotypeAlleles() == false);
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 3);
   assert(*(sampleInfo->getString("GT", 0)) == "0|1");
   assert(*(sampleInfo->getString("GT", 1)) == "0|.");
   assert(*(sampleInfo->getString("GT", 2)) == "1|1");
   assert(sampleInfo->isPhased(0) == true);
   assert(sampleInfo->isPhased(1) == true);
   assert(sampleInfo->isPhased(2) == true);
   assert(sampleInfo->isUnphased(0) == false);
   assert(sampleInfo->isUnphased(1) == false);
   assert(sampleInfo->isUnphased(2) == false);
   assert(record.passedAllFilters() == true);
   assert(record.getNumAlts() == 1);

   assert(reader.readRecord(record) == false);

   assert(reader.getNumKeptRecords() == 3);
   assert(reader.getNumRecords() == 7);

   reader.close();


   ////////////////////////////////
   // Test Discarding minor allele count < 1 with subsetting.
   VcfSubsetSamples minMinorAlleleSubset;
   minMinorAlleleSubset.init(header, true);
   minMinorAlleleSubset.addExcludeSample("NA00002");
   reader.addDiscardMinMinorAlleleCount(1, &minMinorAlleleSubset);

   reader.open("testFiles/vcfFile.vcf", header);
   
   assert(header.getHeaderLine() == HEADER_LINE);
   assert(header.getNumSamples() == NUM_SAMPLES);
   assert(header.getSampleName(2) == SAMPLES[2]);
   assert(header.getSampleName(0) == SAMPLES[0]);
   assert(header.getSampleName(1) == SAMPLES[1]);
   assert(header.getSampleIndex(SAMPLES[1].c_str()) == 1);
   assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
   assert(header.getSampleIndex(SAMPLES[2].c_str()) == 2);
 
   // Read the records.
   assert(reader.readRecord(record));
   assert(record.allPhased() == false);
   assert(record.allUnphased() == false);
   assert(record.hasAllGenotypeAlleles() == true);
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 3);
   assert(*(sampleInfo->getString("GT", 0)) == "0|0");
   assert(*(sampleInfo->getString("GT", 1)) == "1|0");
   assert(*(sampleInfo->getString("GT", 2)) == "1/1");
   assert(sampleInfo->isPhased(0) == true);
   assert(sampleInfo->isPhased(1) == true);
   assert(sampleInfo->isPhased(2) == false);
   assert(sampleInfo->isUnphased(0) == false);
   assert(sampleInfo->isUnphased(1) == false);
   assert(sampleInfo->isUnphased(2) == true);
   assert(record.passedAllFilters() == true);
   assert(record.getNumAlts() == 1);

   assert(reader.readRecord(record));
   assert(record.allPhased() == false);
   assert(record.allUnphased() == false);
   assert(record.hasAllGenotypeAlleles() == true);
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 3);
   assert(*(sampleInfo->getString("GT", 0)) == "0|0");
   assert(*(sampleInfo->getString("GT", 1)) == "0|0");
   assert(*(sampleInfo->getString("GT", 2)) == "0/0");
   assert(sampleInfo->isPhased(0) == true);
   assert(sampleInfo->isPhased(1) == true);
   assert(sampleInfo->isPhased(2) == false);
   assert(sampleInfo->isUnphased(0) == false);
   assert(sampleInfo->isUnphased(1) == false);
   assert(sampleInfo->isUnphased(2) == true);
   assert(record.passedAllFilters() == true);
   assert(record.getNumAlts() == 0);

   assert(reader.readRecord(record));
   assert(record.allPhased() == true);
   assert(record.allUnphased() == false);
   assert(record.hasAllGenotypeAlleles() == false);
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 3);
   assert(*(sampleInfo->getString("GT", 0)) == "0|1");
   assert(*(sampleInfo->getString("GT", 1)) == "0|.");
   assert(*(sampleInfo->getString("GT", 2)) == "1|1");
   assert(sampleInfo->isPhased(0) == true);
   assert(sampleInfo->isPhased(1) == true);
   assert(sampleInfo->isPhased(2) == true);
   assert(sampleInfo->isUnphased(0) == false);
   assert(sampleInfo->isUnphased(1) == false);
   assert(sampleInfo->isUnphased(2) == false);
   assert(record.passedAllFilters() == true);
   assert(record.getNumAlts() == 1);

   assert(reader.readRecord(record) == false);

   assert(reader.getNumKeptRecords() == 3);
   assert(reader.getNumRecords() == 7);

   reader.close();

   ////////////////////////////////
   // Test Discarding minor allele count < 2 without subsetting.
   reader.addDiscardMinMinorAlleleCount(2, &minMinorAlleleSubset);

   reader.open("testFiles/vcfFile.vcf", header);
   
   assert(header.getHeaderLine() == HEADER_LINE);
   assert(header.getNumSamples() == NUM_SAMPLES);
   assert(header.getSampleName(2) == SAMPLES[2]);
   assert(header.getSampleName(0) == SAMPLES[0]);
   assert(header.getSampleName(1) == SAMPLES[1]);
   assert(header.getSampleIndex(SAMPLES[1].c_str()) == 1);
   assert(header.getSampleIndex(SAMPLES[0].c_str()) == 0);
   assert(header.getSampleIndex(SAMPLES[2].c_str()) == 2);
 
   // Read the records.
   assert(reader.readRecord(record));
   assert(record.allPhased() == false);
   assert(record.allUnphased() == false);
   assert(record.hasAllGenotypeAlleles() == true);
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 3);
   assert(*(sampleInfo->getString("GT", 0)) == "0|0");
   assert(*(sampleInfo->getString("GT", 1)) == "1|0");
   assert(*(sampleInfo->getString("GT", 2)) == "1/1");
   assert(sampleInfo->isPhased(0) == true);
   assert(sampleInfo->isPhased(1) == true);
   assert(sampleInfo->isPhased(2) == false);
   assert(sampleInfo->isUnphased(0) == false);
   assert(sampleInfo->isUnphased(1) == false);
   assert(sampleInfo->isUnphased(2) == true);
   assert(record.passedAllFilters() == true);
   assert(record.getNumAlts() == 1);

   assert(reader.readRecord(record));
   assert(record.allPhased() == false);
   assert(record.allUnphased() == false);
   assert(record.hasAllGenotypeAlleles() == true);
   sampleInfo = &(record.getGenotypeInfo());
   assert(sampleInfo->getNumSamples() == 3);
   assert(*(sampleInfo->getString("GT", 0)) == "0|0");
   assert(*(sampleInfo->getString("GT", 1)) == "0|0");
   assert(*(sampleInfo->getString("GT", 2)) == "0/0");
   assert(sampleInfo->isPhased(0) == true);
   assert(sampleInfo->isPhased(1) == true);
   assert(sampleInfo->isPhased(2) == false);
   assert(sampleInfo->isUnphased(0) == false);
   assert(sampleInfo->isUnphased(1) == false);
   assert(sampleInfo->isUnphased(2) == true);
   assert(record.passedAllFilters() == true);
   assert(record.getNumAlts() == 0);

   assert(reader.readRecord(record) == false);

   assert(reader.getNumKeptRecords() == 2);
   assert(reader.getNumRecords() == 7);

   reader.close();
}


void testVcfWriteFile()
{
    VcfFileWriter writer;
    VcfFileReader reader;
    VcfHeader header;
    VcfRecord record;

    assert(reader.open("testFiles/vcfFile.vcf", header) == true);
    assert(writer.open("results/vcfFile.vcf", header, InputFile::DEFAULT)
           == true);
    while(reader.readRecord(record))
    {
        // Write the record.
        assert(writer.writeRecord(record));
    }
    
    assert(reader.open("testFiles/vcfFile.vcf", header) == true);
    assert(writer.open("results/vcfFileNoInfo.vcf", header, 
                       InputFile::DEFAULT) == true);
    while(reader.readRecord(record))
    {
        // Test Clearing the INFO field.
        record.getInfo().clear();
        // Write the record.
        assert(writer.writeRecord(record));
    }

    assert(reader.open("testFiles/vcfFile.vcf", header) == true);
    assert(writer.open("results/vcfFileNoInfoBGZF.vcf", header) == true);
    while(reader.readRecord(record))
    {
        // Test Clearing the INFO field.
        record.getInfo().clear();
        // Write the record.
        assert(writer.writeRecord(record));
    }

    assert(reader.open("testFiles/vcfFile.vcf", header) == true);
    VcfRecordGenotype::addStoreField("GT");
    assert(writer.open("results/vcfFileNoInfoKeepGT.vcf", header, 
                       InputFile::DEFAULT) == true);
    while(reader.readRecord(record))
    {
        // Test Clearing the INFO field.
        record.getInfo().clear();
        // Write the record.
        assert(writer.writeRecord(record));
    }

    assert(reader.open("testFiles/vcfFile.vcf", header) == true);
    // Undo the storing of GT.
    VcfRecordGenotype::storeAllFields();
    VcfRecordGenotype::addStoreField("GQ");
    VcfRecordGenotype::addStoreField("XX");
    VcfRecordGenotype::addStoreField("HQ");
    assert(writer.open("results/vcfFileNoInfoKeepGQHQ.vcf", header, 
                       InputFile::DEFAULT) == true);

    assert(strcmp(writer.getFileName(), "results/vcfFileNoInfoKeepGQHQ.vcf") == 0);

    while(reader.readRecord(record))
    {
        // Test Clearing the INFO field.
        record.getInfo().clear();
        // Write the record.
        assert(writer.writeRecord(record));
    }
}


void testVcfReadSection()
{
    // Test open for read via the constructor with return.
    VcfFileReader reader;
    VcfHeader header;
    VcfRecord record;
    
    const Tabix* tabixPtr = NULL;

    tabixPtr = reader.getVcfIndex();
    assert(tabixPtr == NULL);

    ////////////////////////////////
    // Test the read section logic.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    tabixPtr = reader.getVcfIndex();
    assert(tabixPtr == NULL);
    reader.readVcfIndex();

    //////////////////
    // Test index accessors.
    tabixPtr = reader.getVcfIndex();
    assert(tabixPtr != NULL);
    assert(tabixPtr->getFormat() == Tabix::FORMAT_VCF);
    assert(tabixPtr->getNumRefs() == 2);
    assert(strcmp(tabixPtr->getRefName(0), "1") == 0);
    assert(strcmp(tabixPtr->getRefName(1), "3") == 0);
    bool caughtException = false;
    try
    {
        tabixPtr->getRefName(2);
    }
    catch(std::exception& e)
    {
        caughtException = true;
    }
    assert(caughtException);
    caughtException = false;
    try
    {
        tabixPtr->getRefName(-1);
    }
    catch(std::exception& e)
    {
        caughtException = true;
    }
    assert(caughtException == true);
    caughtException = false;

    reader.set1BasedReadSection("10", 16384, 32767);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("1", 16384, 32767);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("1", 16384, 32768);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("1", 16384, 32769);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("1", 32768, 32769);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("1", 32769, 32767);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("1", 32769, 65537);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("1", 32769, 65537);
    assert(reader.readRecord(record) == false);

    assert(reader.set1BasedReadSection("1", 32768, 65538));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 65537);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    assert(reader.set1BasedReadSection("1", 32769, 65538));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 65537);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    assert(reader.set1BasedReadSection("1", 0, 65538));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 65537);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);


    ////////////////////////////////////////
    // Test selecting whole chroms

    assert(reader.setReadSection("10"));
    assert(reader.readRecord(record) == false);

    assert(reader.setReadSection("1"));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 65537);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    assert(reader.setReadSection("3"));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    ////////////////////////////////////////
    // Test selecting sections with deletions
    reader.set1BasedReadSection("3", 16384, 32767);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 16384, 32768);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 16384, 32769);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 32768, 32769);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 32769, 32767);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 32769, 65537);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 32769, 65537);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 32770, 65537);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 32771, 65537);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 32780, 65537);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 32781, 65537);
    assert(reader.readRecord(record) == false);

    assert(reader.set1BasedReadSection("3", 32768, 65538));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    assert(reader.set1BasedReadSection("3", 32769, 65538));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    assert(reader.set1BasedReadSection("3", 32770, 65538));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    assert(reader.set1BasedReadSection("3", 32771, 65538));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    assert(reader.set1BasedReadSection("3", 0, 65538));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    ////////////////////////////////////////
    // Test selecting sections with deletions for overlapping
    reader.set1BasedReadSection("3", 16384, 32767, true);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 16384, 32768, true);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 16384, 32769, true);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 32768, 32769, true);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 32769, 32767, true);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 32769, 65537, true);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 32769, 65537, true);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 32770, 65537, true);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 32771, 65537, true);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 32780, 65537, true);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);

    reader.set1BasedReadSection("3", 32781, 65537, true);
    assert(reader.readRecord(record) == false);

    assert(reader.set1BasedReadSection("3", 32768, 65538, true));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    assert(reader.set1BasedReadSection("3", 32769, 65538, true));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    assert(reader.set1BasedReadSection("3", 32770, 65538, true));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    assert(reader.set1BasedReadSection("3", 32771, 65538, true));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    assert(reader.set1BasedReadSection("3", 0, 65538, true));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

   reader.close();
}


void testVcfReadSectionNoIndex()
{
    // Test open for read via the constructor with return.
    VcfFileReader reader;
    VcfHeader header;
    VcfRecord record;
    bool caughtException = false;

    ////////////////////////////////
    // Test the read section logic.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    
    reader.set1BasedReadSection("10", 16384, 32767);
    assert(reader.readRecord(record) == false);

     // Can't call setReadSection after reading a record a second time with no index.
    try
    {
        caughtException = false;
        assert(reader.set1BasedReadSection("1", 16384, 32769));
        assert(reader.readRecord(record) == true);
    }
    catch (std::exception& e)
    {
        caughtException = true;
    }
    assert(caughtException == true);

    // Reopen to begining to check chrom 1.
    reader.open("testFiles/testTabix.vcf.bgzf", header);

    reader.set1BasedReadSection("1", 16384, 32768);
    assert(reader.readRecord(record) == false);

    // Reopen to begining to check valid range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    reader.set1BasedReadSection("1", 16384, 32769);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == false);

    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    reader.set1BasedReadSection("1", 32769, 32767);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    reader.set1BasedReadSection("1", 32769, 65537);
    assert(reader.readRecord(record) == false);

    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    assert(reader.set1BasedReadSection("1", 32768, 65538));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 65537);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    assert(reader.set1BasedReadSection("1", 32769, 65538));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 65537);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    assert(reader.set1BasedReadSection("1", 0, 65538));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 65537);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);


    ////////////////////////////////////////
    // Test selecting whole chroms
    try
    {
        caughtException = false;
        assert(reader.setReadSection("1"));
        assert(reader.readRecord(record) == true);
    }
    catch (std::exception& e)
    {
        caughtException = true;
    }
    assert(caughtException);

    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    assert(reader.setReadSection("10"));
    assert(reader.readRecord(record) == false);

    assert(strcmp(reader.getFileName(), "testFiles/testTabix.vcf.bgzf") == 0);

    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    assert(reader.setReadSection("1"));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 65537);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    assert(reader.setReadSection("3"));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    ////////////////////////////////////////
    // Test selecting sections with deletions
    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    reader.set1BasedReadSection("3", 16384, 32767);
    assert(reader.readRecord(record) == false);

    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    reader.set1BasedReadSection("3", 16384, 32768);
    assert(reader.readRecord(record) == false);

    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    reader.set1BasedReadSection("3", 16384, 32769);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == false);

    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    reader.set1BasedReadSection("3", 32768, 32769);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == false);

    ////////////////////////////////////////
    // Test selecting sections with deletions for overlapping
    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    reader.set1BasedReadSection("3", 16384, 32767, true);
    assert(reader.readRecord(record) == false);

    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    reader.set1BasedReadSection("3", 16384, 32768, true);
    assert(reader.readRecord(record) == false);

    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    reader.set1BasedReadSection("3", 16384, 32769, true);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == false);

    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    assert(reader.set1BasedReadSection("3", 32771, 65538, true));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

    // Reopen to begining to check another range.
    reader.open("testFiles/testTabix.vcf.bgzf", header);
    assert(reader.set1BasedReadSection("3", 0, 65538, true));
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32768);
    assert(reader.readRecord(record) == true);
    assert(record.get1BasedPosition() == 32780);
    assert(reader.readRecord(record) == false);
    assert(reader.readRecord(record) == false);

   reader.close();
}


void testVcfReadSectionBadIndex()
{
    // Test open for read via the constructor with return.
    VcfFileReader reader;
    VcfFileWriter writer;
    VcfHeader header;
    VcfRecord record;
    
    ////////////////////////////////
    // Test the read section logic.
    reader.open("testFiles/testTabixBadIndex.vcf.bgzf", header);
    bool hitError = false;
    try
    {
        reader.readVcfIndex();
    }
    catch(std::exception& e)
    {
        hitError = true;
        std::string expectedError = "FAIL_PARSE: ERROR: Tabix file not in VCF format: testFiles/testTabixBadIndex.vcf.bgzf.tbi\nFAIL_IO: Failed to read the vcf Index file: testFiles/testTabixBadIndex.bgzf.tbi";
       assert(expectedError == e.what());
    }
    assert(hitError);
}
