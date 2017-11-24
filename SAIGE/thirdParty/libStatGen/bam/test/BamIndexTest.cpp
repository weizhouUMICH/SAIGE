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

#include "BamIndex.h"
#include "TestValidate.h"
#include "BamIndexTest.h"

#include <assert.h>

void testIndex(BamIndex& bamIndex)
{
#ifdef __ZLIB_AVAILABLE__
    assert(bamIndex.getNumMappedReads(1) == 2);
    assert(bamIndex.getNumUnMappedReads(1) == 0);
    assert(bamIndex.getNumMappedReads(0) == 4);
    assert(bamIndex.getNumUnMappedReads(0) == 1);
    assert(bamIndex.getNumMappedReads(23) == -1);
    assert(bamIndex.getNumUnMappedReads(23) == -1);
    assert(bamIndex.getNumMappedReads(-1) == 0);
    assert(bamIndex.getNumUnMappedReads(-1) == 2);
    assert(bamIndex.getNumMappedReads(-2) == -1);
    assert(bamIndex.getNumUnMappedReads(-2) == -1);
    assert(bamIndex.getNumMappedReads(22) == 0);
    assert(bamIndex.getNumUnMappedReads(22) == 0);

    // Get the chunks for reference id 1.
    Chunk testChunk;
    SortedChunkList chunkList;
    assert(bamIndex.getChunksForRegion(1, -1, -1, chunkList) == true);
    assert(!chunkList.empty());
    testChunk = chunkList.pop();
    assert(chunkList.empty());
    assert(testChunk.chunk_beg == 0x4e7);
    assert(testChunk.chunk_end == 0x599);

    // Get the chunks for reference id 0.
    assert(bamIndex.getChunksForRegion(0, -1, -1, chunkList) == true);
    assert(!chunkList.empty());
    testChunk = chunkList.pop();
    assert(chunkList.empty());
    assert(testChunk.chunk_beg == 0x360);
    assert(testChunk.chunk_end == 0x4e7);


    // Get the chunks for reference id 2.
    assert(bamIndex.getChunksForRegion(2, -1, -1, chunkList) == true);
    assert(!chunkList.empty());
    testChunk = chunkList.pop();
    assert(chunkList.empty());
    assert(testChunk.chunk_beg == 0x599);
    assert(testChunk.chunk_end == 0x5ea);

    // Get the chunks for reference id 3.
    // There isn't one for this ref id, but still successfully read the file,
    // so it should return true, but the list should be empty.
    assert(bamIndex.getChunksForRegion(3, -1, -1, chunkList) == true);
    assert(chunkList.empty());

    // Test reading an indexed bam file.
    SamFile inFile;
    assert(inFile.OpenForRead("testFiles/sortedBam.bam"));
    inFile.setSortedValidation(SamFile::COORDINATE);
    assert(inFile.ReadBamIndex("testFiles/sortedBam.bam.bai"));
    SamFileHeader samHeader;
    assert(inFile.ReadHeader(samHeader));
    SamRecord samRecord;

    // Test getting num mapped/unmapped reads.
    assert(inFile.getNumMappedReadsFromIndex(1) == 2);
    assert(inFile.getNumUnMappedReadsFromIndex(1) == 0);
    assert(inFile.getNumMappedReadsFromIndex(0) == 4);
    assert(inFile.getNumUnMappedReadsFromIndex(0) == 1);
    assert(inFile.getNumMappedReadsFromIndex(23) == -1);
    assert(inFile.getNumUnMappedReadsFromIndex(23) == -1);
    assert(inFile.getNumMappedReadsFromIndex(-1) == 0);
    assert(inFile.getNumUnMappedReadsFromIndex(-1) == 2);
    assert(inFile.getNumMappedReadsFromIndex(-2) == -1);
    assert(inFile.getNumUnMappedReadsFromIndex(-2) == -1);
    assert(inFile.getNumMappedReadsFromIndex(22) == 0);
    assert(inFile.getNumUnMappedReadsFromIndex(22) == 0);

    assert(inFile.getNumMappedReadsFromIndex("2", samHeader) == 2);
    assert(inFile.getNumUnMappedReadsFromIndex("2", samHeader) == 0);
    assert(inFile.getNumMappedReadsFromIndex("1", samHeader) == 4);
    assert(inFile.getNumUnMappedReadsFromIndex("1", samHeader) == 1);
    assert(inFile.getNumMappedReadsFromIndex("22", samHeader) == 0);
    assert(inFile.getNumUnMappedReadsFromIndex("22", samHeader) == 0);
    assert(inFile.getNumMappedReadsFromIndex("", samHeader) == 0);
    assert(inFile.getNumUnMappedReadsFromIndex("*", samHeader) == 2);
    assert(inFile.getNumMappedReadsFromIndex("unknown", samHeader) == -1);
    assert(inFile.getNumUnMappedReadsFromIndex("unknown", samHeader) == -1);
    assert(inFile.getNumMappedReadsFromIndex("X", samHeader) == 0);
    assert(inFile.getNumUnMappedReadsFromIndex("X", samHeader) == 0);

   // Set the read section saying the reads must be fully enclosed in the section.
    assert(inFile.SetReadSection("1", 1010, 1011, false));
    assert(inFile.ReadRecord(samHeader, samRecord) == false);

    assert(inFile.SetReadSection("1", 1011, 1012, false));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead2(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);

    // Section -1 = Ref *: 2 records (8 & 10 from testSam.sam that is reflected
    // in the validation.
    assert(inFile.SetReadSection(-1));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead8(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead10(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);

    // Section 2 = Ref 3: 1 records (9 from testSam.sam that is reflected
    // in the validation.
    assert(inFile.SetReadSection(2));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead9(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);

    // Section 0 = Ref 1: 5 records (3, 4, 1, 2, & 6 from testSam.sam that is
    // reflected in the validation.
    assert(inFile.SetReadSection(0));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead3(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead4(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead1(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead2(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead6(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);

    // Section 1 = Ref 2: 2 records (5 & 7 from testSam.sam that is reflected
    // in the validation.
    assert(inFile.SetReadSection(1));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead5(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead7(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);

    // Section 3 to 22 (ref 4 - 23): 0 records.
    for(int i = 3; i < 23; i++)
    {
        assert(inFile.SetReadSection(i));
        assert(inFile.ReadRecord(samHeader, samRecord) == false);
    }


    // Set the read section.
    assert(inFile.SetReadSection("1", 1010, 1012));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead1(samRecord);
    assert(inFile.GetNumOverlaps(samRecord) == 2);
    assert(samRecord.getNumOverlaps(1010, 1012) == 2);
    assert(samRecord.getNumOverlaps(1010, 1020) == 5);
    assert(samRecord.getNumOverlaps(1010, 1011) == 1);
    assert(samRecord.getNumOverlaps(1011, 1012) == 1);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead2(samRecord);
    assert(inFile.GetNumOverlaps(samRecord) == 0);
    assert(samRecord.getNumOverlaps(1010, 1012) == 0);
    assert(samRecord.getNumOverlaps(1010, 1020) == 0);
    assert(samRecord.getNumOverlaps(1010, 1011) == 0);
    assert(samRecord.getNumOverlaps(1011, 1012) == 0);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);
           
    assert(inFile.SetReadSection("1", 1010, 1020));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead1(samRecord);
    assert(inFile.GetNumOverlaps(samRecord) == 5);
    assert(samRecord.getNumOverlaps(1010, 1012) == 2);
    assert(samRecord.getNumOverlaps(1010, 1020) == 5);
    assert(samRecord.getNumOverlaps(1010, 1011) == 1);
    assert(samRecord.getNumOverlaps(1011, 1012) == 1);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead2(samRecord);
    assert(inFile.GetNumOverlaps(samRecord) == 0);
    assert(samRecord.getNumOverlaps(1010, 1012) == 0);
    assert(samRecord.getNumOverlaps(1010, 1020) == 0);
    assert(samRecord.getNumOverlaps(1010, 1011) == 0);
    assert(samRecord.getNumOverlaps(1011, 1012) == 0);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);
           
    assert(inFile.SetReadSection("1", 1010, 1011));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead1(samRecord);
    assert(inFile.GetNumOverlaps(samRecord) == 1);
    assert(samRecord.getNumOverlaps(1010, 1012) == 2);
    assert(samRecord.getNumOverlaps(1010, 1020) == 5);
    assert(samRecord.getNumOverlaps(1010, 1011) == 1);
    assert(samRecord.getNumOverlaps(1011, 1012) == 1);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);
           
    assert(inFile.SetReadSection("1", 1011, 1012));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead1(samRecord);
    assert(inFile.GetNumOverlaps(samRecord) == 1);
    assert(samRecord.getNumOverlaps(1010, 1012) == 2);
    assert(samRecord.getNumOverlaps(1010, 1020) == 5);
    assert(samRecord.getNumOverlaps(1010, 1011) == 1);
    assert(samRecord.getNumOverlaps(1011, 1012) == 1);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead2(samRecord);
    assert(inFile.GetNumOverlaps(samRecord) == 0);
    assert(samRecord.getNumOverlaps(1010, 1012) == 0);
    assert(samRecord.getNumOverlaps(1010, 1020) == 0);
    assert(samRecord.getNumOverlaps(1010, 1011) == 0);
    assert(samRecord.getNumOverlaps(1011, 1012) == 0);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);
#endif
}


void testBamIndex()
{
    // BAM indexes are compressed, so can't be tested without zlib.
#ifdef __ZLIB_AVAILABLE__
    // Create a bam index.
    BamIndex bamIndex;
    bamIndex.readIndex("testFiles/sortedBam.bam.bai");
    testIndex(bamIndex);

    BamIndexFileTest test1;
    bool caughtException = false;
    try
    {
        // Try reading the bam index without specifying a
        // filename and before opening a bam file.
        assert(test1.ReadBamIndex() == false);
    }
    catch (std::exception& e)
    {
        caughtException = true;
        assert(strcmp(e.what(), "FAIL_ORDER: Failed to read the bam Index file - the BAM file needs to be read first in order to determine the index filename.") == 0);
    }
    // Should have failed and thrown an exception.
    assert(caughtException);

    // Read the bam index with a specified name.
    assert(test1.ReadBamIndex("testFiles/sortedBam.bam.bai"));
    BamIndex* index = test1.getBamIndex();
    assert(index != NULL);
    testIndex(*index);

    // Open the bam file so the index can be opened.
    assert(test1.OpenForRead("testFiles/sortedBam.bam"));
    // Try reading the bam index without specifying a
    // filename after opening a bam file.
    assert(test1.ReadBamIndex() == true);
    index = test1.getBamIndex();
    assert(index != NULL);
    testIndex(*index);

    // Open the bam file so the index can be opened.
    // This time the index file does not have .bam in it.
    assert(test1.OpenForRead("testFiles/sortedBam2.bam"));
    // Try reading the bam index without specifying a
    // filename after opening a bam file.
    assert(test1.ReadBamIndex() == true);
    index = test1.getBamIndex();
    assert(index != NULL);
    testIndex(*index);
#endif
}


