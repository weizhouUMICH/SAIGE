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

#include "TestSamCoordOutput.h"
#include "TestValidate.h"
#include "SamCoordOutput.h"
#include "SamRecordPool.h"
#include <assert.h>

void testSamCoordOutput()
{
    // Call generic test.
    SamCoordOutputTest::testSamCoordOutput();
}


void SamCoordOutputTest::testSamCoordOutput()
{
    SamRecordPool pool(3);

    SamCoordOutput outputBuffer(pool);

    SamFile inSam;
    SamFile outSam;
    SamFileHeader samHeader;
    SamRecord* rec1 = NULL;
    SamRecord* rec2 = NULL;
    SamRecord* rec3 = NULL;

    // Open input file and read the header.
#ifdef __ZLIB_AVAILABLE__
    assert(inSam.OpenForRead("testFiles/testBam.bam"));
#else
    assert(inSam.OpenForRead("testFiles/testSam.sam"));
#endif
    assert(inSam.ReadHeader(samHeader));
    validateHeader(samHeader);

    // Check failed to add empty record.
    assert(!outputBuffer.add(rec1));

    // Read the first 3 records from the input file.
    rec1 = pool.getRecord();
    assert(inSam.ReadRecord(samHeader, *rec1) == true);
    validateRead1(*rec1);
    rec2 = pool.getRecord();
    assert(inSam.ReadRecord(samHeader, *rec2) == true);
    validateRead2(*rec2);
    rec3 = pool.getRecord();
    assert(inSam.ReadRecord(samHeader, *rec3) == true);
    validateRead3(*rec3);
    assert(pool.getRecord() == NULL);

    // Add the first 3 records to the output buffer.
    // Sorted order is rec 3, 1, 2
    assert(outputBuffer.add(rec1));
    assert(outputBuffer.add(rec2));
    assert(outputBuffer.add(rec3));

    // Test writing to the output buffer without having set it.
    // Should flush just rec3 out.
    assert(!outputBuffer.flush(0, 100));
    
    // Set the output buffer.
    outputBuffer.setOutputFile(&outSam, &samHeader);

    // Open output file and write the header.
    assert(outSam.OpenForWrite("results/TestSamCoordOutput.sam"));
    assert(outSam.WriteHeader(samHeader));

    // Read another 1 record (reuse record pointers).
    rec1 = pool.getRecord();
    assert(inSam.ReadRecord(samHeader, *rec1) == true);
    validateRead4(*rec1);
    assert(outputBuffer.add(rec1));

    rec1 = pool.getRecord();
    assert(rec1 == NULL);


    // Flush out just the reads before this position.
    assert(outputBuffer.flush(0, 1011));

    // Read 2 more records.
    rec1 = pool.getRecord();
    assert(inSam.ReadRecord(samHeader, *rec1) == true);
    validateRead5(*rec1);
    assert(outputBuffer.add(rec1));

    rec1 = pool.getRecord();
    assert(inSam.ReadRecord(samHeader, *rec1) == true);
    validateRead6(*rec1);
    assert(outputBuffer.add(rec1));

    // Can get another record (tests to make sure flushes up to and
    // including the specified position).  If it did not
    // flush the specified position, there would not be
    // another record available.
    rec1 = pool.getRecord();
    assert(rec1 != NULL);

    // Flush out just the reads before this position.
    assert(outputBuffer.flush(0, 1012));

    // Read another record.
    assert(inSam.ReadRecord(samHeader, *rec1) == true);
    validateRead7(*rec1);
    assert(outputBuffer.add(rec1));
    assert(pool.getRecord() == NULL);

    // Flush out just the reads on chrom 1 (chrom id 0).
    assert(outputBuffer.flush(0, -1));

    // Read another record.
    rec1 = pool.getRecord();
    assert(inSam.ReadRecord(samHeader, *rec1) == true);
    validateRead8(*rec1);
    assert(outputBuffer.add(rec1));
    assert(pool.getRecord() == NULL);

    // Flush out the chrom 2 (chrom id 1) reads.
    assert(outputBuffer.flush(2, 0));

    // Read the rest of the records.
    rec1 = pool.getRecord();
    assert(inSam.ReadRecord(samHeader, *rec1) == true);
    validateRead9(*rec1);
    assert(outputBuffer.add(rec1));
    rec1 = pool.getRecord();
    assert(inSam.ReadRecord(samHeader, *rec1) == true);
    validateRead10(*rec1);
    assert(outputBuffer.add(rec1));
    assert(pool.getRecord() == NULL);

    // Flush the rest by passing in -1, -1
    assert(outputBuffer.flush(-1, -1));
}
