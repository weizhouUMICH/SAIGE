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

#include "TestSamRecordPool.h"
#include "SamRecordPool.h"
#include <assert.h>

void testSamRecordPool()
{
    // Call generic test.
    SamRecordPoolTest::testSamRecordPool();
}


void SamRecordPoolTest::testSamRecordPool()
{

    // Attempt to allocate with max size 0,
    // fails to get a record.
    SamRecordPool pool(0);
    assert(pool.getRecord() == NULL);

    // Up the max size to 3.
    pool.setMaxAllocatedRecs(3);

    // Successfully get 1st record.
    SamRecord* rec1 = pool.getRecord();
    assert(rec1 != NULL);

    // Successfully get 2nd record.
    SamRecord* rec2 = pool.getRecord();
    assert(rec2 != NULL);
    assert(rec2 != rec1);

    // Successfully get 3rd record.
    SamRecord* rec3 = pool.getRecord();
    assert(rec3 != NULL);
    assert((rec3 != rec1) && (rec3 != rec2));

    // Fail to get a 4th record.
    assert(pool.getRecord() == NULL);

    // Release a record and confirm its reuse.
    pool.releaseRecord(rec2);
    SamRecord* rec = pool.getRecord();
    assert(rec == rec2);

    // Release multiple records and check reuse.
    pool.releaseRecord(rec3);
    pool.releaseRecord(rec1);
    pool.releaseRecord(rec);
    SamRecord* release1 = pool.getRecord();
    SamRecord* release2 = pool.getRecord();
    SamRecord* release3 = pool.getRecord();
    assert(release1 == rec3);
    assert(release2 == rec1);
    assert(release3 == rec);
    assert(pool.getRecord() == NULL);

    // Up the max allocated size but don't allocate any, then
    // reduce the max allocated size and release all the records
    // but the already allocated records will still be used.
    pool.setMaxAllocatedRecs(4);
    pool.setMaxAllocatedRecs(0);
    pool.releaseRecord(release3);
    pool.releaseRecord(release1);
    pool.releaseRecord(release2);
    rec1 = pool.getRecord();
    rec2 = pool.getRecord();
    rec3 = pool.getRecord();
    assert(rec1 == release3);
    assert(rec2 == release1);
    assert(rec3 == release2);
    assert(pool.getRecord() == NULL);


    // Up the max allocated size and allocate another record.
    pool.setMaxAllocatedRecs(4);
    rec = pool.getRecord();
    assert(rec != NULL);
    assert(rec != rec1);
    assert(rec != rec2);
    assert(rec != rec3);
    assert(pool.getRecord() == NULL);
}
