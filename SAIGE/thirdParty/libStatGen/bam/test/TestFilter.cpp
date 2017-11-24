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

#include "TestFilter.h"
#include "TestValidate.h"
#include "SamFilter.h"
#include <assert.h>

void testFilter()
{
    // Call generic test which since the sam and bam are identical, should
    // contain the same results.
    FilterTest::testFilter(FilterTest::SAM);
#ifdef __ZLIB_AVAILABLE__
    FilterTest::testFilter(FilterTest::BAM);
#endif
}


void FilterTest::testFilter(FileType inputType)
{
    SamFile inSam;

    if(inputType == SAM)
    {
        assert(inSam.OpenForRead("testFiles/testSam.sam"));
    }
    else
    {
        assert(inSam.OpenForRead("testFiles/testBam.bam"));
    }

   // Read the SAM Header.
    SamFileHeader samHeader;
    assert(inSam.ReadHeader(samHeader));
    validateHeader(samHeader);

    SamRecord samRecord;
    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead1(samRecord);

    // Clip the read, 2 from the front and 2 from the back, which causes 2D to
    // be dropped.
    assert(SamFilter::softClip(samRecord, 2, 2) == SamFilter::CLIPPED);
    assert(samRecord.get0BasedPosition() == TestValidate::READ1_POS + 2);
    std::string expectedCigar = "2S1M2S";
    assert(samRecord.getCigar() == expectedCigar);
    assert(samRecord.getSequence() == TestValidate::READ1_SEQ);
    assert(samRecord.getQuality() == TestValidate::READ1_QUAL);
    // Only 1 base, so the end is the same as start
    assert(samRecord.get0BasedAlignmentEnd() == TestValidate::READ1_POS + 2);
    assert(samRecord.getAlignmentLength() == 1);
    assert(samRecord.get0BasedUnclippedStart() == TestValidate::READ1_UNCLIP_START);
    // The new unclipped end is not the same as the original end because the
    // 2 deletions are lost.
    assert(samRecord.get0BasedUnclippedEnd() == TestValidate::READ1_UNCLIP_END - 2);


    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead2(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead3(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead4(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead5(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead6(samRecord);

    // Clip the read 2 more from the front and 2 from the back.
    assert(SamFilter::softClip(samRecord, 5, 2) == SamFilter::CLIPPED);
    assert(samRecord.get0BasedPosition() == TestValidate::READ6_POS + 2);
    expectedCigar = "2H5S1M2S";
    assert(samRecord.getCigar() == expectedCigar);
    assert(samRecord.getSequence() == TestValidate::READ6_SEQ);
    assert(samRecord.getQuality() == TestValidate::READ6_QUAL);
    // Only 1 base, so the end is the same as start
    assert(samRecord.get0BasedAlignmentEnd() == TestValidate::READ6_POS + 2);
    assert(samRecord.getAlignmentLength() == 1);
    assert(samRecord.get0BasedUnclippedStart() == TestValidate::READ6_UNCLIP_START);
    assert(samRecord.get0BasedUnclippedEnd() == TestValidate::READ6_UNCLIP_END);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead7(samRecord);

    // Clip the read 2 more from the front and 2 morefrom the back.
    assert(SamFilter::softClip(samRecord, 5, 3) == SamFilter::CLIPPED);
    assert(samRecord.get0BasedPosition() == TestValidate::READ7_POS + 2);
    expectedCigar = "5S1M3S3H";
    assert(samRecord.getCigar() == expectedCigar);
    assert(samRecord.getSequence() == TestValidate::READ7_SEQ);
    assert(samRecord.getQuality() == TestValidate::READ7_QUAL);
    // Only 1 base, so the end is the same as start
    assert(samRecord.get0BasedAlignmentEnd() == TestValidate::READ7_POS + 2);
    assert(samRecord.getAlignmentLength() == 1);
    assert(samRecord.get0BasedUnclippedStart() == TestValidate::READ7_UNCLIP_START);
    assert(samRecord.get0BasedUnclippedEnd() == TestValidate::READ7_UNCLIP_END);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead8(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead9(samRecord);

    assert(inSam.ReadRecord(samHeader, samRecord) == true);
    validateRead10(samRecord);
}

