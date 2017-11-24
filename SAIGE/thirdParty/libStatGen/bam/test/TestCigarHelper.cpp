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

#include "TestCigarHelper.h"
#include "TestValidate.h"
#include "CigarHelper.h"
#include <assert.h>

void testCigarHelper()
{
    // Call generic test.
    CigarHelperTest::testCigarHelper();
}


void CigarHelperTest::testCigarHelper()
{
    testSoftClipBeginByRefPos();
    testSoftClipEndByRefPos();
}


void CigarHelperTest::testSoftClipBeginByRefPos()
{
    SamRecord record;
    CigarRoller newCigar;
    std::string newCigarString;
    int32_t newPos = 0;

    // Setup the current Cigar.
    // Cigar:   HHHSSSMMMDDDMMMIIIMMMPPPMMMDDDMMMSSSHHH
    // ReadPos:    000000   000011111   111   112222
    // ReadPos:    012345   678901234   567   890123
    // RefPos:        111111111   122   222222223
    // RefPos:        012345678   901   234567890
    const char* origCigar = "3H3S3M3D3M3I3M3P3M3D3M3S3H";
    record.setCigar(origCigar);
    record.set0BasedPosition(10);
    record.setSequence("gggAAATTTCCCTTTGGGAAAggg");

    ////////////////////////////////////////////////////////
    // Clip outside of the range (after).  Everything should be clipped.
    assert(CigarHelper::softClipBeginByRefPos(record, 10000, newCigar, newPos) == 23);
    newCigar.getCigarString(newCigarString);
    assert(strcmp(newCigarString.c_str(), "3H24S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip outside of the range (before).  Nothing should change.
    assert(CigarHelper::softClipBeginByRefPos(record, 1, newCigar, newPos) == 
           CigarHelper::NO_CLIP);
    newCigar.getCigarString(newCigarString);
    assert(strcmp(newCigarString.c_str(), origCigar) == 0);


    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    // Test clipping at every position of the read.

    ////////////////////////////////////////////////////////
    // Clip at the first position.
    assert(CigarHelper::softClipBeginByRefPos(record, 10, newCigar, newPos) == 3);
    assert(newPos == 11);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H4S2M3D3M3I3M3P3M3D3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip in the middle of the first Match.
    assert(CigarHelper::softClipBeginByRefPos(record, 11, newCigar, newPos) == 4);
    assert(newPos == 12);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H5S1M3D3M3I3M3P3M3D3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 12, newCigar, newPos) == 5);
    assert(newPos == 16);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H6S3M3I3M3P3M3D3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 13, newCigar, newPos) == 5);
    assert(newPos == 16);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H6S3M3I3M3P3M3D3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 14, newCigar, newPos) == 5);
    assert(newPos == 16);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H6S3M3I3M3P3M3D3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 15, newCigar, newPos) == 5);
    assert(newPos == 16);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H6S3M3I3M3P3M3D3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 16, newCigar, newPos) == 6);
    assert(newPos == 17);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H7S2M3I3M3P3M3D3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 17, newCigar, newPos) == 7);
    assert(newPos == 18);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H8S1M3I3M3P3M3D3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 18, newCigar, newPos) == 11);
    assert(newPos == 19);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H12S3M3P3M3D3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 19, newCigar, newPos) == 12);
    assert(newPos == 20);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H13S2M3P3M3D3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 20, newCigar, newPos) == 13);
    assert(newPos == 21);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H14S1M3P3M3D3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 21, newCigar, newPos) == 14);
    assert(newPos == 22);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H15S3M3D3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 22, newCigar, newPos) == 15);
    assert(newPos == 23);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H16S2M3D3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 23, newCigar, newPos) == 16);
    assert(newPos == 24);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H17S1M3D3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 24, newCigar, newPos) == 17);
    assert(newPos == 28);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H18S3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 25, newCigar, newPos) == 17);
    assert(newPos == 28);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H18S3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 26, newCigar, newPos) == 17);
    assert(newPos == 28);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H18S3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 27, newCigar, newPos) == 17);
    assert(newPos == 28);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H18S3M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 28, newCigar, newPos) == 18);
    assert(newPos == 29);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H19S2M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 29, newCigar, newPos) == 19);
    assert(newPos == 30);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H20S1M3S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 30, newCigar, newPos) == 23);
    assert(newPos == 10);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H24S3H") == 0);

    ////////////////////////////////////////////////////////
    assert(CigarHelper::softClipBeginByRefPos(record, 31, newCigar, newPos) == 23);
    assert(newPos == 10);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H24S3H") == 0);

    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    // Test clipping at every position when insertions & deletions
    // are next to each other.
    origCigar = "3M3D3I3M";
    record.setCigar(origCigar);
    record.setSequence("GGGAAAGGG");
    // Cigar:   MMMDDDIIIMMM
    // ReadPos: 000   000000
    // ReadPos: 012   345678
    // RefPos:  111111   111
    // RefPos:  012345   678
    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 9, newCigar, newPos) == 
           CigarHelper::NO_CLIP);
    assert(newPos == 10);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), origCigar) == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 10, newCigar, newPos) == 0);
    assert(newPos == 11);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "1S2M3D3I3M") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 11, newCigar, newPos) == 1);
    assert(newPos == 12);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "2S1M3D3I3M") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 12, newCigar, newPos) == 5);
    assert(newPos == 16);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "6S3M") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 13, newCigar, newPos) == 5);
    assert(newPos == 16);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "6S3M") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 14, newCigar, newPos) == 5);
    assert(newPos == 16);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "6S3M") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 15, newCigar, newPos) == 5);
    assert(newPos == 16);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "6S3M") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 16, newCigar, newPos) == 6);
    assert(newPos == 17);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "7S2M") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 17, newCigar, newPos) == 7);
    assert(newPos == 18);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "8S1M") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 18, newCigar, newPos) == 8);
    assert(newPos == 10);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "9S") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 19, newCigar, newPos) == 8);
    assert(newPos == 10);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "9S") == 0);

    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    // Test clipping at every position when first non-clip instruction is delete.
    origCigar = "3H3S3D3M3S3H";
    record.setCigar(origCigar);
    record.setSequence("gggAAAggg");
    // Cigar:   HHHSSSDDDMMMSSSHHH
    // ReadPos:    000   000000
    // ReadPos:    012   345678
    // RefPos:        111111
    // RefPos:        012345
    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 9, newCigar, newPos) == 
           CigarHelper::NO_CLIP);
    assert(newPos == 10);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), origCigar) == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 10, newCigar, newPos) == 2);
    assert(newPos == 13);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 11, newCigar, newPos) == 2);
    assert(newPos == 13);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 12, newCigar, newPos) == 2);
    assert(newPos == 13);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 13, newCigar, newPos) == 3);
    assert(newPos == 14);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H4S2M3S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 14, newCigar, newPos) == 4);
    assert(newPos == 15);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H5S1M3S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 15, newCigar, newPos) == 8);
    assert(newPos == 10);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H9S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 16, newCigar, newPos) == 8);
    assert(newPos == 10);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H9S3H") == 0);

    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    // Test clipping at every position when first non-clip instruction is insert.
    origCigar = "3H3S3I3M3S3H";
    record.setCigar(origCigar);
    record.setSequence("gggAAATTTggg");
    // Cigar:   HHHSSSIIIMMMSSSHHH
    // ReadPos:    000000000011
    // ReadPos:    012345678901
    // RefPos:           111
    // RefPos:           012
    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 9, newCigar, newPos) == 
           CigarHelper::NO_CLIP);
    assert(newPos == 10);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), origCigar) == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 10, newCigar, newPos) == 6);
    assert(newPos == 11);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H7S2M3S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 11, newCigar, newPos) == 7);
    assert(newPos == 12);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H8S1M3S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 12, newCigar, newPos) == 11);
    assert(newPos == 10);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H12S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipBeginByRefPos(record, 13, newCigar, newPos) == 11);
    assert(newPos == 10);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H12S3H") == 0);
}


void CigarHelperTest::testSoftClipEndByRefPos()
{
    SamRecord record;
    CigarRoller newCigar;
    std::string newCigarString;

    // Setup the current Cigar.
    // Cigar:   HHHSSSMMMDDDMMMIIIMMMPPPMMMDDDMMMSSSHHH
    // ReadPos:    000000   000011111   111   112222
    // ReadPos:    012345   678901234   567   890123
    // RefPos:        111111111   122   222222223
    // RefPos:        012345678   901   234567890
    const char* origCigar = "3H3S3M3D3M3I3M3P3M3D3M3S3H";
    record.setCigar(origCigar);
    record.set0BasedPosition(10);
    record.setSequence("gggAAATTTCCCTTTGGGAAAggg");

    ////////////////////////////////////////////////////////
    // Clip outside of the range (after).  Nothing should change.
    assert(CigarHelper::softClipEndByRefPos(record, 10000, newCigar) == 
           CigarHelper::NO_CLIP);
    newCigar.getCigarString(newCigarString);
    assert(strcmp(newCigarString.c_str(), origCigar) == 0);

    ////////////////////////////////////////////////////////
    // Clip outside of the range (before).  Everything should be clipped.
    assert(CigarHelper::softClipEndByRefPos(record, 1, newCigar) == 0);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H24S3H") == 0);


    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    // Test clipping at every position of the read.

    ////////////////////////////////////////////////////////
    // Clip at the first position.
    assert(CigarHelper::softClipEndByRefPos(record, 10, newCigar) == 0);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H24S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip in the middle of the first Match.
    assert(CigarHelper::softClipEndByRefPos(record, 11, newCigar) == 4);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S1M20S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip just before the first deletion.
    assert(CigarHelper::softClipEndByRefPos(record, 12, newCigar) == 5);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S2M19S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip at the first deletion.
    assert(CigarHelper::softClipEndByRefPos(record, 13, newCigar) == 6);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M18S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip in the middle of the first deletion.
    assert(CigarHelper::softClipEndByRefPos(record, 14, newCigar) == 6);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M18S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip in the end of the first deletion.
    assert(CigarHelper::softClipEndByRefPos(record, 15, newCigar) == 6);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M18S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip just after the first deletion (should remove the deletion).
    assert(CigarHelper::softClipEndByRefPos(record, 16, newCigar) == 6);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M18S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip in middle of read after 1st deletion.
    assert(CigarHelper::softClipEndByRefPos(record, 17, newCigar) == 7);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3D1M17S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip in middle of read after 1st deletion.
    assert(CigarHelper::softClipEndByRefPos(record, 18, newCigar) == 8);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3D2M16S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip just after the 1st insertion.
    assert(CigarHelper::softClipEndByRefPos(record, 19, newCigar) == 12);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3D3M3I12S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip in middle of the match after 1st insertion.
    assert(CigarHelper::softClipEndByRefPos(record, 20, newCigar) == 13);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3D3M3I1M11S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip in middle of the match after 1st insertion.
    assert(CigarHelper::softClipEndByRefPos(record, 21, newCigar) == 14);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3D3M3I2M10S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip right after the pad
    assert(CigarHelper::softClipEndByRefPos(record, 22, newCigar) == 15);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3D3M3I3M9S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip middle of read after the pad
    assert(CigarHelper::softClipEndByRefPos(record, 23, newCigar) == 16);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3D3M3I3M3P1M8S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip end of read after the pad before deletion
    assert(CigarHelper::softClipEndByRefPos(record, 24, newCigar) == 17);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3D3M3I3M3P2M7S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip at start of 2nd deletion.
    assert(CigarHelper::softClipEndByRefPos(record, 25, newCigar) == 18);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3D3M3I3M3P3M6S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip in 2nd deletion.
    assert(CigarHelper::softClipEndByRefPos(record, 26, newCigar) == 18);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3D3M3I3M3P3M6S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip in 2nd deletion.
    assert(CigarHelper::softClipEndByRefPos(record, 27, newCigar) == 18);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3D3M3I3M3P3M6S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip right after 2nd deletion.
    assert(CigarHelper::softClipEndByRefPos(record, 28, newCigar) == 18);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3D3M3I3M3P3M6S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip in middle of last match.
    assert(CigarHelper::softClipEndByRefPos(record, 29, newCigar) == 19);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3D3M3I3M3P3M3D1M5S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip in middle of last match.
    assert(CigarHelper::softClipEndByRefPos(record, 30, newCigar) == 20);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3M3D3M3I3M3P3M3D2M4S3H") == 0);

    ////////////////////////////////////////////////////////
    // Clip right after the read (no change).
    assert(CigarHelper::softClipEndByRefPos(record, 31, newCigar) == 
           CigarHelper::NO_CLIP);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), origCigar) == 0);

    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    // Test clipping at every position when insertions & deletions
    // are next to each other.
    origCigar = "3M3D3I3M";
    record.setCigar(origCigar);
    record.setSequence("GGGAAAGGG");
    // Cigar:   MMMDDDIIIMMM
    // ReadPos: 000   000000
    // ReadPos: 012   345678
    // RefPos:  111111   111
    // RefPos:  012345   678
    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 9, newCigar) == 0);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "9S") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 10, newCigar) == 0);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "9S") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 11, newCigar) == 1);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "1M8S") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 12, newCigar) == 2);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "2M7S") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 13, newCigar) == 3);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3M6S") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 14, newCigar) == 3);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3M6S") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 15, newCigar) == 3);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3M6S") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 16, newCigar) == 6);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3M3D3I3S") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 17, newCigar) == 7);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3M3D3I1M2S") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 18, newCigar) == 8);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3M3D3I2M1S") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 19, newCigar) == 
           CigarHelper::NO_CLIP);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), origCigar) == 0);

    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    // Test clipping at every position when first non-clip instruction is delete.
    origCigar = "3H3S3D3M3S3H";
    record.setCigar(origCigar);
    record.setSequence("gggAAAggg");
    // Cigar:   HHHSSSDDDMMMSSSHHH
    // ReadPos:    000   000000
    // ReadPos:    012   345678
    // RefPos:        111111
    // RefPos:        012345
    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 9, newCigar) == 0);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H9S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 10, newCigar) == 0);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H9S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 11, newCigar) == 0);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H9S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 12, newCigar) == 0);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H9S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 13, newCigar) == 0);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H9S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 14, newCigar) == 4);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3D1M5S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 15, newCigar) == 5);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3D2M4S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 16, newCigar) == 
           CigarHelper::NO_CLIP);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), origCigar) == 0);

    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    // Test clipping at every position when first non-clip instruction is insert.
    origCigar = "3H3S3I3M3S3H";
    record.setCigar(origCigar);
    record.setSequence("gggAAATTTggg");
    // Cigar:   HHHSSSIIIMMMSSSHHH
    // ReadPos:    000000000011
    // ReadPos:    012345678901
    // RefPos:           111
    // RefPos:           012
    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 9, newCigar) == 0);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H12S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 10, newCigar) == 6);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3I6S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 11, newCigar) == 7);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3I1M5S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 12, newCigar) == 8);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), "3H3S3I2M4S3H") == 0);

    record.setCigar(origCigar);
    assert(CigarHelper::softClipEndByRefPos(record, 13, newCigar) == 
           CigarHelper::NO_CLIP);
    newCigar.getCigarString(newCigarString);
    //std::cout << newCigarString.c_str() << std::endl;
    assert(strcmp(newCigarString.c_str(), origCigar) == 0);
}
