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


// put TEST below here, so that makedepend will see the .h, so that we
// can get a clean dependency for SmithWaterman.o, so that we can at least
// compile the header when we change it.


#include <getopt.h>
#include "Generic.h"
#include "CigarRollerTest.h"

//
// Test the obvious cases.
// Add non-obvious ones as bugs come up.
//
int CigarRollerTest::test(void)
{
    int failures = 0, testNum = 0;
    Cigar::CigarOperator op;

    //   const char *str;
    ////////////////////////////////////
    // Test foundInReference static methods.
    check(failures, ++testNum, "foundInReference(none)", false,
          Cigar::foundInReference(Cigar::none));
    check(failures, ++testNum, "foundInReference(match)", true, 
          Cigar::foundInReference(Cigar::match));
    check(failures, ++testNum, "foundInReference(mismatch)", true, 
          Cigar::foundInReference(Cigar::mismatch));
    check(failures, ++testNum, "foundInReference(insert)", false, 
          Cigar::foundInReference(Cigar::insert));
    check(failures, ++testNum, "foundInReference(del)", true, 
          Cigar::foundInReference(Cigar::del));
    check(failures, ++testNum, "foundInReference(skip)", true, 
          Cigar::foundInReference(Cigar::skip));
    check(failures, ++testNum, "foundInReference(softClip)", false, 
          Cigar::foundInReference(Cigar::softClip));
    check(failures, ++testNum, "foundInReference(hardClip)", false, 
          Cigar::foundInReference(Cigar::hardClip));
    check(failures, ++testNum, "foundInReference(pad)", false, 
          Cigar::foundInReference(Cigar::pad));

    check(failures, ++testNum, 
          "foundInReference('?')", false, Cigar::foundInReference('?'));
    check(failures, ++testNum, 
          "foundInReference('z')", false, Cigar::foundInReference('z'));
    check(failures, ++testNum, 
          "foundInReference('M')", true, Cigar::foundInReference('M'));
    check(failures, ++testNum, 
          "foundInReference('X')", true, Cigar::foundInReference('X'));
    check(failures, ++testNum, 
          "foundInReference('=')", true, Cigar::foundInReference('='));
    check(failures, ++testNum, 
          "foundInReference('I')", false, Cigar::foundInReference('I'));
    check(failures, ++testNum, 
          "foundInReference('D')", true, Cigar::foundInReference('D'));
    check(failures, ++testNum, 
          "foundInReference('N')", true, Cigar::foundInReference('N'));
    check(failures, ++testNum, 
          "foundInReference('S')", false, Cigar::foundInReference('S'));
    check(failures, ++testNum, 
          "foundInReference('H')", false, Cigar::foundInReference('H'));
    check(failures, ++testNum, 
          "foundInReference('P')", false, Cigar::foundInReference('P'));

    op.count = 3;
    op.operation = Cigar::none;
    check(failures, ++testNum, "foundInReference(none)", false,
          Cigar::foundInReference(op));
    op.operation = Cigar::match;
    check(failures, ++testNum, "foundInReference(match)", true,
          Cigar::foundInReference(op));
    op.operation = Cigar::mismatch;
    check(failures, ++testNum, "foundInReference(mismatch)", true,
          Cigar::foundInReference(op));
    op.operation = Cigar::insert;
    check(failures, ++testNum, "foundInReference(insert)", false,
          Cigar::foundInReference(op));
    op.operation = Cigar::del;
    check(failures, ++testNum, "foundInReference(del)", true,
          Cigar::foundInReference(op));
    op.operation = Cigar::skip;
    check(failures, ++testNum, "foundInReference(skip)", true,
          Cigar::foundInReference(op));
    op.operation = Cigar::softClip;
    check(failures, ++testNum, "foundInReference(softClip)", false,
          Cigar::foundInReference(op));
    op.operation = Cigar::hardClip;
    check(failures, ++testNum, "foundInReference(hardClip)", false,
          Cigar::foundInReference(op));
    op.operation = Cigar::pad;
    check(failures, ++testNum, "foundInReference(pad)", false,
          Cigar::foundInReference(op));

    ////////////////////////////////////
    // Test foundInQuery static methods.
    check(failures, ++testNum, "foundInQuery(none)", false,
          Cigar::foundInQuery(Cigar::none));
    check(failures, ++testNum, "foundInQuery(match)", true, 
          Cigar::foundInQuery(Cigar::match));
    check(failures, ++testNum, "foundInQuery(mismatch)", true, 
          Cigar::foundInQuery(Cigar::mismatch));
    check(failures, ++testNum, "foundInQuery(insert)", true, 
          Cigar::foundInQuery(Cigar::insert));
    check(failures, ++testNum, "foundInQuery(del)", false, 
          Cigar::foundInQuery(Cigar::del));
    check(failures, ++testNum, "foundInQuery(skip)", false, 
          Cigar::foundInQuery(Cigar::skip));
    check(failures, ++testNum, "foundInQuery(softClip)", true, 
          Cigar::foundInQuery(Cigar::softClip));
    check(failures, ++testNum, "foundInQuery(hardClip)", false, 
          Cigar::foundInQuery(Cigar::hardClip));
    check(failures, ++testNum, "foundInQuery(pad)", false, 
          Cigar::foundInQuery(Cigar::pad));

    check(failures, ++testNum, 
          "foundInQuery('?')", false, Cigar::foundInQuery('?'));
    check(failures, ++testNum, 
          "foundInQuery('z')", false, Cigar::foundInQuery('z'));
    check(failures, ++testNum, 
          "foundInQuery('M')", true, Cigar::foundInQuery('M'));
    check(failures, ++testNum, 
          "foundInQuery('X')", true, Cigar::foundInQuery('X'));
    check(failures, ++testNum, 
          "foundInQuery('=')", true, Cigar::foundInQuery('='));
    check(failures, ++testNum, 
          "foundInQuery('I')", true, Cigar::foundInQuery('I'));
    check(failures, ++testNum, 
          "foundInQuery('D')", false, Cigar::foundInQuery('D'));
    check(failures, ++testNum, 
          "foundInQuery('N')", false, Cigar::foundInQuery('N'));
    check(failures, ++testNum, 
          "foundInQuery('S')", true, Cigar::foundInQuery('S'));
    check(failures, ++testNum, 
          "foundInQuery('H')", false, Cigar::foundInQuery('H'));
    check(failures, ++testNum, 
          "foundInQuery('P')", false, Cigar::foundInQuery('P'));

    op.count = 3;
    op.operation = Cigar::none;
    check(failures, ++testNum, "foundInQuery(none)", false,
          Cigar::foundInQuery(op));
    op.operation = Cigar::match;
    check(failures, ++testNum, "foundInQuery(match)", true,
          Cigar::foundInQuery(op));
    op.operation = Cigar::mismatch;
    check(failures, ++testNum, "foundInQuery(mismatch)", true,
          Cigar::foundInQuery(op));
    op.operation = Cigar::insert;
    check(failures, ++testNum, "foundInQuery(insert)", true,
          Cigar::foundInQuery(op));
    op.operation = Cigar::del;
    check(failures, ++testNum, "foundInQuery(del)", false,
          Cigar::foundInQuery(op));
    op.operation = Cigar::skip;
    check(failures, ++testNum, "foundInQuery(skip)", false,
          Cigar::foundInQuery(op));
    op.operation = Cigar::softClip;
    check(failures, ++testNum, "foundInQuery(softClip)", true,
          Cigar::foundInQuery(op));
    op.operation = Cigar::hardClip;
    check(failures, ++testNum, "foundInQuery(hardClip)", false,
          Cigar::foundInQuery(op));
    op.operation = Cigar::pad;
    check(failures, ++testNum, "foundInQuery(pad)", false,
          Cigar::foundInQuery(op));


    ////////////////////////////////////
    // Test isClip static methods.
    check(failures, ++testNum, "isClip(none)", false,
          Cigar::isClip(Cigar::none));
    check(failures, ++testNum, "isClip(match)", false, 
          Cigar::isClip(Cigar::match));
    check(failures, ++testNum, "isClip(mismatch)", false, 
          Cigar::isClip(Cigar::mismatch));
    check(failures, ++testNum, "isClip(insert)", false, 
          Cigar::isClip(Cigar::insert));
    check(failures, ++testNum, "isClip(del)", false, 
          Cigar::isClip(Cigar::del));
    check(failures, ++testNum, "isClip(skip)", false, 
          Cigar::isClip(Cigar::skip));
    check(failures, ++testNum, "isClip(softClip)", true, 
          Cigar::isClip(Cigar::softClip));
    check(failures, ++testNum, "isClip(hardClip)", true, 
          Cigar::isClip(Cigar::hardClip));
    check(failures, ++testNum, "isClip(pad)", false, 
          Cigar::isClip(Cigar::pad));

    check(failures, ++testNum, 
          "isClip('?')", false, Cigar::isClip('?'));
    check(failures, ++testNum, 
          "isClip('z')", false, Cigar::isClip('z'));
    check(failures, ++testNum, 
          "isClip('M')", false, Cigar::isClip('M'));
    check(failures, ++testNum, 
          "isClip('X')", false, Cigar::isClip('X'));
    check(failures, ++testNum, 
          "isClip('=')", false, Cigar::isClip('='));
    check(failures, ++testNum, 
          "isClip('I')", false, Cigar::isClip('I'));
    check(failures, ++testNum, 
          "isClip('D')", false, Cigar::isClip('D'));
    check(failures, ++testNum, 
          "isClip('N')", false, Cigar::isClip('N'));
    check(failures, ++testNum, 
          "isClip('S')", true, Cigar::isClip('S'));
    check(failures, ++testNum, 
          "isClip('H')", true, Cigar::isClip('H'));
    check(failures, ++testNum, 
          "isClip('P')", false, Cigar::isClip('P'));

    op.count = 3;
    op.operation = Cigar::none;
    check(failures, ++testNum, "isClip(none)", false,
          Cigar::isClip(op));
    op.operation = Cigar::match;
    check(failures, ++testNum, "isClip(match)", false,
          Cigar::isClip(op));
    op.operation = Cigar::mismatch;
    check(failures, ++testNum, "isClip(mismatch)", false,
          Cigar::isClip(op));
    op.operation = Cigar::insert;
    check(failures, ++testNum, "isClip(insert)", false,
          Cigar::isClip(op));
    op.operation = Cigar::del;
    check(failures, ++testNum, "isClip(del)", false,
          Cigar::isClip(op));
    op.operation = Cigar::skip;
    check(failures, ++testNum, "isClip(skip)", false,
          Cigar::isClip(op));
    op.operation = Cigar::softClip;
    check(failures, ++testNum, "isClip(softClip)", true,
          Cigar::isClip(op));
    op.operation = Cigar::hardClip;
    check(failures, ++testNum, "isClip(hardClip)", true,
          Cigar::isClip(op));
    op.operation = Cigar::pad;
    check(failures, ++testNum, "isClip(pad)", false,
          Cigar::isClip(op));






    ///////////////////////////////

    // Create the CigarRoller.
    CigarRoller cigar;

    //   const char *str;
    String str;
    std::string result;
    
    cigar.getCigarString(str);
    result = str.c_str();
    //    result = str = cigar.getString(); delete str;

    check(failures, ++testNum, "constructor", result, "");    // test empty case

    op.operation = CigarRoller::match;
    op.count = 10;

    cigar += op;

    op.count=5;

    cigar += op;

    op.count=5; op.operation = CigarRoller::mismatch;   // test that match/mismatch get combined

    cigar += op;

    op.count=5; op.operation = CigarRoller::insert;

    cigar += op;

    op.count=5; op.operation = CigarRoller::insert;

    cigar += op;

    op.count=5; op.operation = CigarRoller::del;

    cigar += op;

    op.count=5; op.operation = CigarRoller::mismatch;

    cigar += op;

    op.count=5; op.operation = CigarRoller::match;

    cigar += op;

    op.count=5; op.operation = CigarRoller::skip;

    cigar += op;

    op.count=5; op.operation = CigarRoller::match;

    cigar += op;

    op.count=2; op.operation = CigarRoller::pad;

    cigar += op;

    op.count=3; op.operation = CigarRoller::match;

    cigar += op;

    cigar.getCigarString(str);
    result = str.c_str();
    //    result = str = cigar.getString();  delete str;

    check(failures, ++testNum, "match combining", "20M10I5D10M5N5M2P3M", result);
    check(failures, ++testNum, "length check", 8, cigar.size());

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Test getRefOffset, getQueryIndex, getRefPosition, & getQueryIndex(that takes ref position)

    // Validate the reference offsets when passing in a query index,
    // and the query offsets when passing in a query index.
    // Spot check the offsets out of order - to verify order doesn't matter.
    check(failures, ++testNum, "getRefOffset(20)", -1, cigar.getRefOffset(20));
    check(failures, ++testNum, "getRefOffset(30)", 25, cigar.getRefOffset(30));
    check(failures, ++testNum, "getRefOffset(46)", 46, cigar.getRefOffset(46));
    check(failures, ++testNum, "getRefOffset(0)", 0, cigar.getRefOffset(0));
    check(failures, ++testNum, "getRefPosition(20, 5)", -1, cigar.getRefPosition(20, 5));
    check(failures, ++testNum, "getRefPosition(30, 5)", 30, cigar.getRefPosition(30, 5));
    check(failures, ++testNum, "getRefPosition(46, 5)", 51, cigar.getRefPosition(46, 5));
    check(failures, ++testNum, "getRefPosition(0, 5)", 5, cigar.getRefPosition(0, 5));
    check(failures, ++testNum, "getQueryIndex(20)", CigarRoller::INDEX_NA, cigar.getQueryIndex(20));
    check(failures, ++testNum, "getQueryIndex(30)", 35, cigar.getQueryIndex(30));
    check(failures, ++testNum, "getQueryIndex(35)", Cigar::INDEX_NA, cigar.getQueryIndex(35));
    check(failures, ++testNum, "getQueryIndex(46)", 46, cigar.getQueryIndex(46));
    check(failures, ++testNum, "getQueryIndex(0)", 0, cigar.getQueryIndex(0));
    check(failures, ++testNum, "getQueryIndex(25, 5)", -1, cigar.getQueryIndex(20));
    check(failures, ++testNum, "getQueryIndex(35, 5)", 35, cigar.getQueryIndex(30));
    check(failures, ++testNum, "getQueryIndex(40, 5)", -1, cigar.getQueryIndex(35));
    check(failures, ++testNum, "getQueryIndex(51, 5)", 46, cigar.getQueryIndex(46));
    check(failures, ++testNum, "getQueryIndex(5, 5)", 0, cigar.getQueryIndex(0));

    int i = 0;
    int queryIndex = 0;
    int refOffset = 0;
    // Validate the 20 matches.
    for(i = 0; i < 20; i++)
    {
        check(failures, ++testNum, "getRefOffset(queryIndex)", refOffset, cigar.getRefOffset(queryIndex));
        check(failures, ++testNum, "getQueryIndex(refOffset)", queryIndex, cigar.getQueryIndex(refOffset));
        check(failures, ++testNum, "getRefPosition(queryIndex, 5)", refOffset + 5, cigar.getRefPosition(queryIndex, 5));
        check(failures, ++testNum, "getQueryIndex(refPosition, 5)", queryIndex, cigar.getQueryIndex(refOffset + 5, 5));
        ++queryIndex;
        ++refOffset;
    }
    // Validate the 10 Inserts - exist in query, but not in the reference.
    for(i = 0; i < 10; i++)
    {
        check(failures, ++testNum, "getRefOffset(queryIndex)", -1, cigar.getRefOffset(queryIndex));
        check(failures, ++testNum, "getRefPosition(queryIndex, 5)", -1, cigar.getRefPosition(queryIndex, 5));
        ++queryIndex;
    }
    // Validate the 5 Deletions - exist in the reference, but not the query.
    for(i = 0; i < 5; i++)
    {
        check(failures, ++testNum, "getQueryIndex(refOffset)", -1, cigar.getQueryIndex(refOffset));
        check(failures, ++testNum, "getQueryIndex(refPosition, 5)", -1, cigar.getQueryIndex(refOffset + 5, 5));
        refOffset++;
    }
    // Validate the 10 matches.
    for(i = 0; i < 10; i++)
    {
        check(failures, ++testNum, "getRefOffset(queryIndex)", refOffset, cigar.getRefOffset(queryIndex));
        check(failures, ++testNum, "getQueryIndex(refOffset)", queryIndex, cigar.getQueryIndex(refOffset));
        check(failures, ++testNum, "getRefPosition(queryIndex, 5)", refOffset + 5, cigar.getRefPosition(queryIndex, 5));
        check(failures, ++testNum, "getQueryIndex(refPosition, 5)", queryIndex, cigar.getQueryIndex(refOffset + 5, 5));
        ++queryIndex;
        ++refOffset;
    }
    // Validate the 5 Skips - exist in the reference, but not the query.
    for(i = 0; i < 5; i++)
    {
        check(failures, ++testNum, "getQueryIndex(refOffset)", -1, cigar.getQueryIndex(refOffset));
        check(failures, ++testNum, "getQueryIndex(refPosition, 5)", -1, cigar.getQueryIndex(refOffset + 5, 5));
        refOffset++;
    }
    // Validate the 5 matches.
    for(i = 0; i < 5; i++)
    {
        check(failures, ++testNum, "getRefOffset(queryIndex)", refOffset, cigar.getRefOffset(queryIndex));
        check(failures, ++testNum, "getQueryIndex(refOffset)", queryIndex, cigar.getQueryIndex(refOffset));
        check(failures, ++testNum, "getRefPosition(queryIndex, 5)", refOffset + 5, cigar.getRefPosition(queryIndex, 5));
        check(failures, ++testNum, "getQueryIndex(refPosition, 5)", queryIndex, cigar.getQueryIndex(refOffset + 5, 5));
        ++queryIndex;
        ++refOffset;
    }
    // Nothing to validate for the 2 pads since they do not show up in either the reference or the query.
    // Validate the 3 matches.
    for(i = 0; i < 3; i++)
    {
        check(failures, ++testNum, "getRefOffset(queryIndex)", refOffset, cigar.getRefOffset(queryIndex));
        check(failures, ++testNum, "getQueryIndex(refOffset)", queryIndex, cigar.getQueryIndex(refOffset));
        check(failures, ++testNum, "getRefPosition(queryIndex, 5)", refOffset + 5, cigar.getRefPosition(queryIndex, 5));
        check(failures, ++testNum, "getQueryIndex(refPosition, 5)", queryIndex, cigar.getQueryIndex(refOffset + 5, 5));
        ++queryIndex;
        ++refOffset;
    }

    // Validate that if you go beyond the end, -1 is returned.
    check(failures, ++testNum, "getRefOffset(queryIndex)", -1, cigar.getRefOffset(queryIndex));
    check(failures, ++testNum, "getQueryIndex(refOffset)", -1, cigar.getQueryIndex(refOffset));
    check(failures, ++testNum, "getRefPosition(queryIndex, 5)", -1, cigar.getRefPosition(queryIndex, 5));
    check(failures, ++testNum, "getQueryIndex(refPosition, 5)", -1, cigar.getQueryIndex(refOffset + 5, 5));
    ++queryIndex;
    ++refOffset;
    check(failures, ++testNum, "getRefOffset(queryIndex)", -1, cigar.getRefOffset(queryIndex));
    check(failures, ++testNum, "getQueryIndex(refOffset)", -1, cigar.getQueryIndex(refOffset));
    check(failures, ++testNum, "getRefPosition(queryIndex, 5)", -1, cigar.getRefPosition(queryIndex, 5));
    check(failures, ++testNum, "getQueryIndex(refPosition, 5)", -1, cigar.getQueryIndex(refOffset + 5, 5));

    ////////////////////////////////////////////////////////////////////////
    // Test getNumOverlaps
    
    // When query starts at position 5:
    // Overlaps are at 5-24, 30-39, 45-49, 50-52

    // Test with region [1-5) where query starts at position 5 = 0 overlaps.
    check(failures, ++testNum, "getNumOverlaps(1, 5, 5)", (uint32_t)0, cigar.getNumOverlaps(1, 5, 5));

    // Test with region [53-101) where query starts at position 5 = 0 overlaps.
    check(failures, ++testNum, "getNumOverlaps(53, 101, 5)", (uint32_t)0, cigar.getNumOverlaps(53, 101, 5));

    // Test with region [53-10) where query starts at position 5 = 0 overlaps.
    check(failures, ++testNum, "getNumOverlaps(53, 10, 5)", (uint32_t)0, cigar.getNumOverlaps(53, 10, 5));

    // Test with region [35-10) where query starts at position 5 = 0 overlaps.
    check(failures, ++testNum, "getNumOverlaps(35, 10, 5)", (uint32_t)0, cigar.getNumOverlaps(35, 10, 5));

    // Test with region [35-1) where query starts at position 5 = 0 overlaps.
    check(failures, ++testNum, "getNumOverlaps(35, 1, 5)", (uint32_t)0, cigar.getNumOverlaps(35, 1, 5));

    // Test with region [1-6) where query starts at position 5 - 1 overlap = pos 5.
    check(failures, ++testNum, "getNumOverlaps(1, 6, 5)", (uint32_t)1, cigar.getNumOverlaps(1, 6, 5));

    // Test with region [25-30) where query has only DELETIONS from the reference = 0 overlaps.
    check(failures, ++testNum, "getNumOverlaps(25, 30, 5)", (uint32_t)0, cigar.getNumOverlaps(25, 30, 5));

    // Test with region [24-30) where query has only a match at position 24 = 1 overlap.
    check(failures, ++testNum, "getNumOverlaps(24, 30, 5)", (uint32_t)1, cigar.getNumOverlaps(24, 30, 5));

    // Test with region [25-31) where query has only a match at position 30 = 1 overlap.
    check(failures, ++testNum, "getNumOverlaps(25, 31, 5)", (uint32_t)1, cigar.getNumOverlaps(25, 31, 5));

    // Test with region [1-31), match pos 5-24 & 30 = 21 overlaps
    check(failures, ++testNum, "getNumOverlaps(1, 31, 5)", (uint32_t)21, cigar.getNumOverlaps(1, 31, 5));

    // Test a region that covers the entire range [1-101), match pos 5-24, 30-39, 45-49, & 50-52 = 38 overlaps
    check(failures, ++testNum, "getNumOverlaps(1, 101, 5)", (uint32_t)38, cigar.getNumOverlaps(1, 101, 5));

    // Test a region that covers the entire range [-1--1), (whole region) match pos 5-24, 30-39, 45-49, & 50-52 = 38 overlaps
    check(failures, ++testNum, "getNumOverlaps(-1, -1, 5)", (uint32_t)38, cigar.getNumOverlaps(-1, -1, 5));

    // Test a region that covers the entire range [6-52), match pos 6-24, 30-39, 45-49, & 50-51 = 36 overlaps
    check(failures, ++testNum, "getNumOverlaps(6, 52, 5)", (uint32_t)36, cigar.getNumOverlaps(6, 52, 5));

    // Test with region [40-45) where query has only SKIPS from the reference = 0 overlaps.
    check(failures, ++testNum, "getNumOverlaps(40, 45, 5)", (uint32_t)0, cigar.getNumOverlaps(40, 45, 5));

    // Test a region that covers the range [-1-10), (whole region) match pos 5-9 = 5 overlaps
    check(failures, ++testNum, "getNumOverlaps(-1, 10, 5)", (uint32_t)5, cigar.getNumOverlaps(-1, 10, 5));

    // Test a region that covers the range [50--1), (whole region) match pos 50-52 = 3 overlaps
    check(failures, ++testNum, "getNumOverlaps(50, -1, 5)", (uint32_t)3, cigar.getNumOverlaps(50, -1, 5));

    ////////////////////////////////////////////////////////////////////////////
    // Test a new CIGAR.
    cigar.Set("4M10N4M3I2M4D3M");
    String expectedResult = "4M10N4M3I2M4D3M";
    String cigarString = "HI";
    cigar.getCigarString(cigarString);
    check(failures, ++testNum, "getCigarString", expectedResult, cigarString);

    // 4M10N4M3I2M4D3M
    //                     11111111112222222222
    // ExtIndex: 012345678901234567890123456789
    // ExtCigar: MMMMNNNNNNNNNNMMMMIIIMMDDDDMMM
    //                               111    111
    // QueryIndx:0123          456789012    345
    //                     11111111   112222222
    // RefOffset:012345678901234567   890123456
    //                1111111111222   222222233
    // RefPos:   567890123456789012   345678901

    // Test getExpandedCigarIndex & getCigar op
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(0, 5)", -1, cigar.getExpandedCigarIndexFromRefPos(0,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(1, 5)", -1, cigar.getExpandedCigarIndexFromRefPos(1,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(2, 5)", -1, cigar.getExpandedCigarIndexFromRefPos(2,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(3, 5)", -1, cigar.getExpandedCigarIndexFromRefPos(3,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(4, 5)", -1, cigar.getExpandedCigarIndexFromRefPos(4,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(5, 5)", 0, cigar.getExpandedCigarIndexFromRefPos(5,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(6, 5)", 1, cigar.getExpandedCigarIndexFromRefPos(6,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(7, 5)", 2, cigar.getExpandedCigarIndexFromRefPos(7,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(8, 5)", 3, cigar.getExpandedCigarIndexFromRefPos(8,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(9, 5)", 4, cigar.getExpandedCigarIndexFromRefPos(9,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(10, 5)", 5, cigar.getExpandedCigarIndexFromRefPos(10,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(11, 5)", 6, cigar.getExpandedCigarIndexFromRefPos(11,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(12, 5)", 7, cigar.getExpandedCigarIndexFromRefPos(12,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(13, 5)", 8, cigar.getExpandedCigarIndexFromRefPos(13,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(14, 5)", 9, cigar.getExpandedCigarIndexFromRefPos(14,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(15, 5)", 10, cigar.getExpandedCigarIndexFromRefPos(15,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(16, 5)", 11, cigar.getExpandedCigarIndexFromRefPos(16,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(17, 5)", 12, cigar.getExpandedCigarIndexFromRefPos(17,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(18, 5)", 13, cigar.getExpandedCigarIndexFromRefPos(18,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(19, 5)", 14, cigar.getExpandedCigarIndexFromRefPos(19,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(20, 5)", 15, cigar.getExpandedCigarIndexFromRefPos(20,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(21, 5)", 16, cigar.getExpandedCigarIndexFromRefPos(21,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(22, 5)", 17, cigar.getExpandedCigarIndexFromRefPos(22,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(23, 5)", 21, cigar.getExpandedCigarIndexFromRefPos(23,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(24, 5)", 22, cigar.getExpandedCigarIndexFromRefPos(24,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(25, 5)", 23, cigar.getExpandedCigarIndexFromRefPos(25,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(26, 5)", 24, cigar.getExpandedCigarIndexFromRefPos(26,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(27, 5)", 25, cigar.getExpandedCigarIndexFromRefPos(27,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(28, 5)", 26, cigar.getExpandedCigarIndexFromRefPos(28,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(29, 5)", 27, cigar.getExpandedCigarIndexFromRefPos(29,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(30, 5)", 28, cigar.getExpandedCigarIndexFromRefPos(30,5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefPos(31, 5)", 29, cigar.getExpandedCigarIndexFromRefPos(31,5));

    // 
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(0)", 0, cigar.getExpandedCigarIndexFromRefOffset(0));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(1)", 1, cigar.getExpandedCigarIndexFromRefOffset(1));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(2)", 2, cigar.getExpandedCigarIndexFromRefOffset(2));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(3)", 3, cigar.getExpandedCigarIndexFromRefOffset(3));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(4)", 4, cigar.getExpandedCigarIndexFromRefOffset(4));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(5)", 5, cigar.getExpandedCigarIndexFromRefOffset(5));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(6)", 6, cigar.getExpandedCigarIndexFromRefOffset(6));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(7)", 7, cigar.getExpandedCigarIndexFromRefOffset(7));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(8)", 8, cigar.getExpandedCigarIndexFromRefOffset(8));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(9)", 9, cigar.getExpandedCigarIndexFromRefOffset(9));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(10)", 10, cigar.getExpandedCigarIndexFromRefOffset(10));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(11)", 11, cigar.getExpandedCigarIndexFromRefOffset(11));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(12)", 12, cigar.getExpandedCigarIndexFromRefOffset(12));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(13)", 13, cigar.getExpandedCigarIndexFromRefOffset(13));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(14)", 14, cigar.getExpandedCigarIndexFromRefOffset(14));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(15)", 15, cigar.getExpandedCigarIndexFromRefOffset(15));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(16)", 16, cigar.getExpandedCigarIndexFromRefOffset(16));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(17)", 17, cigar.getExpandedCigarIndexFromRefOffset(17));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(18)", 21, cigar.getExpandedCigarIndexFromRefOffset(18));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(19)", 22, cigar.getExpandedCigarIndexFromRefOffset(19));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(20)", 23, cigar.getExpandedCigarIndexFromRefOffset(20));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(21)", 24, cigar.getExpandedCigarIndexFromRefOffset(21));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(22)", 25, cigar.getExpandedCigarIndexFromRefOffset(22));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(23)", 26, cigar.getExpandedCigarIndexFromRefOffset(23));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(24)", 27, cigar.getExpandedCigarIndexFromRefOffset(24));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(25)", 28, cigar.getExpandedCigarIndexFromRefOffset(25));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(26)", 29, cigar.getExpandedCigarIndexFromRefOffset(26));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(27)", -1, cigar.getExpandedCigarIndexFromRefOffset(27));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(28)", -1, cigar.getExpandedCigarIndexFromRefOffset(28));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(29)", -1, cigar.getExpandedCigarIndexFromRefOffset(29));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(30)", -1, cigar.getExpandedCigarIndexFromRefOffset(30));
    check(failures, ++testNum, "getExpandedCigarIndexFromRefOffset(31)", -1, cigar.getExpandedCigarIndexFromRefOffset(31));

    // 
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(0)", 0, cigar.getExpandedCigarIndexFromQueryIndex(0));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(1)", 1, cigar.getExpandedCigarIndexFromQueryIndex(1));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(2)", 2, cigar.getExpandedCigarIndexFromQueryIndex(2));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(3)", 3, cigar.getExpandedCigarIndexFromQueryIndex(3));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(4)", 14, cigar.getExpandedCigarIndexFromQueryIndex(4));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(5)", 15, cigar.getExpandedCigarIndexFromQueryIndex(5));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(6)", 16, cigar.getExpandedCigarIndexFromQueryIndex(6));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(7)", 17, cigar.getExpandedCigarIndexFromQueryIndex(7));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(8)", 18, cigar.getExpandedCigarIndexFromQueryIndex(8));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(9)", 19, cigar.getExpandedCigarIndexFromQueryIndex(9));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(10)", 20, cigar.getExpandedCigarIndexFromQueryIndex(10));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(11)", 21, cigar.getExpandedCigarIndexFromQueryIndex(11));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(12)", 22, cigar.getExpandedCigarIndexFromQueryIndex(12));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(13)", 27, cigar.getExpandedCigarIndexFromQueryIndex(13));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(14)", 28, cigar.getExpandedCigarIndexFromQueryIndex(14));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(15)", 29, cigar.getExpandedCigarIndexFromQueryIndex(15));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(16)", -1, cigar.getExpandedCigarIndexFromQueryIndex(16));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(17)", -1, cigar.getExpandedCigarIndexFromQueryIndex(17));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(18)", -1, cigar.getExpandedCigarIndexFromQueryIndex(18));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(19)", -1, cigar.getExpandedCigarIndexFromQueryIndex(19));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(20)", -1, cigar.getExpandedCigarIndexFromQueryIndex(20));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(21)", -1, cigar.getExpandedCigarIndexFromQueryIndex(21));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(22)", -1, cigar.getExpandedCigarIndexFromQueryIndex(22));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(23)", -1, cigar.getExpandedCigarIndexFromQueryIndex(23));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(24)", -1, cigar.getExpandedCigarIndexFromQueryIndex(24));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(25)", -1, cigar.getExpandedCigarIndexFromQueryIndex(25));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(26)", -1, cigar.getExpandedCigarIndexFromQueryIndex(26));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(27)", -1, cigar.getExpandedCigarIndexFromQueryIndex(27));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(28)", -1, cigar.getExpandedCigarIndexFromQueryIndex(28));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(29)", -1, cigar.getExpandedCigarIndexFromQueryIndex(29));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(30)", -1, cigar.getExpandedCigarIndexFromQueryIndex(30));
    check(failures, ++testNum, "getExpandedCigarIndexFromQueryIndex(31)", -1, cigar.getExpandedCigarIndexFromQueryIndex(31));

    // Test getCigarCharOp.
    check(failures, ++testNum, "getCigarCharOp(-1)", '?', cigar.getCigarCharOp(-1));
    check(failures, ++testNum, "getCigarCharOp(0)", 'M', cigar.getCigarCharOp(0));
    check(failures, ++testNum, "getCigarCharOp(1)", 'M', cigar.getCigarCharOp(1));
    check(failures, ++testNum, "getCigarCharOp(2)", 'M', cigar.getCigarCharOp(2));
    check(failures, ++testNum, "getCigarCharOp(3)", 'M', cigar.getCigarCharOp(3));
    check(failures, ++testNum, "getCigarCharOp(4)", 'N', cigar.getCigarCharOp(4));
    check(failures, ++testNum, "getCigarCharOp(5)", 'N', cigar.getCigarCharOp(5));
    check(failures, ++testNum, "getCigarCharOp(6)", 'N', cigar.getCigarCharOp(6));
    check(failures, ++testNum, "getCigarCharOp(7)", 'N', cigar.getCigarCharOp(7));
    check(failures, ++testNum, "getCigarCharOp(8)", 'N', cigar.getCigarCharOp(8));
    check(failures, ++testNum, "getCigarCharOp(9)", 'N', cigar.getCigarCharOp(9));
    check(failures, ++testNum, "getCigarCharOp(10)", 'N', cigar.getCigarCharOp(10));
    check(failures, ++testNum, "getCigarCharOp(11)", 'N', cigar.getCigarCharOp(11));
    check(failures, ++testNum, "getCigarCharOp(12)", 'N', cigar.getCigarCharOp(12));
    check(failures, ++testNum, "getCigarCharOp(13)", 'N', cigar.getCigarCharOp(13));
    check(failures, ++testNum, "getCigarCharOp(14)", 'M', cigar.getCigarCharOp(14));
    check(failures, ++testNum, "getCigarCharOp(15)", 'M', cigar.getCigarCharOp(15));
    check(failures, ++testNum, "getCigarCharOp(16)", 'M', cigar.getCigarCharOp(16));
    check(failures, ++testNum, "getCigarCharOp(17)", 'M', cigar.getCigarCharOp(17));
    check(failures, ++testNum, "getCigarCharOp(18)", 'I', cigar.getCigarCharOp(18));
    check(failures, ++testNum, "getCigarCharOp(19)", 'I', cigar.getCigarCharOp(19));
    check(failures, ++testNum, "getCigarCharOp(20)", 'I', cigar.getCigarCharOp(20));
    check(failures, ++testNum, "getCigarCharOp(21)", 'M', cigar.getCigarCharOp(21));
    check(failures, ++testNum, "getCigarCharOp(22)", 'M', cigar.getCigarCharOp(22));
    check(failures, ++testNum, "getCigarCharOp(23)", 'D', cigar.getCigarCharOp(23));
    check(failures, ++testNum, "getCigarCharOp(24)", 'D', cigar.getCigarCharOp(24));
    check(failures, ++testNum, "getCigarCharOp(25)", 'D', cigar.getCigarCharOp(25));
    check(failures, ++testNum, "getCigarCharOp(26)", 'D', cigar.getCigarCharOp(26));
    check(failures, ++testNum, "getCigarCharOp(27)", 'M', cigar.getCigarCharOp(27));
    check(failures, ++testNum, "getCigarCharOp(28)", 'M', cigar.getCigarCharOp(28));
    check(failures, ++testNum, "getCigarCharOp(29)", 'M', cigar.getCigarCharOp(29));
    check(failures, ++testNum, "getCigarCharOp(30)", '?', cigar.getCigarCharOp(30));


    // Test getCigarCharOpFromQueryIndex.
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(-1)", '?', cigar.getCigarCharOpFromQueryIndex(-1));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(0)", 'M', cigar.getCigarCharOpFromQueryIndex(0));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(1)", 'M', cigar.getCigarCharOpFromQueryIndex(1));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(2)", 'M', cigar.getCigarCharOpFromQueryIndex(2));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(3)", 'M', cigar.getCigarCharOpFromQueryIndex(3));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(4)", 'M', cigar.getCigarCharOpFromQueryIndex(4));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(5)", 'M', cigar.getCigarCharOpFromQueryIndex(5));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(6)", 'M', cigar.getCigarCharOpFromQueryIndex(6));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(7)", 'M', cigar.getCigarCharOpFromQueryIndex(7));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(8)", 'I', cigar.getCigarCharOpFromQueryIndex(8));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(9)", 'I', cigar.getCigarCharOpFromQueryIndex(9));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(10)", 'I', cigar.getCigarCharOpFromQueryIndex(10));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(11)", 'M', cigar.getCigarCharOpFromQueryIndex(11));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(12)", 'M', cigar.getCigarCharOpFromQueryIndex(12));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(13)", 'M', cigar.getCigarCharOpFromQueryIndex(13));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(14)", 'M', cigar.getCigarCharOpFromQueryIndex(14));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(15)", 'M', cigar.getCigarCharOpFromQueryIndex(15));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(16)", '?', cigar.getCigarCharOpFromQueryIndex(16));
    check(failures, ++testNum, "getCigarCharOpFromQueryIndex(17)", '?', cigar.getCigarCharOpFromQueryIndex(17));

    // Test getCigarCharOpFromRefOffset.
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(-1)", '?', cigar.getCigarCharOpFromRefOffset(-1));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(0)", 'M', cigar.getCigarCharOpFromRefOffset(0));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(1)", 'M', cigar.getCigarCharOpFromRefOffset(1));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(2)", 'M', cigar.getCigarCharOpFromRefOffset(2));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(3)", 'M', cigar.getCigarCharOpFromRefOffset(3));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(4)", 'N', cigar.getCigarCharOpFromRefOffset(4));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(5)", 'N', cigar.getCigarCharOpFromRefOffset(5));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(6)", 'N', cigar.getCigarCharOpFromRefOffset(6));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(7)", 'N', cigar.getCigarCharOpFromRefOffset(7));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(8)", 'N', cigar.getCigarCharOpFromRefOffset(8));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(9)", 'N', cigar.getCigarCharOpFromRefOffset(9));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(10)", 'N', cigar.getCigarCharOpFromRefOffset(10));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(11)", 'N', cigar.getCigarCharOpFromRefOffset(11));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(12)", 'N', cigar.getCigarCharOpFromRefOffset(12));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(13)", 'N', cigar.getCigarCharOpFromRefOffset(13));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(14)", 'M', cigar.getCigarCharOpFromRefOffset(14));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(15)", 'M', cigar.getCigarCharOpFromRefOffset(15));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(16)", 'M', cigar.getCigarCharOpFromRefOffset(16));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(17)", 'M', cigar.getCigarCharOpFromRefOffset(17));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(18)", 'M', cigar.getCigarCharOpFromRefOffset(18));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(19)", 'M', cigar.getCigarCharOpFromRefOffset(19));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(20)", 'D', cigar.getCigarCharOpFromRefOffset(20));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(21)", 'D', cigar.getCigarCharOpFromRefOffset(21));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(22)", 'D', cigar.getCigarCharOpFromRefOffset(22));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(23)", 'D', cigar.getCigarCharOpFromRefOffset(23));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(24)", 'M', cigar.getCigarCharOpFromRefOffset(24));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(25)", 'M', cigar.getCigarCharOpFromRefOffset(25));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(26)", 'M', cigar.getCigarCharOpFromRefOffset(26));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(27)", '?', cigar.getCigarCharOpFromRefOffset(27));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(28)", '?', cigar.getCigarCharOpFromRefOffset(28));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(29)", '?', cigar.getCigarCharOpFromRefOffset(29));
    check(failures, ++testNum, "getCigarCharOpFromRefOffset(30)", '?', cigar.getCigarCharOpFromRefOffset(30));


    // Test getCigarCharOpFromRefPos.
    check(failures, ++testNum, "getCigarCharOpFromRefPos(-1, 5)", '?', cigar.getCigarCharOpFromRefPos(-1,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(0, 5)", '?', cigar.getCigarCharOpFromRefPos(0,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(1, 5)", '?', cigar.getCigarCharOpFromRefPos(1,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(2, 5)", '?', cigar.getCigarCharOpFromRefPos(2,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(3, 5)", '?', cigar.getCigarCharOpFromRefPos(3,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(4, 5)", '?', cigar.getCigarCharOpFromRefPos(4,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(5, 5)", 'M', cigar.getCigarCharOpFromRefPos(5,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(6, 5)", 'M', cigar.getCigarCharOpFromRefPos(6,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(7, 5)", 'M', cigar.getCigarCharOpFromRefPos(7,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(8, 5)", 'M', cigar.getCigarCharOpFromRefPos(8,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(9, 5)", 'N', cigar.getCigarCharOpFromRefPos(9,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(10, 5)", 'N', cigar.getCigarCharOpFromRefPos(10,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(11, 5)", 'N', cigar.getCigarCharOpFromRefPos(11,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(12, 5)", 'N', cigar.getCigarCharOpFromRefPos(12,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(13, 5)", 'N', cigar.getCigarCharOpFromRefPos(13,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(14, 5)", 'N', cigar.getCigarCharOpFromRefPos(14,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(15, 5)", 'N', cigar.getCigarCharOpFromRefPos(15,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(16, 5)", 'N', cigar.getCigarCharOpFromRefPos(16,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(17, 5)", 'N', cigar.getCigarCharOpFromRefPos(17,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(18, 5)", 'N', cigar.getCigarCharOpFromRefPos(18,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(19, 5)", 'M', cigar.getCigarCharOpFromRefPos(19,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(20, 5)", 'M', cigar.getCigarCharOpFromRefPos(20,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(21, 5)", 'M', cigar.getCigarCharOpFromRefPos(21,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(22, 5)", 'M', cigar.getCigarCharOpFromRefPos(22,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(23, 5)", 'M', cigar.getCigarCharOpFromRefPos(23,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(24, 5)", 'M', cigar.getCigarCharOpFromRefPos(24,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(25, 5)", 'D', cigar.getCigarCharOpFromRefPos(25,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(26, 5)", 'D', cigar.getCigarCharOpFromRefPos(26,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(27, 5)", 'D', cigar.getCigarCharOpFromRefPos(27,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(28, 5)", 'D', cigar.getCigarCharOpFromRefPos(28,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(29, 5)", 'M', cigar.getCigarCharOpFromRefPos(29,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(30, 5)", 'M', cigar.getCigarCharOpFromRefPos(30,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(31, 5)", 'M', cigar.getCigarCharOpFromRefPos(31,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(32, 5)", '?', cigar.getCigarCharOpFromRefPos(32,5));
    check(failures, ++testNum, "getCigarCharOpFromRefPos(33, 5)", '?', cigar.getCigarCharOpFromRefPos(33,5));




    ////////////////////////
    // Test getNumOverlaps.
    check(failures, ++testNum, "getNumOverlaps(5,32,5)", (uint32_t)13, cigar.getNumOverlaps(5,32,5));
    check(failures, ++testNum, "getNumOverlaps(5,31,5)", (uint32_t)12, cigar.getNumOverlaps(5,31,5));
    check(failures, ++testNum, "getNumOverlaps(0,100,5)", (uint32_t)13, cigar.getNumOverlaps(0,100,5));
    check(failures, ++testNum, "getNumOverlaps(-1, -1,5)", (uint32_t)13, cigar.getNumOverlaps(-1, -1,5));
    check(failures, ++testNum, "getNumOverlaps(-1,10,5)", (uint32_t)4, cigar.getNumOverlaps(-1,10,5));
    check(failures, ++testNum, "getNumOverlaps(10,-1,5)", (uint32_t)9, cigar.getNumOverlaps(10,-1,5));
    check(failures, ++testNum, "getNumOverlaps(9,19,5)", (uint32_t)0, cigar.getNumOverlaps(9,19,5));
    check(failures, ++testNum, "getNumOverlaps(9,20,5)", (uint32_t)1, cigar.getNumOverlaps(9,20,5));
    check(failures, ++testNum, "getNumOverlaps(9,6,5)", (uint32_t)0, cigar.getNumOverlaps(9,6,5));
    check(failures, ++testNum, "getNumOverlaps(0,5,5)", (uint32_t)0, cigar.getNumOverlaps(0,5,5));
    check(failures, ++testNum, "getNumOverlaps(32,40,5)", (uint32_t)0, cigar.getNumOverlaps(32,40,5));
    check(failures, ++testNum, "getNumOverlaps(0,5,1)", (uint32_t)4, cigar.getNumOverlaps(0,5,1));
    check(failures, ++testNum, "getNumOverlaps(32,40,32)", (uint32_t)4, cigar.getNumOverlaps(32,40,32));

    // Get Query Index for reference offset 0 - 27
    // 4M
    check(failures, ++testNum, "getQueryIndex(0)", 0, cigar.getQueryIndex(0));
    check(failures, ++testNum, "getQueryIndex(1)", 1, cigar.getQueryIndex(1));
    check(failures, ++testNum, "getQueryIndex(2)", 2, cigar.getQueryIndex(2));
    check(failures, ++testNum, "getQueryIndex(3)", 3, cigar.getQueryIndex(3));
    // 10N
    check(failures, ++testNum, "getQueryIndex(4)", -1, cigar.getQueryIndex(4));
    check(failures, ++testNum, "getQueryIndex(5)", -1, cigar.getQueryIndex(5));
    check(failures, ++testNum, "getQueryIndex(6)", -1, cigar.getQueryIndex(6));
    check(failures, ++testNum, "getQueryIndex(7)", -1, cigar.getQueryIndex(7));
    check(failures, ++testNum, "getQueryIndex(8)", -1, cigar.getQueryIndex(8));
    check(failures, ++testNum, "getQueryIndex(9)", -1, cigar.getQueryIndex(9));
    check(failures, ++testNum, "getQueryIndex(10)", -1, cigar.getQueryIndex(10));
    check(failures, ++testNum, "getQueryIndex(11)", -1, cigar.getQueryIndex(11));
    check(failures, ++testNum, "getQueryIndex(12)", -1, cigar.getQueryIndex(12));
    check(failures, ++testNum, "getQueryIndex(13)", -1, cigar.getQueryIndex(13));
    // 4M
    check(failures, ++testNum, "getQueryIndex(14)", 4, cigar.getQueryIndex(14));
    check(failures, ++testNum, "getQueryIndex(15)", 5, cigar.getQueryIndex(15));
    check(failures, ++testNum, "getQueryIndex(16)", 6, cigar.getQueryIndex(16));
    check(failures, ++testNum, "getQueryIndex(17)", 7, cigar.getQueryIndex(17));
    // 3I - nothing to check - not in reference - covers query indices  8-10
    // 2M
    check(failures, ++testNum, "getQueryIndex(18)", 11, cigar.getQueryIndex(18));
    check(failures, ++testNum, "getQueryIndex(19)", 12, cigar.getQueryIndex(19));
    // 4D
    check(failures, ++testNum, "getQueryIndex(20)", -1, cigar.getQueryIndex(20));
    check(failures, ++testNum, "getQueryIndex(21)", -1, cigar.getQueryIndex(21));
    check(failures, ++testNum, "getQueryIndex(22)", -1, cigar.getQueryIndex(22));
    check(failures, ++testNum, "getQueryIndex(23)", -1, cigar.getQueryIndex(23));
    // 3M
    check(failures, ++testNum, "getQueryIndex(24)", 13, cigar.getQueryIndex(24));
    check(failures, ++testNum, "getQueryIndex(25)", 14, cigar.getQueryIndex(25));
    check(failures, ++testNum, "getQueryIndex(26)", 15, cigar.getQueryIndex(26));

    //  Get Query Index for reference positions 0-33
    // N/A
    check(failures, ++testNum, "getQueryIndex(0, 5)", -1, cigar.getQueryIndex(0, 5));
    check(failures, ++testNum, "getQueryIndex(1, 5)", -1, cigar.getQueryIndex(1, 5));
    check(failures, ++testNum, "getQueryIndex(2, 5)", -1, cigar.getQueryIndex(2, 5));
    check(failures, ++testNum, "getQueryIndex(3, 5)", -1, cigar.getQueryIndex(3, 5));
    check(failures, ++testNum, "getQueryIndex(4, 5)", -1, cigar.getQueryIndex(4, 5));
    // 4M
    check(failures, ++testNum, "getQueryIndex(5, 5)", 0, cigar.getQueryIndex(5, 5));
    check(failures, ++testNum, "getQueryIndex(6, 5)", 1, cigar.getQueryIndex(6, 5));
    check(failures, ++testNum, "getQueryIndex(7, 5)", 2, cigar.getQueryIndex(7, 5));
    check(failures, ++testNum, "getQueryIndex(8, 5)", 3, cigar.getQueryIndex(8, 5));
    // 10N
    check(failures, ++testNum, "getQueryIndex(9, 5)", -1, cigar.getQueryIndex(9, 5));
    check(failures, ++testNum, "getQueryIndex(10, 5)", -1, cigar.getQueryIndex(10, 5));
    check(failures, ++testNum, "getQueryIndex(11, 5)", -1, cigar.getQueryIndex(11, 5));
    check(failures, ++testNum, "getQueryIndex(12, 5)", -1, cigar.getQueryIndex(12, 5));
    check(failures, ++testNum, "getQueryIndex(13, 5)", -1, cigar.getQueryIndex(13, 5));
    check(failures, ++testNum, "getQueryIndex(14, 5)", -1, cigar.getQueryIndex(14, 5));
    check(failures, ++testNum, "getQueryIndex(15, 5)", -1, cigar.getQueryIndex(15, 5));
    check(failures, ++testNum, "getQueryIndex(16, 5)", -1, cigar.getQueryIndex(16, 5));
    check(failures, ++testNum, "getQueryIndex(17, 5)", -1, cigar.getQueryIndex(17, 5));
    check(failures, ++testNum, "getQueryIndex(18, 5)", -1, cigar.getQueryIndex(18, 5));
    // 4M
    check(failures, ++testNum, "getQueryIndex(19, 5)", 4, cigar.getQueryIndex(19, 5));
    check(failures, ++testNum, "getQueryIndex(20, 5)", 5, cigar.getQueryIndex(20, 5));
    check(failures, ++testNum, "getQueryIndex(21, 5)", 6, cigar.getQueryIndex(21, 5));
    check(failures, ++testNum, "getQueryIndex(22, 5)", 7, cigar.getQueryIndex(22, 5));
    // 3I - nothing to check - not in reference - covers query indices  8-10
    // 2M
    check(failures, ++testNum, "getQueryIndex(23, 5)", 11, cigar.getQueryIndex(23, 5));
    check(failures, ++testNum, "getQueryIndex(24, 5)", 12, cigar.getQueryIndex(24, 5));
    // 4D
    check(failures, ++testNum, "getQueryIndex(25, 5)", -1, cigar.getQueryIndex(25, 5));
    check(failures, ++testNum, "getQueryIndex(26, 5)", -1, cigar.getQueryIndex(26, 5));
    check(failures, ++testNum, "getQueryIndex(27, 5)", -1, cigar.getQueryIndex(27, 5));
    check(failures, ++testNum, "getQueryIndex(28, 5)", -1, cigar.getQueryIndex(28, 5));
    // 3M
    check(failures, ++testNum, "getQueryIndex(29, 5)", 13, cigar.getQueryIndex(29, 5));
    check(failures, ++testNum, "getQueryIndex(30, 5)", 14, cigar.getQueryIndex(30, 5));
    check(failures, ++testNum, "getQueryIndex(31, 5)", 15, cigar.getQueryIndex(31, 5));

    // Get reference offset for query index 0 - 17
    // 4M
    check(failures, ++testNum, "getRefOffset(0)", 0, cigar.getRefOffset(0));
    check(failures, ++testNum, "getRefOffset(1)", 1, cigar.getRefOffset(1));
    check(failures, ++testNum, "getRefOffset(2)", 2, cigar.getRefOffset(2));
    check(failures, ++testNum, "getRefOffset(3)", 3, cigar.getRefOffset(3));
    // 10N - nothing to check - not in query - covers ref offsets 4-13
    // 4M
    check(failures, ++testNum, "getRefOffset(4)", 14, cigar.getRefOffset(4));
    check(failures, ++testNum, "getRefOffset(5)", 15, cigar.getRefOffset(5));
    check(failures, ++testNum, "getRefOffset(6)", 16, cigar.getRefOffset(6));
    check(failures, ++testNum, "getRefOffset(7)", 17, cigar.getRefOffset(7));
    // 3I
    check(failures, ++testNum, "getRefOffset(8)", -1, cigar.getRefOffset(8));
    check(failures, ++testNum, "getRefOffset(9)", -1, cigar.getRefOffset(9));
    check(failures, ++testNum, "getRefOffset(10)", -1, cigar.getRefOffset(10));
    // 2M
    check(failures, ++testNum, "getRefOffset(11)", 18, cigar.getRefOffset(11));
    check(failures, ++testNum, "getRefOffset(12)", 19, cigar.getRefOffset(12));
    // 4D - nothing to check - not in query - covers ref offsets 20-23
    // 3M
    check(failures, ++testNum, "getRefOffset(13)", 24, cigar.getRefOffset(13));
    check(failures, ++testNum, "getRefOffset(14)", 25, cigar.getRefOffset(14));
    check(failures, ++testNum, "getRefOffset(15)", 26, cigar.getRefOffset(15));


    // Get reference position for query index 0 - 17
    // 4M
    check(failures, ++testNum, "getRefPosition(0, 5)", 5, cigar.getRefPosition(0, 5));
    check(failures, ++testNum, "getRefPosition(1, 5)", 6, cigar.getRefPosition(1, 5));
    check(failures, ++testNum, "getRefPosition(2, 5)", 7, cigar.getRefPosition(2, 5));
    check(failures, ++testNum, "getRefPosition(3, 5)", 8, cigar.getRefPosition(3, 5));
    // 10N - nothing to check - not in query - covers ref offsets 4-13
    // 4M
    check(failures, ++testNum, "getRefPosition(4, 5)", 19, cigar.getRefPosition(4, 5));
    check(failures, ++testNum, "getRefPosition(5, 5)", 20, cigar.getRefPosition(5, 5));
    check(failures, ++testNum, "getRefPosition(6, 5)", 21, cigar.getRefPosition(6, 5));
    check(failures, ++testNum, "getRefPosition(7, 5)", 22, cigar.getRefPosition(7, 5));
    // 3I
    check(failures, ++testNum, "getRefPosition(8, 5)", -1, cigar.getRefPosition(8, 5));
    check(failures, ++testNum, "getRefPosition(9, 5)", -1, cigar.getRefPosition(9, 5));
    check(failures, ++testNum, "getRefPosition(10, 5)", -1, cigar.getRefPosition(10, 5));
    // 2M
    check(failures, ++testNum, "getRefPosition(11, 5)", 23, cigar.getRefPosition(11, 5));
    check(failures, ++testNum, "getRefPosition(12, 5)", 24, cigar.getRefPosition(12, 5));
    // 4D - nothing to check - not in query - covers ref pos 25-28
    // 3M
    check(failures, ++testNum, "getRefPosition(13, 5)", 29, cigar.getRefPosition(13, 5));
    check(failures, ++testNum, "getRefPosition(14, 5)", 30, cigar.getRefPosition(14, 5));
    check(failures, ++testNum, "getRefPosition(15, 5)", 31, cigar.getRefPosition(15, 5));



    ////////////////////////////////////////////////////////////////////////////
    // Test a new CIGAR set by buffer.
    // 2S 3M 1I 2M 1D 1M 2P 1M 3N 1M 3H
    uint32_t cigarBuf[] = {0x24,  // 2S = 2 << 4 | 4
                           0x30,  // 3M = 3 << 4 | 0
                           0x11,  // 1I = 1 << 4 | 1
                           0x20,  // 2M = 2 << 4 | 0
                           0x12,  // 1D = 1 << 4 | 2
                           0x10,  // 1M = 1 << 4 | 0
                           0x26,  // 2P = 2 << 4 | 6
                           0x10,  // 1m = 1 << 4 | 0
                           0x33,  // 3N = 3 << 4 | 3
                           0x10,  // 1M = 1 << 4 | 0
                           0x35}; // 3H = 3 << 4 | 5
    cigar.Set(cigarBuf, 11);
    cigarString = "HI";
    cigar.getCigarString(cigarString);
    expectedResult = "2S3M1I2M1D1M2P1M3N1M3H";
    check(failures, ++testNum, "getCigarString", expectedResult, cigarString);
    check(failures, ++testNum, "getNumEndClips", 3, cigar.getNumEndClips());
    check(failures, ++testNum, "getNumBeginClips", 2, cigar.getNumBeginClips());

    std::cout << "\nCigarRoller PASS: " << testNum - failures << "  FAIL: " << failures << std::endl;
    // return the number of failures.
    return(failures);
}


int main(int argc, const char **argv)
{
    CigarRollerTest roller;
    
    bool showAllCasesFlag = false;
    int opt;

    while(( opt = getopt(argc, (char **) argv, "v")) != -1) {
        switch(opt) {
            case 'v':
                showAllCasesFlag = true;
                break;
            default:
                std::cerr << "usage: testSW [-v]" << std::endl;
                exit(1);
        }
    }
    if(showAllCasesFlag)
    {
    }

    //
    // do cigar roller tests first
    //
    return(roller.test());

    // CIGAR explanation - for backward SW runs, the corresponding
    // CIGAR string is generated from the back of the string to the
    // front.  Recall that the soft clipping is only done at the
    // "end" of the string, taking direction into account.

    // Comment out this result since it doesn't refelct the results of test.
    //    cout << endl << "Total Errors found: " << errors << endl;
}
