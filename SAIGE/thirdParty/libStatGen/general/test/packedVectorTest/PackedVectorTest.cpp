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

#include <getopt.h>
#include "Generic.h"
#include <stdio.h>
#include "PackedVector.h"
#include "UnitTest.h"

#include <assert.h>
#include <stdlib.h>

class PackedArrayTest : public UnitTest
{
public:
    PackedArrayTest(const char *title) : UnitTest(title) {;}
    void testBool();
    void test2Bit();
    void test4Bit();
    void testResize();

    void test() {
        testBool();
        test2Bit();
        test4Bit();
        testResize();
    }
};

void PackedArrayTest::testBool(void)
{
    PackedVectorBool_t   testVector;

    testVector.resize(11);
    testVector.set(0,0);
    testVector.set(1,1);
    testVector.set(2,0);
    testVector.set(3,1);
    testVector.set(4,1);
    testVector.set(5,0);
    testVector.set(6,1);
    testVector.set(7,0);
    testVector.set(8,0);
    testVector.set(9,0);
    testVector.set(10,1);
    check(m_failures, ++m_testNum, "Access 1 bit element 0", 0U, testVector[0]);
    check(m_failures, ++m_testNum, "Access 1 bit element 1", 1U, testVector[1]);
    check(m_failures, ++m_testNum, "Access 1 bit element 2", 0U, testVector[2]);
    check(m_failures, ++m_testNum, "Access 1 bit element 3", 1U, testVector[3]);
    check(m_failures, ++m_testNum, "Access 1 bit element 4", 1U, testVector[4]);
    check(m_failures, ++m_testNum, "Access 1 bit element 5", 0U, testVector[5]);
    check(m_failures, ++m_testNum, "Access 1 bit element 6", 1U, testVector[6]);
    check(m_failures, ++m_testNum, "Access 1 bit element 7", 0U, testVector[7]);
    check(m_failures, ++m_testNum, "Access 1 bit element 8", 0U, testVector[8]);
    check(m_failures, ++m_testNum, "Access 1 bit element 9", 0U, testVector[9]);
    check(m_failures, ++m_testNum, "Access 1 bit element 10", 1U, testVector[10]);
}

void PackedArrayTest::test2Bit(void)
{
    PackedVector2Bit_t   testVector;

    testVector.resize(11);
    testVector.set(0,0);
    testVector.set(1,1);
    testVector.set(2,2);
    testVector.set(3,3);
    testVector.set(4,3);
    testVector.set(5,2);
    testVector.set(6,1);
    testVector.set(7,0);
    testVector.set(8,2);
    testVector.set(9,1);
    testVector.set(10,3);
    check(m_failures, ++m_testNum, "Access 2 bit element 0", 0U, testVector[0]);
    check(m_failures, ++m_testNum, "Access 2 bit element 1", 1U, testVector[1]);
    check(m_failures, ++m_testNum, "Access 2 bit element 2", 2U, testVector[2]);
    check(m_failures, ++m_testNum, "Access 2 bit element 3", 3U, testVector[3]);
    check(m_failures, ++m_testNum, "Access 2 bit element 4", 3U, testVector[4]);
    check(m_failures, ++m_testNum, "Access 2 bit element 5", 2U, testVector[5]);
    check(m_failures, ++m_testNum, "Access 2 bit element 6", 1U, testVector[6]);
    check(m_failures, ++m_testNum, "Access 2 bit element 7", 0U, testVector[7]);
    check(m_failures, ++m_testNum, "Access 2 bit element 8", 2U, testVector[8]);
    check(m_failures, ++m_testNum, "Access 2 bit element 9", 1U, testVector[9]);
    check(m_failures, ++m_testNum, "Access 2 bit element 10", 3U, testVector[10]);
}

void PackedArrayTest::test4Bit(void)
{
    PackedVector4Bit_t   testVector;

    testVector.resize(11);
    testVector.set(0,0);
    testVector.set(1,1);
    testVector.set(2,2);
    testVector.set(3,3);
    testVector.set(4,4);
    testVector.set(5,5);
    testVector.set(6,6);
    testVector.set(7,7);
    testVector.set(8,8);
    testVector.set(9,9);
    testVector.set(10,10);

    check(m_failures, ++m_testNum, "Access 4 bit element 0", 0U, testVector[0]);
    check(m_failures, ++m_testNum, "Access 4 bit element 1", 1U, testVector[1]);
    check(m_failures, ++m_testNum, "Access 4 bit element 2", 2U, testVector[2]);
    check(m_failures, ++m_testNum, "Access 4 bit element 3", 3U, testVector[3]);
    check(m_failures, ++m_testNum, "Access 4 bit element 4", 4U, testVector[4]);
    check(m_failures, ++m_testNum, "Access 4 bit element 5", 5U, testVector[5]);
    check(m_failures, ++m_testNum, "Access 4 bit element 6", 6U, testVector[6]);
    check(m_failures, ++m_testNum, "Access 4 bit element 7", 7U, testVector[7]);
    check(m_failures, ++m_testNum, "Access 4 bit element 8", 8U, testVector[8]);
    check(m_failures, ++m_testNum, "Access 4 bit element 9", 9U, testVector[9]);
    check(m_failures, ++m_testNum, "Access 4 bit element 10", 10U, testVector[10]);
}

void PackedArrayTest::testResize(void)
{
    PackedVector4Bit_t   testVector;

    testVector.resize(0);
    check(m_failures, ++m_testNum, "New size is 0", 0U, testVector.size());

    testVector.push_back(0);
    testVector.push_back(1);
    testVector.push_back(2);
    testVector.push_back(3);
    testVector.push_back(4);
    testVector.push_back(5);
    testVector.push_back(6);
    testVector.push_back(7);
    testVector.push_back(8);
    testVector.push_back(9);
    testVector.push_back(10);

    check(m_failures, ++m_testNum, "New size is 11", 11U, testVector.size());

    check(m_failures, ++m_testNum, "Access 4 bit element 0", 0U, testVector[0]);
    check(m_failures, ++m_testNum, "Access 4 bit element 1", 1U, testVector[1]);
    check(m_failures, ++m_testNum, "Access 4 bit element 2", 2U, testVector[2]);
    check(m_failures, ++m_testNum, "Access 4 bit element 3", 3U, testVector[3]);
    check(m_failures, ++m_testNum, "Access 4 bit element 4", 4U, testVector[4]);
    check(m_failures, ++m_testNum, "Access 4 bit element 5", 5U, testVector[5]);
    check(m_failures, ++m_testNum, "Access 4 bit element 6", 6U, testVector[6]);
    check(m_failures, ++m_testNum, "Access 4 bit element 7", 7U, testVector[7]);
    check(m_failures, ++m_testNum, "Access 4 bit element 8", 8U, testVector[8]);
    check(m_failures, ++m_testNum, "Access 4 bit element 9", 9U, testVector[9]);
    check(m_failures, ++m_testNum, "Access 4 bit element 10", 10U, testVector[10]);
}

int main(int argc, char **argv)
{
    PackedArrayTest test("PackedArrayTest");
#if 0
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

#endif

    test.test();

    std::cout << test;

    exit(test.getFailureCount());
}
