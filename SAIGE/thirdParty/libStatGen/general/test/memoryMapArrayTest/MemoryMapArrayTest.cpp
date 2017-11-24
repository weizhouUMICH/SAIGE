/*
 *  Copyright (C) 2010-2012  Regents of the University of Michigan
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
#include "MemoryMapArray.h"
#include "MemoryMapArrayTest.h"

#include <assert.h>
#include <stdlib.h>

#define TEST_FILE_NAME "results/testMemoryMapArray.vector"

class MemoryMapArrayTest : public UnitTest
{
public:
    MemoryMapArrayTest(const char *title) : UnitTest(title) {;}
    void testBool();
    void test2Bit();
    void test4Bit();
    void test32Bit();

    void test() {
        testBool();
        test2Bit();
        test4Bit();
        test32Bit();
    }
};

void MemoryMapArrayTest::testBool(void)
{
    mmapArrayBool_t   testVector;

    // ignore return code here
    if(unlink(TEST_FILE_NAME) == 0)
    {
        // Nothing to do, just deleting previous test file
    }
    

    check(m_failures, ++m_testNum, "Create 1 bit vector file", 0,
            testVector.create(TEST_FILE_NAME, 11));
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
    check(m_failures, ++m_testNum, "Close vector file", false, testVector.close());
    check(m_failures, ++m_testNum, "Re-open vector file", false, testVector.open(TEST_FILE_NAME));
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

    check(m_failures, ++m_testNum, "Close vector file", false, testVector.close());
    check(m_failures, ++m_testNum, "Unlink vector file", 0, unlink(TEST_FILE_NAME));


}

void MemoryMapArrayTest::test2Bit(void)
{
    mmapArray2Bit_t   testVector;

    // ignore return code here
    if(unlink(TEST_FILE_NAME) == 0)
    {
        // Nothing to do, just deleting previous test file
    }

    check(m_failures, ++m_testNum, "Create 2 bit vector file", 0,
            testVector.create(TEST_FILE_NAME, 11));

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
    check(m_failures, ++m_testNum, "Close vector file", false, testVector.close());
    check(m_failures, ++m_testNum, "Re-open vector file", false, testVector.open(TEST_FILE_NAME));
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

    check(m_failures, ++m_testNum, "Close vector file", false, testVector.close());
    check(m_failures, ++m_testNum, "Unlink vector file", 0, unlink(TEST_FILE_NAME));

}

void MemoryMapArrayTest::test4Bit(void)
{
    mmapArray4Bit_t   testVector;

    // ignore return code here
    if(unlink(TEST_FILE_NAME) == 0)
    {
        // Nothing to do, just deleting previous test file
    }

    check(m_failures, ++m_testNum, "Create 4 bit vector file", 0,
            testVector.create(TEST_FILE_NAME, 11));

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

    check(m_failures, ++m_testNum, "Close vector file", false, testVector.close());
    check(m_failures, ++m_testNum, "Re-open vector file", false, testVector.open(TEST_FILE_NAME));
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

    check(m_failures, ++m_testNum, "Close vector file", false, testVector.close());
    check(m_failures, ++m_testNum, "Unlink vector file", 0, unlink(TEST_FILE_NAME));
}

void MemoryMapArrayTest::test32Bit(void)
{

    mmapArrayUint32_t   testVector;

    // ignore return code here
    if(unlink(TEST_FILE_NAME) == 0)
    {
        // Nothing to do, just deleting previous test file
    }

    check(m_failures, ++m_testNum, "Create 32 bit vector file", 0,
            testVector.create(TEST_FILE_NAME, 11));

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
    check(m_failures, ++m_testNum, "Access 32 bit element 0", 0U, testVector[0]);
    check(m_failures, ++m_testNum, "Access 32 bit element 1", 1U, testVector[1]);
    check(m_failures, ++m_testNum, "Access 32 bit element 2", 2U, testVector[2]);
    check(m_failures, ++m_testNum, "Access 32 bit element 3", 3U, testVector[3]);
    check(m_failures, ++m_testNum, "Access 32 bit element 4", 4U, testVector[4]);
    check(m_failures, ++m_testNum, "Access 32 bit element 5", 5U, testVector[5]);
    check(m_failures, ++m_testNum, "Access 32 bit element 6", 6U, testVector[6]);
    check(m_failures, ++m_testNum, "Access 32 bit element 7", 7U, testVector[7]);
    check(m_failures, ++m_testNum, "Access 32 bit element 8", 8U, testVector[8]);
    check(m_failures, ++m_testNum, "Access 32 bit element 9", 9U, testVector[9]);
    check(m_failures, ++m_testNum, "Access 32 bit element 10", 10U, testVector[10]);


    check(m_failures, ++m_testNum, "Close vector file", false, testVector.close());
    check(m_failures, ++m_testNum, "Re-open vector file", false, testVector.open(TEST_FILE_NAME));
    check(m_failures, ++m_testNum, "Access 32 bit element 0", 0U, testVector[0]);
    check(m_failures, ++m_testNum, "Access 32 bit element 1", 1U, testVector[1]);
    check(m_failures, ++m_testNum, "Access 32 bit element 2", 2U, testVector[2]);
    check(m_failures, ++m_testNum, "Access 32 bit element 3", 3U, testVector[3]);
    check(m_failures, ++m_testNum, "Access 32 bit element 4", 4U, testVector[4]);
    check(m_failures, ++m_testNum, "Access 32 bit element 5", 5U, testVector[5]);
    check(m_failures, ++m_testNum, "Access 32 bit element 6", 6U, testVector[6]);
    check(m_failures, ++m_testNum, "Access 32 bit element 7", 7U, testVector[7]);
    check(m_failures, ++m_testNum, "Access 32 bit element 8", 8U, testVector[8]);
    check(m_failures, ++m_testNum, "Access 32 bit element 9", 9U, testVector[9]);
    check(m_failures, ++m_testNum, "Access 32 bit element 10", 10U, testVector[10]);

    check(m_failures, ++m_testNum, "Close vector file", false, testVector.close());
    check(m_failures, ++m_testNum, "Unlink vector file", 0, unlink(TEST_FILE_NAME));
}

int main(int argc, char **argv)
{
    MemoryMapArrayTest test("MemoryMapArrayTest");
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
