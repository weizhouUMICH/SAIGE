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
#include "ReusableVector.h"
#include "ReusableVectorTest.h"
#include <assert.h>
#include <iostream>
#include <string.h>

int ReusableVectorTestDataType::ourValue = 0;
int ReusableVectorTestDataType::ourNumDestructs = 0;

int main(int argc, char ** argv)
{

   ReusableVectorTest myTest;

   myTest.test();
}

void ReusableVectorTest::test()
{
    assert(ReusableVectorTestDataType::ourNumDestructs == 0);
    testReuse();
    assert(ReusableVectorTestDataType::ourNumDestructs == 8);
}


void ReusableVectorTest::testReuse()
{
    ReusableVector<ReusableVectorTestDataType> testVector;
    ReusableVector<ReusableVectorTestDataType> testVector2;
    ReusableVectorTestDataType* dataPtr = NULL;

    assert(testVector.size() == 0);
    assert(testInvalidGetIndex(testVector, 0));
    assert(testInvalidGetIndex(testVector, 1));
    testVector.reset();
    assert(testVector.size() == 0);
    assert(testInvalidGetIndex(testVector, 0));
    assert(testInvalidGetIndex(testVector, 1));

    // Get three data pointers and check they are each new.
    dataPtr = &(testVector.getNextEmpty());
    assert(dataPtr->myValue == 0);
    assert(dataPtr->ourValue == 1);
    dataPtr = &(testVector.getNextEmpty());
    assert(dataPtr->myValue == 1);
    assert(dataPtr->ourValue == 2);
    dataPtr = &(testVector.getNextEmpty());
    assert(dataPtr->myValue == 2);
    assert(dataPtr->ourValue == 3);
    assert(testVector.size() == 3);

    // Check a 2nd test vector.
    assert(testVector2.size() == 0);
    assert(testInvalidGetIndex(testVector2, 0));
    assert(testInvalidGetIndex(testVector2, 1));
    testVector2.reset();
    assert(testVector2.size() == 0);
    assert(testInvalidGetIndex(testVector2, 0));
    assert(testInvalidGetIndex(testVector2, 1));

    // Get data pointers and check they are each new.
    dataPtr = &(testVector2.getNextEmpty());
    assert(dataPtr->myValue == 3);
    assert(dataPtr->ourValue == 4);
    dataPtr = &(testVector2.getNextEmpty());
    assert(dataPtr->myValue == 4);
    assert(dataPtr->ourValue == 5);
    assert(testVector2.size() == 2);
    
    // Test the get accessor.
    assert(testVector2.get(1).myValue == 4);
    assert(testVector2.get(0).myValue == 3);
    assert(testInvalidGetIndex(testVector2, 2));
   // Test the get accessor with the first vector.
    assert(testVector.get(1).myValue == 1);
    assert(testVector.get(0).myValue == 0);
    assert(testVector.get(2).myValue == 2);
    assert(testInvalidGetIndex(testVector, 3));

    // Clear the 1st vector.
    testVector.clear();
    assert(testVector.size() == 0);
    assert(testInvalidGetIndex(testVector, 0));
    assert(testInvalidGetIndex(testVector, 1));

    // Check the data values are reused.
    dataPtr = &(testVector.getNextEmpty());
    assert(dataPtr->myValue == 0);
    assert(dataPtr->ourValue == 5);
    assert(testVector.size() == 1);
    dataPtr = &(testVector.getNextEmpty());
    assert(dataPtr->myValue == 1);
    assert(dataPtr->ourValue == 5);
    assert(testVector.size() == 2);
    dataPtr = &(testVector.getNextEmpty());
    assert(dataPtr->myValue == 2);
    assert(dataPtr->ourValue == 5);
    assert(testVector.size() == 3);
    // Test allocating a new value.
    dataPtr = &(testVector.getNextEmpty());
    assert(dataPtr->myValue == 5);
    assert(dataPtr->ourValue == 6);
    assert(testVector.size() == 4);

    // Clear both vectors.
    testVector2.clear();
    testVector.reset();
    assert(testVector.size() == 0);
    assert(testInvalidGetIndex(testVector, 0));
    assert(testInvalidGetIndex(testVector, 1));
    assert(testVector2.size() == 0);
    assert(testInvalidGetIndex(testVector2, 0));
    assert(testInvalidGetIndex(testVector2, 1));

    // Get values for the vectors and verify they are reused.
    dataPtr = &(testVector2.getNextEmpty());
    assert(dataPtr->myValue == 3);
    assert(dataPtr->ourValue == 6);
    assert(testVector2.size() == 1);
    dataPtr = &(testVector.getNextEmpty());
    assert(dataPtr->myValue == 0);
    assert(dataPtr->ourValue == 6);
    assert(testVector.size() == 1);
    dataPtr = &(testVector2.getNextEmpty());
    assert(dataPtr->myValue == 4);
    assert(dataPtr->ourValue == 6);
    assert(testVector2.size() == 2);
    dataPtr = &(testVector2.getNextEmpty());
    assert(dataPtr->myValue == 6);
    assert(dataPtr->ourValue == 7);
    assert(testVector2.size() == 3);
    dataPtr = &(testVector.getNextEmpty());
    assert(dataPtr->myValue == 1);
    assert(dataPtr->ourValue == 7);
    assert(testVector.size() == 2);
    dataPtr = &(testVector.getNextEmpty());
    assert(dataPtr->myValue == 2);
    assert(dataPtr->ourValue == 7);
    assert(testVector.size() == 3);
    dataPtr = &(testVector.getNextEmpty());
    assert(dataPtr->myValue == 5);
    assert(dataPtr->ourValue == 7);
    assert(testVector.size() == 4);
    dataPtr = &(testVector.getNextEmpty());
    assert(dataPtr->myValue == 7);
    assert(dataPtr->ourValue == 8);
    assert(testVector.size() == 5);
}


bool ReusableVectorTest::testInvalidGetIndex(ReusableVector<ReusableVectorTestDataType>& testVector, int index)
{
    bool caught = false;
    try
    {
        testVector.get(index);
    }
    catch(std::exception& e)
    {
        caught = true;
        assert(strcmp(e.what(), "ReusableVector::get called with out of range index.") == 0);
    }
    return(caught);
}


ReusableVectorTestDataType::ReusableVectorTestDataType()
{
    myValue = ourValue++;
}


ReusableVectorTestDataType::~ReusableVectorTestDataType()
{
    ++ourNumDestructs;
}

