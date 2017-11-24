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
#include "StringTest.h"
#include <assert.h>


int main(int argc, char ** argv)
{

    testAsInteger();

    testReadLine();
}

void testAsInteger()
{
    // Test AsInteger with ints & negative ints.
    String intString = "123";
    String negIntString = "-123";
    assert(intString.AsInteger() == 123);
    assert(negIntString.AsInteger() == -123);

    // Run the same tests with AsInteger that returns a bool and takes
    // in a long to set.
    long retValue;
    assert(intString.AsInteger(retValue));
    assert(retValue == 123);
    assert(negIntString.AsInteger(retValue));
    assert(retValue == -123);


    // Strings that are not integers
    // For AsInteger, it returns just the starting integer portion.
    // For AsInteger that returns a bool and a long set, it returns false
    // and sets the long to the starting int.
    String nonIntString = "abd";
    assert(nonIntString.AsInteger() == 0);
    assert(!nonIntString.AsInteger(retValue));

    nonIntString = "12ab33";
    assert(nonIntString.AsInteger() == 12);
    assert(!nonIntString.AsInteger(retValue));
    assert(retValue == 12);
    nonIntString = "as12ab3a4sd";
    assert(nonIntString.AsInteger() == 0);
    assert(!nonIntString.AsInteger(retValue));
    assert(retValue == 0);
    // Negatives are only recognized as the first characer.
    nonIntString = "-12ab3a4sd";
    assert(nonIntString.AsInteger() == -12);
    assert(!nonIntString.AsInteger(retValue));
    assert(retValue == -12);
    nonIntString = "-as12ab3a4sd";
    assert(nonIntString.AsInteger() == 0);
    assert(!nonIntString.AsInteger(retValue));
    assert(retValue == 0);
    nonIntString = "as-12ab3a4sd";
    assert(nonIntString.AsInteger() == 0);
    assert(!nonIntString.AsInteger(retValue));
    assert(retValue == 0);
    nonIntString = "as12-ab3a4sd";
    assert(nonIntString.AsInteger() == 0);
    assert(!nonIntString.AsInteger(retValue));
    assert(retValue == 0);
}

int temp1 = 0;

void testReadLine()
{
    IFILE filePtr = ifopen("testFiles/testFile.txt", "rb");
    assert(filePtr != NULL);
    
    String line = "";
    line.ReadLine(filePtr);

    assert(line == "  Hello, I am a testFile.  ");

    line.Trim();
    assert(line == "Hello, I am a testFile.");


    // Does not compile in current version, but compiles in old verison.
    // This can be added back in to ensure that it will catch the difference
    // in return value for ReadLine (now: int; used to be: string&)
    //    testMethod(line.ReadLine(filePtr));
    line.ReadLine(filePtr);
    assert(temp1 == 0);
    testMethod(line);
    assert(temp1 == 1);

    //    line.ReadLine(filePtr).Trim();
    line.ReadLine(filePtr);
    line.Trim();

    assert(line == "ThirdLine.");

    ifclose(filePtr);
}


void testMethod(String temp)
{
    temp1 = 1;
}
