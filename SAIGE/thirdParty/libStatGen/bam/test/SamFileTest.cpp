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

#include "SamFileTest.h"
#include "SamFile.h"

void testSamFile()
{
    SamFileHeader header;

    // Test open for read via the constructor with return.
    SamFile samInConstructorReadDefault("testFiles/testSam.sam", 
                                        SamFile::READ,
                                        ErrorHandler::RETURN);
    assert(samInConstructorReadDefault.WriteHeader(header) == false);
    assert(samInConstructorReadDefault.ReadHeader(header) == true);

    // Test open for write via the constructor.
    SamFile samInConstructorWrite("results/newWrite.sam", SamFile::WRITE,
                                  ErrorHandler::RETURN);
    assert(samInConstructorWrite.ReadHeader(header) == false);
    assert(samInConstructorWrite.WriteHeader(header) == true);

    // Test open for read via the constructor
    SamFile samInConstructorRead("testFiles/testSam.sam", SamFile::READ);
    bool caughtException = false;
    try
    {
        assert(samInConstructorRead.WriteHeader(header) == false);
    }
    catch (std::exception& e) 
    {
        caughtException = true;
    }
    assert(caughtException);
    assert(samInConstructorRead.ReadHeader(header) == true);

    // Test open for write via child class.
    SamFileWriter samWriteConstructor("results/newWrite1.sam");
    caughtException = false;
    try
    {
        assert(samWriteConstructor.ReadHeader(header) == false);
    }
    catch (std::exception& e) 
    {
        caughtException = true;
    }
    assert(caughtException);
    assert(samWriteConstructor.WriteHeader(header) == true);

    // Test open for read via child class.
    SamFileReader samReadConstructor("testFiles/testSam.sam");
    caughtException = false;
    try
    {
        assert(samReadConstructor.WriteHeader(header) == false);
    }
    catch (std::exception& e) 
    {
        caughtException = true;
    }
    assert(caughtException);
    assert(samReadConstructor.ReadHeader(header) == true);
}
