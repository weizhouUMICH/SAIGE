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

#include "ReadFiles.h"
#include "Validate.h"
#include "GlfException.h"
#include <assert.h>

void testReadGlf()
{
    GlfFile inGlf;
    assert(inGlf.openForRead("testFiles/testGlf.glf"));

    // Read the GLF Header.
    GlfHeader glfHeader;
    assert(inGlf.readHeader(glfHeader));

    validateHeader(glfHeader);

    // TODO, validate the rest of the file.
//    GlfRecord glfRecord;
//     assert(inGlf.ReadRecord(glfHeader, glfRecord) == true);
//     validateRead1(glfRecord);

    // Try opening a file that doesn't exist.
    bool exceptionCaught = false;
    try
    {
        inGlf.openForRead("testFiles/unknown");
    }
    catch(GlfException e)
    {
        exceptionCaught = true;
    }
    assert(exceptionCaught);


}

