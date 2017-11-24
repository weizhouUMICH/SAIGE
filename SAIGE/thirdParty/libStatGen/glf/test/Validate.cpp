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

#include "Validate.h"
#include <assert.h>

//const std::string GLF_HEADER_TEXT = "";

void validateRead1(GlfRecord& glfRecord)
{
    //////////////////////////////////////////
    // Validate Record 1
    // Create record structure for validating.
}


void validateHeader(GlfHeader& glfHeader)
{
    ////////////////////////////////////////////////////////
    // Get the text from the header and verify it is the expected value.
    std::string textString = "DUMMY";
    assert(glfHeader.getHeaderTextString(textString));
    assert(textString == GLF_HEADER_TEXT);
}
