
/*
 *  Copyright (C) 2011  Regents of the University of Michigan,
 *                      Hyun Min Kang, Matthew Flickenger, Matthew Snyder,
 *                      and Goncalo Abecasis
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

#include "VcfHelper.h"

void VcfHelper::parseString(const std::string& inputString, 
                            char delim,
                            ReusableVector<std::string>& outputVector)
{
    if(inputString.empty())
    {
        // Nothing to parse, so just return.
        return;
    }
    std::string* outputStringPtr = &(outputVector.getNextEmpty());
    for(unsigned int i = 0; i < inputString.size(); i++)
    {
        if(inputString[i] == delim)
        {
            // Get a new string to write into and continue
            // to the next character.
            outputStringPtr = &(outputVector.getNextEmpty());
        }
        else
        {
            // Append the character.
            outputStringPtr->push_back(inputString[i]);
        }
    }
}

