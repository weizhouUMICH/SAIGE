/*
 *  Copyright (C) 2012  Regents of the University of Michigan
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

#include "VcfGenotypeField.h"

bool VcfGenotypeField::write(IFILE filePtr)
{
    int numWritten = 0;
    int numExpected = 0;
    std::string* subFieldPtr = NULL;

    // Write the tab before this field.
    numWritten += ifprintf(filePtr, "\t");
    ++numExpected;

    // Loop through and write each subfield.
    for(int i = 0; i < myGenotypeSubFields.size(); i++)
    {
        subFieldPtr = &(myGenotypeSubFields.get(i));
        if(i == 0)
        {
            // First entry, so no ':'
            numWritten += ifprintf(filePtr, "%s", subFieldPtr->c_str());
            numExpected += subFieldPtr->size();
        }
        else
        {
            // Not first entry, so use a ':'
            numWritten += ifprintf(filePtr, ":%s", subFieldPtr->c_str());
            numExpected += 1 + subFieldPtr->size();
        }

    } // End loop through entries.
    return(numWritten == numExpected);
}

 
VcfGenotypeField::SUBFIELD_READ_STATUS 
    VcfGenotypeField::readGenotypeSubField(IFILE filePtr, 
                                           std::string* stringDest)
{
    if(ifeof(filePtr))
    {
        // End of file, so just return END_OF_RECORD.
        return(END_OF_RECORD);
    }
    
    static const std::string fieldStopChars = "\n\t:";
    // Read/parse the field.
    int pos = 0;
    if(stringDest == NULL)
    {
        pos = filePtr->readTilChar(fieldStopChars);
    }
    else
    {
        pos = filePtr->readTilChar(fieldStopChars, *stringDest);
    }
    if(pos == 2)
    {
        // ':' 
        return(MORE_SUBFIELDS);
    }
    else if(pos == 1)
    {
        // '\t'
        return(END_OF_FIELD);
    }
    // '\n' or EOF
    return(END_OF_RECORD);
}


VcfGenotypeField::SUBFIELD_READ_STATUS 
    VcfGenotypeField::getReadStatus(int stopChar)
{
    if(stopChar >= 2)
    {
        // ':' 
        return(MORE_SUBFIELDS);
    }
    else if(stopChar == 1)
    {
        // '\t'
        return(END_OF_FIELD);
    }
    // '\n' or EOF
    return(END_OF_RECORD);
}
