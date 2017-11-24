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

#include "SamHeaderHD.h"

// Constructor
SamHeaderHD::SamHeaderHD()
{
    // Add required tags for this type.
    myType = SamHeaderRecord::HD;
    myTypeString = "HD";
    addRequiredTag("VN");
    myKeyTag.clear();   
}

   
// Destructor
SamHeaderHD::~SamHeaderHD()
{
}


const char* SamHeaderHD::getSortOrder()
{
    return(getTagValue("SO"));
}


SamHeaderRecord* SamHeaderHD ::createCopy() const
{
    SamHeaderHD* newHD = new SamHeaderHD();
    if(newHD == NULL)
    {
        std::cerr << "Failed to create a copy of an HD Header Record\n" ;
        return(NULL);
    }
    internalCopy(*newHD);

    return(newHD);
}
