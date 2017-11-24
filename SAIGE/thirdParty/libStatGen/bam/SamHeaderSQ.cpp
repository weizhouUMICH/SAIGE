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

#include "SamHeaderSQ.h"

// Constructor
SamHeaderSQ::SamHeaderSQ()
{
    // Add required tags for this type.
    myType = SamHeaderRecord::SQ;
    myTypeString = "SQ";
    addRequiredTag("SN");
    addRequiredTag("LN");
    myKeyTag = "SN";   
}

   
// Destructor
SamHeaderSQ::~SamHeaderSQ()
{
}


SamHeaderRecord* SamHeaderSQ::createCopy() const
{
    SamHeaderSQ* newSQ = new SamHeaderSQ();
    if(newSQ == NULL)
    {
        std::cerr << "Failed to create a copy of an SQ Header Record\n" ;
        return(NULL);
    }
    internalCopy(*newSQ);

    return(newSQ);
}
