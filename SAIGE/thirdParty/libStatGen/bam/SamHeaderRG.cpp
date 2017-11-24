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

#include "SamHeaderRG.h"

// Constructor
SamHeaderRG::SamHeaderRG()
{
    // Add required tags for this type.
    myType = SamHeaderRecord::RG;
    myTypeString = "RG";
    addRequiredTag("ID");
    myKeyTag = "ID";   
}

   
// Destructor
SamHeaderRG::~SamHeaderRG()
{
}


SamHeaderRecord* SamHeaderRG::createCopy() const
{
    SamHeaderRG* newRG = new SamHeaderRG();
    if(newRG == NULL)
    {
        std::cerr << "Failed to create a copy of an RG Header Record\n" ;
        return(NULL);
    }
    internalCopy(*newRG);

    return(newRG);
}
