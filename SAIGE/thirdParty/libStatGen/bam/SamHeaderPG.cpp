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

#include "SamHeaderPG.h"

// Constructor
SamHeaderPG::SamHeaderPG()
{
    // Add required tags for this type.
    myType = SamHeaderRecord::PG;
    myTypeString = "PG";
    addRequiredTag("ID");
    myKeyTag = "ID";   
}

   
// Destructor
SamHeaderPG::~SamHeaderPG()
{
}


SamHeaderRecord* SamHeaderPG::createCopy() const
{
    SamHeaderPG* newPG = new SamHeaderPG();
    if(newPG == NULL)
    {
        std::cerr << "Failed to create a copy of an PG Header Record\n" ;
        return(NULL);
    }
    internalCopy(*newPG);

    return(newPG);
}
