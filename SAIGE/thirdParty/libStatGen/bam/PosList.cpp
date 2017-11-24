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

#include "PosList.h"
#include <stdexcept>

PosList::PosList()
    : myNumRefs(24),
      myNumPos(100)
{
    initVars();
}


PosList::PosList(int numRefs, int numPositions)
    : myNumRefs(numRefs),
      myNumPos(numPositions)
{
    initVars();
}

PosList::~PosList()
{
    myPosList.clear();
}


void PosList::addPosition(int refID, int refPosition)
{
    // Check for negative numbers, if so, just return.
    if((refID < 0) || (refPosition < 0))
    {
        return;
    }

    // If the position list is smaller or equal to refID, it cannot handle an index,
    // so increase the size.
    if(myPosList.size() <= (unsigned int)refID)
    {
        // The position list does not currently have space for this reference id,
        // so add it.
        myPosList.resize(refID+1, std::vector<bool>(myNumPos, false));
        myNumRefs = refID + 1;
    }

    // The matrix is now sized for this reference id.
    // Check to see if this id holds this position.
    if((myPosList[refID]).size() <= (unsigned int)refPosition)
    {
        // The index for this position has not yet been created,
        // so increase the size for it.
        if(myNumPos <= refPosition)
        {
            // Our number of positions is smaller than
            // the current reference id, so reset
            // myNumPos for future use to be this position +1.
            myNumPos = refPosition + 1;
        }
        // Increase the size for this reference id to hold at least myNumPos.
        (myPosList[refID]).resize(myNumPos, false);
    }

    // It now holds this position, so set it to true.
    myPosList[refID][refPosition] = true;
}

bool PosList::hasPosition(int refID, int refPosition)
{
    // Check for negative numbers, if so, just return false, not found.
    if((refID < 0) || (refPosition < 0))
    {
        return(false);
    }
    bool found = false;
    try
    {
        if((myPosList.at(refID)).at(refPosition))
        {
            found = true;
        }
    }
    catch (std::out_of_range& oor)
    {
            // Nothing to do here, if it was out of range, then
            // the position was not found (already set to false).
    }
    return(found);
}


void PosList::initVars()
{
    myPosList.clear();
    myPosList.resize(myNumRefs, std::vector<bool>(myNumPos, false));
}
