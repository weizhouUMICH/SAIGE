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

#ifndef __BASE_COMPOSITION_H__
#define __BASE_COMPOSITION_H__

#include <map>

#include "BaseAsciiMap.h"
#include "BaseCount.h"

/// Class that tracks the composition of base by read location.
class BaseComposition
{
public:
    /// Constructor.
    BaseComposition();

    /// Update the composition for the specified index with the specified
    /// character.
    /// \return false if the character is not a valid raw sequence character,
    /// true if it is valid.
    bool updateComposition(unsigned int rawSequenceCharIndex, char baseChar);

    /// Get the space type for this composition.
    BaseAsciiMap::SPACE_TYPE getSpaceType()
    {
        return(myBaseAsciiMap.getSpaceType());
    }

    /// Reset the base map type for this composition.
    void resetBaseMapType()
    {
        myBaseAsciiMap.resetBaseMapType();
    };

    /// Set the base map type for this composition.
    void setBaseMapType(BaseAsciiMap::SPACE_TYPE spaceType)
    {
        myBaseAsciiMap.setBaseMapType(spaceType);
    }

    /// Print the composition.
    void print();

    /// Clear the composition stored in the base count vector.
    void clear();

private:
    // Map of bases used to determine if a character is valid and if so
    // maps it to a number.
    BaseAsciiMap myBaseAsciiMap;

    // Vector used to store the occurrence of each base type at a given 
    // read location.
    vector<BaseCount> myBaseCountVector;
};
#endif
