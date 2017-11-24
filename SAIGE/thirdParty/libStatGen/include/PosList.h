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

#ifndef __POSLIST_H__
#define __POSLIST_H__

#include <vector>

/// Store refID/position, but does not store values < 0.
class PosList
{
public:
    /// Constructor
    PosList();

    /// Reserves space for numRefs reference ids and numPositions for each id.
    PosList(int numRefs, int numPositions);

    /// Destructor
    virtual ~PosList();

    /// Add the specified reference id/position (negative values will not be
    /// added).
    void addPosition(int refID, int refPosition);

    /// Return whether or not this list contains the specified reference ID
    /// and position (negative values will automatically return false).
    bool hasPosition(int refID, int refPosition);

protected:
    PosList(const PosList& p);

    void initVars();

    // 2-D vector.
    // indexed by [referenceID][position].
    std::vector < std::vector<bool> > myPosList;

    int myNumRefs;
    int myNumPos;
};


#endif
