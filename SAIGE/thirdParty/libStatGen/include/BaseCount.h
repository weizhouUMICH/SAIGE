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

#ifndef __BASE_COUNT_H__
#define __BASE_COUNT_H__


/// This class is a wrapper around an array that has one index per base and an
/// extra index for a total count of all bases.  This class is used to keep
/// a count of the number of times each index has occurred.
/// It can print a percentage of the occurrence of each base against the total 
/// number of bases.
class BaseCount
{
public:
    /// Constructor, initializes the array to be all 0s.
    BaseCount();

    /// Update the count for the specified index as well as the overall count 
    /// (The last index).
    /// \return false if the specified index is < 0 or >= myBaseSize-1, otherwise
    /// returns true.  The reason it returns false if it is equal to the size-1
    /// is because the last index is used to track an overall count.
    bool incrementCount(int baseIndex);

    // Print the percentage for each index, 0 to myBaseSize-2, also print
    // the total number of entries (index myBaseSize-1).
    void printPercent();

 private:
    // Constant to size the array and implement the logic for loops as well
    // as tracking the last index for keeping an overall count.
    static const int myBaseSize = 6;

    // Array used to track the occurences of each index.  The last index
    // tracks the total number of occurrences of all the other indexes.
    int myBaseCount[myBaseSize];
};
#endif
