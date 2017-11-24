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

#include <iostream>
#include <iomanip>
#include "BaseCount.h"

// Constructor.  Initializes the array to be all 0s.
BaseCount::BaseCount()
{
   // Init each element of the array to 0.
   for(int i = 0; i < myBaseSize; i++)
   {
      myBaseCount[i] = 0;
   }
}


// Update the count for the specified index as well as the overall count 
// (The last index).
// Returns false if the specified index is < 0 or >= myBaseSize-1.  The
// reason returns false if it is equal to the size-1 is because the last
// index is used to track an overall count.
bool BaseCount::incrementCount(int baseIndex)
{
   // Check to see if the index is within range (>=0 & < myBaseSize-1)
   // The last entry of the array is invalid since it is used to track
   // total occurrence of all other entries.
   if((baseIndex < myBaseSize-1) && (baseIndex >= 0))
   {
      // Valid index, so increment that index as well as the overall
      // count (index myBaseSize-1) and return true.
      myBaseCount[baseIndex]++;
      myBaseCount[myBaseSize-1]++;
      return true;
   }
   else
   {
      // Invalid index, return false
      return false;
   }
}



// Prints the percentage for each index 0 to myBaseSize-2.  Also prints
// the total number of entries (index myBaseSize-1).
void BaseCount::printPercent()
{
   // Do not divide by 0, so check to see if there are any bases by checking
   // the last index of the array.
   if(myBaseCount[myBaseSize-1] == 0)
   {
      // No entries for any index.
      std::cout << "No Valid Bases found.";
   }
   else
   {
      // Print the percentage for each index.
      for(int i = 0; i < myBaseSize -1; i++)
      {
         double percentage = 
            (myBaseCount[i]/(double)myBaseCount[myBaseSize-1]) * 100;
         std::cout << " " << std::setw(7) << percentage;
      }
      // Print the total number of bases.
      std::cout << "\t" << myBaseCount[myBaseSize-1];
   }
   std::cout << std::endl;
}
