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

#include <iomanip>
#include "BaseComposition.h"

// Constructor
// Initialize the base to ascii map based on the specified maptype.
BaseComposition::BaseComposition():
   myBaseAsciiMap()
{
}


// Update the composition for the specified index with the specified character.
// Return false if the character is not a valid raw sequence character.
// Return true if it is valid.
bool BaseComposition::updateComposition(unsigned int rawSequenceCharIndex, 
                                        char baseChar)
{
   bool validIndex  = true;

   // Each time we return to index 0, reset the primer count in the base/ascii
   // map.
   if(rawSequenceCharIndex == 0)
   {
      myBaseAsciiMap.resetPrimerCount();
   }

   // Check to see if the vector size is already sized to include this 
   // index.  If it is not sized appropriately, add entries until it contains
   // the rawSequenceCharIndex.
   while(rawSequenceCharIndex >= myBaseCountVector.size())
   {
      // Add an entry of the base count array object to the vector.
      BaseCount baseCountEntry;
      myBaseCountVector.push_back(baseCountEntry);
   }

   // Get the base info for the specified character. 
   int baseIndex = myBaseAsciiMap.getBaseIndex(baseChar);
   
   // Increment the count for the given character.  This method returns false
   // if the character's index falls outside the range of the base array.
   // This relies on the myBaseAsciiMap indexes and the BaseCOunt object array
   // to use the same indexing values for valid bases.
   validIndex = 
      myBaseCountVector[rawSequenceCharIndex].incrementCount(baseIndex);
   
   // Return whether or not the specified character was valid.
   return(validIndex);
}


// Print the composition.
void BaseComposition::print()
{
   std::cout << std::endl << "Base Composition Statistics:" << std::endl;
   std::cout.precision(2);
   // This assumes the relationship between indexes that are printed
   // by a BaseCount object to be in a specific order based on ATGCN.
   std::cout << std::fixed << "Read Index" 
             << "\t%A" << "\t%C" << "\t%G" << "\t%T" << "\t%N" << "\tTotal Reads At Index" 
             << std::endl;
   for(unsigned int i = 0; i < myBaseCountVector.size(); i++)
   {
      std::cout << std::setw(10) << i << " ";
      myBaseCountVector[i].printPercent();
   }
   std::cout << std::endl;
}


// Clear the vector.
void BaseComposition::clear()
{
   myBaseCountVector.clear();
}
