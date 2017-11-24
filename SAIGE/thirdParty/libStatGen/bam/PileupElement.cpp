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

#include "PileupElement.h"


GenomeSequence* PileupElement::myRefPtr = NULL;


PileupElement::PileupElement()
    : myRefPosition(UNSET_POSITION),
      myChromosome("")
{
}


// NOTE that this method does not actually copy, it just resets.
PileupElement::PileupElement(const PileupElement& q)
    : myRefPosition(UNSET_POSITION),
      myChromosome("")
{
}

PileupElement::~PileupElement()
{
}


// Add an entry to this pileup element.  
void PileupElement::addEntry(SamRecord& record)
{
    if(myChromosome.empty())
    {
        // First entry, save chromosme name.
        myChromosome = record.getReferenceName();
    }
}


// Perform the alalysis associated with this class.  May be a simple print, 
// a calculation, or something else.  Typically performed when this element
// has been fully populated by all records that cover the reference position.
void PileupElement::analyze()
{
    if(myRefPosition != UNSET_POSITION)
    {
        std::cout << myChromosome << "\t" << myRefPosition << "\n";
    }
}


// Resets the entry, setting the new position associated with this element.
void PileupElement::reset(int32_t refPosition)
{
    myChromosome.clear();
    myRefPosition = refPosition;
}


char PileupElement::getRefBase()
{
    if(myRefPtr != NULL)
    {
        // Add 1 to pos because getBase expects 1-based index.
        return(myRefPtr->getBase(myChromosome.c_str(), myRefPosition+1));
    }
    return('N');
}


// Resets the entry, setting the new position associated with this element.
void PileupElement::setReference(GenomeSequence* reference)
{
    myRefPtr = reference;
}
