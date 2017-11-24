/*
 *  Copyright (C) 2012  Regents of the University of Michigan
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

#include "VcfGenotypeFormat.h"
#include "VcfRecordGenotype.h"

VcfGenotypeFormat::VcfGenotypeFormat()
    : VcfGenotypeField(),
      myGTIndex(GENOTYPE_INDEX_NA)
{
}


VcfGenotypeFormat::~VcfGenotypeFormat()
{
}


bool VcfGenotypeFormat::read(IFILE filePtr)
{
    // Clear out any previously set values.
    reset();
    SUBFIELD_READ_STATUS readStatus = MORE_SUBFIELDS;  
    std::string* nextType = &(myGenotypeSubFields.getNextEmpty());
    int subFieldIndex = 0;
    while(readStatus == MORE_SUBFIELDS)
    {
        // more subfields to read.
        readStatus = readGenotypeSubField(filePtr, nextType);
        // Check if this field should be read/stored.
        if(!VcfRecordGenotype::storeField(*nextType))
        {
            // Do not read/store this field.
            myStoreIndices.push_back(false);
            nextType->clear();
        }
        else
        {
            // Check if this is GT.
            if(*nextType == "GT")
            {
                myGTIndex = subFieldIndex;
            }
            myStoreIndices.push_back(true);
            nextType = &(myGenotypeSubFields.getNextEmpty());
        }
        ++subFieldIndex;
    }
    // Since we didn't use the last type that was retrieved, remove it.
    myGenotypeSubFields.rmLast();

    // Return true if there is a tab - it is just END_OF_FIELD.
    return(readStatus == END_OF_FIELD);
}


int VcfGenotypeFormat::getIndex(const std::string& key)
{
    //  Search for this field of the genotypeFormat.
    for(int i = 0; i < myGenotypeSubFields.size(); i++)
    {
        if(myGenotypeSubFields.get(i) == key)
        {
            // Found the type.
            return(i);
        }
    }
    // field was not found, so return null.
    return(GENOTYPE_INDEX_NA);
}


bool VcfGenotypeFormat::storeIndex(unsigned int index)
{
    if(index > myStoreIndices.size())
    {
        return(false);
    }
    return(myStoreIndices[index]);
}


void VcfGenotypeFormat::internal_reset()
{
    myGTIndex = GENOTYPE_INDEX_NA;
    myStoreIndices.clear();
}
