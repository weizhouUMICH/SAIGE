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

#include "VcfGenotypeSample.h"
#include <stdlib.h>
#include <sstream>

const int VcfGenotypeSample::INVALID_GT = -1;
const int VcfGenotypeSample::MISSING_GT = -2;
const std::string VcfGenotypeSample::MISSING_FIELD = ".";

VcfGenotypeSample::VcfGenotypeSample()
    : VcfGenotypeField(),
      myFormatPtr(NULL),
      myPhased(false),
      myUnphased(false),
      myHasAllGenotypeAlleles(false),
      myNewGT(false),
      myGTs()
{
}


VcfGenotypeSample::~VcfGenotypeSample()
{
    myFormatPtr = NULL;
}


bool VcfGenotypeSample::read(IFILE filePtr, VcfGenotypeFormat& format)
{
    static const char* GT_DELIM = "\n\t:|/.";
    static const int END_GT = 2; // Ends at index 2 or less
    static const int PHASED_CHAR_POS = 3;
    static const int UNPHASED_CHAR_POS = 4;
    static const int MISSING_GT_POS = 5;

    // Clear out any previously set values.
    reset();
    
    myFormatPtr = &format;

    int gtIndex = format.getGTIndex();

    // Read the subfields.
    SUBFIELD_READ_STATUS readStatus = MORE_SUBFIELDS;  
    std::string* nextType = NULL;
    int subFieldIndex = 0;
    while(readStatus == MORE_SUBFIELDS)
    {
        // Get the field to write into.
        if(format.storeIndex(subFieldIndex))
        {
            nextType = &(myGenotypeSubFields.getNextEmpty());
            // Check if this is the GT field.
            if(subFieldIndex == gtIndex)
            {
                // There is a GT field, so set that all GT fields are there.
                // if any are missing it will be turned back to false.
                myHasAllGenotypeAlleles = true;
                // This is the GT field, so parse manually looking to see if it
                // is phased and store the genotypes.
                int stopChar = END_GT + 1;
                // Read until a new subfield is found.
                while(stopChar > END_GT)
                {
                    // TODO  have an option to autoparse the genotypes?
                    // todo - store the previous nextType len in order to
                    // do string conversion to ints...
                    stopChar = filePtr->readTilChar(GT_DELIM, *nextType);
                    if(stopChar == PHASED_CHAR_POS)
                    {
                        nextType->push_back('|');
                        myPhased = true;
                    }
                    else if(stopChar == UNPHASED_CHAR_POS)
                    {
                        nextType->push_back('/');
                        myUnphased = true;
                    }
                    else if(stopChar == MISSING_GT_POS)
                    {
                        nextType->push_back('.');
                        myHasAllGenotypeAlleles = false;
                    }
                }
                // Check if this is the END_GT signal.
                readStatus = getReadStatus(stopChar);
            }
            else
            {
                // more subfields to read.
                readStatus = readGenotypeSubField(filePtr, nextType);
            }
        }
        else
        {
            readStatus = readGenotypeSubField(filePtr, NULL);
        }
        ++subFieldIndex;
    }

    // subFieldIndex contains the number of fields in this sample.
    if(subFieldIndex > format.getOrigNumFields())
    {
        throw(std::runtime_error("VCF Number of Fields in a Sample does not match the Format."));
    }
    else if(subFieldIndex < format.getOrigNumFields())
    {
        // If there are no fields for this sample, enter the missing value.
        if(myGenotypeSubFields.size() == 0)
        {
            myGenotypeSubFields.getNextEmpty() = MISSING_FIELD;
        }
    }

    // Return true if there is a tab - it is just END_OF_FIELD.
    return(readStatus == END_OF_FIELD);
}


bool VcfGenotypeSample::write(IFILE filePtr)
{
    if(myNewGT)
    {
        updateGTString();
    }
    return(VcfGenotypeField::write(filePtr));
}

const std::string* VcfGenotypeSample::getString(const std::string& key)
{
    if(myFormatPtr == NULL)
    {
        return(NULL);
    }

    int index = myFormatPtr->getIndex(key);
    if(index != VcfGenotypeFormat::GENOTYPE_INDEX_NA)
    {
        // Check if it is out of range for this sample - means it 
        // is missing for this sample.
        if(index >= myGenotypeSubFields.size())
        {
            // missing for this sample.
            return(&MISSING_FIELD);
        }

        if((key == "GT") && myNewGT)
        {
            updateGTString();
        }
        return(&(myGenotypeSubFields.get(index)));
    }
    // key was not found, so return NULL.
    return(NULL);
}


bool VcfGenotypeSample::setString(const std::string& key, const std::string& value)
{
    if(myFormatPtr == NULL)
    {
        return(false);
    }

    int index = myFormatPtr->getIndex(key);
    if(index != VcfGenotypeFormat::GENOTYPE_INDEX_NA)
    {
        // Found the type, so set it.
        myGenotypeSubFields.get(index) = value;

        if(key == "GT")
        {
            myGTs.clear();
            myNewGT = false;
            if(value.find('|') != std::string::npos)
            {
                myPhased = true;
            }
            if(value.find('/') != std::string::npos)
            {
                myUnphased = true;
            }
            if(value.find('.') != std::string::npos)
            {
                myHasAllGenotypeAlleles = false;
            }
            else
            {
                myHasAllGenotypeAlleles = true;
            }
        }
        return(true);
    }
    // field was not found, so return false.
    return(false);
}


int VcfGenotypeSample::getGT(unsigned int index)
{
    if(myGTs.empty())
    {
        if(!parseGT())
        {
            // Failed to parse GT, so return INVALID_GT.
            return(INVALID_GT);
        }
    }
    
    if(index < myGTs.size())
    {
        return(myGTs[index]);
    }
    // Out of range index.
    return(INVALID_GT);
}


void VcfGenotypeSample::setGT(unsigned int index, int newGt)
{
    if(myGTs.empty())
    {
        if(!parseGT())
        {
            // Failed to parse GT, so return INVALID_GT.
            throw(std::runtime_error("VCF failed to parse GT."));
        }
    }
    
    if(index < myGTs.size())
    {
        if(myGTs[index] != newGt)
        {
            myNewGT = true;
            myGTs[index] = newGt;
        }
    }
    else
    {
        // Out of range index.
        throw(std::runtime_error("VCF setGT called with out of range GT index."));
    }
}


int VcfGenotypeSample::getNumGTs()
{
    if(myGTs.empty())
    {
        if(!parseGT())
        {
            return(0);
        }
    }
    return(myGTs.size());
}


void VcfGenotypeSample::internal_reset()
{
    myFormatPtr = NULL;
    myPhased = false;
    myUnphased = false;
    myHasAllGenotypeAlleles = false;
    myGTs.clear();
    myNewGT = false;
}


bool VcfGenotypeSample::parseGT()
{
    // Parse the GT.
    const std::string* gtStr = getString("GT");
    myGTs.clear();
    myNewGT = false;
    if(gtStr == NULL)
    {
        // GT field not found.
        return(false);
    }

    // Parse til the end of the GT string
    char* startPos = NULL;
    char* endPos = (char*)gtStr->c_str();
    while((endPos != NULL) && (*endPos != '\0'))
    {
        startPos = endPos;
        if(*startPos == '.')
        {
            endPos = startPos + 1;
            // unknown, so set this index to be MISSING_GT.
            myGTs.push_back(MISSING_GT);
            continue;
        }
        if(*startPos == '|')
        {
            endPos = startPos + 1;
            continue;
        }
        if(*startPos == '/')
        {
            endPos = startPos + 1;
            continue;
        }
        // Should be an int, so parse it.
        unsigned long gtLong = strtoul(startPos, &endPos, 10);
        myGTs.push_back((int)gtLong);
    }
    return(true);
}


void VcfGenotypeSample::updateGTString()
{
    if(myNewGT)
    {
        int index = myFormatPtr->getIndex("GT");
        if(index != VcfGenotypeFormat::GENOTYPE_INDEX_NA)
        {
            // Check if it is out of range for this sample - means it 
            // is missing for this sample.
            if(index < myGenotypeSubFields.size())
            {
                std::stringstream gtSS;
                char phaseChar = '/';
                if(myPhased)
                {
                    phaseChar = '|';
                }
                gtSS << myGTs[0];
                for(unsigned int i = 1; i < myGTs.size(); i++)
                {
                    gtSS << phaseChar << myGTs[i];
                }
                myGenotypeSubFields.get(index) = gtSS.str();
                myNewGT = false;
            }
        }
    }
}
