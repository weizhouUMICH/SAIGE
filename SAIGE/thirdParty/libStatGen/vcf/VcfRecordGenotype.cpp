/*
 *  Copyright (C) 2011-2012  Regents of the University of Michigan,
 *                           Hyun Min Kang, Matthew Flickenger, Matthew Snyder,
 *                           and Goncalo Abecasis
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
#include "VcfRecordGenotype.h"
#include <stdlib.h>

std::set <std::string> VcfRecordGenotype::ourStoreFields;


void VcfRecordGenotype::storeAllFields()
{
    ourStoreFields.clear();
}


void VcfRecordGenotype::addStoreField(const char* field)
{
    ourStoreFields.insert(field);
}

bool VcfRecordGenotype::storeField(std::string& field)
{
    if(ourStoreFields.size() == 0)
    {
        // No fields were set so read all fields.
        return(true);
    }
    return(ourStoreFields.find(field) != ourStoreFields.end());
}

VcfRecordGenotype::VcfRecordGenotype()
{
    reset();
}


VcfRecordGenotype::~VcfRecordGenotype()
{
}


bool VcfRecordGenotype::read(IFILE filePtr)
{
    return(read(filePtr, NULL));
}


bool VcfRecordGenotype::read(IFILE filePtr, VcfSubsetSamples* subsetInfo)
{
    // Needed for skipping samples.
    static const std::string fieldEndChars = "\n\t";
    static const int tabPos = 1;

    // Clear out any previously set values.
    reset();
    
    if(ifeof(filePtr))
    {
        // End of file, just return false.
        return(false);
    }
    
    // Read the format.
    if(!myFormat.read(filePtr))
    {
        // No more fields
        return(false);
    }

   // Read all the samples until the end of the line.
    VcfGenotypeSample* nextSample = NULL;
    bool moreSamples = true;
    int sampleIndex = 0;
    while(moreSamples)
    {
        // Done reading the format field, so read the samples.
        // Check if this sample should be kept.
        if(subsetInfo != NULL)
        {
            // Check if this sample should be kept.
            if(!subsetInfo->keep(sampleIndex))
            {
                // this sample should not be kept.
                if(filePtr->readTilChar(fieldEndChars) != tabPos)
                {
                    // Stopped on new line or end of file instead of
                    // a tab, so no more samples to read.
                    moreSamples = false;
                }
                ++sampleIndex;
                continue;
            }
        }

        // Read this sample.
        nextSample = &(mySamples.getNextEmpty());
        if(nextSample == NULL)
        {
            throw(std::runtime_error("VCF failed to get another sample."));
        }
        if(!nextSample->read(filePtr, myFormat))
        {
            // No more fields.
            moreSamples = false;
        }
       ++sampleIndex;
    }
    
    // Return whether or not a tab was found at the end of the field.
    return(false);
}


bool VcfRecordGenotype::write(IFILE filePtr)
{
    bool status = true;

    // Check if there are any fields to write.
    if(myFormat.getNumFields() == 0)
    {
        // Nothing to write.
        return(true);
    }

    // Write the format.
    status &= myFormat.write(filePtr);

    // Loop through and write each sample.
    for(int i = 0; i < mySamples.size(); i++)
    {
        status &= mySamples.get(i).write(filePtr);
    }
    return(status);
}


void VcfRecordGenotype::reset()
{
    myFormat.reset();
    mySamples.reset();
}


const std::string* VcfRecordGenotype::getString(const std::string& key, 
                                                int sampleNum)
{
    if(sampleNum >= mySamples.size())
    {
        // Out of range sample index.
        return(NULL);
    }
    // Get the field from the sample.
    return(mySamples.get(sampleNum).getString(key));
}


bool VcfRecordGenotype::setString(const std::string& key, 
                                  int sampleNum, 
                                  const std::string& value)
{
    if(sampleNum >= mySamples.size())
    {
        // Out of range sample index.
        return(NULL);
    }
    // Set the field in the sample.
    return(mySamples.get(sampleNum).setString(key, value));
}


int VcfRecordGenotype::getGT(int sampleNum, unsigned int gtIndex)
{
    if(sampleNum >= mySamples.size())
    {
        // Out of range sample index.
        return(VcfGenotypeSample::INVALID_GT);
    }
    // Get the field from the sample.
    return(mySamples.get(sampleNum).getGT(gtIndex));

}


void VcfRecordGenotype::setGT(int sampleNum, unsigned int gtIndex, int newGt)
{
    if(sampleNum >= mySamples.size())
    {
        // Out of range sample index.
        throw(std::runtime_error("setGT called with out of range sample."));
    }
    // Set the field for the sample.
    mySamples.get(sampleNum).setGT(gtIndex, newGt);

}


int VcfRecordGenotype::getNumGTs(int sampleNum)
{
    if(sampleNum >= mySamples.size())
    {
        // Out of range sample index, no GTs.
        return(0);
    }
    // Get the field from the sample.
    return(mySamples.get(sampleNum).getNumGTs());

}


bool VcfRecordGenotype::allPhased()
{
    for(int i = 0; i < mySamples.size(); i++)
    {
        if(!mySamples.get(i).isPhased() || mySamples.get(i).isUnphased())
        {
            // found a sample that is not phased or is unphased, so
            // return false.
            return(false);
        }
    }
    // All phased.
    return(true);
}


bool VcfRecordGenotype::allUnphased()
{
    for(int i = 0; i < mySamples.size(); i++)
    {
        if(!mySamples.get(i).isUnphased() || mySamples.get(i).isPhased())
        {
            // found a sample that is not unphased or is phased, so
            // return false.
            return(false);
        }
    }
    // All unphased.
    return(true);
}


bool VcfRecordGenotype::hasAllGenotypeAlleles()
{
    for(int i = 0; i < mySamples.size(); i++)
    {
        if(!mySamples.get(i).hasAllGenotypeAlleles())
        {
            // found a sample that does not have all genotype alleles, so
            // return false.
            return(false);
        }
    }
    // All have all genotype alleles.
    return(true);
}


bool VcfRecordGenotype::isPhased(int sampleNum)
{
    if(sampleNum >= mySamples.size())
    {
        // Out of range sample index.
        return(false);
    }
    return(mySamples.get(sampleNum).isPhased());
}


bool VcfRecordGenotype::isUnphased(int sampleNum)
{
    if(sampleNum >= mySamples.size())
    {
        // Out of range sample index.
        return(false);
    }
    return(mySamples.get(sampleNum).isUnphased());
}


bool VcfRecordGenotype::hasAllGenotypeAlleles(int sampleNum)
{
    if(sampleNum >= mySamples.size())
    {
        // Out of range sample index.
        return(false);
    }
    return(mySamples.get(sampleNum).hasAllGenotypeAlleles());
}
