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

#include "VcfSubsetSamples.h"

void VcfSubsetSamples::reset()
{
    mySampleSubsetIndicator.clear();
    mySampleNames.clear();
}


void VcfSubsetSamples::init(const VcfHeader& header, bool include)
{
    // Get the number of samples from the header.
    unsigned int origNumSamples = header.getNumSamples();

    // Resize the sampleSubsetIndicator to nothing to clear it out.
    mySampleSubsetIndicator.resize(0);

    // Now resize sampleSubsetIndicator to indicate that all of the original
    // samples are to be kept or not kept based on the include parameter.
    // mySampleSubsetIndicator is sized to the original number of samples
    // so it can be used when reading records to determine which ones should
    // be removed/kept.
    mySampleSubsetIndicator.resize(origNumSamples, include);

    // Copy the vector of original sample names.
    mySampleNames.clear();
    mySampleNames.resize(origNumSamples);
    for(unsigned int i = 0; i < origNumSamples; i++)
    {
        mySampleNames[i] = header.getSampleName(i);
    }
}


bool VcfSubsetSamples::addIncludeSample(const char* sampleName)
{
    // Look for the sample name.
    for(unsigned int i = 0; i < mySampleNames.size(); i++)
    {
        if(mySampleNames[i] == sampleName)
        {
            // Found the sample index.
            if(mySampleSubsetIndicator.size() <= i)
            {
                // SampleSubsetIndicator not setup properly.
                return(false);
            }
            mySampleSubsetIndicator[i] = true;
            return(true);
        }
    }
    // Did not find the sample, so can't include it.
    return(false);
}


bool VcfSubsetSamples::addExcludeSample(const char* sampleName)
{
    // Look for the sample name.
    for(unsigned int i = 0; i < mySampleNames.size(); i++)
    {
        if(mySampleNames[i] == sampleName)
        {
            // Found the sample index.
            if(mySampleSubsetIndicator.size() <= i)
            {
                // SampleSubsetIndicator not setup properly.
                return(false);
            }
            mySampleSubsetIndicator[i] = false;
            return(true);
        }
    }
    // Did not find the sample, so can't include it.
    return(false);
}


bool VcfSubsetSamples::init(VcfHeader& header, 
                            const char* includeFileName, 
                            const char* excludeSample, 
                            const char* excludeFileName,
                            const char* delims)
{
    // Setup the sample lists to include/exclude.
    std::set<std::string> includeList;
    std::set<std::string> excludeList;
    if(includeFileName != NULL)
    {
        if(!readSamplesFromFile(includeFileName, includeList, delims))
        {
            // Failed, so return.
            return(false);
        }
    }

    if(excludeFileName != NULL)
    {
        if(!readSamplesFromFile(excludeFileName, excludeList, delims))
        {
            // Failed, so return.
            return(false);
        }
    }
    if(excludeSample != NULL)
    {
        excludeList.insert(excludeSample);
    }

    int origNumSamples = header.getNumSamples();

    // Resize the sampleSubsetIndicator to nothing to clear it out.
    mySampleSubsetIndicator.resize(0);

    // Now resize sampleSubsetIndicator to indicate that all of the original
    // samples are to be kept.  The ones that are not to be kept will be 
    // modified to be unkept (false).
    // mySampleSubsetIndicator is sized to the original number of samples
    // so it can be used when reading records to determine which ones should
    // be removed/kept.
    mySampleSubsetIndicator.resize(origNumSamples, true);

    // if no samples, return.
    if(origNumSamples == 0)
    {
        return(true);
    }

    // Now that the sample lists to include/exclude are setup and the
    // indicator vector is setup, subset the header removing samples that 
    // should not be kept (not in the include list if set or in the exclude 
    // list). Loop from the back of the samples to the beginning since
    // removing samples changes the index of all following samples.
    for(int i = (origNumSamples-1); i >= 0; i--)
    {
        // Check if the sample should be kept.
        const char* sampleName = header.getSampleName(i);
        // Remove the sample if the includeList was specified and the sample
        // was not in it or if the excludeList was specified and the sample 
        // was in it.
        if((!includeList.empty() && 
            (includeList.count(sampleName) == 0)) ||
           (!excludeList.empty() && 
            (excludeList.count(sampleName) != 0)))
        {
            // This sample should be removed.
            header.removeSample(i);
            mySampleSubsetIndicator[i] = false;
        }
    }
    return(true);
}


bool VcfSubsetSamples::keep(unsigned int sampleIndex)
{
    if(sampleIndex >= mySampleSubsetIndicator.size())
    {
        // index out of range.
        return(false);
    }
    return(mySampleSubsetIndicator[sampleIndex]);
}


bool VcfSubsetSamples::readSamplesFromFile(const char* fileName, 
                                           std::set<std::string>& sampleList,
                                           const char* delims)
{
    // Open the file.
    IFILE sampleFile = ifopen(fileName, "r");

    if(sampleFile == NULL)
    {
        // Failed to open.
        return(false);
    }

    // read the file.
    std::string tempString;

    std::string delimString = delims;
    delimString += '\n';

    int readResult = 0;
    while(readResult != -1)
    {
        readResult = sampleFile->readTilChar(delimString, tempString);

        // Check to see if something was read (tempString is not empty).
        if(!tempString.empty())
        {
            // sample name found, so add it to the container.
            sampleList.insert(tempString);
        }
        // Clear the string being read into.
        tempString.clear();
    }
    return(true);
}
