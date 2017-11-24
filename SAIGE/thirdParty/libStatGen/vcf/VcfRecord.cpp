/*
 *  Copyright (C) 2010-2012  Regents of the University of Michigan,
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

#include "VcfRecord.h"

VcfRecord::VcfRecord()
{
    reset();
}


VcfRecord::~VcfRecord()
{
}


bool VcfRecord::read(IFILE filePtr, bool siteOnly,
                     VcfRecordDiscardRules& discardRules,
                     VcfSubsetSamples* sampleSubset)
{
    // Clear out any previously set values.
    reset();
    
    if(filePtr == NULL)
    {
        myStatus.setStatus(StatGenStatus::FAIL_ORDER,
                           "Error reading VCF record before opening the file.");
        return(false);
    }

    if(ifeof(filePtr))
    {
        // End of file, just return false.
        return(false);
    }
    
    // Read the chromosome.
    if(!readTilTab(filePtr, myChrom))
    {
        if(myChrom.empty())
        {
            // EOF.
            return(false);
        }
        // Not an empty line.
        myStatus.setStatus(StatGenStatus::FAIL_PARSE, 
                           "Error reading VCF Record CHROM.");
        return(false);
    }
    // Read the 1-based Position
    std::string strPos;
    if(!readTilTab(filePtr, strPos))
    {
        myStatus.setStatus(StatGenStatus::FAIL_PARSE, 
                           "Error reading VCF Record POS.");
        return(false);
    }
    else
    {
        // Read the position, so convert to an integer.
        my1BasedPosNum = atoi(strPos.c_str());
    }
    // Read the ID.
    if(!readTilTab(filePtr, myID))
    {
        myStatus.setStatus(StatGenStatus::FAIL_PARSE, 
                           "Error reading VCF Record ID.");
        return(false);
    }

    if(discardRules.discardForID(myID))
    {
        // Do not keep this id, so consume the rest of the record and
        // return the next record.
        filePtr->discardLine();
        return(read(filePtr, siteOnly, discardRules, sampleSubset));
    }

    // Read the Ref.
    if(!readTilTab(filePtr, myRef))
    {
        myStatus.setStatus(StatGenStatus::FAIL_PARSE, 
                           "Error reading VCF Record REF.");
        return(false);
    }
    // Read the Alt.
    myAltArray.clear();
    if(!readTilTab(filePtr, myAlt))
    {
        myStatus.setStatus(StatGenStatus::FAIL_PARSE, 
                           "Error reading VCF Record ALT.");
        return(false);
    }
    // Read the Qual.
    if(!readTilTab(filePtr, myQual))
    {
        myStatus.setStatus(StatGenStatus::FAIL_PARSE, 
                           "Error reading VCF Record QUAL.");
        return(false);
    }
    else
    {
        if(myQual != ".")
        {
            // Read the quality, so convert to an integer.
            myQualNum = atof(myQual.c_str());
        }
        else
        {
            myQualNum = -1;
        }
    }
    // Read the Filter.
    if(!myFilter.read(filePtr))
    {
        myStatus.setStatus(StatGenStatus::FAIL_PARSE, 
                           "Error reading VCF Record FILTER.");
        return(false);
    }
    // Read the Info (could be the last word in the line or file).
    if(!myInfo.read(filePtr))
    {
        // Found the end of the line after the info field, so return true,
        // successfully read the record.
        return(true);
    }

    if(siteOnly)
    {
        // Do not store genotypes, so just consume the rest of the line.
        filePtr->readTilChar("\n");
    }
    else
    {
        // Not yet at the end of the line, so read the genotype fields
        // (format & samples)
        try
        {
            myGenotype.read(filePtr, sampleSubset);
        }
        catch(std::exception& e)
        {
            myDummyString = "Failed parsing the Genotype Fields of " + myChrom + ":" + 
                std::to_string(my1BasedPosNum) + " (chr:pos) - " + e.what();
            myStatus.setStatus(StatGenStatus::FAIL_PARSE, myDummyString.c_str());
            return(false);
        }
    }
    // Found the end of the line, return true since all required fields
    // were read.
    return(true);
}


bool VcfRecord::write(IFILE filePtr, bool siteOnly)
{
    if(filePtr == NULL)
    {
        myStatus.setStatus(StatGenStatus::FAIL_ORDER,
                           "Error writing VCF record before opening the file.");
        return(false);
    }

    int numWritten = 0;
    int numExpected = 0;
    if(myChrom.length() == 0)
    {
        numWritten += ifprintf(filePtr, ".\t");
        numExpected += 2;
    }
    else
    {
        numWritten += ifprintf(filePtr, "%s\t", myChrom.c_str());
        numExpected += myChrom.length() + 1;
    }
    if(false) //my1BasedPos.length() == 0)
    {
        numWritten += ifprintf(filePtr, ".\t");
        numExpected += 2;
    }
    else
    {
        std::string strPos = std::to_string(my1BasedPosNum);
        numWritten += ifprintf(filePtr, "%s\t", strPos.c_str());
        numExpected += strPos.length() + 1;
    }
    if(myID.length() == 0)
    {
        numWritten += ifprintf(filePtr, ".\t");
        numExpected += 2;
    }
    else
    {
        numWritten += ifprintf(filePtr, "%s\t", myID.c_str());
        numExpected += myID.length() + 1;
    }
    if(myRef.length() == 0)
    {
        numWritten += ifprintf(filePtr, ".\t");
        numExpected += 2;
    }
    else
    {
        numWritten += ifprintf(filePtr, "%s\t", myRef.c_str());
        numExpected += myRef.length() + 1;
    }
    if(myAlt.length() == 0)
    {
        numWritten += ifprintf(filePtr, ".\t");
        numExpected += 2;
    }
    else
    {
        numWritten += ifprintf(filePtr, "%s\t", myAlt.c_str());
        numExpected += myAlt.length() + 1;
    }
    if(myQual.length() == 0)
    {
        numWritten += ifprintf(filePtr, ".\t");
        numExpected += 2;
    }
    else
    {
        numWritten += ifprintf(filePtr, "%s\t", myQual.c_str());
        numExpected += myQual.length() + 1;
    }
    const std::string& filterString = myFilter.getString();
    if(filterString.length() == 0)
    {
        numWritten += ifprintf(filePtr, ".\t");
        numExpected += 2;
    }
    else
    {
        numWritten += ifprintf(filePtr, "%s\t", filterString.c_str());
        numExpected += filterString.length() + 1;
    }

    // Write the info.
    bool writeSuccess = myInfo.write(filePtr);

    // Only write the format & genotype if we are not just writing siteOnly
    // data and there is at least one sample 
    if((!siteOnly) && (myGenotype.getNumSamples() != 0))
    {
        writeSuccess &= myGenotype.write(filePtr);
    }

    // Write the new line.
    numWritten += ifprintf(filePtr, "\n");
    numExpected += 1;

    return((numWritten == numExpected) && writeSuccess);
}



void VcfRecord::reset()
{
    myChrom.clear();
    my1BasedPosNum = 0;
    //my1BasedPos.clear();
    myID.clear();
    myRef.clear();
    myAlt.clear();
    myAltArray.clear();
    myQualNum = 0;;
    myQual.clear();
    myFilter.clear();
    myInfo.clear();
    myGenotype.clear();
    myAlleleCount.clear();
    myStatus = StatGenStatus::SUCCESS;
}


// Return the error after a failed call.
const StatGenStatus& VcfRecord::getStatus()
{
    return(myStatus);
}


const char* VcfRecord::getAlleles(unsigned int index)
{
    if(index == 0)
    {
        return(myRef.c_str());
    }
    if(index > getNumAlts())
    {
        // Index out of range.
        // Throw an exception.
        throw(std::runtime_error("VcfRecord::getAlleles called with an index that is greater than the number of alternates."));
        return(NULL);
    }
    // Alternate allele, so return the alternate.
    return(myAltArray.get(index-1).c_str());
}


int VcfRecord::getIntAllele(unsigned int index)
{
    const char* alleles = getAlleles(index);
    switch(alleles[0])
    {
        case 'A':
            return(1);
            break;
        case 'C':
            return(2);
            break;
        case 'G':
            return(3);
            break;
        case 'T':
            return(4);
            break;
        default:
            std::cerr << "VcfRecord::getIntAllele, unknown allele, " 
                      << alleles[0] << std::endl;
    }
    return(0);
}


unsigned int VcfRecord::getNumAlts()
{
    int numAlts = myAltArray.size();
    if(numAlts != 0)
    {
        // Already parsed so just return the number of alternates.
        return(numAlts);
    }
    // Check if it is just '.'.
    if((myAlt.length() == 1) && (myAlt == "."))
    {
        // No alternates.
        return(0);
    }

    // Parse the alternates by looping looking for commas.
    std::string* altStr = &(myAltArray.getNextEmpty());
    for(std::string::iterator iter = myAlt.begin(); iter != myAlt.end(); iter++)
    {
        if(*iter == ',')
        {
            altStr = &(myAltArray.getNextEmpty());
        }
        else
        {
            altStr->push_back(*iter);
        }
    }
    return(myAltArray.size());
}


int VcfRecord::getAlleleCount(unsigned int index,
                              VcfSubsetSamples* sampleSubset)
{
    unsigned int numAlts = getNumAlts();
    if(index > numAlts)
    {
        // Index out of range.
        // Throw an exception.
        throw(std::runtime_error("VcfRecord::getAlleles called with an index that is greater than the number of alternates."));
        return(-1);
    }

    if(myAlleleCount.size() == 0)
    {
        unsigned int gt = 0;
        myAlleleCount.resize(numAlts+1, 0);

        // Loop through the samples, counting the number of each allele.
        for(int sampleNum = 0; sampleNum < getNumSamples(); 
            sampleNum++)
        {
            if((sampleSubset != NULL) &&
               !(sampleSubset->keep(sampleNum)))
            {
                // Skip this sample.
                continue;
            }
            for(int gtNum = 0; gtNum < getNumGTs(sampleNum); gtNum++)
            {
                gt = getGT(sampleNum, gtNum);
                if((gt < 0) || (gt > numAlts))
                {
                    // Out of range GT, so continue to the next gt
                    continue;
                }
                // Increment the minor allele count
                ++myAlleleCount[gt];
            }
        }
    }

    // Alternate allele, so return the alternate.
    return(myAlleleCount[index]);
}


bool VcfRecord::readTilTab(IFILE filePtr, std::string& stringRef)
{
    int charRead = 0;
    while(1)
    {
        charRead = ifgetc(filePtr);
        
        if((charRead == '\n') || (charRead == EOF))
        {
            // Didn't find a tab, found a '\n' or eof
            // It still populated the string with values up
            // until the tab.
            return(false);
        }
        if(charRead == '\t')
        {
            // hit the tab character, so exit the loop.
            break;
        }
        stringRef += charRead;
    }
    return(true);
}
