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

#include "VcfFileReader.h"

VcfFileReader::VcfFileReader()
    : VcfFile(),
      myVcfIndex(NULL),
      myNewSection(false),
      mySectionChrom(""),
      mySection1BasedStartPos(-1),
      mySection1BasedEndPos(-1),
      mySectionOverlap(false),
      myRecordDiscardRules(),
      mySampleSubset(),
      myUseSubset(false),
      myMinAltAlleleCount(UNSET_MIN_ALT_ALLELE_COUNT),
      myAltAlleleCountSubset(NULL),
      myMinMinorAlleleCount(UNSET_MIN_MINOR_ALLELE_COUNT),
      myMinorAlleleCountSubset(NULL),
      myDiscardRules(0),
      myNumKeptRecords(0),
      myTotalRead(0)
{
  myFilePtr = NULL;
}


VcfFileReader::~VcfFileReader() 
{
    resetFile();
}


bool VcfFileReader::open(const char* filename, VcfHeader& header)
{
    // Close an already open file.
    close();

    myStatus = StatGenStatus::SUCCESS;
    if(VcfFile::open(filename, "r"))
    {
        // Successfully opened, so read the header.
        if(!header.read(myFilePtr))
        {
            // Failed, so copy the status.
            myStatus = header.getStatus();
            return(false);
        }
    }
    else
    {
        // Failed, status set by VcfFile::open.
        return(false);
    }

    // Successfully opened and read the header.
    return(true);
}


bool VcfFileReader::open(const char* filename, VcfHeader& header,
                         const char* includeFileName, const char* excludeSample,
                         const char* excludeFileName, const char* delims)
{
    if(!open(filename, header))
    {
        // Failed to open & read header, so return.
        return(false);
    }

    // Successfully opened and read the header, so setup the sample subset
    // object based on the specified sample files and the header.
    if(!mySampleSubset.init(header, includeFileName, excludeSample,
                            excludeFileName, delims))
    {
        // Failed to setup the subsetting.
        std::cerr << "VcfFileReader - failed to setup sample subsetting\n";
    }

    myUseSubset = true;

    // Successfully opened and read the header.
    return(true);
}


// Read VCF Index file.
bool VcfFileReader::readVcfIndex(const char* vcfIndexFilename)
{
    // Cleanup a previously setup index.
    if(myVcfIndex != NULL)
    {
        delete myVcfIndex;
        myVcfIndex = NULL;
    }

    // Create a new vcf index.
    myVcfIndex = new Tabix();
    StatGenStatus::Status indexStat = myVcfIndex->readIndex(vcfIndexFilename);

    if(indexStat != StatGenStatus::SUCCESS)
    {
        std::string errorMessage = "Failed to read the vcf Index file: ";
        errorMessage += vcfIndexFilename;
        myStatus.setStatus(indexStat, errorMessage.c_str());
        delete myVcfIndex;
        myVcfIndex = NULL;
        return(false);
    }

    if(myVcfIndex->getFormat() != Tabix::FORMAT_VCF)
    {
        std::string errorMessage = "ERROR: Tabix file not in VCF format: ";
        errorMessage += vcfIndexFilename;
        myStatus.setStatus(StatGenStatus::FAIL_PARSE, errorMessage.c_str());
        delete myVcfIndex;
        myVcfIndex = NULL;
        return(false);
    }

    myStatus = StatGenStatus::SUCCESS;
    return(true);
}


// Read VCF Index file.
bool VcfFileReader::readVcfIndex()
{
    if(myFilePtr == NULL)
    {
        // Can't read the vcf index file because the VCF file has not yet been
        // opened, so we don't know the base filename for the index file.
        std::string errorMessage = "Failed to read the vcf Index file -"
            " the VCF file needs to be read first in order to determine"
            " the index filename.";
        myStatus.setStatus(StatGenStatus::FAIL_ORDER, errorMessage.c_str());
        return(false);
    }

    const char* vcfBaseName = myFilePtr->getFileName();
    
    std::string indexName = vcfBaseName;
    indexName += ".tbi";

    bool foundFile = true;
    std::string failMessage = "";
    try
    {
        if(readVcfIndex(indexName.c_str()) == false)
        {
            foundFile = false;
        }
    }
    catch (std::exception& e)
    {
        foundFile = false;
        failMessage = e.what();
    }

    // Check to see if the index file was found.
    if(!foundFile)
    {
        // Not found - try without the vcf extension.
        // Locate the start of the vcf extension
        size_t startExt = indexName.find(".vcf");
        if(startExt == std::string::npos)
        {
            // Could not find the .vcf extension, so just return false since the
            // call to readVcfIndex set the status.
            return(false);
        }
        // Remove ".vcf" and try reading the index again.
        indexName.erase(startExt,  4);
        try
        {
            return(readVcfIndex(indexName.c_str()));
        }
        catch (std::exception& e)
        {
            failMessage += "\n";
            failMessage += e.what();
            throw(std::runtime_error(failMessage));
            return(false);
        }
    }
    return(true);
}


// return a pointer to the VCF Index file.
const Tabix* VcfFileReader::getVcfIndex()
{
    return(myVcfIndex);
}


bool VcfFileReader::readRecord(VcfRecord& record, VcfSubsetSamples* subset)
{
    myStatus = StatGenStatus::SUCCESS;
    // Subset the read if there are subsets specified.
    VcfSubsetSamples* subsetPtr = subset;
    if((subsetPtr == NULL) && myUseSubset)
    {
        subsetPtr = &mySampleSubset;
    }

    // Check to see if a new region has been set.  If so, setup for that region.
    bool searchChrom = false;
    if(myNewSection)
    {
        if(myVcfIndex != NULL)
        {
            // Have an index file so use
            if(!processNewSection())
            {
                // processNewSection sets the status appropriately on failure.
                return(false);
            }
        }
        else if(myTotalRead == 0)
        {
            // ReadSection without an index only works if no records
            // have been read.
            searchChrom = true;
            myNewSection = false;
        }
        else
        {
            myNewSection = false;
            myStatus.setStatus(StatGenStatus::FAIL_ORDER, 
                               "Cannot set read section with no index after reading records");
            return(false);
        }
    }

    // Keep looping until a desired record is found.
    bool recordFound = false;
    while(!recordFound)
    {
        if(!record.read(myFilePtr, mySiteOnly, myRecordDiscardRules, subsetPtr))
        {
            myStatus = record.getStatus();
            myTotalRead += myRecordDiscardRules.getNumDiscarded();
            myNumRecords += myRecordDiscardRules.getNumDiscarded();
            myRecordDiscardRules.clearNumDiscarded();
            return(false);
        }

        ++myTotalRead;
        myTotalRead += myRecordDiscardRules.getNumDiscarded();

        // Check to see if the record is in the section.
        // First check the chromosome.
        if(!mySectionChrom.empty() && (mySectionChrom != record.getChromStr()))
        {
            if(searchChrom)
            {
                // Still searching for the chromosome, so continue
                // to the next record.
                continue;
            }

            // Record is not within the correct chromosome, so return failure.
            myStatus = StatGenStatus::NO_MORE_RECS;
           return(false);
        }
        searchChrom = false;

        // Check if the record is after the section end if applicable.
        if((mySection1BasedEndPos != -1) && 
           (record.get1BasedPosition() >= mySection1BasedEndPos))
        {
            myStatus = StatGenStatus::NO_MORE_RECS;
            return(false);
        }
        
        // Check if the record is prior to the section start if applicable.
        // Determinine the VCF record end position.
        // If we are not requiring overlap, then we only need to check
        // the start position, but if overlap is required, then it needs
        // to incrment the start by the length-1.
        int numIncBases = 0;
        if(mySectionOverlap)
        {
            // The VCF record end position is the start position + length of the
            // reference string - 1.
            numIncBases = record.getNumRefBases() - 1;
        }
        if((mySection1BasedStartPos != -1) &&
           ((record.get1BasedPosition() + numIncBases)
            < mySection1BasedStartPos))
        {
            // This record is prior to the section, so keep reading.
            continue;
        }

        ++myNumRecords;
        myNumRecords += myRecordDiscardRules.getNumDiscarded();
        myRecordDiscardRules.clearNumDiscarded();
        
        // Record successfully read, so check to see if it is discarded.
        if((myDiscardRules & DISCARD_NON_PHASED) && !record.allPhased())
        {
            // Not all samples are phased, so discard this record.
            continue;
        }
        if((myDiscardRules & DISCARD_MISSING_GT) &&
           !record.hasAllGenotypeAlleles())
        {
            // discard missing GTs and this record had missing alleles,
            // so keep reading.
            continue;
        }
        if((myDiscardRules & DISCARD_FILTERED) && 
           !(record.getFilter().passedAllFilters()))
        {
            // Record was filtered, so discard it.
            continue;
        }
        if((myDiscardRules & DISCARD_MULTIPLE_ALTS) &&
           (record.getNumAlts() > 1))
        {
            // Record had multiple alternates, so discard.
            continue;
        }

        // Check allele counts for discarding.
        if(myMinAltAlleleCount != UNSET_MIN_ALT_ALLELE_COUNT)
        {
            // Count the number of alternates.
            int32_t altCount = 0;
            for(int sampleNum = 0; sampleNum < record.getNumSamples(); 
                sampleNum++)
            {
                if((myAltAlleleCountSubset != NULL) &&
                   !(myAltAlleleCountSubset->keep(sampleNum)))
                {
                    // Skip this sample.
                    continue;
                }
                for(int gtNum = 0; gtNum < record.getNumGTs(sampleNum); gtNum++)
                {
                    if(record.getGT(sampleNum, gtNum) > 0)
                    {
                        // Alternate, so increment the count.
                        ++altCount;
                    }
                }
            }
            if(altCount < myMinAltAlleleCount)
            {
                // Not enough alternates so continue to the next sample.
                continue;
            }
        }

        // Check to see if the minimum alternate allele count is met.
        if(myMinMinorAlleleCount != UNSET_MIN_MINOR_ALLELE_COUNT)
        {
            // Get the number of possible alternates.
            unsigned int numAlts = record.getNumAlts();

            // Verify that each allele has the min count.
            bool failMinorAlleleCount = false;
            for(unsigned int i = 0; i <= numAlts; i++)
            {
                if(record.getAlleleCount(i, myMinorAlleleCountSubset) 
                   < myMinMinorAlleleCount)
                {
                    // Not enough of one gt, so not ok.
                    failMinorAlleleCount = true;
                    break;
                }
            }
            if(failMinorAlleleCount)
            {
                // not enough alleles, so continue to the next record.
                continue;
            }
        }

        // Record was not discarded.
        recordFound = true;
    }

    // Increment the number of kept records.
    ++myNumKeptRecords;
    return(true);
}


bool VcfFileReader::setReadSection(const char* chromName)
{
    return(set1BasedReadSection(chromName, -1, -1));
}


bool VcfFileReader::set1BasedReadSection(const char* chromName, 
                                         int32_t start, int32_t end,
                                         bool overlap)
{
    myNewSection = true;
    mySectionChrom = chromName;
    mySection1BasedStartPos = start;
    mySection1BasedEndPos = end;
    mySectionOverlap = overlap;
    return(true);
}


// Returns whether or not the end of the file has been reached.
// return: int - true = EOF; false = not eof.
bool VcfFileReader::isEOF()
{
    if (myFilePtr != NULL)
    {
        // File Pointer is set, so return if eof.
        return(ifeof(myFilePtr));
    }
    // File pointer is not set, so return true, eof.
    return true;
}


bool VcfFileReader::setExcludeIDs(const char* filename)
{
    return(myRecordDiscardRules.setExcludeIDs(filename));
}


bool VcfFileReader::setIncludeIDs(const char* filename)
{
    return(myRecordDiscardRules.setIncludeIDs(filename));
}


void VcfFileReader::addDiscardMinAltAlleleCount(int32_t minAltAlleleCount, 
                                                VcfSubsetSamples* subset)
{
    myMinAltAlleleCount = minAltAlleleCount;
    myAltAlleleCountSubset = subset;
}


void VcfFileReader::rmDiscardMinAltAlleleCount()
{
    myMinAltAlleleCount = UNSET_MIN_ALT_ALLELE_COUNT;
    myAltAlleleCountSubset = NULL;
}


void VcfFileReader::addDiscardMinMinorAlleleCount(int32_t minMinorAlleleCount, 
                                                  VcfSubsetSamples* subset)
{
    myMinMinorAlleleCount = minMinorAlleleCount;
    myMinorAlleleCountSubset = subset;
}


void VcfFileReader::rmDiscardMinMinorAlleleCount()
{
    myMinMinorAlleleCount = UNSET_MIN_ALT_ALLELE_COUNT;
    myMinorAlleleCountSubset = NULL;
}


void VcfFileReader::resetFile()
{
    myRecordDiscardRules.reset(),
    mySampleSubset.reset();
    myUseSubset = false;
    myNumKeptRecords = 0;
    myTotalRead = 0;
    myNewSection = false;
    mySectionChrom = "";
    mySection1BasedStartPos = -1;
    mySection1BasedEndPos = -1;
    mySectionOverlap = false;

    if(myVcfIndex != NULL)
    {
        delete myVcfIndex;
        myVcfIndex = NULL;
    }
}


bool VcfFileReader::processNewSection()
{
    myNewSection = false;
    
    // Check to see if the index file has been read.
    if(myVcfIndex == NULL)
    {
        myStatus.setStatus(StatGenStatus::FAIL_ORDER, 
                           "Cannot read section since there is no index file open");
        throw(std::runtime_error("SOFTWARE BUG: trying to read a VCF record by section prior to opening the VCF Index file."));
        return(false);
    }

    if(myFilePtr == NULL)
    {
        myStatus.setStatus(StatGenStatus::FAIL_ORDER, 
                           "Cannot read section without first opening the VCF file.");
        throw(std::runtime_error("SOFTWARE BUG: trying to read a VCF record by section prior to opening the VCF file."));
        return(false);
    }

    // Using random access, so can't buffer
    myFilePtr->disableBuffering();

    uint64_t startPos = 0;
    // Find where this section starts in the file.
    if(!myVcfIndex->getStartPos(mySectionChrom.c_str(),
                                mySection1BasedStartPos, 
                                startPos))
    {
        // Didn't find the position.
        myStatus = StatGenStatus::NO_MORE_RECS;
        return(false);
    }
    if(startPos != (uint64_t)iftell(myFilePtr))
    {
        // Seek to the start position.
        if(ifseek(myFilePtr, startPos, SEEK_SET) != true)
        {
            // seek failed, return failure.
            myStatus.setStatus(StatGenStatus::FAIL_IO, 
                               "Failed to seek to the specified section");
            return(false);
        }
    }
    return(true);
}
