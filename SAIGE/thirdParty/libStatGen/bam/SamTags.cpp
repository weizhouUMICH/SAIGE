/*
 *  Copyright (C) 2011  Regents of the University of Michigan
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

#include "SamTags.h"
#include "BaseUtilities.h"

const char* SamTags::BQ_TAG = "BQ";
const char SamTags::BQ_TAG_TYPE = 'Z';
const char* SamTags::MD_TAG = "MD";
const char SamTags::MD_TAG_TYPE = 'Z';
const char* SamTags::ORIG_POS_TAG = "OP";
const char SamTags::ORIG_POS_TAG_TYPE = 'i';
const char* SamTags::ORIG_CIGAR_TAG = "OC";
const char SamTags::ORIG_CIGAR_TAG_TYPE = 'Z';
const char* SamTags::ORIG_QUAL_TAG = "OQ";
const char SamTags::ORIG_QUAL_TAG_TYPE = 'Z';


// Create the MD tag for the specified input record and the genome.
bool SamTags::createMDTag(String& outputMDtag, SamRecord& inputRec, 
                          GenomeSequence& genome)
{
    outputMDtag.Clear();
    // Get the cigar to use for determing alignment.
    Cigar* cigarInfo = inputRec.getCigarInfo();
    if(cigarInfo == NULL)
    {
        throw(std::runtime_error("Cannot createMDTag - failed to get the cigar"));
        return(false);
    }
    int32_t queryIndex = Cigar::INDEX_NA;

    // get where this read starts on the reference.
    uint32_t startOfReadOnRefIndex = 
        genome.getGenomePosition(inputRec.getReferenceName());
    if(startOfReadOnRefIndex == (uint32_t)INVALID_CHROMOSOME_INDEX)
    {
        // Failed to find the reference for this chromosome, so return false.
        return(false);
    }
    startOfReadOnRefIndex += inputRec.get0BasedPosition();

    // Track the number of consecutive matches.
    int32_t matchCount = 0;
    // Track if it is currently in a deletion so it knows when not to add
    // a '^'.
    bool currentDeletion = false;

    // Loop through the Reference indices (ignores insertions/pads/clips).
    for(int refOffset = 0; 
        refOffset < cigarInfo->getExpectedReferenceBaseCount();
        ++refOffset)
    {
        // Get the query index for this reference position..
        queryIndex = cigarInfo->getQueryIndex(refOffset);

        char refBase = genome[startOfReadOnRefIndex + refOffset];

        if(queryIndex != Cigar::INDEX_NA)
        {
            // Both the reference and the read have a base, so get the bases.
            char readBase = inputRec.getSequence(queryIndex);
            currentDeletion = false;

            // If neither base is unknown and they are the same, count it
            // as a match.
            if(!BaseUtilities::isAmbiguous(readBase) && 
               !BaseUtilities::isAmbiguous(refBase) && 
               (BaseUtilities::areEqual(readBase, refBase)))
            {
                // Match, so update counter.
                ++matchCount;
            }
            else
            {
                // Mismatch, so output the number of matches if any.
                if(matchCount != 0)
                {
                    outputMDtag += matchCount;
                    matchCount = 0;
                }
                outputMDtag += refBase;
            }
        }
        else
        {
            // This reference position is not in the query, so it is a deletion.
            // Deletion, so output the number of matches if any.
            if(matchCount != 0)
            {
                outputMDtag += matchCount;
                matchCount = 0;
            }

            if(!currentDeletion)
            {
                // Not currently in a deletion, so add the ^
                outputMDtag += '^';
            }
            // Add the deleted base.
            outputMDtag += refBase;
            currentDeletion = true;
        }
    }

    // output the match count at the end.
    outputMDtag += matchCount;
    return(true);
}

// Check to see if the MD tag in the record is accurate.
bool SamTags::isMDTagCorrect(SamRecord& inputRec, GenomeSequence& genome)
{
    String calcMDtag;
    if(!createMDTag(calcMDtag, inputRec, genome))
    {
        // Could not generate the MD tag, so just return that it is incorrect.
        return(false);
    }
    
    const String* origMDtag = inputRec.getStringTag(MD_TAG);
    
    if(origMDtag == NULL)
    {
        // There was no tag.
        // if there is not a new tag, then they are the same and true
        // should be returned.  If there is a new tag, then the old one was
        // wrong so false should be returned.  So just return the result of
        // IsEmpty.
        return(calcMDtag.IsEmpty());
    }
    else
    {
        // origMDtag is not NULL, so just compare the two tags.
        return(calcMDtag == *origMDtag);
    }
}


// Update/Add the MD tag in the inputRec.
bool SamTags::updateMDTag(SamRecord& inputRec, GenomeSequence& genome)
{
    // Calculate the new MD tag.
    String calcMDtag;
    createMDTag(calcMDtag, inputRec, genome);
    
    // Add the MD tag.  If it is already there and is different it will
    // replace it.  If it is already there and it is the same, it won't
    // do anything.
    return(inputRec.addTag(MD_TAG, MD_TAG_TYPE, calcMDtag.c_str()));
}
