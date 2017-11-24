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

//////////////////////////////////////////////////////////////////////////


#include "SamFilter.h"

#include "SamQuerySeqWithRefHelper.h"
#include "BaseUtilities.h"
#include "SamFlag.h"

SamFilter::FilterStatus SamFilter::clipOnMismatchThreshold(SamRecord& record, 
                                                           GenomeSequence& refSequence,
                                                           double mismatchThreshold)
{
    // Read & clip from the left & right.    
    SamQuerySeqWithRefIter iterFromFront(record, refSequence, true);
    SamQuerySeqWithRefIter iterFromBack(record, refSequence, false);

    SamSingleBaseMatchInfo baseMatchInfo;

    int32_t readLength = record.getReadLength();
    // Init last front clip to be prior to the lastFront index (0).
    const int32_t initialLastFrontClipPos = -1;
    int32_t lastFrontClipPos = initialLastFrontClipPos;
    // Init first back clip to be past the last index (readLength).
    int32_t firstBackClipPos = readLength;

    bool fromFrontComplete = false;
    bool fromBackComplete = false;
    int32_t numBasesFromFront = 0;
    int32_t numBasesFromBack = 0;
    int32_t numMismatchFromFront = 0;
    int32_t numMismatchFromBack = 0;

    //////////////////////////////////////////////////////////
    // Determining the clip positions.
    while(!fromFrontComplete || !fromBackComplete)
    {
        // Read from the front (left to right) of the read until
        // more have been read from that direction than the opposite direction.
        while(!fromFrontComplete && 
              ((numBasesFromFront <= numBasesFromBack) ||
               (fromBackComplete)))
        {
            if(iterFromFront.getNextMatchMismatch(baseMatchInfo) == false)
            {
                // Nothing more to read in this direction.
                fromFrontComplete = true;
                break;
            }
            // Got a read.  Check to see if it is to or past the last clip.
            if(baseMatchInfo.getQueryIndex() >= firstBackClipPos)
            {
                // This base is past where we are clipping, so we
                // are done reading in this direction.
                fromFrontComplete = true;
                break;
            }
            // This is an actual base read from the left to the
            // right, so up the counter and determine if it was a mismatch.
            ++numBasesFromFront;

            if(baseMatchInfo.getType() == SamSingleBaseMatchInfo::MISMATCH)
            {
                // Mismatch
                ++numMismatchFromFront;
                // Check to see if it is over the threshold.
                double mismatchPercent = 
                    (double)numMismatchFromFront / numBasesFromFront;
                if(mismatchPercent > mismatchThreshold)
                {
                    // Need to clip.
                    lastFrontClipPos = baseMatchInfo.getQueryIndex();
                    // Reset the counters.
                    numBasesFromFront = 0;
                    numMismatchFromFront = 0;
                }
            }
        }

        // Now, read from right to left until more have been read
        // from the back than from the front.
        while(!fromBackComplete && 
              ((numBasesFromBack <= numBasesFromFront) ||
               (fromFrontComplete)))
        {
            if(iterFromBack.getNextMatchMismatch(baseMatchInfo) == false)
            {
                // Nothing more to read in this direction.
                fromBackComplete = true;
                break;
            }
            // Got a read.  Check to see if it is to or past the first clip.
            if(baseMatchInfo.getQueryIndex() <= lastFrontClipPos)
            {
                // This base is past where we are clipping, so we
                // are done reading in this direction.
                fromBackComplete = true;
                break;
            }
            // This is an actual base read from the right to the
            // left, so up the counter and determine if it was a mismatch.
            ++numBasesFromBack;

            if(baseMatchInfo.getType() == SamSingleBaseMatchInfo::MISMATCH)
            {
                // Mismatch
                ++numMismatchFromBack;
                // Check to see if it is over the threshold.
                double mismatchPercent = 
                    (double)numMismatchFromBack / numBasesFromBack;
                if(mismatchPercent > mismatchThreshold)
                {
                    // Need to clip.
                    firstBackClipPos = baseMatchInfo.getQueryIndex();
                    // Reset the counters.
                    numBasesFromBack = 0;
                    numMismatchFromBack = 0;
                }
            }
        }
    }

    //////////////////////////////////////////////////////////
    // Done determining the clip positions, so clip.
    // To determine the number of clips from the front, add 1 to the
    // lastFrontClipPos since the index starts at 0.
    // To determine the number of clips from the back, subtract the
    // firstBackClipPos from the readLength.
    // Example:
    // Pos:  012345
    // Read: AAAAAA
    // Read Length = 6.  If lastFrontClipPos = 2 and firstBackClipPos = 4, numFrontClips = 3 & numBack = 2.
    return(softClip(record, lastFrontClipPos + 1, readLength - firstBackClipPos));
}


// Soft clip the record from the front and/or the back.
SamFilter::FilterStatus SamFilter::softClip(SamRecord& record, 
                                            int32_t numFrontClips,
                                            int32_t numBackClips)
{
    //////////////////////////////////////////////////////////
    Cigar* cigar = record.getCigarInfo();
    FilterStatus status = NONE;
    int32_t startPos = record.get0BasedPosition();
    CigarRoller updatedCigar;

    status = softClip(*cigar, numFrontClips, numBackClips, 
                      startPos, updatedCigar);

    if(status == FILTERED)
    {
        /////////////////////////////
        // The entire read is clipped, so rather than clipping it,
        // filter it out.
        filterRead(record);
        return(FILTERED);
    }
    else if(status == CLIPPED)
    {
        // Part of the read was clipped, and now that we have
        // an updated cigar, update the read.
        record.setCigar(updatedCigar);

        // Update the starting position.
        record.set0BasedPosition(startPos);
    }
    return(status);
}


// Soft clip the cigar from the front and/or the back, writing the value
// into the new cigar.
SamFilter::FilterStatus SamFilter::softClip(Cigar& oldCigar, 
                                            int32_t numFrontClips,
                                            int32_t numBackClips,
                                            int32_t& startPos,
                                            CigarRoller& updatedCigar)
{
    int32_t readLength = oldCigar.getExpectedQueryBaseCount();
    int32_t endClipPos = readLength - numBackClips;
    FilterStatus status = NONE;

    if((numFrontClips != 0) || (numBackClips != 0))
    {
        // Clipping from front and/or from the back.

        // Check to see if the entire read was clipped.
        int32_t totalClips = numFrontClips + numBackClips;
        if(totalClips >= readLength)
        {
            /////////////////////////////
            // The entire read is clipped, so rather than clipping it,
            // filter it out.
            return(FILTERED);
        }
         
        // Part of the read was clipped.
        status = CLIPPED;
            
        // Loop through, creating an updated cigar.
        int origCigarOpIndex = 0;
        
        // Track how many read positions are covered up to this
        // point by the cigar to determine up to up to what
        // point in the cigar is affected by this clipping.
        int32_t numPositions = 0;
        
        // Track if any non-clips are in the new cigar.
        bool onlyClips = true;

        const Cigar::CigarOperator* op = NULL;

        //////////////////
        // Clip from front
        while((origCigarOpIndex < oldCigar.size()) &&
              (numPositions < numFrontClips))
        {
            op = &(oldCigar.getOperator(origCigarOpIndex));
            switch(op->operation)
            {
                case Cigar::hardClip:
                    // Keep this operation as the new clips do not
                    // affect other clips.
                    updatedCigar += *op;
                    break;
                case Cigar::del:
                case Cigar::skip:
                    // Skip and delete are going to be dropped, and
                    // are not in the read, so the read index doesn't
                    // need to be updated
                    break;
                case Cigar::insert:
                case Cigar::match:
                case Cigar::mismatch:
                case Cigar::softClip:
                    // Update the read index as these types
                    // are found in the read.
                    numPositions += op->count;
                    break;
                case Cigar::none:
                default:
                    // Nothing to do for none.
                    break;
            };
            ++origCigarOpIndex;
        }
    
        // If bases were clipped from the front, add the clip and
        // any partial cigar operation as necessary.
        if(numFrontClips != 0)
        {
            // Add the softclip to the front of the read.
            updatedCigar.Add(Cigar::softClip, numFrontClips);
        
            // Add the rest of the last Cigar operation if
            // it is not entirely clipped.
            int32_t newCount = numPositions - numFrontClips;
            if(newCount > 0)
            {
                // Before adding it, check to see if the same
                // operation is clipped from the end.
                // numPositions greater than the endClipPos
                // means that it is equal or past that position,
                // so shorten the number of positions.
                if(numPositions > endClipPos)
                {
                    newCount -= (numPositions - endClipPos);
                }
                if(newCount > 0)
                {
                    updatedCigar.Add(op->operation, newCount);
                    if(!Cigar::isClip(op->operation))
                    {
                        onlyClips = false;
                    }
                }
            }
        }
    
        // Add operations until the point of the end clip is reached.
        // For example...
        //   2M1D3M = MMDMMM  readLength = 5
        // readIndex: 01 234
        //   at cigarOpIndex 0 (2M), numPositions = 2.
        //   at cigarOpIndex 1 (1D), numPositions = 2.
        //   at cigarOpIndex 2 (3M), numPositions = 5.
        // if endClipPos = 2, we still want to consume the 1D, so
        // need to keep looping until numPositions > endClipPos
        while((origCigarOpIndex < oldCigar.size()) &&
              (numPositions <= endClipPos))
        {
            op = &(oldCigar.getOperator(origCigarOpIndex));
            
            // Update the numPositions count if the operations indicates
            // bases within the read.
            if(!Cigar::foundInQuery(op->operation))
            {
                // This operation is not in the query read sequence,
                // so it is not yet to the endClipPos, just add the
                // operation do not increment the number of positions.
                updatedCigar += *op;
                if(!Cigar::isClip(op->operation))
                {
                    onlyClips = false;
                }
            }
            else
            {
                // This operation appears in the query sequence, so
                // check to see if the clip occurs in this operation.
                
                // endClipPos is 0 based & numPositions is a count.
                // If endClipPos is 4, then it is the 5th position.
                // If 4 positions are covered so far (numPositions = 4), 
                // then we are right at endCLipPos: 4-4 = 0, none of 
                // this operation should be kept. 
                // If only 3 positions were covered, then we are at offset
                // 3, so offset 3 should be added: 4-3 = 1.
                uint32_t numPosTilClip = endClipPos - numPositions;
                
                if(numPosTilClip < op->count)
                {
                    // this operation is partially clipped, write the part
                    // that was not clipped if it is not all clipped.
                    if(numPosTilClip != 0)
                    {
                        updatedCigar.Add(op->operation,
                                     numPosTilClip);
                        if(!Cigar::isClip(op->operation))
                        {
                            onlyClips = false;
                        }
                    }
                }
                else
                {
                    // This operation is not clipped, so add it
                    updatedCigar += *op;
                    if(!Cigar::isClip(op->operation))
                    {
                        onlyClips = false;
                    }
                }
                // This operation occurs in the query sequence, so 
                // increment the number of positions covered.
                numPositions += op->count;
            }

            // Move to the next cigar position.
            ++origCigarOpIndex;
        }
            
        //////////////////
        // Add the softclip to the back.
        if(numBackClips != 0)
        {
            // Add the softclip to the end
            updatedCigar.Add(Cigar::softClip, numBackClips);
        }
        
        //////////////////
        // Add any hardclips remaining in the original cigar to the back.
        while(origCigarOpIndex < oldCigar.size())
        {
            op = &(oldCigar.getOperator(origCigarOpIndex));
            if(op->operation == Cigar::hardClip)
            {
                // Keep this operation as the new clips do not
                // affect other clips.
                updatedCigar += *op;
            }
            ++origCigarOpIndex;
        }
        
        // Check to see if the new cigar is only clips.
        if(onlyClips)
        {
            // Only clips in the new cigar, so mark the read as filtered
            // instead of updating the cigar.
            /////////////////////////////
            // The entire read was clipped.
            status = FILTERED;
        }
        else
        {
            // Part of the read was clipped.
            // Update the starting position if a clip was added to
            // the front.
            if(numFrontClips > 0)
            {
                // Convert from query index to reference position (from the
                // old cigar)
                // Get the position for the last front clipped position by
                // getting the position associated with the clipped base on
                // the reference.  Then add one to get to the first
                // non-clipped position.
                int32_t lastFrontClipPos = numFrontClips - 1;
                int32_t newStartPos = oldCigar.getRefPosition(lastFrontClipPos, 
                                                              startPos);
                if(newStartPos != Cigar::INDEX_NA)
                {
                    // Add one to get first non-clipped position.
                    startPos = newStartPos + 1;
                }
            }
        }
    }
    return(status);
}


SamFilter::FilterStatus SamFilter::filterOnMismatchQuality(SamRecord& record, 
                                                           GenomeSequence& refSequence,
                                                           uint32_t qualityThreshold, 
                                                           uint8_t defaultQualityInt)
{
    uint32_t totalMismatchQuality = 
        sumMismatchQuality(record, refSequence, defaultQualityInt); 
    
    // If the total mismatch quality is over the threshold, 
    // filter the read.
    if(totalMismatchQuality > qualityThreshold)
    {
        filterRead(record);
        return(FILTERED);
    }
    return(NONE);
}


// NOTE: Only positions where the reference and read both have bases that
//       are different and not 'N' are considered mismatches.
uint32_t SamFilter::sumMismatchQuality(SamRecord& record, 
                                       GenomeSequence& refSequence,
                                       uint8_t defaultQualityInt)
{
    // Track the mismatch info.
    int mismatchQual = 0;
    int numMismatch = 0;

    SamQuerySeqWithRefIter sequenceIter(record, refSequence);

    SamSingleBaseMatchInfo baseMatchInfo;
    while(sequenceIter.getNextMatchMismatch(baseMatchInfo))
    {
        if(baseMatchInfo.getType() == SamSingleBaseMatchInfo::MISMATCH)
        {
            // Got a mismatch, get the associated quality.
            char readQualityChar = 
                record.getQuality(baseMatchInfo.getQueryIndex());
            uint8_t readQualityInt = 
                BaseUtilities::getPhredBaseQuality(readQualityChar);
            
            if(readQualityInt == BaseUtilities::UNKNOWN_QUALITY_INT)
            {
                // Quality was not specified, so use the configured setting.
                readQualityInt = defaultQualityInt;
            }
            mismatchQual += readQualityInt;
            ++numMismatch;
        }
    }

    return(mismatchQual);
}


void SamFilter::filterRead(SamRecord& record)
{
    // Filter the read by marking it as unmapped.
    uint16_t flag = record.getFlag(); 
    SamFlag::setUnmapped(flag);
    // Clear N/A flags.
    flag &= ~SamFlag::PROPER_PAIR;
    flag &= ~SamFlag::SECONDARY_ALIGNMENT;
    flag &= ~SamFlag::SUPPLEMENTARY_ALIGNMENT;
    record.setFlag(flag);
    // Clear Cigar
    record.setCigar("*");
    // Clear mapping quality
    record.setMapQuality(0);
}
