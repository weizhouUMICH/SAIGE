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

//////////////////////////////////////////////////////////////////////////

#include "CigarHelper.h"

// Soft Clip from the beginning of the read to the specified reference position.
int32_t CigarHelper::softClipBeginByRefPos(SamRecord& record, 
                                           int32_t refPosition0Based,
                                           CigarRoller& newCigar,
                                           int32_t &new0BasedPosition)
{
    newCigar.clear();
    Cigar* cigar = record.getCigarInfo();
    if(cigar == NULL)
    {
        // Failed to get the cigar.
        ErrorHandler::handleError("Soft clipping, but failed to read the cigar");
        return(NO_CLIP);
    }

    // No cigar or position in the record, so return no clip.
    if((cigar->size() == 0) || (record.get0BasedPosition() == -1))
    {
        return(NO_CLIP);
    }

    // Check to see if the reference position occurs before the record starts,
    // if it does, do no clipping.
    if(refPosition0Based < record.get0BasedPosition())
    {
        // Not within this read, so nothing to clip.
        newCigar.Set(record.getCigar());
        return(NO_CLIP);
    }

    // The position falls after the read starts, so loop through until the
    // position or the end of the read is found.
    int32_t readClipPosition = 0;
    bool clipWritten = false;
    new0BasedPosition = record.get0BasedPosition();
    for(int i = 0; i < cigar->size(); i++)
    {
        const Cigar::CigarOperator* op = &(cigar->getOperator(i));

        if(clipWritten)
        {
            // Clip point has been found, so just add everything.
            newCigar += *op;
            // Go to the next operation.
            continue;
        }

        // The clip point has not yet been found, so check to see if we found 
        // it now.

        // Not a clip, check to see if the operation is found in the
        // reference.
        if(Cigar::foundInReference(*op))
        {
            // match, mismatch, deletion, skip

            // increment the current reference position to just past this
            // operation.
            new0BasedPosition += op->count;

            // Check to see if this is also in the query, because otherwise
            // the operation is still being consumed.
            if(Cigar::foundInQuery(*op))
            {
                // Also in the query, determine if the entire thing should
                // be clipped or just part of it.

                uint32_t numKeep = 0;
                // Check to see if we have hit our clip position.
                if(refPosition0Based < new0BasedPosition)
                {
                    // The specified clip position is in this cigar operation.
                    numKeep = new0BasedPosition - refPosition0Based - 1;
                    
                    if(numKeep > op->count)
                    {
                        // Keep the entire read.  This happens because
                        // we keep reading until the first match/mismatch
                        // after the clip.
                        numKeep = op->count;
                    }
                }

                // Add the part of this operation that is being clipped
                // to the clip count.
                readClipPosition += (op->count - numKeep);

                // Only write the clip if we found a match/mismatch
                // to write.  Otherwise we will keep accumulating clips
                // for the case of insertions.
                if(numKeep > 0)
                {
                    new0BasedPosition -= numKeep;

                    newCigar.Add(Cigar::softClip, readClipPosition);
                    
                    // Add the clipped part of this cigar to the clip
                    // position.
                    newCigar.Add(op->operation, numKeep);
                    
                    // Found a match after the clip point, so stop
                    // consuming cigar operations.
                    clipWritten = true;
                    continue;
                }
            }
        }
        else
        {
            // Only add hard clips.  The softclips will be added in
            // when the total number is found.
            if(op->operation == Cigar::hardClip)
            {
                // Check if this is the first operation, if so, just write it.
                if(i == 0)
                {
                    newCigar += *op;
                }
                // Check if it is the last operation (otherwise skip it).
                else if(i == (cigar->size() - 1))
                {
                    // Check whether or not the clip was ever written, and if
                    // not, write it.
                    if(clipWritten == false)
                    {
                        newCigar.Add(Cigar::softClip, readClipPosition);
                        // Since no match/mismatch was ever found, set
                        // the new ref position to the original one.
                        new0BasedPosition = record.get0BasedPosition();
                        clipWritten = true;
                    }
                    // Add the hard clip.
                    newCigar += *op;
                }
            }
            // Not yet to the clip position, so do not add this operation.
            if(Cigar::foundInQuery(*op))
            {
                // Found in the query, so update the read clip position.
                readClipPosition += op->count;
            }
        }
    } // End loop through cigar.


    // Check whether or not the clip was ever written, and if
    // not, write it.
    if(clipWritten == false)
    {
        newCigar.Add(Cigar::softClip, readClipPosition);
        // Since no match/mismatch was ever found, set
        // the new ref position to the original one.
        new0BasedPosition = record.get0BasedPosition();
    }

    // Subtract 1 since readClipPosition atually contains the first 0based 
    // position that is not clipped.
    return(readClipPosition - 1);
}


// Soft Clip from the end of the read at the specified reference position.
int32_t CigarHelper::softClipEndByRefPos(SamRecord& record, 
                                         int32_t refPosition0Based,
                                         CigarRoller& newCigar)
{
    newCigar.clear();
    Cigar* cigar = record.getCigarInfo();
    if(cigar == NULL)
    {
        // Failed to get the cigar.
        ErrorHandler::handleError("Soft clipping, but failed to read the cigar");
        return(NO_CLIP);
    }

    // No cigar or position in the record, so return no clip.
    if((cigar->size() == 0) || (record.get0BasedPosition() == -1))
    {
        return(NO_CLIP);
    }

    // Check to see if the reference position occurs after the record ends,
    // if so, do no clipping.
    if(refPosition0Based > record.get0BasedAlignmentEnd())
    {
        // Not within this read, so nothing to clip.
        newCigar.Set(record.getCigar());
        return(NO_CLIP);
    }

    // The position falls before the read ends, so loop through until the
    // position is found.
    int32_t currentRefPosition = record.get0BasedPosition();
    int32_t readClipPosition = 0;
    for(int i = 0; i < cigar->size(); i++)
    {
        const Cigar::CigarOperator* op = &(cigar->getOperator(i));

        // If the operation is found in the reference, increase the
        // reference position.
        if(Cigar::foundInReference(*op))
        {
            // match, mismatch, deletion, skip
            // increment the current reference position to just past
            // this operation.
            currentRefPosition += op->count;
        }
         
        // Check to see if we have hit our clip position.
        if(refPosition0Based < currentRefPosition)
        {
            // If this read is also in the query (match/mismatch), 
            // write the partial op to the new cigar.
            int32_t numKeep = 0;
            if(Cigar::foundInQuery(*op))
            {
                numKeep = op->count - (currentRefPosition - refPosition0Based);
                if(numKeep > 0)
                {
                    newCigar.Add(op->operation, numKeep);
                    readClipPosition += numKeep;
                }
            }
            else if(Cigar::isClip(*op))
            {
                // This is a hard clip, so write it.
                newCigar.Add(op->operation, op->count);
            }
            else
            {

                // Not found in the query (skip/deletion),
                // so don't write any of the operation.
            }
            // Found the clip point, so break.
            break;
        }
        else if(refPosition0Based == currentRefPosition)
        {
            newCigar += *op;
            if(Cigar::foundInQuery(*op))
            {
                readClipPosition += op->count;
            }
        }
        else
        {
            // Not yet to the clip position, so add this operation/size to
            // the new cigar.
            newCigar += *op;
            if(Cigar::foundInQuery(*op))
            {
                // Found in the query, so update the read clip position.
                readClipPosition += op->count;
            }
        }
    } // End loop through cigar.

    // Before adding the softclip, read from the end of the cigar checking to
    // see if the operations are in the query, removing operations that are
    // not (pad/delete/skip) until a hardclip or an operation in the query is
    // found. We do not want a pad/delete/skip right before a softclip.
    for(int j = newCigar.size() - 1; j >= 0; j--)
    {
        const Cigar::CigarOperator* op = &(newCigar.getOperator(j));
        if(!Cigar::foundInQuery(*op) && !Cigar::isClip(*op))
        {
            // pad/delete/skip
            newCigar.Remove(j);
        }
        else if(Cigar::foundInQuery(*op) & Cigar::isClip(*op))
        {
            // Soft clip, so increment the clip position for the return value.
            // Remove the softclip since the readClipPosition is used to
            // calculate teh size of the soft clip added.
            readClipPosition -= op->count;
            newCigar.Remove(j);
        }
        else
        {
            // Found a cigar operation that should not be deleted, so stop deleting.
            break;
        }
    } 

    // Determine the number of soft clips.
    int32_t numSoftClips = record.getReadLength() - readClipPosition;
    // NOTE that if the previous operation is a softclip, the CigarRoller logic
    // will merge this with that one.
    newCigar.Add(Cigar::softClip, numSoftClips);

    // Check if an ending hard clip needs to be added.
    if(cigar->size() != 0)
    {
        const Cigar::CigarOperator* lastOp = 
            &(cigar->getOperator(cigar->size() - 1));
        if(lastOp->operation == Cigar::hardClip)
        {
            newCigar += *lastOp;
        }
    }

    return(readClipPosition);
}
