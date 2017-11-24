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

#ifndef __CIGAR_HELPER_H__
#define __CIGAR_HELPER_H__

#include "SamRecord.h"

/// Class for helping to filter a SAM/BAM record.
class CigarHelper
{
public:
    static const int32_t NO_CLIP = -1;

    /// Soft clip the cigar from the beginning of the read at the specified
    /// reference position.  If the clip position is deleted/skipped
    /// or is immediately followed by a deletion/skip/pad/insert, that entire 
    /// CIGAR operation is also removed.
    /// Nothing is clipped if the reference position is before the read starts,
    /// everything is clipped if the reference position is after the read ends.
    /// \param record record to calculate the clip for.
    /// \param refPosition0Based 0-based reference position to end the clip at
    /// \param newCigar cigar object to set with the updated cigar.
    /// \param new0BasedPosition new 0-based reference position of the read.
    /// \param read position where the clip ends (last clipped position) or 
    //         NO_CLIP if nothing is clipped.
    static int32_t softClipBeginByRefPos(SamRecord& record, 
                                         int32_t refPosition0Based,
                                         CigarRoller& newCigar,
                                         int32_t &new0BasedPosition);

    /// Soft clip the cigar from the back of the read at the specified
    /// reference position.  If the clip position is deleted/skipped
    /// or is immediately preceded by a deletion/skip/pad, that entire CIGAR
    /// operation is also removed.  If the clip position is immediately
    /// preceded by an insertion, the insertion is left in the CIGAR.
    /// Nothing is clipped if the reference position is after the read ends,
    /// everything is clipped if the reference position is before the read
    /// starts (including insertions).
    /// \param record record to calculate the clip for.
    /// \param refPosition0Based 0-based reference position to start clip at
    /// \param newCigar cigar object to set with the updated cigar.
    /// \param read position where the clip starts or 
    //         NO_CLIP if nothing is clipped.
    static int32_t softClipEndByRefPos(SamRecord& record, 
                                       int32_t refPosition0Based,
                                       CigarRoller& newCigar);
};

#endif

