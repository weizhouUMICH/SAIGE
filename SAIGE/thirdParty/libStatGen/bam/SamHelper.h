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

#ifndef __SAM_HELPER_H__
#define __SAM_HELPER_H__

#include <stdint.h>

#ifdef DUPLICATE
#undef DUPLICATE
#endif

/// Class for extracting information from a SAM Flag.
class SamHelper
{
public:

    /// Helper method that combines the chromosome ID and position into a 
    /// 64bit number by shifting the chromosome ID to the upper bits.
    static inline uint64_t combineChromPos(int32_t chromID, int32_t position)
    {
        return(((uint64_t)chromID << 32) | (position & 0xFFFFFFFF));
    }

private:
    SamHelper();
};


#endif
