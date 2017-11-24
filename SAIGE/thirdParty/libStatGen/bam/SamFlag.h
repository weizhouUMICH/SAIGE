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

#ifndef __SAM_FLAG_H__
#define __SAM_FLAG_H__

#include <stdint.h>

#ifdef DUPLICATE
#undef DUPLICATE
#endif

/// Class for extracting information from a SAM Flag.
class SamFlag
{
public:
    ///////////////////////
    /// @name  Constants for parsing a flag.
    //@{
    static const int16_t PAIRED              = 0x0001;
    static const int16_t PROPER_PAIR         = 0x0002;
    static const int16_t UNMAPPED            = 0x0004;
    static const int16_t MATE_UNMAPPED       = 0x0008;
    static const int16_t REVERSE             = 0x0010;
    static const int16_t MATE_REVERSED       = 0x0020;
    static const int16_t FIRST_READ          = 0x0040;
    static const int16_t SECOND_READ         = 0x0080;
    static const int16_t SECONDARY_ALIGNMENT = 0x0100;
    static const int16_t FAILED_QUALITY      = 0x0200;
    static const int16_t DUPLICATE           = 0x0400;
    static const int16_t SUPPLEMENTARY_ALIGNMENT = 0x0800;
    static const int16_t FRAGMENT_INFO       = 0x00C0;
    static const int16_t FRAGMENT_SHIFT      = 6;
    //@}

    ///////////////////////
    /// @name  Static methods for determining the contents of a flag.
    //@{
    static inline bool isMapped(uint16_t flag) {return(!(flag & UNMAPPED));}
    static inline bool isMateMapped(uint16_t flag) {return(!(flag & MATE_UNMAPPED));}

    static inline bool isPaired(uint16_t flag) {return(flag & PAIRED);}
    static inline bool isReverse(uint16_t flag) {return(flag & REVERSE);}
    static inline bool isMateReverse(uint16_t flag) {return(flag & MATE_REVERSED);}
    static inline bool isProperPair(uint16_t flag) 
    {
        // Proper pair is only applicable if also paired.
        return(isPaired(flag) && (flag & PROPER_PAIR));
    }
    static inline bool isDuplicate(uint16_t flag) {return(flag & DUPLICATE);}
    static inline bool isQCFailure(uint16_t flag) {return(flag & FAILED_QUALITY);}

    static inline bool isSecondary(uint16_t flag) {return(flag & SECONDARY_ALIGNMENT);}

    /// Return if it is the first fragment or not
    /// (if FIRST_READ is set and SECOND_READ is not).
    static inline bool isFirstFragment(uint16_t flag) 
    {
        // first fragment if FIRST_READ is set and SECOND_READ is not.
        return((flag & FIRST_READ) && !(flag & SECOND_READ));
    }
    /// Return if it is the last fragment or not
    /// (if FIRST_READ is not set and SECOND_READ is).
    static inline bool isLastFragment(uint16_t flag) 
    {
        // last fragment if FIRST_READ is not set and SECOND_READ is set.
        return(!(flag & FIRST_READ) && (flag & SECOND_READ));
    }
    /// Return if it is a middle fragment or not
    /// (if FIRST_READ is set and SECOND_READ is also set).
    static inline bool isMidFragment(uint16_t flag) 
    {
        // mid fragment if both FIRST_READ and SECOND_READ are set.
        return((flag & FIRST_READ) && (flag & SECOND_READ));
    }
    /// Return if it is an unknown fragment fragment or not
    /// (if FIRST_READ is not set and SECOND_READ is also not set).
    static inline bool isUnknownFragment(uint16_t flag) 
    {
        // unknown fragment index if neither FIRST_READ nor SECOND_READ are not.
        return(!(flag & FIRST_READ) && !(flag & SECOND_READ));
    }

    static inline uint8_t getFragmentType(uint16_t flag)
    {
        return((flag & FRAGMENT_INFO) >> FRAGMENT_SHIFT);
    }

    /// Mark the passed in flag as unmapped.
    static inline void setUnmapped(uint16_t& flag) { flag |= UNMAPPED;}
    /// Mark the passed in flag as not duplicate.
    static inline void setNotDuplicate(uint16_t& flag) { flag ^= DUPLICATE;}
    /// Mark the passed in flag as not duplicate.
    static inline void setDuplicate(uint16_t& flag) { flag |= DUPLICATE;}
    //@}

private:
    SamFlag();
};


#endif
