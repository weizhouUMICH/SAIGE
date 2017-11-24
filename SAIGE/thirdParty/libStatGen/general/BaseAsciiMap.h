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

#ifndef _BASE_ASCII_MAP_H
#define _BASE_ASCII_MAP_H

#include "StringBasics.h"

/// Map between characters and the associated base type.
class BaseAsciiMap
{
public:
    /// Value associated with 'N' in the ascii to base map (bad read).
    static const int baseNIndex = 004;
    /// Value associated with any non-base character in the ascii to base
    /// map (unknown, bad data).
    static const int baseXIndex = 005;

    // Two arrays for converting back and forth between base pair character
    // value (ASCII) to a base integer in the range 0..3.  Note there is actually
    // a value 4 and 5, for 'N' (indelible) and 'M' (unknown to me).
    //
    /// Convert from int representation to the base.
    static const char int2base[];
    /// Convert from int representation to colorspace representation.
    static const char int2colorSpace[];
    static unsigned char base2complement[];

    /// The type of space (color or base) to use in the mapping.
    enum SPACE_TYPE {
        /// Base decision on the first raw seq character/type has yet 
        /// to be determined.
        UNKNOWN,
        BASE_SPACE, ///< Bases only (A,C,G,T,N).
        COLOR_SPACE ///< Color space only (0,1,2,3,.).
    };

    /// Map ASCII values to a 2 (or 3) bit encoding for the base pair value for
    /// both base and color space.
    /// 'A'/'a'/'0' -> 0; 'C'/'c'/'1' -> 1; 'G'/'g'/'2' -> 2; 'T'/'t'/'3' -> 3;
    /// 'N'/'n'/'4' -> 4; anything else -> 5.
    static unsigned char baseColor2int[256+1];   // base space read (ATCG)
    /// Map ASCII values to a 2 (or 3) bit encoding for the base pair value for
    /// just base space (ACTGNactgn).
    /// 'A'/'a' -> 0;  'C'/'c' -> 1;  'G'/'g' -> 2;  'T'/'t' -> 3;
    /// 'N'/'n' -> 4; anything else -> 5.
    static unsigned char base2int[256+1];        // base space read (ATCG)
    /// Map ASCII values to a 2 (or 3) bit encoding for the base pair value for
    /// just color space (0123).
    /// '0' -> 0; '1' -> 1; '2' -> 2; '3' -> 3; '4' -> 4; anything else -> 5.
    static unsigned char color2int[256+1];       // base space read (ATCG)

public:
    BaseAsciiMap();
    ~BaseAsciiMap();

    /// Set the base type based on the passed in option.
    inline void setBaseMapType(SPACE_TYPE spaceType)
    {
        resetPrimerCount();
        //First check to see if it is in base space.
        switch (spaceType)
        {
            case BASE_SPACE:
                // base space.
                myBase2IntMapPtr = base2int;
                break;
            case COLOR_SPACE:
                // color space.
                myBase2IntMapPtr = color2int;
                break;
            default:
                // Unknown map type, zero the pointer.
                myBase2IntMapPtr = NULL;
                break;
        }
    };

    /// Returns the baseIndex value for the character passed in.
    inline int getBaseIndex(const char& letter)
    {
        if (myBase2IntMapPtr == NULL)
        {
            // Check to see if we have hit the number of primer bases.
            if (myPrimerCount < myNumPrimerBases)
            {
                // Still expecting primer bases, so lookup
                // the letter in the base map.
                ++myPrimerCount;
                return(base2int[(int)letter]);
            }

            // Have already processed all the primers, so determine
            // whether this is base or color space.

            // Need to determime the base type.
            setBaseMapType(letter);

            // If it is still null, return invalid.  Will be set when the first
            // letter is either color or base.
            if (myBase2IntMapPtr == NULL)
            {
                return(baseXIndex);
            }
        }

        // Also check if configured as color space that the primers are correct.
        if ((myBase2IntMapPtr == color2int) && (myPrimerCount < myNumPrimerBases))
        {
            // Still expecting primer bases, so lookup
            // the letter in the base map.
            ++myPrimerCount;
            return(base2int[(int)letter]);
        }

        return myBase2IntMapPtr[(int)letter];
    }

    /// Return the space type that is currently set.
    inline SPACE_TYPE getSpaceType()
    {
        if (myBase2IntMapPtr == base2int)
        {
            return(BASE_SPACE);
        }
        else if (myBase2IntMapPtr == color2int)
        {
            return(COLOR_SPACE);
        }
        else
        {
            return(UNKNOWN);
        }
    }

    /// Set the number of primer bases expected before the actual
    /// base/color space type occurs for the rest of the entries.
    void setNumPrimerBases(int numPrimerBases)
    {
        myNumPrimerBases = numPrimerBases;
    }

    /// Reset the number of primers to 0.
    void resetPrimerCount()
    {
        myPrimerCount = 0;
    };

    /// Reset the base mapping type to UNKNOWN.
    void resetBaseMapType()
    {
        myBase2IntMapPtr = NULL;
        resetPrimerCount();
    };

private:
    // Set the base type based on the passed in letter.
    // If the letter is in neither the color space or the base space, both
    // will be allowed.
    inline void setBaseMapType(const char& letter)
    {
        //First check to see if it is in base space.
        if (base2int[(int)letter] != baseXIndex)
        {
            // This is a valid base space index, so it is base space.
            myBase2IntMapPtr = base2int;
        }
        else if (color2int[(int)letter] != baseXIndex)
        {
            // This is a valid color space index, so it is base space.
            myBase2IntMapPtr = color2int;
        }
        else
        {
            // Unknown map type, zero the pointer.
            myBase2IntMapPtr = NULL;
        }
    };


    // The number of primer bases to expect for a color-space file.
    unsigned int myNumPrimerBases;

    // This is the number of primer bases that have been seen since
    // the map type was set/reset.
    unsigned int myPrimerCount;

    unsigned char* myBase2IntMapPtr;
};

#endif
