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

#include "BaseAsciiMap.h"

//
// Map ASCII values to a 2 (or 3) bit encoding for the base pair value for
// both base and color space.
//  class 0 -> 'A' (Adenine - 0x41 and 0x61)
//  class 1 -> 'C' (Cytosine - 0x43 and 0x63)
//  class 2 -> 'G' (Guanine - 0x47 and 0x67)
//  class 3 -> 'T' (Thymine - 0x54 and 0x74)
//  class 4 -> 'N' (Unknown - read error or incomplete data - 0x4E and 0x6E)
//  class 5 -> not a valid DNA base pair character
//
// Note: The +1 array size is for the terminating NUL character
//
// NB: This table also maps 0, 1, 2, and 3 to the corresponding integers,
// and '.' to class 4.  This allows ABI SOLiD reads to be converted
// to integers via ReadIndexer::Word2Integer.
//
unsigned char BaseAsciiMap::baseColor2int[256+1] =
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x00-0x0F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x10-0x1F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\004\005"  // 0x20-0x2F
    "\000\001\002\003\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x30-0x3F
    "\005\000\005\001\005\005\005\002\005\005\005\005\005\005\004\005"  // 0x40-0x4F
    "\005\005\005\005\003\005\005\005\005\005\005\005\005\005\005\005"  // 0x50-0x5F
    "\005\000\005\001\005\005\005\002\005\005\005\005\005\005\004\005"  // 0x60-0x6F
    "\005\005\005\005\003\005\005\005\005\005\005\005\005\005\005\005"  // 0x70-0x7F
// not used, but included for completeness:
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x80-0x8F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x90-0x9F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xA0-0xAF
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xB0-0xBF
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xC0-0xCF
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xD0-0xDF
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xE0-0xEF
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xF0-0xFF
    ;

// Map ASCII values to a 2 (or 3) bit encoding for the base pair value for
// just base space (ACTGNactgn).
unsigned char BaseAsciiMap::base2int[256+1] =
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x00-0x0F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x10-0x1F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x20-0x2F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x30-0x3F
    "\005\000\005\001\005\005\005\002\005\005\005\005\005\005\004\005"  // 0x40-0x4F
    "\005\005\005\005\003\005\005\005\005\005\005\005\005\005\005\005"  // 0x50-0x5F
    "\005\000\005\001\005\005\005\002\005\005\005\005\005\005\004\005"  // 0x60-0x6F
    "\005\005\005\005\003\005\005\005\005\005\005\005\005\005\005\005"  // 0x70-0x7F
// not used, but included for completeness:
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x80-0x8F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x90-0x9F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xA0-0xAF
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xB0-0xBF
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xC0-0xCF
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xD0-0xDF
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xE0-0xEF
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xF0-0xFF
    ;

// Map ASCII values to a 2 (or 3) bit encoding for the base pair value for
// just color space (0123).
unsigned char BaseAsciiMap::color2int[256+1] =
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x00-0x0F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x10-0x1F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\004\005"  // 0x20-0x2F
    "\000\001\002\003\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x30-0x3F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x40-0x4F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x50-0x5F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x60-0x6F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x70-0x7F
// not used, but included for completeness:
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x80-0x8F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0x90-0x9F
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xA0-0xAF
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xB0-0xBF
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xC0-0xCF
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xD0-0xDF
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xE0-0xEF
    "\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005\005"  // 0xF0-0xFF
    ;


//
// This is obviously for base space use only:
//
const char BaseAsciiMap::int2base[] = "ACGTNMXXXXXXXXXX";
//
// convert int to color space value
//
const char BaseAsciiMap::int2colorSpace[] = "0123NXXXXXXXXXXX";

/// This table maps 5' base space to the 3' complement base space
/// values, as well as 5' color space values to the corresponding
/// 3' complement color space values.
///
/// In both cases, invalids are mapped to 'N', which isn't accurate
/// for ABI SOLiD, but internally it shouldn't matter (on output it
/// will).
unsigned char BaseAsciiMap::base2complement[256+1 /* for NUL char */] =
    "NNNNNNNNNNNNNNNN"  // 0x00-0x0F
    "NNNNNNNNNNNNNNNN"  // 0x10-0x1F
    "NNNNNNNNNNNNNNNN"  // 0x20-0x2F
    "0123NNNNNNNNNNNN"  // 0x30-0x3F
    "NTNGNNNCNNNNNNNN"  // 0x40-0x4F
    "NNNNANNNNNNNNNNN"  // 0x50-0x5F
    "NTNGNNNCNNNNNNNN"  // 0x60-0x6F
    "NNNNANNNNNNNNNNN"  // 0x70-0x7F
// not used, but included for completeness:
    "NNNNNNNNNNNNNNNN"  // 0x80-0x8F
    "NNNNNNNNNNNNNNNN"  // 0x90-0x9F
    "NNNNNNNNNNNNNNNN"  // 0xA0-0xAF
    "NNNNNNNNNNNNNNNN"  // 0xB0-0xBF
    "NNNNNNNNNNNNNNNN"  // 0xC0-0xCF
    "NNNNNNNNNNNNNNNN"  // 0xD0-0xDF
    "NNNNNNNNNNNNNNNN"  // 0xE0-0xEF
    "NNNNNNNNNNNNNNNN"  // 0xF0-0xFF
    ;

BaseAsciiMap::BaseAsciiMap()
        : myNumPrimerBases(1)
{
    myBase2IntMapPtr = NULL;
}

BaseAsciiMap::~BaseAsciiMap()
{
}
