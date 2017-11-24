/*
 *  Copyright (C) 2010-2012  Regents of the University of Michigan
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

#ifndef __BASE_UTILITIES_H__
#define __BASE_UTILITIES_H__

#include <stdint.h>
#include <string>


/// This class contains static utilities for performing
/// basic operations on bases.
class BaseUtilities
{
public:
    /// Returns whether or not the specified bases is
    /// an indicator for ambiguity.
    /// \return true if base = 'n' or 'N' or '.'
    static bool isAmbiguous(char base);

    /// Returns whether or not two bases are equal (case insensitive),
    /// if one of the bases is '=', the bases are consided
    /// to be equal.
    static bool areEqual(char base1, char base2);

    /// Get phred base quality from the specified ascii quality.
    static uint8_t getPhredBaseQuality(char charQuality);

    /// Get ascii quality from the specified phred quality.
    static char getAsciiQuality(uint8_t phredQuality);

    static void reverseComplement(std::string& sequence);

    /// Character used when the quality is unknown.
    static const char UNKNOWN_QUALITY_CHAR = ' ';
    /// Int value used when the quality is unknown.
    static const uint8_t UNKNOWN_QUALITY_INT = 0xFF;
};


#endif
