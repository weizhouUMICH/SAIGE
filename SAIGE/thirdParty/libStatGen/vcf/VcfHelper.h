/*
 *  Copyright (C) 2011  Regents of the University of Michigan,
 *                      Hyun Min Kang, Matthew Flickenger, Matthew Snyder,
 *                      and Goncalo Abecasis
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


#ifndef __VCF_HELPER_H__
#define __VCF_HELPER_H__

#include <string>
#include "ReusableVector.h"

/// This header file provides helper methods for dealing with VCF Files.
class  VcfHelper
{
public:
    /// Parse the string at the specified delimiters into
    /// the specified reusable vector.
    static void parseString(const std::string& inputString, 
                            char delim,
                            ReusableVector<std::string>& outputVector);
};

#endif
