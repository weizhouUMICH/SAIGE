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

#include "BaseQualityHelper.h"

#include <math.h>

baseQualityConvertor bQualityConvertor;

baseQualityConvertor::baseQualityConvertor()
{
    // Create a quick lookup table to speed up conversion of
    // base quality values stored as log10 (error rates) into
    // fractional error rates
    for (int i = 0; i <= 255; i++)
        doubleLookup[i] = pow(0.1, i * 0.1);
    // doubleLookup[255] = 0.0;
}

double baseQualityConvertor::toDouble(unsigned char bq)
{
    return doubleLookup[bq];
}

