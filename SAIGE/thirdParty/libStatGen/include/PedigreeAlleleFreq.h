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

#ifndef __ALLELEFREQUENCIES_H__
#define __ALLELEFREQUENCIES_H__

#include "Pedigree.h"

int  CountAlleles(Pedigree & ped, int marker);
void LumpAlleles(Pedigree & ped, int marker, double threshold, bool reorder);

#define FREQ_ALL        0
#define FREQ_FOUNDERS   1
#define FREQ_EQUAL      2

// Returns true if frequencies estimated, false if previous information okay
bool EstimateFrequencies(Pedigree & ped, int marker, int estimator);

#endif


