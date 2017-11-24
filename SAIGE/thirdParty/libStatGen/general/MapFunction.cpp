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

#include "MapFunction.h"
#include "MathConstant.h"

#include <math.h>

double DistanceToRecombination(double distance)
{
    return (1.0 - exp(-2.0 * distance)) * 0.5;
}

double RecombinationToDistance(double recombination)
{
    return (log(max(1.0 - 2 * recombination, 1e-7)) * -0.5);
}

double KosambiDistanceToRecombination(double distance)
{
    double e_to_4x = exp(4.0 * distance);

    return (0.5 *(e_to_4x - 1.0) / (e_to_4x + 1.0));
}

double RecombinationToKosambiDistance(double theta)
{
    return 0.25 * log((1.0 + 2*theta) / max(1.0 - 2.0*theta, 1e-7));
}
