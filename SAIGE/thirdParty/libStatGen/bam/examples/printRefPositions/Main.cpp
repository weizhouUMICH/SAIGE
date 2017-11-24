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

#include <iostream>
#include <string.h>
#include <stdlib.h>

#include "PrintRefPositions.h"

int main(int argc, char ** argv)
{
    std::string inFile = "../../test/testFiles/sortedBam.bam";
    std::string indexFile = "../../test/testFiles/sortedBam.bam.bai";
    std::string rname = "1";
    int startPosition = 1013;
    int endPosition = 1751;
    if(argc == 6)
    {
        inFile = argv[1];
        indexFile = argv[2];
        rname = argv[3];
        startPosition = atoi(argv[4]);
        endPosition = atoi(argv[5]);
    }
    printRefPositions(inFile, indexFile, rname, startPosition, endPosition);
    return(0);
}



