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

#ifndef __PEDFAMILY_H__
#define __PEDFAMILY_H__

#include "PedigreeAlleles.h"
#include "PedigreePerson.h"
#include "StringBasics.h"

class Pedigree;

class Family
{
public:
    Pedigree & ped;
    String   famid;
    int      serial;
    int      first, last;    // sentinel family members
    int      count;          // number of individuals in pedigree
    int      founders;       // number of founders in pedigree
    int      nonFounders;    // number of non-founders in pedigree
    int      mzTwins;        // number of MZ twins, excluding 1st twin in set
    int      * path;         // traverses the pedigree so that ancestors
    // preceed their descendants

    int      generations;    // Rough classification as:
    //  1 -- all individuals are unrelated
    //  2 -- two generations (inc. multiple couples)
    //  3 -- three or more generations

    bool   isNuclear()
    {
        return (generations == 2) && (founders == 2);
    }

    Family(Pedigree & ped, int top, int bottom, int serial = 0);
    ~Family();

    int  ConnectedGroups(IntArray * groupMembership = NULL);

private:
    void ShowInvalidCycles();

    Family & operator = (Family & rhs);
//      void Mark(int who, int group, IntArray * stack, IntArray & group_id );
};

#endif

