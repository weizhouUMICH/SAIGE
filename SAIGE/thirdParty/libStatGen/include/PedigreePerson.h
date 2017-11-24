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

#ifndef __PEDPERSON_H__
#define __PEDPERSON_H__

#include "Constant.h"
#include "PedigreeAlleles.h"
#include "PedigreeGlobals.h"
#include "StringArray.h"
#include "IntArray.h"

#define  SEX_MALE       1
#define  SEX_FEMALE     2
#define  SEX_UNKNOWN    0

class Person : public PedigreeGlobals
{
public:
    String      famid;
    String      pid;
    String      motid;
    String      fatid;
    int         sex;
    int         zygosity;
    int         serial, traverse;

    Alleles *   markers;
    double *    traits;
    char *      affections;
    double *    covariates;
    String *    strings;

    Person *    father;
    Person *    mother;

    int         sibCount;
    Person **   sibs;

    int         ngeno;

    bool        filter;

    Person();
    ~Person();

    bool isHalfSib(Person & sib)
    {
        return hasBothParents &&
               ((sib.father == father) ^(sib.mother == mother));
    }

    bool isSib(Person & sib)
    {
        return hasBothParents &&
               (sib.father == father) && (sib.mother == mother);
    }

    bool isTwin(Person & twin)
    {
        return (zygosity != 0) && (zygosity == twin.zygosity) && isSib(twin);
    }

    bool isMzTwin(Person & mzTwin)
    {
        return (zygosity & 1) && (zygosity == mzTwin.zygosity) && isSib(mzTwin);
    }

    // Check that both parents or none are available
    // Verify that fathers are male and mothers are female
    bool CheckParents();

    // Assess status before using quick diagnostics functions
    void AssessStatus();

    // Quick diagnostics
    bool isFounder()
    {
        return !hasBothParents;
    }
    bool isSexed()
    {
        return sex != 0;
    }
    bool isGenotyped(int m)
    {
        return markers[m].isKnown();
    }
    bool isFullyGenotyped()
    {
        return ngeno == markerCount;
    }
    bool isControlled(int c)
    {
        return covariates[c] != _NAN_;
    }
    bool isFullyControlled()
    {
        return hasAllCovariates;
    }
    bool isPhenotyped(int t)
    {
        return traits[t] != _NAN_;
    }
    bool isFullyPhenotyped()
    {
        return hasAllTraits;
    }
    bool isDiagnosed(int a)
    {
        return affections[a] != 0;
    }
    bool isFullyDiagnosed()
    {
        return hasAllAffections;
    }
    bool haveData();
    bool isAncestor(Person * descendant);

    int GenotypedMarkers();

    static void Order(Person * & p1, Person * & p2);

    void Copy(Person & rhs);
    void CopyIDs(Person & rhs);
    void CopyPhenotypes(Person & rhs);
    void WipePhenotypes(bool remove_genotypes = true);

private:

    bool hasAllCovariates, hasAllTraits,
    hasAllAffections, hasBothParents;
};

#endif




