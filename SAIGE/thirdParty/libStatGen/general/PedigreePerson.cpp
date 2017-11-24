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

#include "PedigreePerson.h"
#include "Constant.h"
#include "StringArray.h"
#include "Error.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

Person::Person()
{
    zygosity = sex = 0;
    serial = traverse = -1;

    markers = new Alleles [markerCount];
    traits = new double [traitCount];
    covariates = new double [covariateCount];
    affections = new char [affectionCount];
    strings = new String [stringCount];

    for (int i = 0; i < traitCount; i++) traits[i] = _NAN_;
    for (int i = 0; i < covariateCount; i++) covariates[i] = _NAN_;
    for (int i = 0; i < affectionCount; i++) affections[i] = 0;

    filter = false;

    father = mother = NULL;
    sibs = NULL;
    sibCount = 0;

    ngeno = 0;
    hasBothParents = hasAllTraits = hasAllAffections = hasAllCovariates = false;
}

Person::~Person()
{
    delete [] markers;
    delete [] traits;
    delete [] affections;
    delete [] covariates;
    delete [] strings;

    if (sibCount) delete [] sibs;
}

void Person::Copy(Person & rhs)
{
    CopyIDs(rhs);
    CopyPhenotypes(rhs);
}

void Person::CopyPhenotypes(Person & rhs)
{
    for (int i = 0; i < Person::traitCount; i++)
        traits[i] = rhs.traits[i];
    for (int i = 0; i < Person::affectionCount; i++)
        affections[i] = rhs.affections[i];
    for (int i = 0; i < Person::covariateCount; i++)
        covariates[i] = rhs.covariates[i];
    for (int i = 0; i < Person::markerCount; i++)
        markers[i] = rhs.markers[i];
    ngeno = rhs.ngeno;
}

void Person::WipePhenotypes(bool remove_genotypes)
{
    for (int i = 0; i < traitCount; i++) traits[i] = _NAN_;
    for (int i = 0; i < covariateCount; i++) covariates[i] = _NAN_;
    for (int i = 0; i < affectionCount; i++) affections[i] = 0;

    if (remove_genotypes)
    {
        for (int i = 0; i < markerCount; i++)
            markers[i][0] = markers[i][1] = 0;
        ngeno = 0;
    }
}

void Person::CopyIDs(Person & rhs)
{
    famid = rhs.famid;
    pid = rhs.pid;
    fatid = rhs.fatid;
    motid = rhs.motid;
    sex = rhs.sex;
    zygosity = rhs.zygosity;
}

bool Person::CheckParents()
{
    hasBothParents = father != NULL && mother != NULL;

    if (!hasBothParents)
    {
        if (father != NULL || mother != NULL)
        {
            printf("Parent named %s for Person %s in Family %s is missing\n",
                   (father == NULL) ? (const char *) fatid : (const char *) motid,
                   (const char *) pid, (const char *) famid);
            return false;
        }
        else
            return true;
    }

    if (father->sex == SEX_FEMALE || mother->sex == SEX_MALE)
        // If parents are switched around, we can fix it...
    {
        Person * swap = father;
        father = mother;
        mother = swap;

        String temp = fatid;
        fatid = motid;
        motid = temp;
    }

    if (father->sex == SEX_FEMALE || mother->sex == SEX_MALE)
        // If things still don't make sense then the problem is more serious ...
    {
        printf("Parental sex codes don't make sense for Person %s in Family %s\n",
               (const char *) pid, (const char *) famid);
        return false;
    }

    return true;
}

void Person::AssessStatus()
{
    hasBothParents = father != NULL && mother != NULL;

    hasAllTraits = hasAllAffections = hasAllCovariates = true;

    ngeno = 0;
    for (int m = 0; m < markerCount; m++)
        if (isGenotyped(m))
            ngeno++;

    for (int t = 0; t < traitCount; t++)
        if (!isPhenotyped(t))
        {
            hasAllTraits = false;
            break;
        }

    for (int c = 0; c < covariateCount; c++)
        if (!isControlled(c))
        {
            hasAllCovariates = false;
            break;
        }

    for (int a = 0; a < affectionCount; a++)
        if (!isDiagnosed(a))
        {
            hasAllAffections = false;
            break;
        }
}

void Person::Order(Person * & p1, Person * & p2)
{
    if (p1->traverse > p2->traverse)
    {
        Person * temp = p1;
        p1 = p2;
        p2 = temp;
    }
}

int Person::GenotypedMarkers()
{
    int count = 0;

    for (int m = 0; m < Person::markerCount; m++)
        if (markers[m].isKnown())
            count++;

    return count;
}

bool Person::haveData()
{
    if (ngeno)
        return true;

    for (int i = 0; i < affectionCount; i++)
        if (affections[i] != 0)
            return true;

    for (int i = 0; i < traitCount; i++)
        if (traits[i] != _NAN_)
            return true;

    return false;
}

bool Person::isAncestor(Person * descendant)
{
    if (traverse > descendant->traverse)
        return false;

    if (serial == descendant->serial)
        return true;

    if (descendant->isFounder())
        return false;

    return (isAncestor(descendant->mother) ||
            isAncestor(descendant->father));
}
