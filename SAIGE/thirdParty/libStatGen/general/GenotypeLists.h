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

#ifndef __GENOTYPE_ELIMINATION__
#define __GENOTYPE_ELIMINATION__

#include "Pedigree.h"

class GenotypeList
{
public:

    IntArray allele1, allele2;
    IntArray alleles;

    bool ignore;
    int  checked;

    GenotypeList();

    static bool EliminateGenotypes(Pedigree & ped, Family * family, int marker);

    void   Dimension(int genotypes);
    void   Delete(int genotype);

    bool   Matches(int genotype, int allele);
    bool   Matches(int allele);

    int    SaveGenotype(int genotype);
    void   SetGenotype(int genotype, int al1, int al2);

private:
    static void InitializeList(GenotypeList * list, Pedigree & p, Family * f, int marker);
    static bool PairwiseCheck(GenotypeList * list, Pedigree & p, Family * f);
    static bool FamilyCheck(GenotypeList * list, Pedigree & p, Family * f);

    static bool CheckTrio(GenotypeList * list, int fatid, int motid, int child, int i, int j, int k);
    static bool TrimParent(GenotypeList * list, Person & person, int fatid, int motid);
    static bool Cleanup(GenotypeList * list, Person & person, int fatid, int motid, int child, int geno);

    static void Print(GenotypeList * List, Pedigree & p, Family * f, int marker);
};



#endif
