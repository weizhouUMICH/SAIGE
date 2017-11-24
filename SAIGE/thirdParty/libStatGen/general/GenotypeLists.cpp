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

#include "GenotypeLists.h"

// When the next line is uncommented, the genotype elimination routines
// produce a lot of output useful for debugging
// #define DEBUG_ELIMINATOR

GenotypeList::GenotypeList()
{
    ignore = false;
}

bool GenotypeList::EliminateGenotypes(Pedigree & ped, Family * family, int marker)
{
    // First, allocate a genotype list for the family
    GenotypeList * list = new GenotypeList [family->count];

    // Next, update the possible allele lists for each individual
    InitializeList(list, ped, family, marker);

    // Then, do multiple rounds of elimination until a problem is found
    // or no changes are made

#ifdef DEBUG_ELIMINATOR
    Print(list, ped, family, marker);
#endif

    while (PairwiseCheck(list, ped, family) || FamilyCheck(list, ped, family))
#ifdef DEBUG_ELIMINATOR
        Print(list, ped, family, marker)
#endif
        ;

    for (int i = 0; i < family->count; i++)
        if (!list[i].ignore && list[i].allele1.Length() == 0)
        {
            printf("%s - Family %s has a subtle genotype inconsistency\n",
                   (const char *) ped.markerNames[marker], (const char *) family->famid);

            delete [] list;
            return false;
        }

    delete [] list;
    return true;
}

void GenotypeList::InitializeList(GenotypeList * list, Pedigree & ped, Family * family, int marker)
{
    for (int i = family->count - 1; i >= 0; i--)
    {
        Person & person = ped[family->path[i]];
        int id = person.traverse;
        bool maleX = person.sex == SEX_MALE && ped.chromosomeX;

#ifdef DEBUG_ELIMINATOR
        printf("Initializing genotype list for %s ...\n", (const char *) person.pid);
#endif

        // If an individual is genotyped ...
        if (person.markers[marker].isKnown())
        {
            // Their genotype list starts with just one entry!
            list[id].Dimension(1);
            list[id].SetGenotype(0, person.markers[marker][0], person.markers[marker][1]);
            list[id].alleles.Clear();
            list[id].alleles.Push(person.markers[marker][0]);
            list[id].alleles.PushIfNew(person.markers[marker][1]);
            list[id].ignore = false;

            // "Heterozygous" males have no possible genotypes
            if (maleX && person.markers[marker].isHeterozygous())
                list[id].Dimension(0);
        }
        else if (list[id].alleles.Length())
            if (person.sex == SEX_MALE && ped.chromosomeX)
            {
                // Males only carry one X chromosome
                list[id].Dimension(list[id].alleles.Length() + 1);

                for (int i = 0, out = 0; i < list[id].alleles.Length(); i++)
                    list[id].SetGenotype(out++, list[id].alleles[i], list[id].alleles[i]);
                list[id].SetGenotype(list[id].alleles.Length(), -1, -1);

                list[id].ignore = false;
            }
            else
            {
                // Build the genotype list based on the available allele lists
                int count = list[id].alleles.Length() * (list[id].alleles.Length() + 3) / 2 + 1;

                list[id].Dimension(count);

                for (int i = 0, out = 0; i < list[id].alleles.Length(); i++)
                {
                    // Allow for all pairs of "transmitted" alleles
                    for (int j = 0; j <= i; j++)
                        list[id].SetGenotype(out++, list[id].alleles[i], list[id].alleles[j]);

                    // Allow for an unstransmitted allele
                    list[id].SetGenotype(out++, list[id].alleles[i], -1);
                }

                // Allow for a pair of untransmitted alleles
                list[id].SetGenotype(count - 1, -1, -1);

                list[id].ignore = false;
            }
        else
            list[id].ignore = true;

        // If the individual is a founder this is all there is to it
        if (i < family->founders) continue;

        // If the individual is not a founder, update the parental genotype lists...
        int fatid = person.father->traverse;
        int motid = person.mother->traverse;

        for (int i = 0; i < list[id].alleles.Length(); i++)
        {
            list[motid].alleles.PushIfNew(list[id].alleles[i]);
            if (!maleX) list[fatid].alleles.PushIfNew(list[id].alleles[i]);
        }
    }
}

bool GenotypeList::PairwiseCheck(GenotypeList * list, Pedigree & ped, Family * family)
{
#ifdef DEBUG_ELIMINATOR
    printf("Checking Relative Pairs ...\n");
#endif

    bool changed = false;

    for (int i = family->count - 1; i >= family->founders; i--)
    {
        Person & person = ped[family->path[i]];

        int id = person.traverse;
        int fatid = person.father->traverse;
        int motid = person.mother->traverse;

        bool maleX = person.sex == SEX_MALE && ped.chromosomeX;

        if (list[id].ignore) continue;

        // Check if genotypes are consistent with paternal genotypes
        for (int i = 0; i < list[id].allele1.Length(); i++)
        {
            int al1 = list[id].allele1[i];
            int al2 = list[id].allele2[i];

            // Remove offspring genotypes incompatible with parental genotypes
            if ((maleX && !list[motid].Matches(al1) && al1 != -1) ||
                    (!maleX && !(al1 == -1 && al2 == -1) &&
                     !(list[fatid].Matches(al1) && (al2 == -1 || list[motid].Matches(al2))) &&
                     !((al2 == -1 || list[fatid].Matches(al2)) && list[motid].Matches(al1))))
            {
                list[id].Delete(i--);
                changed = true;
            }
        }

        // The offspring genotype list allows for a wild-card untransmitted allele
        // so any single parental genotype is possible
        if (list[id].Matches(-1))
            continue;

        // Check if genotypes are consistent with offspring genotypes
        for (int i = 0; i < list[motid].allele1.Length(); i++)
        {
            int al1 = list[motid].allele1[i];
            int al2 = list[motid].allele2[i];

            // Remove genotypes incompatible with offspring genotypes
            if (!list[id].Matches(al1) &&
                    !list[id].Matches(al2))
            {
                list[motid].Delete(i--);
                changed = true;
            }
        }

        // Males don't affect genotype lists for their fathers
        if (maleX) continue;

        // Check if genotypes are consistent with offspring genotypes
        for (int i = 0; i < list[fatid].allele1.Length(); i++)
        {
            int al1 = list[fatid].allele1[i];
            int al2 = list[fatid].allele2[i];

            // Remove genotypes incompatible with offspring genotypes
            if (!list[id].Matches(al1) &&
                    !list[id].Matches(al2))
            {
                list[fatid].Delete(i--);
                changed = true;
            }
        }

#ifdef DEBUG_ELIMINATOR
        printf("Done checking individual %s\n", (const char *) person.pid);
        Print(list, ped, family, 0);
#endif
    }

    return changed;
}


bool GenotypeList::FamilyCheck(GenotypeList * list, Pedigree & ped, Family * family)
{
#ifdef DEBUG_ELIMINATOR
    printf("Checking Nuclear Families ...\n");
#endif

    bool changed = false;

    for (int i = family->count - 1; i >= family->founders; i--)
    {
        Person & person = ped[family->path[i]];

        int fatid = person.father->traverse;
        int motid = person.mother->traverse;

        // Only go through the loop once per sibship
        if (person.sibs[0] != &person || list[fatid].ignore || list[motid].ignore)
            continue;

#ifdef DEBUG_ELIMINATOR
        printf("Checking Sibship with %s ...\n", (const char *) person.pid);
#endif

        // Reset checked genotypes for the mother, father and child
        list[fatid].checked = 0;
        list[motid].checked = 0;

        for (int i = 0; i < person.sibCount; i++)
            list[person.sibs[i]->traverse].checked = 0;

        // Go through each of the paternal genotypes
        changed |= TrimParent(list, person, fatid, motid);

        // Go through each of maternal genotypes
        changed |= TrimParent(list, person, motid, fatid);

        // Sort out the unchecked offspring genotypes ...
        for (int i = 0; i < person.sibCount; i++)
        {
            int sibid = person.sibs[i]->traverse;
            bool maleX = person.sibs[i]->sex == SEX_MALE && ped.chromosomeX;

            // For dealing with male X chromosomes, the pairwise check is all we need
            if (maleX) continue;

            for (int j = list[sibid].checked; j < list[sibid].allele1.Length(); j++)
                changed |= Cleanup(list, person, motid, fatid, sibid, j);
        }

#ifdef DEBUG_ELIMINATOR
//      Print(list, ped, family, 0);
#endif
    }

    return changed;
}

bool GenotypeList::Matches(int genotype, int allele)
{
    return allele1[genotype] == allele || allele2[genotype] == allele;
}

bool GenotypeList::Matches(int allele)
{
    return allele1.Find(allele) != -1 || allele2.Find(allele) != -1;
}

int GenotypeList::SaveGenotype(int genotype)
{
    if (checked > genotype)
        return genotype;

    if (checked != genotype)
    {
        allele1.Swap(genotype, checked);
        allele2.Swap(genotype, checked);
    }

    return checked++;
}

bool GenotypeList::CheckTrio(GenotypeList * list, int fatid, int motid, int child,
                             int i, int j, int k)
{
    // TODO: add tests for this code...
    return   (list[fatid].Matches(i, list[child].allele1[k]) &&
             (list[motid].Matches(j, list[child].allele2[k]) || list[child].allele2[k] == -1)) ||
             ((list[fatid].Matches(i, list[child].allele2[k]) || list[child].allele2[k] == -1) &&
             list[motid].Matches(j, list[child].allele1[k])) ||
             (list[child].allele1[k] == -1 && list[child].allele2[k] == -1);
}

void GenotypeList::Dimension(int genotypes)
{
    allele1.Dimension(genotypes);
    allele2.Dimension(genotypes);
}

void GenotypeList::SetGenotype(int genotype, int al1, int al2)
{
    allele1[genotype] = al1;
    allele2[genotype] = al2;
}

void GenotypeList::Delete(int genotype)
{
    allele1.Delete(genotype);
    allele2.Delete(genotype);
}

bool GenotypeList::TrimParent(GenotypeList * list, Person & person, int motid, int fatid)
{
    bool trimmed = false;

    while (list[motid].checked < list[motid].allele1.Length())
    {
        int  current = list[motid].allele1.Length() - 1;
        bool saved = false;

        // Pair it with each possible paternal genotype
        for (int i = list[fatid].allele1.Length() - 1; i >= 0; i--)
        {
            int matches = 0;

            // Find out if the pairing is compatible with at least one genotype for each child
            for (int j = 0; j < person.sibCount; j++)
            {
                int sibid = person.sibs[j]->traverse;
                int maleX = person.sibs[j]->sex == SEX_MALE && person.chromosomeX;

                // Since we have done the pairwise check, there is nothing more
                // to do for males ...
                if (list[sibid].ignore || maleX)
                {
                    matches++;
                    continue;
                }

                for (int k = list[sibid].allele1.Length() - 1; k >= 0; k--)
                    if (CheckTrio(list, motid, fatid, sibid, current, i, k))
                    {
                        matches++;
                        break;
                    }

                if (matches != j + 1)
                    break;
            }

            // Save maternal and paternal genotypes, mark all compatible sibling genotypes
            if (matches == person.sibCount)
            {
                for (int j = 0; j < person.sibCount; j++)
                {
                    int sibid = person.sibs[j]->traverse;

                    for (int k = list[sibid].checked; k < list[sibid].allele1.Length(); k++)
                        if (CheckTrio(list, motid, fatid, sibid, current, i, k))
                            list[sibid].SaveGenotype(k);
                }

                list[motid].SaveGenotype(current);
                list[fatid].SaveGenotype(i);

                saved = true;

                break;
            }
        }

        if (!saved)
        {
            list[motid].Delete(current);
            trimmed = true;
        }
    }

    return trimmed;
}

bool GenotypeList::Cleanup(GenotypeList * list, Person & person, int motid, int fatid, int child, int geno)
{
    for (int current = 0; current < list[motid].allele1.Length(); current++)
        for (int i = list[fatid].allele1.Length() - 1; i >= 0; i--)
            if (CheckTrio(list, motid, fatid, child, current, i, geno))
            {
                int matches = 0;

                // Find out if the pairing is compatible with at least one genotype for each child
                for (int j = 0; j < person.sibCount; j++)
                {
                    int sibid = person.sibs[j]->traverse;
                    int maleX = person.sibs[j]->sex == SEX_MALE && person.chromosomeX;

                    // After completing the pairwise check, all males are guaranteed
                    // to be compatible with their mothers
                    if (list[sibid].ignore || maleX)
                    {
                        matches++;
                        continue;
                    }

                    for (int k = list[sibid].allele1.Length() - 1; k >= 0; k--)
                        if (CheckTrio(list, motid, fatid, sibid, current, i, k))
                        {
                            matches++;
                            break;
                        }

                    if (matches != j + 1)
                        break;
                }

                // Update list of compatible sibling genotypes
                if (matches == person.sibCount)
                    for (int j = 0; j < person.sibCount; j++)
                    {
                        int sibid = person.sibs[j]->traverse;

                        for (int k = list[sibid].checked; k < list[sibid].allele1.Length(); k++)
                            if (CheckTrio(list, motid, fatid, sibid, current, i, k))
                                list[sibid].SaveGenotype(k);

                        return false;
                    }
            }

    list[child].Delete(geno);

    return true;
}

void GenotypeList::Print(GenotypeList * list, Pedigree & ped, Family * family, int marker)
{
    MarkerInfo * info = ped.GetMarkerInfo(marker);

    for (int i = 0; i < family->count; i++)
    {
        printf("%s - ", (const char *) ped[family->path[i]].pid);

        for (int j = 0; j < list[i].allele1.Length(); j++)
        {
            if (list[i].allele1[j] == -1)
                printf("*/");
            else
                printf("%s/", (const char *) info->GetAlleleLabel(list[i].allele1[j]));

            if (list[i].allele2[j] == -1)
                printf("* ");
            else
                printf("%s ", (const char *) info->GetAlleleLabel(list[i].allele2[j]));
        }

        printf("\n");
    }
    printf("\n");
}

