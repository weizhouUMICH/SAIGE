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

#include "Pedigree.h"
#include "GenotypeLists.h"
#include "MemoryInfo.h"
#include "Constant.h"
#include "Error.h"
#include "Sort.h"

#include <stdlib.h>

bool   Pedigree::sexAsCovariate = false;
String Pedigree::missing("-99.999");

Pedigree::Pedigree() : pd()
{
    haveTwins = count = 0;
    size = 10000;
    persons = new Person *[size];
    familyCount = 0;
    families = new Family * [1];
    multiPd = NULL;
    multiFileCount = 0;
}

Pedigree::~Pedigree()
{
    for (int i = 0; i < count; i++)
        delete persons[i];

    for (int i = 0; i < familyCount; i++)
        delete families[i];

    delete [] families;
    delete [] persons;

    if (multiPd != NULL)
        delete [] multiPd;
}

void Pedigree::Sort()
{
    QuickSort(persons, count, sizeof(Person *),
              COMPAREFUNC Pedigree::ComparePersons);

    haveTwins = 0;

    // Check for structural problems in input pedigree
    bool problem = false;

    // Check that we have no duplicates...
    for (int i = 1; i < count; i++)
        if (ComparePersons((const Person **) &persons[i-1],
                           (const Person **) &persons[i]) == 0)
        {
            printf("Family %s: Person %s is duplicated\n",
                   (const char *) persons[i]->famid,
                   (const char *) persons[i]->pid);
            problem = true;

            do
            {
                i++;
            }
            while (i < count &&
                    ComparePersons((const Person **) &persons[i-1],
                                   (const Person **) &persons[i]) == 0);
        }

    // Assign parents...
    for (int i = 0; i < count; i++)
    {
        persons[i]->serial = i;
        persons[i]->father = FindPerson(persons[i]->famid, persons[i]->fatid);
        persons[i]->mother = FindPerson(persons[i]->famid, persons[i]->motid);

        problem |= !persons[i]->CheckParents();

        persons[i]->AssessStatus();

        // Check if we have any twins...
        haveTwins |= persons[i]->zygosity;
    }

    if (problem)
        error("Please correct problems with pedigree structure\n");

    MakeSibships();
    MakeFamilies();
}

void Pedigree::MakeSibships()
{
    Person ** sibs = new Person * [count];
    for (int i = 0; i < count; i++)
        sibs[i] = persons[i];

    QuickSort(sibs, count, sizeof(Person *),
              COMPAREFUNC Pedigree::CompareParents);

    for (int first = 0; first < count; first++)
        if (!sibs[first]->isFounder())
        {
            int last = first + 1;
            while (last < count)
                if (sibs[first]-> mother != sibs[last]->mother ||
                        sibs[first]-> father != sibs[last]->father)
                    break;
                else last++;
            last --;

            for (int j = first; j <= last; j++)
            {
                if (sibs[j]->sibCount) delete [] sibs[j]->sibs;
                sibs[j]->sibCount = last - first + 1;
                sibs[j]->sibs = new Person * [sibs[j]->sibCount];
                for (int k = first; k <= last; k++)
                    sibs[j]->sibs[k - first] = sibs[k];
            }
            first = last;
        }
    delete [] sibs;
}

void Pedigree::MakeFamilies()
{
    for (int i = 0; i < familyCount; i++)
        delete families[i];
    delete [] families;

    familyCount = 0;
    families = new Family * [count];

    for (int first=0; first < count; first++)
    {
        int last = first;
        while (last < count)
            if (SlowCompare(persons[first]->famid, persons[last]->famid) == 0)
                last++;
            else break;

        families[familyCount] = new Family(*this, first, --last, familyCount);

        first = last;
        familyCount++;
    }
}

// Utility functions for finding a person in a pedigree

struct PedigreeKey
{
    const char * famid;
    const char * pid;
};

int CompareKeyToPerson(PedigreeKey * key, Person ** p)
{
    int result = SlowCompare(key->famid, (**p).famid);

    if (result != 0)
        return result;

    return SlowCompare(key->pid, (**p).pid);
}

int CompareKeyToFamily(PedigreeKey * key, Family ** f)
{
    return SlowCompare(key->famid, (**f).famid);
}

Person * Pedigree::FindPerson(const char * famid, const char * pid)
{
    PedigreeKey key;
    key.famid = famid;
    key.pid   = pid;

    Person ** result = (Person **) BinarySearch
                       (&key, persons, count, sizeof(Person *),
                        COMPAREFUNC CompareKeyToPerson);

    return (result == NULL) ? (Person *) NULL : *result;
}

Person * Pedigree::FindPerson(const char *famid, const char *pid, int universe)
{
    PedigreeKey key;
    key.famid = famid;
    key.pid   = pid;

    Person ** result = (Person **) BinarySearch
                       (&key, persons, universe, sizeof(Person *),
                        COMPAREFUNC CompareKeyToPerson);

    return (result == NULL) ? (Person *) NULL : *result;
}

Family * Pedigree::FindFamily(const char * famid)
{
    PedigreeKey key;
    key.famid = famid;

    Family ** result = (Family **) BinarySearch
                       (&key, families, familyCount, sizeof(Family *),
                        COMPAREFUNC CompareKeyToFamily);

    return (result == NULL) ? (Family *) NULL : *result;
}

int Pedigree::CountAlleles(int marker)
{
    return ::CountAlleles(*this, marker);
}

void Pedigree::LumpAlleles(double min, bool reorder)
{
    if (min > 0.0)
        printf("Lumping alleles with frequencies of %.2f or less...\n\n", min);

    for (int m=0; m < markerCount; m++)
        ::LumpAlleles(*this, m, min, reorder);
}

void Pedigree::EstimateFrequencies(int estimator, bool quiet)
{
    bool estimated = false;
    int  line = 3;

    const char * estimators[] = 
        { "using all genotypes", "using founder genotypes", "assumed equal" };

    bool condensed = markerCount > 100;
    int  grain = markerCount / 50, estimates = 0;

    for (int m=0; m < markerCount; m++)
        if (::EstimateFrequencies(*this, m, estimator))
            if (!quiet)
            {
                if (!estimated)
                    printf("Estimating allele frequencies... [%s]\n   ",
                           estimators[estimator]), estimated = true;

                if (condensed)
                {
                    if (estimates++ % grain == 0)
                    {
                        printf(".");
                        fflush(stdout);
                    }
                    continue;
                }

                if (line + markerNames[m].Length() + 1 > 79)
                    printf("\n   "), line = 3;

                printf("%s ", (const char *) markerNames[m]);
                line += markerNames[m].Length() + 1;
            }

    if (estimated)
        printf(condensed ? "\nDone estimating frequencies for %d markers\n\n" : "\n\n", estimates);
}

int Pedigree::ComparePersons(const Person ** p1, const Person ** p2)
{
    int result = SlowCompare((**p1).famid, (**p2).famid);

    if (result != 0) return result;

    return SlowCompare((**p1).pid, (**p2).pid);
}

int Pedigree::CompareParents(const Person ** p1, const Person ** p2)
{
    int result = SlowCompare((**p1).famid, (**p2).famid);

    if (result) return result;

    result = SlowCompare((**p1).fatid, (**p2).fatid);

    if (result) return result;

    return SlowCompare((**p1).motid, (**p2).motid);
}

void Pedigree::Grow()
{
    size *= 2;

    Person ** temp = new Person * [size];
    if (temp == NULL) error("Out of memory");

    for (int i=0; i<count; i++)
        temp[i] = persons[i];

    delete [] persons;
    persons = temp;
}

void Pedigree::Add(Person & rhs)
{
    if (count == size)
        Grow();

    persons[count] = new Person();
    persons[count++]->Copy(rhs);
}

void Pedigree::WriteDataFile(FILE * output)
{
    // write in the following order:
    // markers, traits, affections, covariates

    if (haveTwins)
        fprintf(output, " Z  Zygosity\n");

    for (int m = 0; m < markerCount; m++)
        fprintf(output, " M  %s\n", (const char *) markerNames[m]);

    for (int t = 0; t < traitCount; t++)
        fprintf(output, " T  %s\n", (const char *) traitNames[t]);

    for (int a = 0; a < affectionCount; a++)
        fprintf(output, " A  %s\n", (const char *) affectionNames[a]);

    for (int c = 0; c < covariateCount; c++)
        fprintf(output, " C  %s\n", (const char *) covariateNames[c]);

    for (int s = 0; s < stringCount; s++)
        fprintf(output, " $  %s\n", (const char *) stringNames[s]);

    fprintf(output, " E  END-OF-DATA \n");
}

void Pedigree::WritePedigreeFile(FILE * output)
{
    MarkerInfo ** info = new MarkerInfo * [markerCount];

    for (int i = 0; i < markerCount; i++)
        info[i] = GetMarkerInfo(i);

    for (int i = 0; i < count; i++)
        WriteRecodedPerson(output, i, info);
    fprintf(output, "end\n");

    delete [] info;
}

void Pedigree::WritePerson(FILE * output, int person, const char * famid,
                           const char * pid, const char * fatid, const char * motid)
{
    WriteRecodedPerson(output, person, NULL, famid, pid, fatid, motid);
}

void Pedigree::WriteRecodedPerson(
    FILE * output, int person, MarkerInfo ** markerInfo,
    const char * famid, const char * pid, const char * fatid,
    const char * motid)
{
    Person * p = persons[person];

    if (famid == NULL) famid = p->famid;
    if (pid == NULL)   pid = p->pid;
    if (fatid == NULL) fatid = p->fatid;
    if (motid == NULL) motid = p->motid;

    // write in the following order:
    // markers, traits, affections, covariates

    fprintf(output, "%s\t%s\t%s\t%s\t%d\t",
            famid, pid, fatid, motid, p->sex);

    const char * twinCodes[] = {"0", "MZ", "DZ"};

    if (haveTwins)
    {
        if (p->zygosity <= 2)
            fprintf(output, "%s\t", twinCodes[p->zygosity]);
        else
            fprintf(output, "%d\t", p->zygosity);
    }

    for (int m = 0; m < markerCount; m++)
        if (markerInfo == NULL)
            fprintf(output, markerCount < 20 ? "%3d/%3d\t" : "%d/%d\t",
                    p->markers[m][0], p->markers[m][1]);
        else
            fprintf(output, markerCount < 20 ? "%3s/%3s\t" : "%s/%s\t",
                    (const char *) markerInfo[m]->GetAlleleLabel(p->markers[m][0]),
                    (const char *) markerInfo[m]->GetAlleleLabel(p->markers[m][1]));

    for (int t = 0; t < traitCount; t++)
        if (p->isPhenotyped(t))
            fprintf(output, "%.3f\t", p->traits[t]);
        else
            fprintf(output, "x\t");

    for (int a = 0; a < affectionCount; a++)
        if (p->isDiagnosed(a))
            fprintf(output, "%d\t", p->affections[a]);
        else
            fprintf(output, "x\t");

    for (int c = 0; c < covariateCount; c++)
        if (p->isControlled(c))
            fprintf(output, "%.3f\t", p->covariates[c]);
        else
            fprintf(output, "x\t");

    for (int s = 0; s < stringCount; s++)
        if (!p->strings[s].IsEmpty())
            fprintf(output, "%s\t", (const char *) p->strings[s]);
        else
            fprintf(output, ".\t");

    fprintf(output, "\n");
}

void Pedigree::WriteDataFile(const char * output)
{
    FILE * f = fopen(output, "wt");
    if (f == NULL) error("Couldn't open data file %s", output);
    WriteDataFile(f);
    fclose(f);
}

void Pedigree::WritePedigreeFile(const char * output)
{
    FILE * f = fopen(output, "wt");
    if (f == NULL) error("Couldn't open pedigree file %s", output);
    WritePedigreeFile(f);
    fclose(f);
}

void Pedigree::PrepareDichotomization()
{

    for (int t = 0; t < traitCount; t++)
    {
        String new_affection = traitNames[t] + "*";
        GetAffectionID(new_affection);
    }
}

int Pedigree::Dichotomize(int t, double mean)
{
    String new_affection = traitNames[t] + "*";

    int af = GetAffectionID(new_affection);

    if (mean == _NAN_)
    {
        mean  = 0.0;
        double dcount = 0;
        for (int i = 0; i < count; i++)
            if (persons[i]->isPhenotyped(t) &&
                    !persons[i]->isFounder())
            {
                mean += persons[i]->traits[t];
                dcount ++;
            }

        if (!dcount) return af;

        mean /= dcount;
    }

    printf("Dichotomizing %s around mean of %.3f ...\n",
           (const char *) traitNames[t], mean);

    for (int i = 0; i < count; i++)
        if (persons[i]->isPhenotyped(t) && !persons[i]->isFounder())
            persons[i]->affections[af] = persons[i]->traits[t] > mean ? 2 : 1;
        else
            persons[i]->affections[af] = 0;

    Sort();

    return af;
}

void Pedigree::DichotomizeAll(double mean)
{
    for (int t = 0; t < traitCount; t++)
        Dichotomize(t, mean);
}

bool Pedigree::InheritanceCheck(bool abortIfInconsistent)
{
    bool fail = false;

    if (haveTwins) fail |= TwinCheck();

    if (chromosomeX)
        fail |= SexLinkedCheck();
    else
        fail |= AutosomalCheck();

    if (fail && abortIfInconsistent)
        error("Mendelian inheritance errors detected\n");

    return !fail;
}

bool Pedigree::AutosomalCheck()
{
    // Arrays indicating which alleles and homozygotes occur
    IntArray haplos, genos, counts, failedFamilies;

    bool fail = false;

    // For each marker ...
    for (int m = 0; m < markerCount; m++)
    {
        MarkerInfo * info = GetMarkerInfo(m);

        // Summary for marker
        int alleleCount = CountAlleles(m);
        int genoCount = alleleCount * (alleleCount + 1) / 2;

        // Initialize arrays
        haplos.Dimension(alleleCount + 1);
        haplos.Set(-1);

        genos.Dimension(genoCount + 1);
        genos.Set(-1);

        failedFamilies.Dimension(familyCount);
        failedFamilies.Zero();

        counts.Dimension(alleleCount + 1);

        for (int f = 0; f < familyCount; f++)
            for (int i = families[f]->first; i <= families[f]->last; i++)
                if (!persons[i]->isFounder() && persons[i]->sibs[0] == persons[i])
                {
                    // This loop runs once per sibship
                    Alleles fat = persons[i]->father->markers[m];
                    Alleles mot = persons[i]->mother->markers[m];
                    bool    fgeno = fat.isKnown();
                    bool    mgeno = mot.isKnown();

                    // Number of alleles, homozygotes and genotypes in this sibship
                    int haplo = 0, homo = 0, diplo = 0;

                    // No. of different genotypes per allele
                    counts.Zero();

                    // In general, there should be no more than 3 genotypes per allele
                    bool too_many_genos = false;

                    for (int j = 0; j < persons[i]->sibCount; j++)
                        if (persons[i]->sibs[j]->isGenotyped(m))
                        {
                            Alleles geno = persons[i]->sibs[j]->markers[m];

                            int fat1 = fat.hasAllele(geno.one);
                            int fat2 = fat.hasAllele(geno.two);
                            int mot1 = mot.hasAllele(geno.one);
                            int mot2 = mot.hasAllele(geno.two);

                            if ((fgeno && mgeno && !((fat1 && mot2) || (fat2 && mot1))) ||
                                    (fgeno && !(fat1 || fat2)) || (mgeno && !(mot1 || mot2)))
                            {
                                printf("%s - Fam %s: Child %s [%s/%s] has ",
                                       (const char *) markerNames[m],
                                       (const char *) persons[i]->sibs[j]->famid,
                                       (const char *) persons[i]->sibs[j]->pid,
                                       (const char *) info->GetAlleleLabel(geno.one),
                                       (const char *) info->GetAlleleLabel(geno.two));

                                if (!fgeno || !mgeno)
                                    printf("%s [%s/%s]\n",
                                           fgeno ? "father" : "mother",
                                           (const char *) info->GetAlleleLabel(fgeno ? fat.one : mot.one),
                                           (const char *) info->GetAlleleLabel(fgeno ? fat.two : mot.two));
                                else
                                    printf("parents [%s/%s]*[%s/%s]\n",
                                           (const char *) info->GetAlleleLabel(fat.one),
                                           (const char *) info->GetAlleleLabel(fat.two),
                                           (const char *) info->GetAlleleLabel(mot.one),
                                           (const char *) info->GetAlleleLabel(mot.two));

                                fail = true;
                                failedFamilies[f] = true;
                            }
                            else
                            {
                                if (haplos[geno.one] != i)
                                {
                                    haplo++;
                                    haplos[geno.one] = i;
                                };
                                if (haplos[geno.two] != i)
                                {
                                    haplo++;
                                    haplos[geno.two] = i;
                                };

                                int index = geno.SequenceCoded();

                                if (genos[index] != i)
                                {
                                    genos[index] = i;
                                    diplo++;
                                    counts[geno.one]++;
                                    if (geno.isHomozygous())
                                        homo++;
                                    else
                                        counts[geno.two]++;
                                    if (counts[geno.one] > 2) too_many_genos = true;
                                    if (counts[geno.two] > 2) too_many_genos = true;
                                }
                            }
                        }

                    if (fgeno)
                    {
                        if (haplos[fat.one] != i)
                        {
                            haplo++;
                            haplos[fat.one] = i;
                        }
                        if (haplos[fat.two] != i)
                        {
                            haplo++;
                            haplos[fat.two] = i;
                        }
                        homo += fat.isHomozygous();
                    }

                    if (mgeno)
                    {
                        if (haplos[mot.one] != i)
                        {
                            haplo++;
                            haplos[mot.one] = i;
                        }
                        if (haplos[mot.two] != i)
                        {
                            haplo++;
                            haplos[mot.two] = i;
                        }
                        homo += mot.isHomozygous();
                    }

                    if (diplo > 4 || haplo + homo > 4 || (haplo == 4 && too_many_genos))
                    {
                        printf("%s - Fam %s: ",
                               (const char *) markerNames[m],
                               (const char *) persons[i]->famid);
                        if (persons[i]->father->markers[m].isKnown())
                            printf("Father %s [%s/%s] has children [",
                                   (const char *) persons[i]->father->pid,
                                   (const char *) info->GetAlleleLabel(fat.one),
                                   (const char *) info->GetAlleleLabel(fat.two));
                        else if (persons[i]->mother->markers[m].isKnown())
                            printf("Mother %s [%s/%s] has children [",
                                   (const char *) persons[i]->mother->pid,
                                   (const char *) info->GetAlleleLabel(mot.one),
                                   (const char *) info->GetAlleleLabel(mot.two));
                        else
                            printf("Couple %s * %s has children [",
                                   (const char *) persons[i]->mother->pid,
                                   (const char *) persons[i]->father->pid);

                        for (int j = 0; j < persons[i]->sibCount; j++)
                            printf("%s%s/%s", j == 0 ? "" : " ",
                                   (const char *) info->GetAlleleLabel(persons[i]->sibs[j]->markers[m].one),
                                   (const char *) info->GetAlleleLabel(persons[i]->sibs[j]->markers[m].two));
                        printf("]\n");

                        fail = true;
                        failedFamilies[f] = true;
                    }
                }

        for (int f = 0; f < familyCount; f++)
            if (!failedFamilies[f] &&
                    (families[f]->count > families[f]->founders + 1) &&
                    !families[f]->isNuclear())
                fail |= !GenotypeList::EliminateGenotypes(*this, families[f], m);
    }

    if (fail)
        printf("\nMendelian inheritance errors detected\n");

    return fail;
}

bool Pedigree::SexLinkedCheck()
{
    bool fail = false;

    // Keep track of what families fail the basic inheritance check,
    // so that we can run later run genotype elimination check on the remainder
    IntArray failedFamilies(familyCount);

    // For each marker ...
    for (int m = 0; m < markerCount; m++)
    {
        MarkerInfo * info = GetMarkerInfo(m);

        failedFamilies.Zero();

        // Check for homozygous males
        for (int f = 0; f < familyCount; f++)
            for (int i = families[f]->first; i <= families[f]->last; i++)
                if (persons[i]->sex == SEX_MALE && persons[i]->markers[m].isKnown() &&
                        !persons[i]->markers[m].isHomozygous())
                {
                    printf("%s - Fam %s: Male %s has two X alleles [%s/%s]\n",
                           (const char *) markerNames[m],
                           (const char *) persons[i]->famid, (const char *) persons[i]->pid,
                           (const char *) info->GetAlleleLabel(persons[i]->markers[m].one),
                           (const char *) info->GetAlleleLabel(persons[i]->markers[m].two));

                    // Wipe this genotype so we don't get cascading errors below
                    persons[i]->markers[m][0] = persons[i]->markers[m][1] = 0;

                    fail = true;
                    failedFamilies[f] = true;
                }

        // Check full sibships for errors
        // TODO -- We could do better by grouping male half-sibs
        for (int f = 0; f < familyCount; f++)
            for (int i = families[f]->first; i <= families[f]->last; i++)
                if (!persons[i]->isFounder() && persons[i]->sibs[0] == persons[i])
                {
                    // This loop runs once per sibship
                    Alleles fat = persons[i]->father->markers[m];
                    Alleles mot = persons[i]->mother->markers[m];

                    bool fgeno = fat.isKnown();
                    bool mgeno = mot.isKnown();

                    Alleles inferred_mother = mot;
                    Alleles first_sister;
                    Alleles inferred_father;

                    bool mother_ok = true;

                    int sisters = 0;

                    for (int j = 0; j < persons[i]->sibCount; j++)
                        if (persons[i]->sibs[j]->isGenotyped(m))
                        {
                            Alleles geno = persons[i]->sibs[j]->markers[m];

                            bool fat1 = fat.hasAllele(geno.one);
                            bool fat2 = fat.hasAllele(geno.two);
                            bool mot1 = mot.hasAllele(geno.one);
                            bool mot2 = mot.hasAllele(geno.two);

                            int sex = persons[i]->sibs[j]->sex;

                            if (sex == SEX_MALE)
                            {
                                if (mgeno && !mot1)
                                {
                                    printf("%s - Fam %s: Child %s [%s/Y] has mother [%s/%s]\n",
                                           (const char *) markerNames[m],
                                           (const char *) persons[i]->famid,
                                           (const char *) persons[i]->sibs[j]->pid,
                                           (const char *) info->GetAlleleLabel(geno.one),
                                           (const char *) info->GetAlleleLabel(mot.one),
                                           (const char *) info->GetAlleleLabel(mot.two));
                                    fail = true;
                                    failedFamilies[f] = true;
                                }
                                else
                                    mother_ok &= inferred_mother.AddAllele(geno.one);
                            }
                            if (sex == SEX_FEMALE)
                            {
                                if ((fgeno && mgeno && !((fat1 && mot2) || (fat2 && mot1))) ||
                                        (fgeno && !(fat1 || fat2)) || (mgeno && !(mot1 || mot2)))
                                {
                                    printf("%s - Fam %s: Child %s [%s/%s] has ",
                                           (const char *) markerNames[m],
                                           (const char *) persons[i]->famid,
                                           (const char *) persons[i]->sibs[j]->pid,
                                           (const char *) info->GetAlleleLabel(geno.one),
                                           (const char *) info->GetAlleleLabel(geno.two));

                                    if (!fgeno)
                                        printf("mother [%s/%s]\n",
                                               (const char *) info->GetAlleleLabel(mot.one),
                                               (const char *) info->GetAlleleLabel(mot.two));
                                    else if (!mgeno)
                                        printf("father [%s/Y]\n",
                                               (const char *) info->GetAlleleLabel(fat.one));
                                    else
                                        printf("parents [%s/Y]*[%s/%s]\n",
                                               (const char *) info->GetAlleleLabel(fat.one),
                                               (const char *) info->GetAlleleLabel(mot.one),
                                               (const char *) info->GetAlleleLabel(mot.two));

                                    fail = true;
                                    failedFamilies[f] = true;
                                }
                                else
                                {
                                    if (!sisters++)
                                        inferred_father = first_sister = geno;
                                    else if (first_sister != geno)
                                    {
                                        inferred_father.Intersect(geno);

                                        mother_ok &= inferred_mother.AddAllele(
                                                         geno.otherAllele(inferred_father.one));
                                        mother_ok &= inferred_mother.AddAllele(
                                                         first_sister.otherAllele(inferred_father.one));
                                    }

                                    if (!fgeno && (mot1 ^ mot2))
                                        inferred_father.Intersect(mot1 ? geno.two : geno.one);

                                    if (!mgeno && (fat1 ^ fat2))
                                        mother_ok &= inferred_mother.AddAllele(fat1 ? geno.two : geno.one);
                                }
                            }
                        }

                    if (!mother_ok || (sisters && !inferred_father.isKnown()))
                    {
                        printf("%s - Fam %s: ",
                               (const char *) markerNames[m],
                               (const char *) persons[i]->famid);
                        if (fgeno)
                            printf("Father %s [%s/Y] has children [",
                                   (const char *) persons[i]->father->pid,
                                   (const char *) info->GetAlleleLabel(fat.one));
                        else if (mgeno)
                            printf("Mother %s [%s/%s] has children [",
                                   (const char *) persons[i]->mother->pid,
                                   (const char *) info->GetAlleleLabel(mot.one),
                                   (const char *) info->GetAlleleLabel(mot.two));
                        else
                            printf("Couple %s * %s has children [",
                                   (const char *) persons[i]->mother->pid,
                                   (const char *) persons[i]->father->pid);

                        for (int j = 0; j < persons[i]->sibCount; j++)
                            printf(
                                persons[i]->sibs[j]->sex == SEX_MALE ? "%s%s/Y" : "%s%s/%s",
                                j == 0 ? "" : " ",
                                (const char *) info->GetAlleleLabel(persons[i]->sibs[j]->markers[m].one),
                                (const char *) info->GetAlleleLabel(persons[i]->sibs[j]->markers[m].two));
                        printf("]\n");
                        fail = true;
                        failedFamilies[f] = true;
                    }
                }

        for (int f = 0; f < familyCount; f++)
            if (!failedFamilies[f] &&
                    (families[f]->count > families[f]->founders + 1) &&
                    !families[f]->isNuclear())
                fail |= !GenotypeList::EliminateGenotypes(*this, families[f], m);
    }

    if (fail)
        printf("\nMendelian inheritance errors detected\n");

    return fail;
}

void Pedigree::ExtractFamily(int id, Pedigree & single_fam_ped)
{
    for (int i = families[id]->first; i <= families[id]->last; i++)
        single_fam_ped.Add(*persons[i]);

    single_fam_ped.Sort();
}

void Pedigree::ExtractOnAffection(int a, Pedigree & new_ped, int target_status)
{
    for (int i = 0; i < count; i++)
        if (persons[i]->affections[a] == target_status)
            new_ped.Add(*persons[i]);
        else
        {
            Person blank_person;
            blank_person.CopyIDs(*persons[i]);
            new_ped.Add(blank_person);
        }

    new_ped.Sort();
}

void Pedigree::Filter(IntArray & filter)
{
    if (filter.Length() != count)
        error("Pedigree:Size of pedigree filter doesn't match number of persons in pedigree");

    for (int i = 0; i < count; i++)
        if (filter[i] == 1)
        {
            persons[i]->WipePhenotypes();
            persons[i]->filter = true;
        }
}

void Pedigree::AddPerson(const char * famid, const char * pid,
                         const char * fatid, const char * motid,
                         int sex, bool delay_sort)
{
    if (count == size) Grow();

    persons[count] = new Person;

    persons[count]->famid = famid;
    persons[count]->pid = pid;
    persons[count]->fatid = fatid;
    persons[count]->motid = motid;
    persons[count]->sex = sex;

    count++;

    if (!delay_sort) Sort();
}

void Pedigree::ShowMemoryInfo()
{
    unsigned int bytes = 0;

    for (int i = 0; i < count; i++)
        bytes += persons[i]->famid.BufferSize() + persons[i]->pid.BufferSize() +
                 persons[i]->fatid.BufferSize() + persons[i]->motid.BufferSize();

    bytes += count * (markerCount * sizeof(Alleles) + traitCount * sizeof(double) +
                      covariateCount * sizeof(double) + affectionCount * sizeof(char) +
                      sizeof(Person));

    printf("   %40s %s\n", "Pedigree file ...", (const char *) MemoryInfo(bytes));
}


