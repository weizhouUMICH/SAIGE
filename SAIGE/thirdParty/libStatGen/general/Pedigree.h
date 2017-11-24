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

#ifndef _PEDIGREE_H_
#define _PEDIGREE_H_

#include "Constant.h"

#include <stdio.h>

#include "PedigreeAlleles.h"
#include "PedigreePerson.h"
#include "PedigreeGlobals.h"
#include "PedigreeFamily.h"
#include "PedigreeDescription.h"
#include "PedigreeAlleleFreq.h"

class Pedigree : public PedigreeGlobals
{
public:
    static bool          sexAsCovariate;
    static String        missing;

    int                  size;
    int                  count;
    Person **            persons;
    int                  familyCount;
    Family **            families;
    int                  haveTwins;

    PedigreeDescription  pd;
    PedigreeDescription *multiPd;
    int                  multiFileCount;

    Pedigree();
    ~Pedigree();

    void Prepare(IFILE & input);       // Read pedigree parameters from data file
    void Load(IFILE & input);          // Read pedigree from pedigree file
    void LoadMendel(IFILE & input);    // Read pedigree in Mendel format
    void Prepare(const char * input);  // Read pedigree parameters from named file

    // Read pedigree parameters from named file, stop program on failure
    // depending on setting of allow failures
    void Load(const char * input, bool allowFailures = false);

    // I/O related utility functions
    int  TranslateSexCode(const char * code, bool & failure);

    void PrepareDichotomization();   // Register dummy affections for each trait
    int  Dichotomize(int trait, double mean = _NAN_);
    void DichotomizeAll(double mean = _NAN_);

    void WriteDataFile(FILE * output);           // Write data file
    void WritePedigreeFile(FILE * output);       // Write pedigree file
    void WriteDataFile(const char * output);     // Write named data file
    void WritePedigreeFile(const char * output); // Write named pedigree file
    void WritePerson(FILE * output, int who,     // Write a single person
                     const char * famid = NULL,              // if supplied, famid, pid,
                     const char * pid = NULL,                // fatid and motid allow a
                     const char * fatid = NULL,              // pedigree or person to
                     const char * motid = NULL);             // be renamed / restructured
    void WriteRecodedPerson(                     // Like write person, but uses
        FILE * output, int who,                 // user supplied markerInfo
        MarkerInfo ** markerInfo,               // array to recode marker
        const char * famid = NULL,              // alleles as they are written
        const char * pid = NULL,
        const char * fatid = NULL,
        const char * motid = NULL);

    void Sort();                              // Sorts the pedigree items
    Family * FindFamily(const char * famid);  // Find a family
    Person * FindPerson(const char * famid,   // Find an individual
                        const char * pid);

    // functions dealing with genetic markers
    // Counts the alleles at a marker
    int  CountAlleles(int marker);

    // Lumps together rare alleles and, depending on reorder flag,
    // sorts alleles so the most common allele has the lowest index
    void LumpAlleles(double treshold, bool reorder = true);

    // Calculate allele frequencies
    void EstimateFrequencies(int estimator, bool quiet = false);

    // shorthand operators
    Person & operator [](int i)
    {
        return *(persons[i]);
    }

    // Perform a basic inheritance check
    bool InheritanceCheck(bool abortIfInconsistent = true);
    bool AutosomalCheck();
    bool SexLinkedCheck();
    bool TwinCheck();

    // Merge twins into a single individual
    void MergeTwins();

    // Remove individuals with no data from pedigree
    void Trim(bool quiet = false, int * informative = NULL);

    // Add a single individual to a pedigree
    void AddPerson(const char * famid, const char * pid,
                   const char * fatid, const char * motid,
                   int sex, bool delay_sort = false);

    // Add all individuals in family with famid = id to new_ped
    void ExtractFamily(int id, Pedigree & new_ped);
    // Add individuals with affection status target_status for affection a to new_ped
    void ExtractOnAffection(int a, Pedigree & new_ped, int target_status = 2);

    // Remove all covariate, affection and genotype information from persons for which filter[i] = 0
    void Filter(IntArray & filter);

    // Reports memory usage for storing the pedigree
    void ShowMemoryInfo();

private:
    void Grow();
    void Add(Person & rhs);

    static int ComparePersons(const Person ** p1, const Person ** p2);
    static int CompareParents(const Person ** p1, const Person ** p2);

    void MakeSibships();
    void MakeFamilies();

    Person * FindPerson(const char * famid, const char * pid, int universe);

    void ShowTrimHeader(bool & flag);
};

#endif




