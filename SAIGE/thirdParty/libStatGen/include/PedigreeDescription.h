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

#ifndef __PEDDESCRIBE_H__
#define __PEDDESCRIBE_H__

#include "PedigreeGlobals.h"
#include "PedigreePerson.h"
#include "StringArray.h"
#include "IntArray.h"

#include <stdio.h>

// Possible pedigree columns
#define  pcSkip      0
#define  pcMarker    1
#define  pcTrait     2
#define  pcAffection 3
#define  pcCovariate 4
#define  pcString    5
#define  pcZygosity  6
#define  pcEnd       7

// Undocumented pedigree column types -- not recommended
#define  pcUndocumentedTraitCovariate   1001

class PedigreeDescription : public PedigreeGlobals
{
public:
    int      columnCount;
    IntArray columns, columnHash;

    PedigreeDescription();
    ~PedigreeDescription();

    void Load(IFILE & Input, bool warnIfLinkage = false);
    void Load(const char * filename, bool warnIfLinkage = false);

    void LoadLinkageDataFile(IFILE & input);
    void LoadLinkageDataFile(const char * filename);

    void LoadMendelDataFile(IFILE & input);
    void LoadMendelDataFile(const char * filename);

    void LoadMap(IFILE & Input);
    void LoadMap(const char * filename);

    PedigreeDescription & operator = (PedigreeDescription & rhs);

    int CountTextColumns();

    // returns a string summarizing column contents
    const char * ColumnSummary(String & string);

    // Flag specifying Mendel format
    bool mendelFormat;

    String filename;

    void AddMarkerColumn(const char * markerName);
    void AddTraitColumn(const char * traitName);
    void AddAffectionColumn(const char * affectionName);
    void AddCovariateColumn(const char * covariateName);
    void AddStringColumn(const char * stringName);
    void AddZygosityColumn();
    void AddSkippedColumn();

private:
    int ReadLineHelper(IFILE & input, String & buffer, StringArray & tokens);

    int CountColumns(int type);
    void UpdateSummary(String & string, int type, const char * label);
};

#endif

