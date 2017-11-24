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
#include "FortranFormat.h"
#include "Error.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

void Pedigree::Prepare(IFILE & input)
{
    pd.Load(input);
}

void Pedigree::Load(IFILE & input)
{
    if (pd.mendelFormat)
    {
        LoadMendel(input);
        return;
    }

    int sexCovariate = sexAsCovariate ? GetCovariateID("sex") : -1;

    int textCols = pd.CountTextColumns() + 5;
    int oldCount = count;
    bool warn    = true;
    int line     = 0;

    String      buffer;
    StringArray tokens;

    while (!ifeof(input))
    {
        int field = 0;

        buffer.ReadLine(input);

        tokens.Clear();
        tokens.AddTokens(buffer, WHITESPACE);

        if (tokens.Length() == 0) continue;
        if (tokens[0].SlowCompare("end") == 0) break;

        line++;

        if (tokens.Length() < textCols)
        {
            if (buffer.Length() > 79)
            {
                buffer.SetLength(75);
                buffer += " ...";
            }

            String description;

            pd.ColumnSummary(description);
            error("Loading Pedigree...\n\n"
                  "Expecting %d columns (%s),\n"
                  "but read only %d columns in line %d.\n\n"
                  "The problem line is transcribed below:\n%s\n",
                  textCols, (const char *) description,
                  tokens.Length(), line, (const char  *) buffer);
        }

        if (tokens.Length() > textCols && warn && textCols > 5)
        {
            pd.ColumnSummary(buffer);
            printf("WARNING -- Trailing columns in pedigree file will be ignored\n"
                   "  Expecting %d data columns (%s)\n"
                   "  However line %d, for example, has %d data columns\n\n",
                   textCols - 5, (const char *) buffer, line, tokens.Length() - 5);
            warn = false;
        }

        Person * p;

        // create a new person if necessary
        if (oldCount==0 || (p = FindPerson(tokens[0], tokens[1], oldCount))==NULL)
        {
            if (count == size) Grow();

            p = persons[count++] = new Person;
        }

        p->famid = tokens[field++];         // famid
        p->pid = tokens[field++];           // pid
        p->fatid = tokens[field++];         // fatid
        p->motid = tokens[field++];         // motid

        bool failure = false;
        p->sex = TranslateSexCode(tokens[field++], failure);
        if (failure)
            error("Can't interpret the sex of individual #%d\n"
                  "Family: %s  Individual: %s  Sex Code: %s", count,
                  (const char *) p->famid, (const char *) p->pid,
                  (const char *) tokens[field-1]);

        if (sexAsCovariate)
        {
            if (p->sex)
                p->covariates[sexCovariate] = p->sex;
            else
                p->covariates[sexCovariate] = _NAN_;
        }

        for (int col = 0; col < pd.columnCount; col++)
            switch (pd.columns[col])
            {
                case pcAffection :
                {
                    int a = pd.columnHash[col];
                    int new_status;

                    const char * affection = tokens[field++];

                    switch (toupper(affection[0]))
                    {
                        case '1' :
                        case 'N' :
                        case 'U' :
                            new_status = 1;
                            break;
                        case '2' :
                        case 'D' :
                        case 'A' :
                        case 'Y' :
                            new_status = 2;
                            break;
                        default :
                            new_status = atoi(affection);
                            if (new_status < 0 || new_status > 2)
                                error("Incorrect formating for affection status "
                                      "Col %d, Affection %s\n"
                                      "Family: %s  Individual: %s  Status: %s",
                                      col, (const char *) affectionNames[a],
                                      (const char *) p->famid, (const char *) p->pid,
                                      affection);
                    }
                    if (new_status != 0 && p->affections[a] != 0 &&
                            new_status != p->affections[a])
                        error("Conflict with previous affection status - "
                              "Col %d, Affection %s\n"
                              "Family: %s  Individual: %s  Old: %d New: %d",
                              col, (const char *) affectionNames[a],
                              (const char *) p->famid, (const char *) p->pid,
                              p->affections[a], new_status);
                    if (new_status) p->affections[a] = new_status;
                    break;
                }
                case pcMarker :
                {
                    int m = pd.columnHash[col];

                    Alleles new_genotype;

                    new_genotype[0] = LoadAllele(m, tokens[field++]);
                    new_genotype[1] = LoadAllele(m, tokens[field++]);

                    if (p->markers[m].isKnown() && new_genotype.isKnown() &&
                            new_genotype != p->markers[m])
                    {
                        MarkerInfo * info = GetMarkerInfo(m);

                        error("Conflict with previous genotype - Col %d, Marker %s\n"
                              "Family: %s  Individual: %s  Old: %s/%s New: %s/%s",
                              col, (const char *) markerNames[m],
                              (const char *) p->famid, (const char *) p->pid,
                              (const char *) info->GetAlleleLabel(p->markers[m][0]),
                              (const char *) info->GetAlleleLabel(p->markers[m][1]),
                              (const char *) info->GetAlleleLabel(new_genotype[0]),
                              (const char *) info->GetAlleleLabel(new_genotype[1]));
                    }

                    if (new_genotype.isKnown()) p->markers[m] = new_genotype;
                    break;
                }
                case pcTrait :
                case pcUndocumentedTraitCovariate :
                {
                    int t = pd.columnHash[col];
                    double new_pheno = _NAN_;

                    if (pd.columns[col] == pcUndocumentedTraitCovariate)
                        t = t / 32768;

                    const char * value = tokens[field++];
                    char * flag = NULL;

                    if (missing == (const char *) NULL || strcmp(value, missing) != 0)
                        new_pheno = strtod(value, &flag);
                    if (flag != NULL && *flag) new_pheno = _NAN_;

                    if (p->traits[t] != _NAN_ && new_pheno != _NAN_ &&
                            new_pheno != p->traits[t])
                        error("Conflict with previous phenotype - Col %d, Trait %s\n"
                              "Family: %s  Individual: %s  Old: %f New: %f",
                              col, (const char *) traitNames[t],
                              (const char *) p->famid, (const char *) p->pid,
                              p->traits[t], new_pheno);

                    if (new_pheno != _NAN_) p->traits[t] = new_pheno;
                    if (pd.columns[col] == pcTrait) break;
                }
                case pcCovariate :
                {
                    int c = pd.columnHash[col];
                    double new_covar = _NAN_;

                    if (pd.columns[col] == pcUndocumentedTraitCovariate)
                    {
                        c = c % 32768;
                        field--;
                    }

                    const char * value = tokens[field++];
                    char * flag = NULL;

                    if (missing == (const char *) NULL || strcmp(value, missing) != 0)
                        new_covar = strtod(value, &flag);
                    if (flag != NULL && *flag) new_covar = _NAN_;

                    if (p->covariates[c] != _NAN_ && new_covar != _NAN_ &&
                            new_covar != p->covariates[c])
                        error("Conflict with previous value - Col %d, Covariate %s\n"
                              "Family: %s  Individual: %s  Old: %f New: %f",
                              col, (const char *) covariateNames[c],
                              (const char *) p->famid, (const char *) p->pid,
                              p->covariates[c], new_covar);

                    if (new_covar != _NAN_) p->covariates[c] = new_covar;
                    break;
                }
                case pcString :
                {
                    int c = pd.columnHash[col];

                    if (!p->strings[c].IsEmpty() && p->strings[c] != tokens[field])
                        error("Conflict with previous value - Col %d, String %s\n"
                              "Family: %s  Individual: %s  Old: %s New: %s",
                              col, (const char *) stringNames[c],
                              (const char *) p->famid, (const char *) p->pid,
                              (const char *) p->strings[c], (const char *) tokens[field]);

                    p->strings[c] = tokens[field++];

                    break;
                }
                case pcSkip :
                    field++;
                    break;
                case pcZygosity :
                {
                    int new_zygosity;

                    const char * zygosity = tokens[field++];

                    switch (zygosity[0])
                    {
                        case 'D' :
                        case 'd' :
                            new_zygosity = 2;
                            break;
                        case 'M' :
                        case 'm' :
                            new_zygosity = 1;
                            break;
                        default :
                            new_zygosity = atoi(zygosity);
                    }
                    if (p->zygosity != 0 && new_zygosity != p->zygosity)
                        error("Conflict with previous zygosity - "
                              "Column %d in pedigree\n"
                              "Family: %s  Individual: %s  Old: %d New: %d\n",
                              col, (const char *) p->famid, (const char *) p->pid,
                              p->zygosity, new_zygosity);
                    p->zygosity = new_zygosity;
                    break;
                }
                case pcEnd :
                    break;
                default :
                    error("Inconsistent Pedigree Description -- Internal Error");
            }
    }

    Sort();
}

void Pedigree::LoadMendel(IFILE & input)
{
    // First, retrieve the two format statements from file
    String familyHeader;
    String individualRecord;

    familyHeader.ReadLine(input);
    individualRecord.ReadLine(input);

    // Then create two FORTRAN input streams...
    // One will be used for retrieving family labels and sizes, the other
    // will be used for individual information
    FortranFormat headers, records;

    headers.SetInputFile(input);
    headers.SetFormat(familyHeader);

    records.SetInputFile(input);
    records.SetFormat(individualRecord);

    // Storage for key pieces of information
    String famid;
    String phenotype;
    String affectionCode;
    String affectionStem;
    int    familySize;

    String allele1, allele2;

    int sexCovariate = sexAsCovariate ? GetCovariateID("sex") : -1;

    while (!ifeof(input))
    {
        if (count == size)
            Grow();

        // Retrieve header for next family
        familySize = headers.GetNextInteger();
        headers.GetNextField(famid);
        headers.Flush();

        if (famid.IsEmpty())
        {
            if (ifeof(input) && familySize == 0)
                break;
            else
                error("Blank family id encountered\n");
        }

        // Retrieve each individual in the family
        for (int i = 0; i < familySize; i++)
        {
            Person * p = persons[count++] = new Person;

            // Retrieve basic pedigree structure
            p->famid = famid;
            records.GetNextField(p->pid);
            records.GetNextField(p->fatid);
            records.GetNextField(p->motid);

            if (p->pid.IsEmpty())
                error("No unique identifier for individual #%d in family %s\n",
                      i + 1, (const char *) famid);

            if (p->pid.Compare(".") == 0)
                error("Family %s has an individual named '.', but this code is\n"
                      "reserved to indicate missing parents\n");

            if (p->fatid.IsEmpty()) p->fatid = ".";
            if (p->motid.IsEmpty()) p->motid = ".";

            // Retrieve and decode sex code
            char sex = records.GetNextCharacter();

            switch (sex)
            {
                case '0' :
                case 'x' :
                case 'X' :
                case '?' :
                case 0 :
                    p->sex = 0;
                    break;
                case '1' :
                case 'm' :
                case 'M' :
                    p->sex = 1;
                    break;
                case '2' :
                case 'f' :
                case 'F' :
                    p->sex = 2;
                    break;
                default :
                    error("Can't interpret the sex of individual #%d\n"
                          "Family: %s  Individual: %s  Sex Code: %s", count,
                          (const char *) p->famid, (const char *) p->pid, sex);
            };

            if (sexAsCovariate)
            {
                if (p->sex)
                    p->covariates[sexCovariate] = p->sex;
                else
                    p->covariates[sexCovariate] = _NAN_;
            }

            // Retrieve and decode zygosity
            char zygosity = records.GetNextCharacter();

            // Mendel uses a unique character to indicate each MZ pair,
            // we use a unique odd number...
            if (zygosity)
                p->zygosity = (zygosity - ' ') * 2 - 1;

            affectionStem.Clear();
            for (int col = 0; col < pd.columnCount; col++)
                switch (pd.columns[col])
                {
                    case pcAffection :
                    {
                        int a = pd.columnHash[col];

                        // We expand each Mendel non-codominant trait into multiple
                        // affection status column... First, if this is  not a
                        // continuation of a previous expansion we first retrieve
                        // and encode the affection status.
                        if (affectionStem.Length() == 0 ||
                                affectionNames[a].CompareToStem(affectionStem) != 0)
                        {
                            affectionStem.Copy(affectionNames[a], 0, affectionNames[a].FindChar('>') + 1);
                            records.GetNextField(phenotype);
                            affectionCode = affectionStem + phenotype;
                        }

                        // Then encode each phenotype appropriately
                        if (phenotype.IsEmpty())
                            p->affections[a] = 0;
                        else
                            p->affections[a] = affectionCode.Compare(affectionNames[a]) == 0 ? 2 : 1;

                        break;
                    }
                    case pcMarker :
                    {
                        int m = pd.columnHash[col];

                        records.GetNextField(phenotype);

                        if (phenotype.IsEmpty())
                        {
                            p->markers[m].one = p->markers[m].two = 0;
                            continue;
                        }

                        int separator = phenotype.FindChar('/');
                        if (separator == -1) separator = phenotype.FindChar('|');

                        if (separator == -1)
                            error("At marker %s, person %s in family %s has genotype %s.\n"
                                  "This genotype is not in the 'al1/al2' format.\n",
                                  (const char *) markerNames[m],
                                  (const char *) p->pid,
                                  (const char *) p->famid,
                                  (const char *) phenotype);

                        allele1.Copy(phenotype, 0, separator);
                        allele1.Trim();

                        allele2.Copy(phenotype, separator + 1, 8);
                        allele2.Trim();

                        MarkerInfo * info = GetMarkerInfo(m);

                        int one = info->alleleNumbers.Integer(allele1);

                        if (one < 0)
                        {
                            if (info->freq.Length() == 0)
                                one = info->NewAllele(allele1);
                            else
                                error("At marker %s, person %s in family %s has genotype %s.\n"
                                      "However, '%s' is not a valid allele for this marker.\n",
                                      (const char *) markerNames[m],
                                      (const char *) p->pid,
                                      (const char *) p->famid,
                                      (const char *) phenotype,
                                      (const char *) allele1);
                        }

                        int two = info->alleleNumbers.Integer(allele2);

                        if (two < 0)
                        {
                            if (info->freq.Length() == 0)
                                two = info->NewAllele(allele2);
                            else
                                error("At marker %s, person %s in family %s has genotype %s.\n"
                                      "However, '%s' is not a valid allele for this marker.\n",
                                      (const char *) markerNames[m],
                                      (const char *) p->pid,
                                      (const char *) p->famid,
                                      (const char *) phenotype,
                                      (const char *) allele2);
                        }

                        p->markers[m].one = one;
                        p->markers[m].two = two;
                        break;
                    }
                    case pcEnd :
                        break;
                    case pcTrait :
                    case pcCovariate :
                    case pcSkip :
                    case pcZygosity :
                    default:
                        error("Inconsistent Pedigree Description -- Internal Error");
                }

            records.Flush();
        }
    }

    Sort();
}

void Pedigree::Prepare(const char * filename)
{
    // Clear any previously loaded pedigree description
    if (multiPd != NULL)
        delete [] multiPd;

    multiFileCount = 1;

    // Enable multifile support
    StringArray filenames;

    filenames.AddColumns(filename, ',');

    if (filenames.Length() <= 1)
        pd.Load(filename);
    else
    {
        printf("AUTOMATIC MERGE ENABLED: Detected multiple datafile names, separated by commas...\n");

        multiPd = new PedigreeDescription[filenames.Length()];

        for (int i = 0; i < filenames.Length(); i++)
        {
            printf("  AUTOMATIC MERGE: Reading data file '%s' ...\n", (const char *) filenames[i]);
            multiPd[i].Load(filenames[i], false);
        }

        multiFileCount = filenames.Length();
    }
}

void Pedigree::Load(const char * filename, bool allowFailures)
{
    if (multiFileCount <= 1)
    {
        IFILE f = ifopen(filename, "rb");

        if (f == NULL && allowFailures)
            return;

        if (f == NULL)
            error(
                "The pedigree file %s cannot be opened\n\n"
                "Common causes for this problem are:\n"
                "  * You might not have used the correct options to specify input file names,\n"
                "    please check the program documentation for information on how to do this\n\n"
                "  * The file doesn't exist or the filename might have been misspelt\n\n"
                "  * The file exists but it is being used by another program which you will need\n"
                "    to close\n\n"
                "  * The file is larger than 2GB and you haven't compiled this application with\n"
                "    large file support.\n\n",
                filename);

        Load(f);
        ifclose(f);
    }
    else
    {
        StringArray filenames;

        filenames.AddColumns(filename, ',');

        if (filenames.Length() != multiFileCount)
            error("Different numbers of comma separated data and pedigree file names provided\n");

        for (int i = 0; i < filenames.Length(); i++)
        {
            printf("  AUTOMATIC MERGE: Datafile '%s' matched to pedigree '%s' ...\n",
                   (const char *) multiPd[i].filename, (const char *) filenames[i]);

            pd = multiPd[i];

            IFILE f = ifopen(filenames[i], "rb");

            if (f == NULL)
                error("The pedigree file '%s' cannot be opened\n\n",
                      (const char *) filenames[i]);

            Load(f);
            ifclose(f);
        }

        printf("\n");
    }
}

int Pedigree::TranslateSexCode(const char * code, bool & failure)
{
    failure = false;

    switch (code[0])
    {
        case 'x' :
        case 'X' :
        case '?' :
            return 0;
        case '1' :
        case 'm' :
        case 'M' :
            return 1;
        case '2' :
        case 'f' :
        case 'F' :
            return 2;
        default :
        {
            int result = atoi(code);

            if (result != 0 && result != 1 && result != 2)
            {
                failure = true;
                result = 0;
            }

            return result;
        }
    };
}
