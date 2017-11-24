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

#include "PedigreeDescription.h"
#include "MapFunction.h"
#include "MathVector.h"
#include "Constant.h"
#include "FortranFormat.h"
#include "Error.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

PedigreeDescription::PedigreeDescription()
{
    columnCount = 0;
    mendelFormat = false;
}

PedigreeDescription::~PedigreeDescription()
{ };

PedigreeDescription & PedigreeDescription::operator = (PedigreeDescription & rhs)
{
    columnCount = rhs.columnCount;

    columns = rhs.columns;
    columnHash = rhs.columnHash;

    return *this;
};

void PedigreeDescription::Load(IFILE & input, bool warnIfLinkage)
{
    // Check if we are dealing with a linkage format data file
    String      buffer;
    StringArray tokens;

    mendelFormat = false;

    ReadLineHelper(input, buffer, tokens);
    ifrewind(input);

    if (tokens.Length() == 4 && isdigit(tokens[0][0]))
    {
        if (warnIfLinkage) printf("Data file looks like a LINKAGE format file...\n\n");
        LoadLinkageDataFile(input);
        return;
    }

    if (buffer.Length() > 18
            && (buffer.SubStr(8,8).SlowCompare("AUTOSOME") == 0 ||
                buffer.SubStr(8,8).SlowCompare("X-LINKED") == 0)
            && (isdigit(buffer[16])  || isdigit(buffer[17]))
            && (isdigit(buffer[18])  || isdigit(buffer[19]) ||
                (buffer.Length() > 19 && isdigit(buffer[20]))))
    {
        printf("Data file looks like a MENDEL format file...\n"
               "   Activating EXPERIMENTAL support for this format\n\n");
        LoadMendelDataFile(input);
        return;
    }

    // Reset things
    ifrewind(input);
    int done = 0;
    int line = 0;

    columns.Clear();
    columnHash.Clear();
    columnCount = 0;

    while (!ifeof(input) && !done)
    {
        int   i;

        buffer.ReadLine(input);
        line++;

        tokens.Clear();
        tokens.AddTokens(buffer, WHITESPACE);

        if (tokens.Length() < 1) continue;

        if (tokens.Length() == 1)
            error("Problem reading data file:\n"
                  "Item #%d (of type %s) has no name.",
                  columnCount+1, (const char *) tokens[0]);

        switch (toupper(tokens[0][0]))
        {
            case 'A' :
                columnHash.Push(GetAffectionID(tokens[1]));
                columns.Push(pcAffection);
                columnCount++;
                break;
            case 'M' :
                columnHash.Push(GetMarkerID(tokens[1]));
                columns.Push(pcMarker);
                columnCount++;
                break;
            case 'T' :
                columnHash.Push(GetTraitID(tokens[1]));
                columns.Push(pcTrait);
                columnCount++;
                break;
            case 'C' :
                columnHash.Push(GetCovariateID(tokens[1]));
                columns.Push(pcCovariate);
                columnCount++;
                break;
            case '$' :
                columnHash.Push(GetStringID(tokens[1]));
                columns.Push(pcString);
                columnCount++;
                break;
            case 'S' :
                i = (int) tokens[0].SubStr(1);
                i = i > 0 ? i : 1;
                while (i--)
                {
                    columns.Push(pcSkip);
                    columnHash.Push(0);
                    columnCount++;
                }
                break;
            case 'Z' :
                columnHash.Push(0);
                columns.Push(pcZygosity);
                columnCount++;
                break;
            case 'V' :
                GetMarkerID(tokens[1]);
                break;
            case 'E' :
                done = 1;
                break;
            case 'U' :
                if (toupper(tokens[0][1]) == 'T' && toupper(tokens[0][2]) == 'C')
                {
                    int c = GetCovariateID(tokens[1]);
                    int t = GetTraitID(tokens[1]);

                    if (c >= 32767 || t >= 32767)
                        error("Internal error processing data file\n");

                    columnHash.Push(t * 32768 + c);
                    columns.Push(pcUndocumentedTraitCovariate);
                    columnCount++;
                    break;
                }
            default :
                error("Problem in data file (line %d):\n%s\n",
                      line, (const char *) buffer);
        }
    }

    columns.Push(pcEnd);
    columnHash.Push(0);
};

void PedigreeDescription::Load(const char * iFilename, bool warnIfLinkage)
{
    IFILE f = ifopen(iFilename, "rb");

    if (f == NULL)
        error(
            "The datafile %s cannot be opened\n\n"
            "Common causes for this problem are:\n"
            "  * You might not have used the correct options to specify input file names,\n"
            "    please check the program documentation for information on how to do this\n\n"
            "  * The file doesn't exist or the filename might have been misspelt\n\n"
            "  * The file exists but it is being used by another program which you will need\n"
            "    to close before continuing\n\n"
            "  * The file is larger than 2GB and you haven't compiled this application with\n"
            "    large file support.\n\n",
            iFilename);

    Load(f, warnIfLinkage);
    ifclose(f);

    filename = iFilename;
};

void PedigreeDescription::LoadMap(const char * iFilename)
{
    IFILE f = ifopen(iFilename, "rb");

    if (f == NULL)
        error(
            "The mapfile %s cannot be opened\n\n"
            "Please check that the file exists and is not being used by another program\n"
            "To find out how to set input filenames, check the documentation\n",
            iFilename);

    LoadMap(f);
    ifclose(f);
};

void PedigreeDescription::LoadMap(IFILE & input)
{
    columns.Clear();
    columnHash.Clear();
    columnCount = 0;

    int         lastposition = 0;
    String      buffer;
    StringArray tokens;

    buffer.ReadLine(input);
    tokens.AddTokens(buffer, WHITESPACE);

    while (tokens.Length() == 0 && !ifeof(input))
    {
        buffer.ReadLine(input);
        tokens.AddTokens(buffer, WHITESPACE);
    }

    if (tokens.Length() != 3)
        error("Error reading map file header, which has %d columns.\n"
              "Three columns were expected, corresponding to\n"
              "MARKER_ID, MARKER_NAME and BASE_PAIR_POSITION\n"
              "The offending header is transcribed below:\n\n"
              "%s", tokens.Length(), (const char *) buffer);
    else
        printf("Map file column labels\n"
               "  -- COLUMN 1, Expecting MARKER_ID, Read %s\n"
               "  -- COLUMN 2, Expecting MARKER_NAME, Read %s\n"
               "  -- COLUMN 3, Expection BASE_PAIR_POSITION, Read %s\n\n",
               (const char *)(tokens[0]), (const char *)(tokens[1]),
               (const char *)(tokens[2]));

    int line = 1;
    while (!ifeof(input))
    {
        int    serial;
        long   position;

        buffer.ReadLine(input);
        line++;

        tokens.Clear();
        tokens.AddTokens(buffer, WHITESPACE);

        if (tokens.Length() < 1) continue;
        if (tokens.Length() != 3)
            error("Each line in the map file should have 3 tokens, corresponding\n"
                  "to MARKER_ID, MARKER_NAME and BASE_PAIR_POSITION respectively\n"
                  "However, there are %d tokens in line %d, transcribed below:\n\n"
                  "%s", tokens.Length(), line, (const char *) buffer);

        serial = (int) tokens[0];
        if (serial != columnCount + 1)
            error("Reading Marker Index from Map File...\n"
                  "Markers should be indexed consecutively starting at 1\n"
                  "Marker %d does not fit this pattern\n", columnCount + 1);

        position = (int) tokens[2];
        if (position < lastposition)
            error("Reading Marker Position from Map File...\n"
                  "Marker position should be in base-pairs\n"
                  "and markers should be in map order\n");

        // TODO -- store marker locations somewhere!
        lastposition = position;

        columnHash.Push(GetMarkerID(tokens[1]));
        columns.Push(pcMarker);
        columnCount++;

        GetMarkerInfo(tokens[1])->position = position * 1e-8;
    }

    columns.Push(pcEnd);
    columnHash.Push(0);
};

int PedigreeDescription::CountTextColumns()
{
    int count = 0;

    for (int i = 0; i < columnCount; i++, count++)
        if (columns[i] == pcMarker)
            count++;

    return count;
}

void PedigreeDescription::LoadLinkageDataFile(const char * iFilename)
{
    IFILE f = ifopen(iFilename, "rb");

    if (f == NULL)
        error(
            "The linkage format datafile %s cannot be opened\n\n"
            "Please check that the file exists and is not being used by another program\n"
            "To find out how to set input filenames, check the documentation\n",
            iFilename);

    LoadLinkageDataFile(f);
    ifclose(f);

    filename = iFilename;
};

void PedigreeDescription::LoadLinkageDataFile(IFILE & input)
{
    columns.Clear();
    columnHash.Clear();
    columnCount = 0;

    String      buffer, label;
    StringArray tokens;

    ReadLineHelper(input, buffer, tokens);

    if (tokens.Length() != 4 || tokens[2].AsInteger() != (int) chromosomeX ||
            tokens[0].AsInteger() < 0)
        error("Cannot handle first line of data file\n\n"
              "Expecting four (4) numeric values, which correspond to:\n"
              "   num-loci   -- number of loci in the pedigree\n"
              "                 this value must be positive\n"
              "   risk-locus -- locus for which risks should be calculated\n"
              "                 this value will be ignored\n"
              "   sex-link   -- are the loci sex linked [0 - No, 1 - Yes]\n"
              "                 %s\n"
              "   program    -- which LINKAGE program do you want to use?\n"
              "                 this value will also be ignored\n\n"
              "The actual input read:\n%s\n",
              chromosomeX ? "expecting X-linked data, so this value must be ONE (1)"
              : "expecting autosomal data, so this must be ZERO (0)",
              (const char *) buffer);

    int numloci = tokens[0];

    ReadLineHelper(input, buffer, tokens);

    if (tokens.Length() != 4 ||
            tokens[0].AsInteger() != 0 ||
            tokens[3].AsInteger() != 0)
        error("Cannot handle second line of data file\n\n"
              "Expecting four (4) numeric values, which correspond to:\n"
              "   mutation-model         -- must be zero, corresponding to no mutation\n"
              "   male-mutation-rate     -- ignored\n"
              "   female-mutation-rate   -- ignored\n"
              "   linkage-disequilibrium -- must be zero, may be used in the future to\n"
              "                             read haplotype frequencies\n\n"
              "The actual input read:\n%s\n", (const char *) buffer);

    StringArray markerOrder;
    int         unknown = 0;

    ReadLineHelper(input, buffer, markerOrder);

    if (markerOrder.Length() > numloci)
        error("The third line of the data file lists marker order\n\n"
              "Although %d loci are defined [in the first line],\n"
              "this line includes %d values:\n%s\n",
              numloci, markerOrder.Length(), (const char *) buffer);

    IntArray    locus;
    bool need_blank_line = false;

    while (!ifeof(input) && numloci--)
    {
        if (ReadLineHelper(input, buffer, tokens) == 0)
            error("Linkage data file ends unexpectedly");

        if (tokens.Length() < 2)
            error("Incomplete locus information in data file\n"
                  "Information for each locus should include 2 or more fiels\n"
                  "The expected fields are:\n"
                  "   field_type  -- indicator of locus type (trait, marker,...)\n"
                  "   alleles     -- number of alleles\n"
                  "   name        -- locus name, preceded by hash (#) sign\n\n"
                  "The actual input read:\n%s\n", (const char *) buffer);

        int locus_type = (int) tokens[0];
        int alleles    = (int) tokens[1];

        String locus_name("LOCUS");
        locus_name += ++unknown;

        if (tokens.Length() > 2 && tokens[2][0] == '#')
        {
            if (tokens[2][1] != 0)
                locus_name = tokens[2].SubStr(1);
            else if (tokens.Length() > 3)
                locus_name = tokens[3];
        }

        if ((locus_type == 4 && alleles == 0) ||
                (locus_type == 4 && alleles == 1))
        {
            columnHash.Push(GetCovariateID(locus_name));
            columns.Push(pcCovariate);
            columnCount++;
            continue;
        }

        if (locus_type == 0 && alleles == 0)
        {
            columnHash.Push(GetTraitID(locus_name));
            columns.Push(pcTrait);
            columnCount++;
            continue;
        }

        if (ReadLineHelper(input, buffer, tokens) != alleles)
            error("Expecting %d allele frequencies, but input has %d columns:\n"
                  "%s\n", alleles, tokens.Length(), (const char *) buffer);

        Vector frequencies(alleles + 1);

        frequencies[0] = 0.0;
        for (int i = 1; i <= alleles; i++)
            frequencies[i] = (double) tokens[i - 1];

        double sum = frequencies.Sum();

        if (sum <= 0.0)
            error("Allele frequencies at %s sum to %f, which doesn't make sense\n",
                  (const char *) locus_name, sum);

        if (fabs(sum - 1.0) > 1.2e-5)
        {
            printf("Allele frequencies at %s sum to %f, adjusted to 1.0\n",
                   (const char *) locus_name, sum);
            need_blank_line = true;
        }

        if (sum != 1.0)
            frequencies *= 1.0 / sum;

        switch (locus_type)
        {
            case 1 :
            {
                // Affection
                columnHash.Push(GetAffectionID(locus_name));
                columns.Push(pcAffection);
                columnCount++;

                // Read number of liability classes
                if (ReadLineHelper(input, buffer, tokens) == 0)
                    error("Linkage data file ends unexpectedly\n");

                // Skip liability class data
                int classes = tokens[0];
                if (classes > 1)
                {
                    columnHash.Push(0);
                    columns.Push(pcSkip);
                    columnCount++;
                }

                // Separate liability class rows for males and females for X-linked data
                if (chromosomeX) classes *= 2;

                while (classes--)
                    if (ReadLineHelper(input, buffer, tokens) == 0)
                        error("Linkage data file ends unexpectedly\n");

                // Ignore map location for quantitative variables
                locus.Push(-1);
            }
            break;
            case 3 :
            {
                columnHash.Push(GetMarkerID(locus_name));
                columns.Push(pcMarker);
                columnCount++;

                // Store allele frequencies
                MarkerInfo * info = GetMarkerInfo(locus_name);

                info->freq = frequencies;

                // Initialize allele labels
                info->alleleLabels.Clear();
                for (int i = 0; i < frequencies.Length(); i++)
                    info->alleleLabels.Push(label = i);
                info->IndexAlleles();

                // Store marker id, so that we can track map location
                locus.Push(GetMarkerID(locus_name));
            }
            break;
            case 0 :
            {
                // Read number of quantitative variables
                if (ReadLineHelper(input, buffer, tokens) == 0)
                    error("Linkage data file ends unexpectedly\n");

                // Add each quantitative variable to pedigree
                // Discard information on means
                for (int vars = tokens[0], i = 0; i < vars; i++)
                {
                    if (ReadLineHelper(input, buffer, tokens) == 0)
                        error("Linkage data file ends unexpectedly\n");

                    String trait_name(locus_name);

                    if (i)
                    {
                        trait_name += ".";
                        trait_name += i + 1;
                    }

                    columnHash.Push(GetTraitID(trait_name));
                    columns.Push(pcTrait);
                    columnCount++;
                }

                // Skip var-covar matrix
                if (ReadLineHelper(input, buffer, tokens) == 0)
                    error("Linkage data file ends unexpectedly\n");

                // Skip heterozygote scaling factor for var-covar matrix
                if (ReadLineHelper(input, buffer, tokens) == 0)
                    error("Linkage data file ends unexpectedly\n");

                // Ignore map location for quantitative variables
                locus.Push(-1);
            }
            break;
            case 2 :
                error("The data file includes binary factors\n"
                      "Regretably, loci of this type are not supported\n\n");
                break;
            default :
                error("Unsupported locus type [%d] in data file", locus_type);
                break;
        }
    }

    if (need_blank_line) printf("\n");

    columns.Push(pcEnd);
    columnHash.Push(0);

    ReadLineHelper(input, buffer, tokens);
    int sexDifference = tokens.Length() ? tokens[0].AsInteger() : -1;
    if (tokens.Length() != 2 ||
            (sexDifference != 0 && sexDifference != 2) ||
            tokens[1].AsInteger() != 0)
        error("Error retrieving recombination information\n\n"
              "Expecting two (2) numeric values, which correspond to:\n"
              "   sex-difference   -- must be zero (no difference) or two (sex specific recombination)\n"
              "   map-function     -- must be zero, that is, no interference\n"
              "The actual input read:\n%s\n", (const char *) buffer);

    Vector distances[2];
    bool distance_in_centimorgans = false;

    for (int r = 0; r <= sexDifference; r += 2)
    {
        ReadLineHelper(input, buffer, tokens);
        if (tokens.Length() != markerOrder.Length() - 1)
            error("Error retrieving recombination information\n\n"
                  "Expecting %d recombination fractions (current map includes %d loci)\n"
                  "Instead the following line was input:\n%s\n",
                  markerOrder.Length() - 1, markerOrder.Length(), (const char *) buffer);

        distances[r >> 1].Dimension(tokens.Length());
        for (int i = 0; i < tokens.Length(); i++)
            distances[r >> 1][i] = (double) tokens[i];

        if (distances[r >> 1].Min() < 0.0)
            error("Linkage datafile specifies negative recombination fractions");

        bool centimorgans = distances[r >> 1].Max() > 0.5;

        if (centimorgans && !distance_in_centimorgans)
            printf("  Some recombination fractions in datafile are greater than 0.5,\n"
                   "  so recombination fractions will be interpreted as cM distances\n\n");

        distance_in_centimorgans |= centimorgans;
    }

    double position = 0.0, positionMale = 0.0;

    for (int i = 0, moving = false; i < markerOrder.Length(); i++)
    {
        int m = markerOrder[i].AsInteger() - 1;

        if (m < 0 || m >= locus.Length())
            error("The marker order in the linkage datafile is invalid\n");

        m = locus[m];

        if (m != -1)
        {
            MarkerInfo * info = GetMarkerInfo(m);
            info->chromosome = chromosomeX ? 9999 : 0;

            if (sexDifference == 2)
                info->position = (position + positionMale) * 0.5,
                                 info->positionFemale = position,
                                                        info->positionMale = positionMale;
            else
                info->position = info->positionMale = info->positionFemale = position;

            moving = true;
        }

        if (i < markerOrder.Length() - 1 && moving)
            position += distance_in_centimorgans ?
                        0.01 * distances[0][i] : RecombinationToDistance(distances[0][i]);

        if (sexDifference == 2 && i < markerOrder.Length() - 1 && moving)
            positionMale += distance_in_centimorgans ?
                            0.01 * distances[1][i] : RecombinationToDistance(distances[1][i]);
    }
}

int PedigreeDescription::ReadLineHelper(IFILE & input,
                                        String & buffer,
                                        StringArray & tokens)
{
    do
    {
        // Read Line
        buffer.ReadLine(input);
        buffer.Trim();

        // Strip comments marked with >>
        int pos = buffer.FastFind(">>");
        if (pos == -1) pos = buffer.FastFind("<<");
        if (pos == -1) pos = buffer.Length() + 1;
        if (buffer[0] == '#') pos = 0;

        // Find space/tab delimited tokens
        tokens.Clear();
        tokens.AddTokens(buffer.Left(pos - 1), WHITESPACE);

    }
    while (tokens.Length() == 0 && !ifeof(input));

    return tokens.Length();
}

void PedigreeDescription::LoadMendelDataFile(const char * iFilename)
{
    IFILE f = ifopen(iFilename, "rb");

    if (f == NULL)
        error(
            "The MENDEL format datafile %s cannot be opened\n\n"
            "Please check that the file exists and is not being used by another program\n"
            "To find out how to set input filenames, check the documentation\n",
            iFilename);

    LoadMendelDataFile(f);
    ifclose(f);
};

void PedigreeDescription::LoadMendelDataFile(IFILE & file)
{
    // Processes mendel format file
    mendelFormat = true;

    // Codominant markers are mapped to markers
    // Non-codominant markers are mapped into multiple "affection status"
    // (Y/N) variables
    columns.Clear();
    columnHash.Clear();
    columnCount = 0;

    FortranFormat parser;

    // Variables for storing parsed input
    String locusName;
    String locusType;
    String alleleLabel;
    String alleleFreq;
    String phenotype;
    String genotype;
    int phenoCount;
    int alleleCount;

    while (!ifeof(file))
    {
        // Cycle through headers for each locus
        parser.SetInputFile(file);
        parser.SetFormat("(2A8,I2,I3)");

        // After retrieving locus name, check that we haven't tried to
        // read past the end-of-file
        parser.GetNextField(locusName);
        parser.GetNextField(locusType);
        alleleCount = parser.GetNextInteger();
        phenoCount = parser.GetNextInteger();

        if (locusName.IsEmpty() && locusType.IsEmpty() && alleleCount == 0 &&
                phenoCount == 0 && ifeof(file))
            break;

        // Only recognize autosomal and x-linked loci
        if (locusType.Compare("AUTOSOME") != 0 && locusType.Compare("X-LINKED"))
            error("Unrecognized locus type '%s' in Mendel data file\n\n"
                  "Recognized locus types are \"AUTOSOME\" and \"X-LINKED\".",
                  (const char *) locusType);

        if (locusType.Compare("AUTOSOME") == 0 && chromosomeX)
            error("The data file indicates that locus %s is AUTOSOMAL, but\n"
                  "X-LINKED loci were expected as input\n",
                  (const char *) locusName);

        if (locusType.Compare("X-LINKED") == 0 && !chromosomeX)
            error("The data file indicates that locus %s is X-LINKED, but\n"
                  "AUTOSOMAL loci were expected as input\n",
                  (const char *) locusName);

        if (locusName.IsEmpty())
            error("Blank locus name encountered in data file\n");

        if (phenoCount == 0)
        {
            // Co-dominant marker
            columns.Push(pcMarker);
            columnHash.Push(GetMarkerID(locusName));
            columnCount++;

            // Update marker info with allele labels and frequencies
            MarkerInfo * info = GetMarkerInfo(locusName);

            info->alleleLabels.Clear();
            info->alleleLabels.Push("");
            info->freq.Clear();

            parser.SetFormat("(2A8)");

            // Mendel allows allele names to be specified with frequencies
            // left blank
            for (int i = 0; i < alleleCount; i++)
            {
                parser.GetNextField(alleleLabel);
                parser.GetNextField(alleleFreq);

                if (alleleLabel.IsEmpty())
                    error("Locus %s is missing allele label for allele #%d\n",
                          (const char *) locusName, i+1);

                info->alleleLabels.Push(alleleLabel);

                if (!alleleFreq.IsEmpty())
                {
                    if (info->freq.Length() == 0)
                        info->freq.Push(0.0);

                    info->freq.Push(alleleFreq.AsDouble());
                }
            }
            info->IndexAlleles();

            if (info->alleleLabels.Length() != info->freq.Length() &&
                    info->freq.Length() != 0)
                error("Locus %s is missing allele frequency information for %d alleles\n",
                      (const char *) locusName,
                      info->alleleLabels.Length() - info->freq.Length());
        }
        else
        {
            // Non-codominant marker, which we decompose into multiple traits...
            parser.SetFormat("(2A8)");

            // First skip allele frequency information
            for (int i = 0; i < alleleCount; i++)
            {
                parser.GetNextField(alleleLabel);
                parser.GetNextField(alleleFreq);
            }

            // Then read in each phenotype
            for (int i = 0; i < alleleCount; i++)
            {
                parser.SetFormat("(A8,I3)");
                parser.GetNextField(phenotype);
                int genoCount = parser.GetNextInteger();

                parser.SetFormat("(A17)");
                for (int j = 0; j < genoCount; j++)
                    parser.GetNextField(genotype);

                columns.Push(pcAffection);
                columnHash.Push(GetAffectionID(locusName + "->" + phenotype));
                columnCount++;
            }
        }
    }

    columns.Push(pcEnd);
    columnHash.Push(0);
}

int PedigreeDescription::CountColumns(int type)
{
    int count = 0;

    for (int i = 0; i < columns.Length(); i++)
        if (columns[i] == type)
            count++;

    return count;
}

const char * PedigreeDescription::ColumnSummary(String & string)
{
    string.Clear();
    UpdateSummary(string, pcMarker, " markers [x2 cols]");
    UpdateSummary(string, pcTrait, " traits");
    UpdateSummary(string, pcAffection, " discrete traits");
    UpdateSummary(string, pcCovariate, " covariates");
    UpdateSummary(string, pcString, " strings");
    UpdateSummary(string, pcZygosity, " zygosity");
    UpdateSummary(string, pcSkip, " skipped");
    return string;
}

void PedigreeDescription::UpdateSummary(String & string, int type, const char * label)
{
    int count = CountColumns(type);

    if (count)
    {
        if (string.Length())
            string += ", ";
        string += count;
        string += label;
    }
}


void PedigreeDescription::AddMarkerColumn(const char * markerName)
{
    if (columns.Last() == pcEnd)
    {
        columns.Pop();
        columnHash.Pop();
    }

    columnHash.Push(GetMarkerID(markerName));
    columns.Push(pcMarker);
    columnCount++;
}

void PedigreeDescription::AddCovariateColumn(const char * covariateName)
{
    if (columns.Last() == pcEnd)
    {
        columns.Pop();
        columnHash.Pop();
    }

    columnHash.Push(GetCovariateID(covariateName));
    columns.Push(pcCovariate);
    columnCount++;
}

void PedigreeDescription::AddTraitColumn(const char * traitName)
{
    if (columns.Last() == pcEnd)
    {
        columns.Pop();
        columnHash.Pop();
    }

    columnHash.Push(GetCovariateID(traitName));
    columns.Push(pcTrait);
    columnCount++;
}

void PedigreeDescription::AddAffectionColumn(const char * affectionName)
{
    if (columns.Last() == pcEnd)
    {
        columns.Pop();
        columnHash.Pop();
    }

    columnHash.Push(GetAffectionID(affectionName));
    columns.Push(pcAffection);
    columnCount++;
}

void PedigreeDescription::AddStringColumn(const char * stringName)
{
    if (columns.Last() == pcEnd)
    {
        columns.Pop();
        columnHash.Pop();
    }

    columnHash.Push(GetStringID(stringName));
    columns.Push(pcString);
    columnCount++;
}

void PedigreeDescription::AddZygosityColumn()
{
    if (columns.Last() == pcEnd)
    {
        columns.Pop();
        columnHash.Pop();
    }

    columnHash.Push(0);
    columns.Push(pcZygosity);
    columnCount++;
}

void PedigreeDescription::AddSkippedColumn()
{
    if (columns.Last() == pcEnd)
    {
        columns.Pop();
        columnHash.Pop();
    }

    columnHash.Push(0);
    columns.Push(pcSkip);
    columnCount++;
}

