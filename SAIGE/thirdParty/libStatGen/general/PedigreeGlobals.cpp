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

#include "PedigreeGlobals.h"
#include "Sort.h"
#include "Error.h"

#include <math.h>
#include <string.h>
#include <ctype.h>

int PedigreeGlobals::traitCount = 0;
int PedigreeGlobals::affectionCount = 0;
int PedigreeGlobals::covariateCount = 0;
int PedigreeGlobals::markerCount = 0;
int PedigreeGlobals::stringCount = 0;

// If this value isn't set, all X chromosome data will be rejected
bool PedigreeGlobals::chromosomeX = false;
bool PedigreeGlobals::sexSpecificMap = false;

StringArray   PedigreeGlobals::traitNames;
StringArray   PedigreeGlobals::markerNames;
StringArray   PedigreeGlobals::covariateNames;
StringArray   PedigreeGlobals::affectionNames;
StringArray   PedigreeGlobals::stringNames;
StringIntHash PedigreeGlobals::markerLookup;
StringIntHash PedigreeGlobals::traitLookup;
StringIntHash PedigreeGlobals::affectionLookup;
StringIntHash PedigreeGlobals::covariateLookup;
StringIntHash PedigreeGlobals::stringLookup;

int PedigreeGlobals::markerInfoCount = 0;
int PedigreeGlobals::markerInfoSize = 0;

MarkerInfo ** PedigreeGlobals::markerInfo = NULL;
MarkerInfo ** PedigreeGlobals::markerInfoByInteger = NULL;
StringHash    PedigreeGlobals::markerInfoByName;

int MarkerInfo::count = 0;

int MarkerInfo::ComparePosition(MarkerInfo ** left, MarkerInfo ** right)
{
    if ((*left)->chromosome != (*right)->chromosome)
        return (*left)->chromosome - (*right)->chromosome;

    double difference = (*left)->position - (*right)->position;

    if (difference >  0.0)
        return 1;
    else if (difference == 0.0)
        return (*left)->serial - (*right)->serial;
    else
        return -1;
}

String MarkerInfo::GetAlleleLabel(int allele)
{
    if (alleleLabels.Length() > allele && alleleLabels[allele].Length())
        return alleleLabels[allele];
    else if (alleleLabels.Length() <= allele)
        alleleLabels.Dimension(allele + 1);
    return alleleLabels[allele] = allele;
}

bool MarkerInfo::AdjustFrequencies()
{
    if (freq.Length() <= 1)
    {
        freq.Clear();
        return false;
    }

    if (freq.Min() < 0.0)
        error("Locus %s has negative allele frequencies\n", (const char *) name);

    double sum = freq.Sum();

    if (sum <= 0.0)
        error("Locus %s frequencies sum to %f, which doesn't make sense\n",
              (const char *) name, sum);

    if (sum != 1.0)
        freq *= 1.0 / sum;

    if (fabs(sum - 1.0) > 1.2e-5)
    {
        printf("Locus %s frequencies sum to %f, adjusted to 1.0\n",
               (const char *) name, sum);
        return true;
    }

    return false;
}

void MarkerInfo::IndexAlleles()
{
    if (alleleLabels.Length() >= 255)
        error("Marker %s has more than 254 distinct alleles\n",
              (const char *) name);

    alleleNumbers.Clear();

    for (int i = 1; i < alleleLabels.Length(); i++)
        alleleNumbers.SetInteger(alleleLabels[i], i);
}

int MarkerInfo::NewAllele(const String & label)
{
    if (alleleLabels.Length() == 0)
        alleleLabels.Push("");

    if (alleleLabels.Length() >= 255)
        error("Marker %s has more than 254 distinct alleles\n",
              (const char *) name);

    alleleNumbers.SetInteger(label, alleleLabels.Length());
    alleleLabels.Push(label);

    return alleleLabels.Length() - 1;
}

int PedigreeGlobals::GetTraitID(const char * name)
{
    int idx = traitLookup.Integer(name);

    if (idx != -1) return idx;

    traitNames.Add(name);
    traitLookup.SetInteger(name, traitCount);
    return traitCount++;
}

int PedigreeGlobals::GetAffectionID(const char * name)
{
    int idx = affectionLookup.Integer(name);

    if (idx != -1) return idx;

    affectionNames.Add(name);
    affectionLookup.SetInteger(name, affectionCount);
    return affectionCount++;
}

int PedigreeGlobals::GetCovariateID(const char * name)
{
    int idx = covariateLookup.Integer(name);

    if (idx != -1) return idx;

    covariateNames.Add(name);
    covariateLookup.SetInteger(name, covariateCount);
    return covariateCount++;
}

int PedigreeGlobals::GetStringID(const char * name)
{
    int idx = stringLookup.Integer(name);

    if (idx != -1) return idx;

    stringNames.Add(name);
    stringLookup.SetInteger(name, stringCount);
    return stringCount++;
}

int PedigreeGlobals::GetMarkerID(const char * name)
{
    int idx = markerLookup.Integer(name);

    if (idx != -1) return idx;

    markerNames.Add(name);
    markerLookup.SetInteger(name, markerCount);

    // Grow the marker info key ...
    if (markerCount == 0)
    {
        markerInfoByInteger = new MarkerInfo * [16];

        for (int i = 0; i < 16; i++)
            markerInfoByInteger[i] = NULL;
    }
    else if ((markerCount & (markerCount - 1)) == 0 && markerCount > 15)
    {
        MarkerInfo ** newKey = new MarkerInfo * [markerCount * 2];

        for (int i = 0; i < markerCount; i++)
            newKey[i] = markerInfoByInteger[i];

        for (int i = markerCount; i < markerCount * 2; i++)
            newKey[i] = NULL;

        delete [] markerInfoByInteger;

        markerInfoByInteger = newKey;
    }

    return markerCount++;
}

MarkerInfo * PedigreeGlobals::GetMarkerInfo(String & name)
{
    MarkerInfo * info = (MarkerInfo *) markerInfoByName.Object(name);

    if (info != NULL) return info;

    info = new MarkerInfo(name);
    markerInfoByName.Add(name, info);

    if (markerInfoCount >= markerInfoSize)
        GrowMarkerInfo();

    markerInfo[markerInfoCount++] = info;

    int markerId = LookupMarker(name);
    if (markerId >= 0) markerInfoByInteger[markerId] = info;

    return info;
}

MarkerInfo * PedigreeGlobals::GetMarkerInfo(int markerId)
{
    if (markerId >= markerCount)
        error("Attempted to retrieve MarkerInfo using out-of-bounds index\n");

    if (markerInfoByInteger[markerId] != NULL)
        return markerInfoByInteger[markerId];
    else
        return GetMarkerInfo(markerNames[markerId]);
}

void PedigreeGlobals::GrowMarkerInfo()
{
    int newSize = markerInfoSize ? 2 * markerInfoSize : 32;

    MarkerInfo ** newArray = new MarkerInfo * [newSize];

    if (markerInfoSize)
    {
        memcpy(newArray, markerInfo, sizeof(MarkerInfo *) * markerInfoSize);
        delete [] markerInfo;
    }

    markerInfo = newArray;
    markerInfoSize = newSize;
}

void PedigreeGlobals::FlagMissingMarkers(IntArray & missingMarkers)
{
    int skipped_markers = 0;

    if (missingMarkers.Length())
    {
        StringArray names;

        printf("These markers couldn't be placed and won't be analysed:");

        for (int i = 0; i < missingMarkers.Length(); i++)
            names.Push(GetMarkerInfo(missingMarkers[i])->name);
        names.Sort();

        for (int i = 0, line = 80, lines = 0; i < missingMarkers.Length(); i++)
        {
            if (line + names[i].Length() + 1 > 79)
                printf("\n   "), line = 3, lines++;

            if (lines < 5)
            {
                printf("%s ", (const char *) names[i]);
                line += names[i].Length() + 1;
            }
            else
                skipped_markers++;
        }

        if (skipped_markers)
            printf("as well as %d other unlisted markers...", skipped_markers);

        printf("\n\n");
    }
}

void PedigreeGlobals::GetOrderedMarkers(IntArray & markers)
{
    if (markers.Length() == 0)
    {
        markers.Dimension(markerCount);
        markers.SetSequence(0, 1);
    }

    MarkerInfo ** subset = new MarkerInfo * [markers.Length()];

    int count = 0;
    IntArray missingMarkers;

    for (int i = 0; i < markers.Length(); i++)
    {
        MarkerInfo * info = GetMarkerInfo(markers[i]);

        if (info->chromosome != -1)
            subset[count++] = info;
        else
            missingMarkers.Push(i);
    }

    FlagMissingMarkers(missingMarkers);

    QuickSort(subset, count, sizeof(MarkerInfo *),
              COMPAREFUNC MarkerInfo::ComparePosition);

    markers.Clear();
    for (int i = 0; i < count; i++)
        markers.Push(GetMarkerID(subset[i]->name));
}

int PedigreeGlobals::SortMarkersInMapOrder(IntArray & markers, int chromosome)
{
    if (markers.Length() == 0)
    {
        markers.Dimension(markerCount);
        markers.SetSequence(0, 1);
    }

    MarkerInfo ** subset = new MarkerInfo * [markers.Length()];

    int count = 0;
    IntArray missingMarkers;

    for (int i = 0; i < markers.Length(); i++)
    {
        MarkerInfo * info = GetMarkerInfo(markers[i]);

        if (info->chromosome != -1)
            subset[count++] = info;
        else if (chromosome == -1)
            missingMarkers.Push(i);
    }

    if (chromosome == -1)
        FlagMissingMarkers(missingMarkers);

    QuickSort(subset, count, sizeof(MarkerInfo *),
              COMPAREFUNC MarkerInfo::ComparePosition);

    markers.Clear();

    int  current_chromosome = -1, next_chromosome = 0;

    for (int i = 0; i < count; i++)
        if (subset[i]->chromosome < chromosome)
            continue;
        else if (current_chromosome == -1 ||
                 subset[i]->chromosome == current_chromosome)
        {
            markers.Push(GetMarkerID(subset[i]->name));
            current_chromosome = subset[i]->chromosome;
        }
        else if (!next_chromosome)
        {
            next_chromosome = subset[i]->chromosome;
            break;
        }

    delete [] subset;

    return next_chromosome;
}

void PedigreeGlobals::VerifySexSpecificOrder()
{
    if (markerCount <= 1)
        return;

    MarkerInfo ** sortedMarkers = new MarkerInfo * [markerCount];

    for (int i = 0; i < markerCount; i++)
        sortedMarkers[i] = GetMarkerInfo(i);

    QuickSort(sortedMarkers, markerCount, sizeof(MarkerInfo *),
              COMPAREFUNC MarkerInfo::ComparePosition);

    double prev_female = sortedMarkers[0]->positionFemale;
    double prev_male = sortedMarkers[0]->positionMale;
    double curr_female, curr_male;

    int    prev_chromosome = sortedMarkers[0]->chromosome;
    int    curr_chromosome;

    for (int i = 1; i < markerCount; i++)
    {
        curr_chromosome = sortedMarkers[i]->chromosome;
        curr_female = sortedMarkers[i]->positionFemale;
        curr_male = sortedMarkers[i]->positionMale;

        if (curr_chromosome == prev_chromosome &&
                (curr_female < prev_female || curr_male < prev_male))
            error("Sex-specific and sex-averaged maps are inconsistent.\n\n"
                  "In the sex-averaged map, marker %s (%.2f cM) follows marker %s (%.2f cM).\n"
                  "In the %smale map, marker %s (%.2f cM) PRECEDES marker %s (%.2f cM).\n",
                  (const char *) sortedMarkers[i]->name,
                  sortedMarkers[i]->position * 100,
                  (const char *) sortedMarkers[i-1]->name,
                  sortedMarkers[i-1]->position * 100,
                  curr_female < prev_female ? "fe" : "",
                  (const char *) sortedMarkers[i]->name,
                  (curr_female < prev_female ? curr_female : curr_male) * 100,
                  (const char *) sortedMarkers[i-1]->name,
                  (curr_female < prev_female ? prev_female : prev_male) * 100);

        prev_chromosome = curr_chromosome;
        prev_female = curr_female;
        prev_male = curr_male;
    }

    delete [] sortedMarkers;
}

void PedigreeGlobals::LoadAlleleFrequencies(const char * filename, bool required)
{
    // This function is often called with an empty string, and not
    // all implementations of the C library like that ...
    if (filename[0] == 0)
    {
        if (required)
            error("No name provided for required allele freuquency file\n");
        else
            return;
    }

    // If we get here, the filename is not empty and things should
    // work as planned
    IFILE f = ifopen(filename, "rb");

    if (f == NULL)
    {
        if (required)
            error("Failed to open required allele frequency file '%s'",
                  (const char *) filename);
        else
            return;
    }

    LoadAlleleFrequencies(f);
    ifclose(f);
}

void PedigreeGlobals::LoadAlleleFrequencies(IFILE & input)
{
    int         done = 0;
    String      buffer;
    StringArray tokens;
    MarkerInfo *info = NULL;

    bool need_blank_line = false;
    int  allele_size, old_max, next_allele = 0; // Initialization avoids compiler warning

    while (!ifeof(input) && !done)
    {
        int   i, j;

        buffer.ReadLine(input);

        tokens.Clear();
        tokens.AddTokens(buffer, WHITESPACE);

        if (tokens.Length() < 1) continue;

        switch (toupper(tokens[0][0]))
        {
            case 'M' :
                if (tokens.Length() == 1)
                    error("Unnamed marker in allele frequency file");
                if (info != NULL)
                    need_blank_line |= info->AdjustFrequencies();
                info = GetMarkerInfo(tokens[1]);
                info->freq.Clear();
                info->freq.Push(0.0);
                next_allele = 1;
                break;
            case 'F' :
                if (info != NULL)
                    for (i = 1; i < tokens.Length(); i++)
                    {
                        buffer = next_allele++;

                        int allele = LoadAllele(info, buffer);

                        if (allele >= info->freq.Length())
                        {
                            old_max = info->freq.Length();
                            info->freq.Dimension(allele + 1);
                            for (j = old_max; j < allele; j++)
                                info->freq[j] = 0.0;
                        }

                        info->freq[allele] = tokens[i].AsDouble();
                    }
                break;
            case 'A' :
                if (info == NULL) continue;

                if (tokens.Length() != 3)
                    error("Error reading named allele frequencies for locus %s\n"
                          "Lines with named alleles should have the format\n"
                          "   A allele_label allele_frequency\n\n"
                          "But the following line was read:\n%s\n",
                          (const char *) info->name, (const char *) buffer);

                allele_size = LoadAllele(info, tokens[1]);
                next_allele = atoi(tokens[1]) + 1;

                if (allele_size < 1)
                    error("Error reading named allele frequencies for locus %s\n"
                          "An invalid allele label was encountered\n",
                          (const char *) info->name);

                if (allele_size >= info->freq.Length())
                {
                    old_max = info->freq.Length();
                    info->freq.Dimension(allele_size + 1);
                    for (i = old_max; i < allele_size; i++)
                        info->freq[i] = 0.0;
                }

                info->freq[allele_size] = tokens[2];
                break;
            case 'E' :
                done = 1;
                break;
            default :
                error("Problem in allele frequency file.\n"
                      "Lines in this file should be of two types:\n"
                      "  -- Marker name lines begin with an M\n"
                      "  -- Frequency lines begin with an F\n\n"
                      "However the following line is different:\n%s\n",
                      (const char *) buffer);
        }
    }

    if (info != NULL)
        need_blank_line |= info->AdjustFrequencies();

    if (need_blank_line) printf("\n");
}

void PedigreeGlobals::LoadMarkerMap(const char * filename, bool filter)
{
    IFILE f = ifopen(filename, "rb");
    if (f == NULL) return;
    LoadMarkerMap(f, filter);
    ifclose(f);
}

void PedigreeGlobals::LoadMarkerMap(IFILE & input, bool filter)
{
    String      buffer;
    StringArray tokens;
    bool first_pass = true;

    while (!ifeof(input))
    {
        buffer.ReadLine(input);

        tokens.Clear();
        tokens.AddTokens(buffer, WHITESPACE);

        if (tokens.Length() < 1) continue;

        if (first_pass)
        {
            sexSpecificMap = (tokens.Length() == 5);

            // if (sexSpecificMap)
            //   printf("\n  Found sex-specific map ...\n\n");

            first_pass = false;
        }

        if (tokens.Length() != 3 && !sexSpecificMap)
            error("Error reading map file\n"
                  "Each line in this file should include 3 fields:\n"
                  "CHROMOSOME, MARKER_NAME, and POSITION\n"
                  "However the following line has %d fields\n%s\n",
                  tokens.Length(), (const char *) buffer);

        if (tokens.Length() != 5 && sexSpecificMap)
            error("Error reading map file\n"
                  "Each line in this file should include 5 fields:\n\n"
                  "CHROMOSOME, MARKER_NAME, SEX_AVERAGED_POS, FEMALE_POS AND MALE_POS\n\n"
                  "However the following line has %d fields\n%s\n",
                  tokens.Length(), (const char *) buffer);

        bool previous_state = String::caseSensitive;
        String::caseSensitive = false;

        if ((tokens[0] == "CHR" || tokens[0] == "CHROMOSOME") &&
                (tokens[1] == "MARKER" || tokens[1] == "MARKER_NAME" || tokens[1] == "MRK") &&
                (tokens[2] == "KOSAMBI" || tokens[2] == "POS" || tokens[2] == "POSITION" ||
                 tokens[2] == "SEX_AVERAGED_POS" || tokens[2] == "CM" || tokens[2] == "HALDANE"))
            continue;

        String::caseSensitive = previous_state;

        if (filter)
            if (LookupMarker(tokens[1]) < 0)
                continue;

        MarkerInfo * info = GetMarkerInfo(tokens[1]);

        int chr = (tokens[0][0] == 'x' || tokens[0][0] == 'X') ? 999 : (int) tokens[0];

        info->chromosome = chr;
        info->position = (double) tokens[2] * 0.01;

        if (sexSpecificMap)
        {
            char * flag;

            double female = strtod(tokens[3], &flag);
            if (*flag)
                error("In the map file, the female cM position for marker\n"
                      "%s is %s. This is not a valid number.",
                      (const char *) tokens[1], (const char *) tokens[3]);

            double male = strtod(tokens[4], &flag);
            if (*flag)
                error("In the map file, the male cM position for marker\n"
                      "%s is %s. This is not a valid number.",
                      (const char *) tokens[1], (const char *) tokens[4]);

            info->positionFemale = (double) female * 0.01;
            info->positionMale = (double) male * 0.01;
        }
        else
            info->positionFemale = info->positionMale = info->position;
    }

    if (sexSpecificMap) VerifySexSpecificOrder();
}

void PedigreeGlobals::LoadBasepairMap(const char * filename)
{
    IFILE f = ifopen(filename, "rb");
    if (f == NULL)
        error("The map file [%s] could not be opened\n\n"
              "Please check that the filename is correct and that the file is\n"
              "not being used by another program", filename);
    LoadBasepairMap(f);
    ifclose(f);
}

void PedigreeGlobals::LoadBasepairMap(IFILE & input)
{
    String      buffer;
    StringArray tokens;

    sexSpecificMap = false;

    while (!ifeof(input))
    {
        buffer.ReadLine(input);

        tokens.Clear();
        tokens.AddTokens(buffer, WHITESPACE);

        if (tokens.Length() < 1) continue;

        if (tokens.Length() != 3)
            error("Error reading map file\n"
                  "Each line in this file should include 3 fields:\n"
                  "CHROMOSOME, MARKER_NAME, and POSITION\n"
                  "However the following line has %d fields\n%s\n",
                  tokens.Length(), (const char *) buffer);

        bool previous_state = String::caseSensitive;
        String::caseSensitive = false;

        if ((tokens[0] == "CHR" || tokens[0] == "CHROMOSOME") &&
                (tokens[1] == "MARKER" || tokens[1] == "MARKER_NAME" || tokens[1] == "MRK") &&
                (tokens[2] == "BASEPAIR" || tokens[2] == "POS" || tokens[2] == "POSITION"))
            continue;

        String::caseSensitive = previous_state;

        MarkerInfo * info = GetMarkerInfo(tokens[1]);

        int chr = (tokens[0][0] == 'x' || tokens[0][0] == 'X') ? 999 : (int) tokens[0];

        info->chromosome = chr;
        info->position = (double) tokens[2];
    }
}

int PedigreeGlobals::instanceCount = 0;

PedigreeGlobals::~PedigreeGlobals()
{
    if (--instanceCount == 0 && markerInfoSize)
    {
        for (int i = 0; i < markerInfoCount; i++)
            delete markerInfo[i];
        delete [] markerInfo;
        delete [] markerInfoByInteger;
    }
}

void PedigreeGlobals::WriteMapFile(const char * filename)
{
    if (!MarkerPositionsAvailable())
        return;

    FILE * output = fopen(filename, "wt");

    if (output == NULL)
        error("Creating map file \"%s\"", filename);

    WriteMapFile(output);
    fclose(output);
}

void PedigreeGlobals::WriteMapFile(FILE * output)
{
    if (!sexSpecificMap)
        fprintf(output, "CHR  MARKER    POS\n");
    else
        fprintf(output, "CHR  MARKER    POS   POSF   POSM\n");

    for (int i = 0; i < markerInfoCount; i++)
    {
        if (markerInfo[i]->chromosome != -1)
        {
            if (!sexSpecificMap)
                fprintf(output, "%3d  %-10s  %g\n",
                        markerInfo[i]->chromosome,
                        (const char *) markerInfo[i]->name,
                        markerInfo[i]->position * 100.0);
            else
                fprintf(output, "%3d %-10s %g  %g  %g\n",
                        markerInfo[i]->chromosome,
                        (const char *) markerInfo[i]->name,
                        markerInfo[i]->position * 100.0,
                        markerInfo[i]->positionFemale * 100.0,
                        markerInfo[i]->positionMale * 100.0);
        }
    }
}

void PedigreeGlobals::WriteFreqFile(const char * filename, bool old_format)
{
    FILE * output = fopen(filename, "wt");

    if (output == NULL)
        error("Creating allele frequency file \"%s\"", filename);

    WriteFreqFile(output, old_format);
    fclose(output);
}

void PedigreeGlobals::WriteFreqFile(FILE * output, bool old_format)
{
    for (int i = 0; i < markerInfoCount; i++)
    {
        MarkerInfo * info = markerInfo[i];

        if (info->freq.Length() == 0) continue;

        fprintf(output, "M %s\n", (const char *) info->name);

        if (old_format && info->alleleLabels.Length() == 0)
            for (int j = 1; j < info->freq.Length(); j++)
                fprintf(output, "%s%.5f%s",
                        j % 7 == 1 ? "F " : "", info->freq[j],
                        j == info->freq.Length() - 1 ? "\n" : j % 7 == 0 ? "\n" : " ");
        else
            for (int j = 1; j < info->freq.Length(); j++)
                if (info->freq[j] > 1e-7)
                    fprintf(output, "A %5s %.5f\n",
                            (const char *) info->GetAlleleLabel(j), info->freq[j]);
    }
}

bool PedigreeGlobals::MarkerPositionsAvailable()
{
    for (int i = 0; i < markerInfoCount; i++)
        if (markerInfo[i]->chromosome != -1)
            return true;

    return false;
}

bool PedigreeGlobals::AlleleFrequenciesAvailable()
{
    for (int i = 0; i < markerInfoCount; i++)
        if (markerInfo[i]->freq.Length() > 1)
            return true;

    return false;
}

int PedigreeGlobals::LoadAllele(int marker, String & token)
{
    return LoadAllele(GetMarkerInfo(marker), token);
}

int PedigreeGlobals::LoadAllele(MarkerInfo * info, String & token)
{
    int allele = info->GetAlleleNumber(token);

    if (allele >= 0) return allele;

    static unsigned char lookup[128];
    static bool init = false;

    if (!init)
    {
        init = true;

        for (int i = 0; i < 128; i++)
            lookup[i] = 0;

        for (int i = '1'; i <= '9'; i++)
            lookup[i] = 1;

        lookup[int('a')] = lookup[int('A')] = lookup[int('c')] = lookup[int('C')] = 2;
        lookup[int('g')] = lookup[int('G')] = lookup[int('t')] = lookup[int('T')] = 2;
    }

    int  first = token[0];
    bool goodstart = first > 0 && first < 128;

    if (token.Length() == 1 && goodstart && lookup[int(token[0])])
        return info->NewAllele(token);

    if (!goodstart || lookup[int(token[0])] != 1)
        return 0;

    int integer = token.AsInteger();
    token = integer;

    allele = info->GetAlleleNumber(token);

    if (allele > 0)
        return allele;

    if (integer <= 0) return 0;

    if (integer > 1000000)
    {
        static bool warn_user = true;

        if (warn_user)
        {
            printf("Some allele numbers for marker %s are > 1000000\n"
                   "All allele numbers >1000000 will be treated as missing\n\n",
                   (const char *) info->name);
            warn_user = false;
        }

        return 0;
    }

    return info->NewAllele(token);
}

std::ostream &operator << (std::ostream &stream, MarkerInfo &m)
{
    stream << "MarkerInfo for marker " << m.name << std::endl;
    stream << " located on chromsome " << m.chromosome << ":" << (int64_t)(100 * m.position) << std::endl;
    stream << " allele count = " << m.freq.Length() << std::endl;
    stream << " label count = " << m.alleleLabels.Length() << std::endl;
    if (m.freq.Length() == m.alleleLabels.Length())
    {
        for (int i=0; i<m.freq.Length(); i++)
        {
            stream << " " << m.alleleLabels[i] << ":" << m.freq[i];
        }
        stream << std::endl;
    }
    else
    {
        stream << " output truncated - counts appear corrupt." << std::endl;
    }
    return stream;
}

