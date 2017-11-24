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

#include "PedigreeAlleleFreq.h"
#include "QuickIndex.h"
#include "Error.h"

#include <math.h>

int CountAlleles(Pedigree & /* ped */, int marker)
{
    // With automatic recoding in the pedigree loader there
    // is no need to iterate through the pedigree ...
    MarkerInfo * info = Pedigree::GetMarkerInfo(marker);

    return info->CountAlleles();
}

void LumpAlleles(Pedigree & ped, int marker, double threshold, bool reorder)
{
    // find out how many alleles there are
    int alleles = ped.CountAlleles(marker);

    if (alleles < 2) return;

    MarkerInfo * info = PedigreeGlobals::GetMarkerInfo(marker);

    if (alleles < info->freq.Length())
        alleles = info->freq.Length() - 1;

    IntArray counts(alleles + 1);
    counts.Zero();

    // Count number of occurrences for each allele
    for (int i = 0; i < ped.count; i++)
    {
        counts[int(ped[i].markers[marker][0])]++;
        counts[int(ped[i].markers[marker][1])]++;
    }

    // Calculate treshold for lumping alleles
    int total = 0;
    for (int i = 1; i <= alleles; i++)
        total += counts[i];
    int thresh = int(total * threshold);

    // If threshold is set at zero, we artificially increase
    // counts for alleles that do not appear in the pedigree
    // but whose frequencies are set > 0.0. This ensures that
    // allele frequency data does not get discarded when simply
    // recoding alleles (vs. lumping)
    if (thresh == 0)
        for (int i = 1; i < info->freq.Length(); i++)
            if (counts[i] == 0 && info->freq[i] > 0.0)
                counts[i] = 1, total++;

    // If allele reordering is disabled, put in dummy allele
    // counts so as to ensure that allele have desired ordering
    if (!reorder)
    {
        QuickIndex index(info->alleleLabels);
        index.Reverse();

        for (int i = 0; i < index.Length(); i++)
            counts[index[i]] = i + 1;

        total = counts.Sum(1, counts.Length() - 1);
    }

    // Order all alleles according to their frequency
    // Zero, which corresponds to missing values, stays put!
    counts[0] = total + 1;
    QuickIndex index(counts);
    index.Reverse();

    // recode alleles
    // all alleles where frequency < thresh are labelled N
    // use counts array to keep track of labels
    int  N = 0;
    bool rare = false;
    for (int i = 0; i <= alleles; i++)
        if (counts[index[i]] > thresh)
        {
            counts[index[i]] = i;
            N++;
        }
        else
        {
            if (counts[index[i]] > 0)
                rare = true;
            counts[index[i]] = N;
        }

    // This loop does the recoding
    for (int i = 0; i < ped.count; i++)
    {
        Alleles & current = ped[i].markers[marker];
        current[0] = counts[current[0]];
        current[1] = counts[current[1]];
    }

    StringArray  oldLabels(info->alleleLabels);
    String       label;

    info->alleleLabels.Clear();
    info->alleleNumbers.Clear();

    for (int i = 0; i < N; i++)
    {
        if (oldLabels.Length() <= index[i])
            info->alleleLabels.Push(label = index[i]);
        else
            info->alleleLabels.Push(oldLabels[index[i]]);

        if (i) info->alleleNumbers.SetInteger(info->alleleLabels.Last(), i);
    }

    // Reorder allele frequencies if necessary
    if (info->freq.Length())
    {
        Vector freq(info->freq);

        info->freq.Dimension(N);
        info->freq[0] = 0.0;

        for (int i = 1; i < N; i++)
        {
            info->freq[i]  = freq[index[i]];
            freq[index[i]] = 0;
        }

        if ((1.0 - info->freq.Sum()) > 1e-10)
            rare = true;

        if (rare)
        {
            info->freq.Dimension(N + 1);
            info->freq[N] = 1.0 - info->freq.Sum();
        }
    }

    if (rare)
    {
        info->alleleLabels.Push("OTHER");
        info->alleleNumbers.SetInteger("OTHER", info->alleleLabels.Length());
    }
}

bool EstimateFrequencies(Pedigree & ped, int marker, int estimator)
{
    int alleleCount = CountAlleles(ped, marker);

    IntArray founder(alleleCount + 1);
    IntArray all(alleleCount + 1);

    founder.Zero();
    all.Zero();

    for (int i = 0; i < ped.count; i++)
    {
        // When counting alleles, note that males only carry one X chromosome
        // and are arbitrarily scored as homozygous.
        all[ped[i].markers[marker][0]]++;
        if (!ped.chromosomeX || ped[i].sex != SEX_MALE)
            all[ped[i].markers[marker][1]]++;
        if (!ped[i].isFounder()) continue;
        founder[ped[i].markers[marker][0]]++;
        if (!ped.chromosomeX || ped[i].sex != SEX_MALE)
            founder[ped[i].markers[marker][1]]++;
    }

    MarkerInfo * info = ped.GetMarkerInfo(marker);

    if (info->freq.dim > 0)
    {
        // previous allele frequency information is available
        if (alleleCount >= info->freq.dim)
            error("For marker %s, input files define %d alleles, but at least\n"
                  "one other allele (named '%s') occurs in the pedigree\n",
                  (const char *) info->name, info->freq.dim - 1,
                  (const char *) info->GetAlleleLabel(alleleCount));

        for (int i = 1; i <= alleleCount; i++)
            if (all[i] > 0 && info->freq[i] <= 0.0)
                error("Although allele %s for marker %s has frequency zero,\n"
                      "it occurs %d times in the pedigree",
                      (const char *) info->GetAlleleLabel(i), (const char *) info->name, all[i]);

        return false;
    }
    else
    {
        if (alleleCount < 1)
        {
            // If no one is genotyped, default to two equifrequent allele
            // since some programs do not like monomorphic markers
            info->freq.Dimension(3);
            info->freq[0] = 0.0;
            info->freq[1] = 0.99999;
            info->freq[2] = 0.00001;
            return true;
        }

        info->freq.Dimension(alleleCount + 1);
        info->freq.Zero();

        if (estimator == FREQ_FOUNDERS && founder.Sum() > founder[0])
        {
            // Make sure the frequency of alleles occuring in the pedigree
            // is never zero
            for (int i = 1; i <= alleleCount; i++)
                if (founder[i] == 0 && all[i] > 0)
                    founder[i] = 1;

            // To get frequencies, just multiply counts by 1 / total_counts
            double factor = 1.0 / (founder.Sum() - founder[0]);

            for (int i = 1; i <= alleleCount; i++)
                info->freq[i] = founder[i] * factor;
        }
        else if (estimator == FREQ_ALL || estimator == FREQ_FOUNDERS)
        {
            // To get frequencies, just multiply counts by 1 / total_counts
            double factor = 1.0 / (all.Sum() - all[0]);

            for (int i = 1; i <= alleleCount; i++)
                info->freq[i] = all[i] * factor;
        }
        else if (estimator == FREQ_EQUAL)
            // Assume all alleles have equal frequency
        {
            // Count the number of observed alleles
            all[0] = 0;
            int alleles = all.CountIfGreater(0);
            double freq = 1.0 / alleles;

            // Set equal frequencies for all occuring alleles
            for (int i = 0; i <= alleleCount; i++)
                info->freq[i] = all[i] ? freq : 0.0;
        }
    }

    return true;
}

