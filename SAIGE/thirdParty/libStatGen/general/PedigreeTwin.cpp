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
#include "Error.h"

#include <stdio.h>

bool Pedigree::TwinCheck()
{
    bool fail = false;
    IntArray mzTwins;

    for (int f = 0; f < familyCount; f++)
    {
        mzTwins.Clear();

        for (int i = families[f]->first, j; i <= families[f]->last; i++)
            // Is this person an identical twin?
            if (persons[i]->isMzTwin(*persons[i]))
            {
                // Have we got another identical sib yet?
                for (j = 0; j < mzTwins.Length(); j++)
                    if (persons[i]->isMzTwin(*persons[mzTwins[j]]))
                        break;

                // If not, add to list of twins
                if (j == mzTwins.Length())
                {
                    mzTwins.Push(i);
                    continue;
                }

                // Check that their genotypes are compatible and
                // merge new twin's genotypes into original twin...
                Person * original = persons[mzTwins[j]];
                Person * twin = persons[i];

                for (int m = 0; m < Person::markerCount; m++)
                {
                    if (!original->markers[m].isKnown())
                        original->markers[m] = twin->markers[m];
                    else if (twin->markers[m].isKnown() &&
                             twin->markers[m] != original->markers[m])
                        printf("MZ Twins %s and %s in family %s have "
                               "different %s genotypes\n",
                               (const char *) original->pid,
                               (const char *) twin->pid,
                               (const char *) original->famid,
                               (const char *) Person::markerNames[m]),
                        fail = true;

                    if (twin->sex != original->sex)
                        printf("MZ Twins %s and %s in family %s have "
                               "different sexes\n",
                               (const char *) original->pid,
                               (const char *) twin->pid,
                               (const char *) original->famid),
                        fail = true;
                }
            }

        if (mzTwins.Length() == 0) continue;

        // In the second pass we copy merged twin genotypes
        // from original twin to other twins
        for (int i = families[f]->first, j; i <= families[f]->last; i++)
            if (persons[i]->isMzTwin(*persons[i]))
            {
                for (j = 0; j < mzTwins.Length(); j++)
                    if (persons[i]->isMzTwin(*persons[mzTwins[j]]))
                        break;

                if (mzTwins[j] == i) continue;

                Person * original = persons[mzTwins[j]];
                Person * twin = persons[i];

                for (int m = 0; m < Person::markerCount; m++)
                    twin->markers[m] = original->markers[m];
            }
    }
    return fail;
}

void Pedigree::MergeTwins()
{
    if (!haveTwins) return;

    IntArray mzTwins, surplus;

    printf("Merging MZ twins into a single individual...\n");

    for (int f = 0; f < familyCount; f++)
    {
        mzTwins.Clear();

        for (int i = families[f]->first, j; i <= families[f]->last; i++)
            if (persons[i]->isMzTwin(*persons[i]))
            {
                // Have we got another identical sib yet?
                for (j = 0; j < mzTwins.Length(); j++)
                    if (persons[i]->isMzTwin(*persons[mzTwins[j]]))
                        break;

                // If not, add to list of twins
                if (j == mzTwins.Length())
                {
                    mzTwins.Push(i);
                    continue;
                }

                // Append name to first twins name
                persons[mzTwins[j]]->pid += ((char) '$') + persons[i]->pid;

                // Set the first twin to affected if at least one of the cotwins is affected
                for (int j = 0; j < affectionCount; j++)
                    if (persons[i]->affections[j] == 2)
                        persons[mzTwins[j]]->affections[j] = 2;

                surplus.Push(i);
            }

        // More than one representative of each twin-pair...
        if (surplus.Length())
        {
            // Reassign parent names for each offspring
            for (int i = families[f]->first, j; i < families[f]->last; i++)
                if (!persons[i]->isFounder())
                {
                    if (persons[i]->father->isMzTwin(*persons[i]->father))
                    {
                        for (j = 0; j < mzTwins.Length(); j++)
                            if (persons[i]->father->isMzTwin(*persons[mzTwins[j]]))
                                break;
                        persons[i]->fatid = persons[mzTwins[j]]->pid;
                    }
                    if (persons[i]->mother->isMzTwin(*persons[i]->mother))
                    {
                        for (j = 0; j < mzTwins.Length(); j++)
                            if (persons[i]->mother->isMzTwin(*persons[mzTwins[j]]))
                                break;
                        persons[i]->motid = persons[mzTwins[j]]->pid;
                    }
                }

            // Delete surplus individuals
            while (surplus.Length())
            {
                int serial = surplus.Pop();

                delete persons[serial];

                for (count--; serial < count; serial++)
                    persons[serial] = persons[serial + 1];
            }

            // Resort pedigree
            Sort();
        }
    }
}




