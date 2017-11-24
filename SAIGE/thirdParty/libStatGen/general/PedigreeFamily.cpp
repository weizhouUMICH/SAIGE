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
#include "Constant.h"
#include "MathConstant.h"
#include "Error.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

Family::Family(Pedigree & pedigree, int _first, int _last, int _serial) :
        ped(pedigree)
{
    serial = _serial;
    first = _first;
    last = _last;
    count = last - first + 1;
    path = new int [count];
    famid = ped[first].famid;

    founders = mzTwins = 0;

    for (int i=first; i<=last; i++)
        if (ped[i].isFounder())
        {
            ped[i].traverse = founders;
            path[founders++] = ped[i].serial;
        }
        else
        {
            ped[i].traverse = -1;
            if (ped[i].isMzTwin(ped[i]))
                for (int j = first; j < i; j++)
                    if (ped[i].isMzTwin(ped[j]))
                    {
                        mzTwins++;
                        break;
                    }
        }

    nonFounders = count - founders;
    generations = nonFounders == 0 ? 1 : 2;

    int next = founders;
    while (next < count)
    {
        bool check = false;

        // Create traversal where path ancestors precede their offspring
        for (int i=first; i<=last; i++)
            if (ped[i].traverse == -1)
            {
                int fatherSerial = ped[i].father->traverse;
                int motherSerial = ped[i].mother->traverse;

                if (fatherSerial >= 0 && motherSerial >= 0)
                {
                    check = true;

                    ped[i].traverse = next;
                    path[next++] = i;

                    if (fatherSerial >= founders || motherSerial >= founders)
                        generations = 3;

                    // If this individual is part of a set of MZ twins
                    if (ped[i].zygosity & 1)
                        for (int j = 0; j < ped[i].sibCount; j++)
                        {
                            Person & sib = *ped[i].sibs[j];

                            // Insert all co-twins at the same position in traversal
                            // order
                            if (sib.traverse == -1 && ped[i].zygosity == sib.zygosity)
                            {
                                sib.traverse = next;
                                path[next++] = sib.serial;
                            }
                        }
                }
            }

        if (!check) ShowInvalidCycles();
    }
}

Family::~Family()
{
    delete [] path;
}

void Family::ShowInvalidCycles()
{
    // Try and identify key individuals responsible for
    // pedigree mess-up ... when this function is called
    // pedigree has been traversed top-down and individuals
    // that are correctly specified have IDs of >= 0.

    // This routine traverses the pedigree bottom up to
    // identify a subset of individuals likely to be causing
    // the problem
    IntArray descendants(ped.count);
    descendants.Zero();

    for (int i = first; i <= last; i++)
        if (ped[i].traverse == -1)
        {
            descendants[ped[i].father->serial]++;
            descendants[ped[i].mother->serial]++;
        }

    IntArray stack;

    for (int i = first; i <= last; i++)
        if (ped[i].traverse == -1 && descendants[i] == 0)
        {
            stack.Push(i);

            do
            {
                int j = stack.Pop();

                if (ped[j].traverse != -1) continue;

                ped[j].traverse = 9999;

                if (--descendants[ped[j].father->serial] == 0)
                    stack.Push(ped[j].father->serial);
                if (--descendants[ped[j].mother->serial] == 0)
                    stack.Push(ped[j].mother->serial);
            }
            while (stack.Length());
        }

    printf("The structure of family %s requires\n"
           "an individual to be his own ancestor.\n\n"
           "To identify the problem(s), examine the\n"
           "following key individuals:\n\n",
           (const char *) famid);

    for (int i = first; i <= last; i++)
        if (ped[i].traverse == -1)
            printf("Problem Person: %s\n", (const char *) ped[i].pid);

    error("Invalid pedigree structure.");
}

int Family::ConnectedGroups(IntArray * groupMembership)
{
    IntArray groups(count);

    // Use the quick union algorithm to identify connected groups
    groups.SetSequence(0, 1);
    for (int i = count - 1; i >= founders; i--)
    {
        // Lookup parents
        int group0 = i;
        int group1 = ped[path[i]].father->traverse;
        int group2 = ped[path[i]].mother->traverse;

        // Identify their corresponding groupings
        while (groups[group0] != group0) group0 = groups[group0];
        while (groups[group1] != group1) group1 = groups[group1];
        while (groups[group2] != group2) group2 = groups[group2];

        int group = group1 < group2 ? group1 : group2;
        if (group0 < group) group = group0;

        groups[group0] = groups[group1] = groups[group2] = group;
    }

    // Count groupings
    int groupCount = 0;
    for (int i = 0; i < founders; i++)
        if (groups[i] == i)
            groupCount++;

    if (groupMembership == NULL)
        return groupCount;

    // Flatten tree so all items point to root
    for (int i = 1; i < count; i++)
        groups[i] = groups[groups[i]];

    // Update group membership info
    int group = 0;
    groupMembership->Dimension(count);
    for (int i = 0; i < count; i++)
        if (groups[i] == i)
            (*groupMembership)[i] = ++group;
        else
            (*groupMembership)[i] = (*groupMembership)[groups[i]];

#if 0
    // This stretch of code outputs family structure and group membership
    // And should usually be commented out!
    for (int j = first; j <= last; j++)
        printf("%s %s %s %s %d %d\n",
               (const char *) famid, (const char *) ped[j].pid,
               (const char *) ped[j].fatid, (const char *) ped[j].motid,
               ped[j].sex, groups[ped[j].traverse]);
#endif

    return groupCount;
}

/*
int Family::ConnectedGroups(IntArray * groupMembership)
   {
   IntArray * stack = new IntArray[count];
   IntArray groups(count);

   groups.Zero();

   int group = 0;
   int seed = count - 1;

   // Search for connected sets of individuals until everyone is accounted for
   while (true)
      {
      while ((seed >= 0) && (groups[seed] != 0))
         seed--;

      if (seed == -1)
         break;

      Mark(seed, ++group, stack, groups);

      for (int j = seed; j >= founders; j--)
         if (groups[j] == 0)
            {
            int fat_j = ped[path[j]].father->traverse;
            int mot_j = ped[path[j]].mother->traverse;

            if (groups[fat_j] == group || groups[mot_j] == group)
               Mark(j, group, stack, groups);
            else
               stack[mot_j].Push(j),
               stack[fat_j].Push(j);
            }

      for (int j = 0; j < count; j++)
         stack[j].Clear();
      }

   if (groupMembership != NULL)
      (*groupMembership) = groups;

   // This stretch of code outputs family structure and group membership
   // And should usually be commented out!
#if 0
   for (int j = first; j <= last; j++)
      printf("%s %s %s %s %d %d\n",
         (const char *) famid, (const char *) ped[j].pid,
         (const char *) ped[j].fatid, (const char *) ped[j].motid,
         ped[j].sex, groups[ped[j].traverse]);
#endif

   delete [] stack;

   return group;
   }

void Family::Mark(int j, int group, IntArray * stack, IntArray & groups)
   {
   if (groups[j] == group) return;

   groups[j] = group;

   while (stack[j].Length())
      Mark(stack[j].Pop(), group, stack, groups);

   if (j < founders) return;

   Mark(ped[path[j]].father->traverse, group, stack, groups);
   Mark(ped[path[j]].mother->traverse, group, stack, groups);
   }
*/
