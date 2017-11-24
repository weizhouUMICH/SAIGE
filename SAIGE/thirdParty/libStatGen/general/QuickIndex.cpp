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

#include "QuickIndex.h"
#include "Error.h"

#define __QI_INVALID          0
#define __QI_VECTOR           1
#define __QI_INTARRAY         2
#define __QI_STRINGARRAY      3

QuickIndex::QuickIndex()
{
    source = NULL;
    datatype = __QI_INVALID;
}

void QuickIndex::Index(const IntArray & source_data)
{
    source = (const void *) &source_data;
    datatype = __QI_INTARRAY;

    Dimension(source_data.Length());
    SetSequence();
    Sort();
}

void QuickIndex::Index(const Vector & source_data)
{
    source = (const void *) &source_data;
    datatype = __QI_VECTOR;

    Dimension(source_data.Length());
    SetSequence();
    Sort();
}

void QuickIndex::Index(const StringArray & source_data)
{
    source = (const void *) &source_data;
    datatype = __QI_STRINGARRAY;

    Dimension(source_data.Length());
    SetSequence();
    Sort();
}

void QuickIndex::IndexCounts(const StringIntMap & source_data)
{
    IntArray counts(source_data.Length());

    for (int i = 0; i < source_data.Length(); i++)
        counts[i] = source_data.GetCount(i);

    Index(counts);
}

void QuickIndex::IndexCounts(const StringIntHash & source_data)
{
    IntArray counts(source_data.Capacity());

    for (int i = 0; i < source_data.Capacity(); i++)
        if (source_data.SlotInUse(i))
            counts[i] = source_data.Integer(i);
        else
            counts[i] = -1;

    Index(counts);

    Reverse();
    Dimension(source_data.Entries());
    Reverse();
}

bool QuickIndex::IsBefore(int i, int j)
{
    i = (*this)[i];
    j = (*this)[j];

    switch (datatype)
    {
        case __QI_VECTOR :
        {
            const Vector & data = * (const Vector *) source;
            return data[i] < data[j];
        }
        case __QI_INTARRAY :
        {
            const IntArray & data = * (const IntArray *) source;
            return data[i] < data[j];
        }
        case __QI_STRINGARRAY :
        {
            const StringArray & data = * (const StringArray *) source;
            return data[i].SlowCompare(data[j]) < 0;
        }
    }
    return 0;
}

void QuickIndex::Sort()
{
    struct __QuickIndexStack
    {
        int left, right;
    };

    if (Length() <= 1)
        return;

    // Create a pseudo-stack to avoid recursion
    __QuickIndexStack stack[32];

    int stackIdx = 0;

    // Size of minimum partition to median of three
    const int Threshold = 7;

    // current partitions
    int lsize, rsize;
    int l, mid, r;
    int scanl, scanr, pivot;

    l = 0;
    r = Length() - 1;

    while (1)
    {
        while (r > l)
        {
            if (r - l > Threshold)
                // QuickSort : median of three partitioning
            {
                mid = (r + l) / 2;

                // sort l, mid, and r
                if (IsBefore(mid, l))
                    Swap(mid, l);

                if (IsBefore(r, l))
                    Swap(r, l);

                if (IsBefore(r, mid))
                    Swap(r, mid);

                // set up for partitioning...
                pivot = r - 1;

                Swap(mid, pivot);

                scanl = l + 1;
                scanr = r - 2;
            }
            else
            {
                // set up random partition -- faster
                pivot = r;
                scanl = l;
                scanr = r - 1;
            }

            while (1)
            {
                // scan from left for element >= pivot
                while ((scanl < r) && IsBefore(scanl, pivot))
                    ++scanl;

                while ((scanr > l) && IsBefore(pivot, scanr))
                    --scanr;

                // if scans have met, we are done
                if (scanl >= scanr)
                    break;

                Swap(scanl, scanr);

                if (scanl < r)
                    ++scanl;

                if (scanr > l)
                    --scanr;
            }

            // Exchange final element
            Swap(pivot, scanl);

            // Place largest partition on stack
            lsize = scanl - l;
            rsize = r - scanl;

            if (lsize > rsize)
            {
                // if size is one we are done
                ++ stackIdx;

                stack[stackIdx].left = l;
                stack[stackIdx].right = scanl - 1;

                if (rsize != 0)
                    l = scanl + 1;
                else
                    break;
            }
            else
            {
                // if size is one we are done
                ++ stackIdx;

                stack[stackIdx].left = scanl + 1;
                stack[stackIdx].right = r;

                if (lsize != 0)
                    r = scanl - 1;
                else
                    break;
            }
        }

        // iterate with values from stack
        if (stackIdx)
        {
            l = stack[stackIdx].left;
            r = stack[stackIdx].right;

            --stackIdx;
        }
        else
            break;
    }
}







