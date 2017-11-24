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

#include "Sort.h"
#include "Error.h"

#include <stddef.h>
#include <string.h>


#define Item(b)         (base_char+(b)*width)
#define IsBefore(x,y)   ((cmp(Item(x),Item(y)))<0)
#define Exchange(x,y)   {\
                        memcpy(tmp,Item(x),width);\
                        memcpy(Item(x),Item(y),width);\
                        memcpy(Item(y),tmp,width);\
                        }
#define TRUE   1

void QuickSort(void *base, size_t nelem, size_t width,
               int (*cmp)(const void *, const void *))
{
    struct __QuickSortStack
    {
        size_t left, right;
    };

    if (nelem <= 1)
        return;

    // Create a pseudo-stack to avoid recursion

    char * base_char = (char *) base;
    const size_t stackSize = 128;

    __QuickSortStack * stack = new __QuickSortStack[stackSize];
    char * tmp = new char [width];

    if ((stack == NULL) || (tmp == NULL))
        error("Out of memory in QuickSort routine");

    size_t stackIdx = 0;

    // Size of minimum partition to median of three
    const size_t Threshold = 7;

    // current partitions

    size_t lsize, rsize;
    size_t l, mid, r;
    size_t scanl, scanr, pivot;

    l = 0;
    r = nelem - 1;

    while (TRUE)
    {
        while (r > l)
        {
            if (r - l > Threshold)
                // QuickSort : median of three partitioning
            {
                mid = (r + l) / 2;

                // sort l, mid, and r
                if (IsBefore(mid, l))
                    Exchange(mid, l);

                if (IsBefore(r, l))
                    Exchange(r, l);

                if (IsBefore(r, mid))
                    Exchange(r, mid);

                // set up for partitioning...
                pivot = r - 1;

                Exchange(mid, pivot);

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

            while (TRUE)
            {
                // scan from left for element >= pivot
                while ((scanl < r) && IsBefore(scanl, pivot))
                    ++scanl;

                while ((scanr > l) && IsBefore(pivot, scanr))
                    --scanr;

                // if scans have met, we are done
                if (scanl >= scanr)
                    break;

                Exchange(scanl, scanr);

                if (scanl < r)
                    ++scanl;

                if (scanr > l)
                    --scanr;
            }

            // Exchange final element
            Exchange(pivot, scanl);

            // Place largest partition on stack
            lsize = scanl - l;
            rsize = r - scanl;

            if (lsize > rsize)
            {
                // if size is one we are done
                ++ stackIdx;

                if (stackIdx == stackSize)
                    error("Out of Stack in QuickSort routine");

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

                if (stackIdx == stackSize)
                    error("Out of Stack in QuickSort routine");

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

    delete [] stack;
    delete [] tmp;
}

#define Item2(b)        (base_char2+(b)*width)
#define Exchange2(x,y)  {\
                        memcpy(tmp,Item(x),width);\
                        memcpy(Item(x),Item(y),width);\
                        memcpy(Item(y),tmp,width);\
                        memcpy(tmp,Item2(x),width);\
                        memcpy(Item2(x),Item2(y),width);\
                        memcpy(Item2(y),tmp,width);\
                        }


void QuickSort2(void *base, void *base2, size_t nelem, size_t width,
                int (*cmp)(const void *, const void *))
{
    struct __QuickSortStack
    {
        size_t left, right;
    };

    if (nelem <= 1)
        return;

    // Create a pseudo-stack to avoid recursion

    char * base_char = (char *) base;
    char * base_char2 = (char *) base2;
    const size_t stackSize = 128;

    __QuickSortStack * stack = new __QuickSortStack[stackSize];
    char * tmp = new char [width];

    if ((stack == NULL) || (tmp == NULL))
        error("Out of memory in QuickSort routine");

    size_t stackIdx = 0;

    // Size of minimum partition to median of three
    const size_t Threshold = 7;

    // current partitions

    size_t lsize, rsize;
    size_t l, mid, r;
    size_t scanl, scanr, pivot;

    l = 0;
    r = nelem - 1;

    while (TRUE)
    {
        while (r > l)
        {
            if (r - l > Threshold)
                // QuickSort : median of three partitioning
            {
                mid = (r + l) / 2;

                // sort l, mid, and r
                if (IsBefore(mid, l))
                    Exchange2(mid, l);

                if (IsBefore(r, l))
                    Exchange2(r, l);

                if (IsBefore(r, mid))
                    Exchange2(r, mid);

                // set up for partitioning...
                pivot = r - 1;

                Exchange2(mid, pivot);

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

            while (TRUE)
            {
                // scan from left for element >= pivot
                while ((scanl < r) && IsBefore(scanl, pivot))
                    ++scanl;

                while ((scanr > l) && IsBefore(pivot, scanr))
                    --scanr;

                // if scans have met, we are done
                if (scanl >= scanr)
                    break;

                Exchange2(scanl, scanr);

                if (scanl < r)
                    ++scanl;

                if (scanr > l)
                    --scanr;
            }

            // Exchange final element
            Exchange2(pivot, scanl);

            // Place largest partition on stack
            lsize = scanl - l;
            rsize = r - scanl;

            if (lsize > rsize)
            {
                // if size is one we are done
                ++ stackIdx;

                if (stackIdx == stackSize)
                    error("Out of Stack in QuickSort routine");

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

                if (stackIdx == stackSize)
                    error("Out of Stack in QuickSort routine");

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

    delete [] stack;
    delete [] tmp;
}

void * BinarySearch(const void *key, const void *base,
                    size_t nelem, size_t width,
                    int (*cmp)(const void *, const void *))
{
    if (nelem == 0)
        return NULL;

    char * base_char = (char *) base;

    int left = 0;
    int right = nelem - 1;

    while (right >= left)
    {
        int probe = (left + right) / 2;
        int test  = cmp(key, Item(probe));

        if (test == 0)
            return (void *) Item(probe);

        if (test < 0)
            right = probe - 1;
        else
            left  = probe + 1;
    }

    return NULL;
}

