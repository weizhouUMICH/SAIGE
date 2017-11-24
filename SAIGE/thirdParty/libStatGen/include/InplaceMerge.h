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

#ifndef _INPLACE_MERGE_H
#define _INPLACE_MERGE_H

#include <algorithm>
#if defined(DEBUG_INPLACE_MERGE)
#include "Generic.h"
#include <iostream>
#endif
#include <stdexcept>
#include <vector>



//
// given a partially ordered vector of values, use
// inplace_merge to merge the ordered subsets together in some
// reasonable fashion.
//
// On output, values is sorted in ascending order.
//
// the counts vector is also modified, the result being
// undefined, except that counts[0] == values.size() at final exit.
//
template<typename T> void inplace_merge(
    std::vector<int> &indeces,
    std::vector<int> &counts,
    int                 first,
    int                 last,
    std::vector<T> &values)
{
    if (first == (last)) return;    // empty set -> no merge
    if (first == (last-1)) return;  // only one set present -> no merge

    // here we see if we have non-adjacent sets to merge,
    // if so, do them independently, then we can do a final
    // merge next
    if (first != (last - 2))
    {
        int middle = (first + last) / 2;
        inplace_merge(indeces, counts, middle, last, values);
#if defined(DEBUG_INPLACE_MERGE)
        std::cout << values;
#endif
        inplace_merge(indeces, counts, first, middle, values);
#if defined(DEBUG_INPLACE_MERGE)
        std::cout << values;
#endif

        // get ready to drop through to below code which will
        // merge our two merged subsets
        last = middle + 1;
    }

    // inplace_merge just two adjacent sets
    typename std::vector<T>::iterator startIterator = values.begin()+indeces[first];
    typename std::vector<T>::iterator middleIterator = values.begin() + indeces[last-1];
    typename std::vector<T>::iterator endIterator = values.begin() + indeces[last-1] + counts[last - 1];
    std::inplace_merge(startIterator, middleIterator, endIterator);
    counts[first] += counts[last - 1];
#if defined(DEBUG_INPLACE_MERGE)
    std::cout << values;
#endif

    return;
}

#endif
