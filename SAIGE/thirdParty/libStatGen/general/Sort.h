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

#ifndef __SORT_H__
#define __SORT_H__

#include "Constant.h"

#include <stddef.h>

void QuickSort(void *base, size_t nelem, size_t width,
               int (*cmp)(const void *, const void *));

void QuickSort2(void *base, void * base2, size_t nelem, size_t width,
                int (*cmp)(const void *, const void *));

void * BinarySearch(const void *key, const void *base,
                    size_t nelem, size_t width,
                    int (*cmp)(const void *, const void *));

#endif
