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

#ifndef __QUICKINDEX_H__
#define __QUICKINDEX_H__

#include "MathVector.h"
#include "StringArray.h"
#include "StringHash.h"
#include "IntArray.h"
#include "StringMap.h"

class QuickIndex : public IntArray
{
public:
    QuickIndex();
    QuickIndex(const IntArray & source_data)
    {
        Index(source_data);
    }
    QuickIndex(const StringArray & source_data)
    {
        Index(source_data);
    }
    QuickIndex(const Vector & source_data)
    {
        Index(source_data);
    }

    void Index(const IntArray & source_data);
    void Index(const StringArray & source_data);
    void Index(const Vector & source_data);
    void IndexCounts(const StringIntMap & source_data);
    void IndexCounts(const StringIntHash & source_data);

private:
    const void * source;
    int    datatype;

    bool IsBefore(int i, int j);
    void Sort();
};

#endif

