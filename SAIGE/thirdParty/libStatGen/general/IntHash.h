/*
 *  Copyright (C) 2000-2007 Goncalo Abecasis
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

//////////////////////////////////////////////////////////////////////
// libsrc/IntHash.h
// (c) 2000-2007 Goncalo Abecasis
//
// This file is distributed as part of the MaCH source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
//
// Permission is granted for you to use this file to compile MaCH.
//
// All computer programs have bugs. Use this file at your own risk.
//
// Monday October 29, 2007
//

#ifndef __INTHASH_H__
#define __INTHASH_H__

#include <stdlib.h>

class IntHash
{
protected:
    bool           * objects;
    unsigned int      * keys;
    unsigned int count, size;
    unsigned int        mask;

public:
    IntHash(int startsize = 32);
    virtual ~IntHash();

    void Grow()
    {
        SetSize(size * 2);
    }
    void Shrink()
    {
        SetSize(size / 2);
    }

    void SetSize(int newsize);

    void Clear();

    int  Capacity() const
    {
        return size;
    }
    int  Entries() const
    {
        return count;
    }

    bool Object(int i) const
    {
        return objects[i];
    }

    void SetObject(int i, bool object)
    {
        objects[i] = object;
    }

    int Add(int key, bool object = true);
    int Find(int key);
    int Rehash(int key, int h);

    IntHash & operator = (const IntHash & rhs);

    bool operator [](int i) const
    {
        return objects[i];
    }

    void Delete(unsigned int index);

    bool SlotInUse(int index)
    {
        return objects[index] != false;
    }

private:
    unsigned int Iterate(unsigned int key) const
    {
        unsigned int h = key & mask;

        while (objects[h] != false && keys[h] != key)
            h = (h + 1) & mask;

        return h;
    }

    unsigned int ReIterate(unsigned int key, unsigned int h) const
    {
        h = (h + 1) & mask;

        while (objects[h] != false && keys[h] != key)
            h = (h + 1) & mask;

        return h;
    }
};

#endif

