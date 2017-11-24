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

#ifndef __BASICHASH_H__
#define __BASICHASH_H__

#include <stdlib.h>

class BasicHash
{
protected:
    void          ** objects;
    unsigned int      * keys;
    unsigned int count, size;
    unsigned int        mask;

public:
    BasicHash(int startsize = 32);
    virtual ~BasicHash();

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

    void * Object(int i) const
    {
        return objects[i];
    }

    void SetObject(int i, void * object)
    {
        objects[i] = object;
    }

    int Add(int key, void * object = NULL);
    int Find(int key);
    int Rehash(int key, int h);

    BasicHash & operator = (const BasicHash & rhs);

    void * operator [](int i) const
    {
        return objects[i];
    }

    void Delete(unsigned int index);

    bool SlotInUse(int index)
    {
        return objects[index] != NULL;
    }

private:
    unsigned int Iterate(unsigned int key) const
    {
        unsigned int h = key & mask;

        while (objects[h] != NULL && keys[h] != key)
            h = (h + 1) & mask;

        return h;
    }

    unsigned int ReIterate(unsigned int key, unsigned int h) const
    {
        h = (h + 1) & mask;

        while (objects[h] != NULL && keys[h] != key)
            h = (h + 1) & mask;

        return h;
    }
};

#endif
