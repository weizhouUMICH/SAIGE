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

#ifndef __LONGHASH_H__
#define __LONGHASH_H__

#include "Error.h"

#include <limits.h>

#ifdef   UINT_MAX
#define  LH_NOTFOUND     (UINT_MAX)
#else
#define  LH_NOTFOUND     0xFFFFFFFF
#endif

template <class ObjectT> class LongHash
{
protected:
    ObjectT        * objects;
    long long      * keys;
    bool           * occupancy;
    unsigned int     count, size;
    unsigned int     mask;
    bool             allowDuplicates;

public:
    LongHash(int startsize = 32)
    {
        count = 0;
        size  = startsize;
        mask  = startsize - 1;

        // In this implementation, the size of hash tables must be a power of two
        if (startsize & mask)
            error("LongHash: Hash table size must be a power of two.\n");

        occupancy = new bool [size];
        objects = new ObjectT [size];
        keys    = new long long [size];

        allowDuplicates = false;

        for (unsigned int i = 0; i < size; i++)
        {
            occupancy[i] = false;
        }
    };

    ~LongHash()
    {
        delete [] occupancy;
        delete [] objects;
        delete [] keys;
    }

    void Grow()
    {
        SetSize(size * 2);
    }
    void Shrink()
    {
        SetSize(size / 2);
    }

    void SetSize(int newsize)
    {
        int newmask = newsize - 1;

        bool      * newoccupancy = new bool [newsize];
        ObjectT   * newobjects = new ObjectT [newsize];
        long long * newkeys = new long long [newsize];

        for (int i = 0; i < newsize; i++)
            newoccupancy[i] = false;

        if (count)
            for (unsigned int i = 0; i < size; i++)
                if (occupancy[i] != false)
                {
                    long long  key = keys[i];
                    unsigned int h = newmask & (unsigned int) key;

                    while (newoccupancy[h] == true && (newkeys[h] != key || allowDuplicates))
                        h = (h + 1) & newmask;

                    if (newoccupancy[h])
                        count--;

                    newkeys[h] = key;
                    newobjects[h] = objects[i];
                    newoccupancy[h] = true;
                }

        delete [] occupancy;
        delete [] objects;
        delete [] keys;

        occupancy = newoccupancy;
        objects = newobjects;
        keys = newkeys;
        size = newsize;
        mask = newmask;
    }

    void Clear()
    {
        count = 0;

        if (size > 32)
            SetSize(32);

        for (unsigned int i = 0; i < size; i++)
            occupancy[i] = false;
    }

    int  Capacity() const
    {
        return size;
    }
    int  Entries() const
    {
        return count;
    }

    ObjectT Object(int i) const
    {
        return objects[i];
    }
    ObjectT & Object(int i)
    {
        return objects[i];
    }

    void SetObject(int i, ObjectT object)
    {
        objects[i] = object;
    }

    unsigned int Add(long long key, ObjectT object)
    {
        if (count * 2 > size)
            Grow();

        unsigned int h = Iterate(key);

        while (allowDuplicates && occupancy[h] && objects[h] != object)
            h = ReIterate(key, h);

        if (!occupancy[h])
        {
            occupancy[h] = true;
            keys[h] = key;
            count++;
        }

        objects[h] = object;

        return h;
    }

    unsigned int Find(long long key)
    {
        unsigned int h = Iterate(key);

        return occupancy[h] ? h : LH_NOTFOUND;
    }

    unsigned int Rehash(long long key, unsigned int h)
    {
        h = ReIterate(key, h);

        return occupancy[h] ? h : LH_NOTFOUND;
    }

    LongHash & operator = (const LongHash & rhs);

    ObjectT operator [](int i) const
    {
        return objects[i];
    }
    ObjectT operator [](unsigned int i) const
    {
        return objects[i];
    }

    void Delete(unsigned int index)
    {
        if (index >= size || !occupancy[index])
            return;

        occupancy[index] = false;
        count--;

        if (count * 8 < size && size > 32)
            Shrink();
        else
        {
            // rehash the next entries until we find empty slot
            index = (index + 1) & mask;

            while (occupancy[index])
            {
                if ((keys[index] & mask) != index)
                {
                    unsigned int h = Iterate(keys[index]);

                    while (occupancy[h] && objects[h] != objects[index])
                        h = ReIterate(keys[index], h);

                    if (h != (unsigned int) index)
                    {
                        keys[h] = keys[index];
                        occupancy[h] = true;
                        objects[h] = objects[index];

                        occupancy[index] = false;
                    }
                }

                index = (index + 1) & mask;
            }
        }
    }


    bool SlotInUse(int index) const
    {
        return occupancy[index] == true;
    }
    bool SlotInUse(unsigned int index) const
    {
        return occupancy[index] == true;
    }

    // Accessor to get a key.
    long long GetKey(int index) const
    {
        return keys[index];
    }

    long long GetKey(const unsigned int index) const
    {
        return keys[index];
    }

    void SetAllowDuplicateKeys(bool toggle)
    {
        allowDuplicates = toggle;

        if (count && !allowDuplicates)
            SetSize(size);
    }

private:
    unsigned int Iterate(long long key) const
    {
        unsigned int h = mask & (unsigned int) key;

        while (occupancy[h] == true && keys[h] != key)
            h = (h + 1) & mask;

        return h;
    }

    unsigned int ReIterate(long long key, unsigned int h) const
    {
        h = (h + 1) & mask;

        while (occupancy[h] == true && keys[h] != key)
            h = (h + 1) & mask;

        return h;
    }
};

#endif
