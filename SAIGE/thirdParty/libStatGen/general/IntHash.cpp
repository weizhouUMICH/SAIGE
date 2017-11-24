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

#include "IntHash.h"
#include "Error.h"

#include <stdio.h>

IntHash::IntHash(int startsize)
{
    count = 0;
    size  = startsize;
    mask  = startsize - 1;

    // In this implementation, the size of hash tables must be a power of two
    if (startsize & mask)
        error("IntHash: Hash table size must be a power of two.\n");

    objects = new bool [size];
    keys    = new unsigned int [size];

    for (unsigned int i = 0; i < size; i++)
    {
        objects[i] = false;
    }
};

IntHash::~IntHash()
{
    delete [] objects;
    delete [] keys;
}

void IntHash::Clear()
{
//   printf("Clearing...\n");

    count = 0;

    if (size > 16)
        SetSize(16);

    for (unsigned int i = 0; i < size; i++)
        objects[i] = false;
}

void IntHash::SetSize(int newsize)
{
    int newmask = newsize - 1;

    bool      * newobjects = new bool [newsize];
    unsigned int * newkeys = new unsigned int [newsize];

    for (int i = 0; i < newsize; i++)
    {
        newobjects[i] = false;
    }

    if (count)
        for (unsigned int i = 0; i < size; i++)
            if (objects[i] != false)
            {
                unsigned int key = keys[i];
                unsigned int h   = key & newmask;

                while (newobjects[h] != false && newkeys[h] != h)
                    h = (h + 1) & newmask;

                newkeys[h] = key;
                newobjects[h] = objects[i];
            }

    delete [] objects;
    delete [] keys;

    objects = newobjects;
    keys = newkeys;
    size = newsize;
    mask = newmask;
}

int IntHash::Add(int key, bool object)
{
    if (count * 2 > size)
        Grow();

    unsigned int h = Iterate(key);

    while ((objects[h] != false) && (objects[h] != object))
        h = ReIterate(key, h);

    if (objects[h] == false)
    {
//      printf("At position %d, inserted %x\n", h, key);
        keys[h] = key;
        count++;
    }

    objects[h] = object;

    return h;
}

int IntHash::Find(int key)
{
    int h = Iterate(key);

    return objects[h] == false ? -1 : h;
}

int IntHash::Rehash(int key, int h)
{
    h = ReIterate(key, h);

    return objects[h] == false ? -1 : h;
}

void IntHash::Delete(unsigned int index)
{
    if (index >= size || objects[index] == false)
        return;

    objects[index] = false;
    count--;

    if (count * 8 < size && size > 32)
        Shrink();
    else
    {
        // rehash the next entries until we find empty slot
        index = (index + 1) & mask;

        while (objects[index] != false)
        {
            if ((keys[index] & mask) != index)
            {
                unsigned int h = Iterate(keys[index]);

                while ((objects[h] != false) && (objects[h] != objects[index]))
                    h = ReIterate(keys[index], h);

                if (h != (unsigned int) index)
                {
                    keys[h] = keys[index];
                    objects[h] = objects[index];
                    objects[index] = false;
                }
            }

            index = (index + 1) & mask;
        }
    }
}

