/*
 *  Copyright (C) 2010-2012  Regents of the University of Michigan
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

#include "StringHash.h"
#include "InputFile.h"
#include "Error.h"

StringHash::StringHash(int startsize)
    : StringHashBase()
{
    count = 0;
    size  = startsize;
    mask  = startsize - 1;

    // In this implementation, the size of hash tables must be a power of two
    if (startsize & mask)
        error("StringHash: Hash table size must be a power of two.\n");

    strings = new String * [size];
    objects = new void * [size];
    keys    = new unsigned int [size];

    for (unsigned int i = 0; i < size; i++)
    {
        strings[i] = NULL;
        objects[i] = NULL;
    }
};

StringHash::~StringHash()
{
    for (unsigned int i = 0; i < size; i++)
        if (strings[i] != NULL)
            delete strings[i];

    if(strings) delete [] strings;
    if(objects) delete [] objects;
    if(keys) delete [] keys;
}

void StringHash::Clear()
{
    for (unsigned int i = 0; i < size; i++)
        if (strings[i] != NULL)
        {
            delete strings[i];
            strings[i] = NULL;
        }

    count = 0;

    if (size > 256)
        SetSize(256);
}

void StringHash::SetSize(int newsize)
{
    int newmask = newsize - 1;

    String   ** newstrings = new String * [newsize];
    void     ** newobjects = new void * [newsize];
    unsigned int * newkeys = new unsigned int [newsize];

    for (int i = 0; i < newsize; i++)
    {
        newstrings[i] = NULL;
        newobjects[i] = NULL;
    }

    if (count)
        for (unsigned int i = 0; i < size; i++)
            if (strings[i] != NULL)
            {
                unsigned int key = keys[i];
                unsigned int h   = key & newmask;

                while (newstrings[h] != NULL &&
                        (newkeys[h] != key || 
                         (!stringsEqual(*(newstrings[h]), *(strings[i])))))
                    h = (h + 1) & newmask;

                newkeys[h] = key;
                newstrings[h] = strings[i];
                newobjects[h] = objects[i];
            }

    if(strings) delete [] strings;
    if(objects) delete [] objects;
    if(keys) delete [] keys;

    strings = newstrings;
    objects = newobjects;
    keys = newkeys;
    size = newsize;
    mask = newmask;
}

int StringHash::Add(const String & string, void * object)
{
    unsigned int key = getKey(string);
    unsigned int h   = Iterate(key, string);

    if (strings[h] == NULL)
        Insert(h, key, string);

    objects[h] = object;

    if (count * 2 > size)
    {
        Grow();
        return Iterate(key, string);
    }

    return h;
}

int StringHash::Find(const String & string,  void *(*create_object)())
{
    unsigned int key = getKey(string);
    unsigned int h   = Iterate(key, string);

    if (strings[h] == NULL && create_object == NULL)
        return -1;

    if (strings[h] == NULL && create_object != NULL)
    {
        Insert(h, key, string);
        objects[h] = create_object();

        if (count * 2 > size)
        {
            Grow();
            return Iterate(key, string);
        }
    }

    return h;
}

int StringHash::Find(const String & string) const
{
    unsigned int key = getKey(string);
    unsigned int h   = Iterate(key, string);

    if (strings[h] == NULL)
        return -1;

    return h;
}
void * StringHash::CreateHash()
{
    return (void *) new StringHash();
}

void StringHash::Delete(unsigned int index)
{
    if (index >= size || strings[index] == NULL)
        return;

    delete strings[index];
    strings[index] = NULL;
    count--;

    if (count * 8 < size && size > 32)
        Shrink();
    else
    {
        // rehash the next strings until we find empty slot
        index = (index + 1) & mask;

        while (strings[index] != NULL)
        {
            if ((keys[index] & mask) != index)
            {
                unsigned int h = Iterate(keys[index], *strings[index]);

                if (h != (unsigned int) index)
                {
                    keys[h] = keys[index];
                    strings[h] = strings[index];
                    objects[h] = objects[index];

                    strings[index] = NULL;
                    objects[index] = NULL;
                }
            }

            index = (index + 1) & mask;
        }
    }
}

void StringHash::ReadLinesFromFile(const char * filename)
{
    IFILE f = ifopen(filename, "rb");
    if (f == NULL) return;
    ReadLinesFromFile(f);
    ifclose(f);
}

void StringHash::ReadLinesFromFile(FILE * f)
{
    String buffer;
    
    while (!feof(f))
    {
        buffer.ReadLine(f);
        Add(buffer.Trim());
    }
}

void StringHash::ReadLinesFromFile(IFILE & f)
{
    String buffer;

    while (!ifeof(f))
    {
        buffer.ReadLine(f);
        Add(buffer.Trim());
    }
}

// StringIntHash implementation

StringIntHash::StringIntHash(int startsize)
    : StringHashBase()
{
    count = 0;
    size  = startsize;
    mask  = startsize - 1;

    // In this implementation, the size of hash tables must be a power of two
    if (startsize & mask)
        error("StringIntHash: Hash table size must be a power of two.\n");

    strings  = new String * [size];
    integers = new int [size];
    keys     = new unsigned int [size];

    for (unsigned int i = 0; i < size; i++)
        strings[i] = NULL;
};

StringIntHash::~StringIntHash()
{
    for (unsigned int i = 0; i < size; i++)
        if (strings[i] != NULL)
            delete strings[i];

    if(strings) delete [] strings;
    if(integers) delete [] integers;
    if(keys) delete [] keys;
}

void StringIntHash::SetSize(int newsize)
{
    int newmask = newsize - 1;

    String   ** newstrings = new String * [newsize];
    int      * newintegers = new int [newsize];
    unsigned int * newkeys = new unsigned int [newsize];

    for (int i = 0; i < newsize; i++)
        newstrings[i] = NULL;

    for (unsigned int i = 0; i < size; i++)
        if (strings[i] != NULL)
        {
            unsigned int key = keys[i];
            unsigned int h   = key & newmask;

            while (newstrings[h] != NULL &&
                   (newkeys[h] != key || (!stringsEqual(*(newstrings[h]), *(strings[i])))))
                h = (h + 1) & newmask;

            newkeys[h] = key;
            newstrings[h] = strings[i];
            newintegers[h] = integers[i];
        }

    if(strings) delete [] strings;
    if(integers) delete [] integers;
    if(keys) delete [] keys;

    strings = newstrings;
    integers = newintegers;
    keys = newkeys;
    size = newsize;
    mask = newmask;
}

void StringIntHash::Clear()
{
    for (unsigned int i = 0; i < size; i++)
        if (strings[i] != NULL)
        {
            delete strings[i];
            strings[i] = NULL;
        }

    count = 0;

    if (size > 256)
        SetSize(256);
}

int StringIntHash::Add(const String & string, int value)
{
    unsigned int key = getKey(string);
    unsigned int h   = Iterate(key, string);

    if (strings[h] == NULL)
        Insert(h, key, string);

    integers[h] = value;

    if (count * 2 > size)
    {
        Grow();
        return Iterate(key, string);
    }

    return h;
}

int StringIntHash::Find(const String & string,  int defaultValue)
{
    unsigned int key = getKey(string);
    unsigned int h   = Iterate(key, string);

    if (strings[h] == NULL)
    {
        Insert(h, key, string);
        integers[h] = defaultValue;

        if (count * 2 > size)
        {
            Grow();
            return Iterate(key, string);
        }
    }

    return h;
}

int StringIntHash::Find(const String & string) const
{
    unsigned int key = getKey(string);
    unsigned int h   = Iterate(key, string);

    if (strings[h] == NULL)
        return -1;

    return h;
}

void StringIntHash::Delete(unsigned int index)
{
    if (index >= size || strings[index] == NULL)
        return;

    delete strings[index];
    strings[index] = NULL;
    count--;

    if (count * 8 < size && size > 32)
        Shrink();
    else
    {
        // rehash the next strings until we find empty slot
        index = (index + 1) & mask;

        while (strings[index] != NULL)
        {
            if ((keys[index] & mask) != index)
            {
                unsigned int h = Iterate(keys[index], *strings[index]);

                if (h != (unsigned int) index)
                {
                    keys[h] = keys[index];
                    strings[h] = strings[index];
                    integers[h] = integers[index];

                    strings[index] = NULL;
                }
            }

            index = (index + 1) & mask;
        }
    }
}

// StringDoubleHash implementation

StringDoubleHash::StringDoubleHash(int startsize)
    : StringHashBase()
{
    count = 0;
    size  = startsize;
    mask  = startsize - 1;

    // In this implementation, the size of hash tables must be a power of two
    if (startsize & mask)
        error("StringDoubleHash: Hash table size must be a power of two.\n");

    strings  = new String * [size];
    doubles  = new double [size];
    keys     = new unsigned int [size];

    for (unsigned int i = 0; i < size; i++)
        strings[i] = NULL;
};

StringDoubleHash::~StringDoubleHash()
{
    for (unsigned int i = 0; i < size; i++)
        if (strings[i] != NULL)
            delete strings[i];

    if(strings) delete [] strings;
    if(doubles) delete [] doubles;
    if(keys) delete [] keys;
}

void StringDoubleHash::SetSize(int newsize)
{
    int newmask = newsize - 1;

    String   ** newstrings = new String * [newsize];
    double    * newdoubles = new double [newsize];
    unsigned int * newkeys = new unsigned int [newsize];

    for (int i = 0; i < newsize; i++)
        newstrings[i] = NULL;

    for (unsigned int i = 0; i < size; i++)
        if (strings[i] != NULL)
        {
            unsigned int key = keys[i];
            unsigned int h   = key & newmask;

            while (newstrings[h] != NULL &&
                    (newkeys[h] != key || (!stringsEqual(*(newstrings[h]), *(strings[i])))))
                h = (h + 1) & newmask;

            newkeys[h] = key;
            newstrings[h] = strings[i];
            newdoubles[h] = doubles[i];
        }

    if(strings) delete [] strings;
    if(doubles) delete [] doubles;
    if(keys) delete [] keys;

    strings = newstrings;
    doubles = newdoubles;
    keys = newkeys;
    size = newsize;
    mask = newmask;
}

int StringDoubleHash::Add(const String & string, double value)
{
    unsigned int key = getKey(string);
    unsigned int h   = Iterate(key, string);

    if (strings[h] == NULL)
        Insert(h, key, string);

    doubles[h] = value;

    if (count * 2 > size)
    {
        Grow();
        return Iterate(key, string);
    }

    return h;
}

int StringDoubleHash::Find(const String & string, double defaultValue)
{
    unsigned int key = getKey(string);
    unsigned int h   = Iterate(key, string);

    if (strings[h] == NULL)
    {
        Insert(h, key, string);
        doubles[h] = defaultValue;

        if (count * 2 > size)
        {
            Grow();
            return Iterate(key, string);
        }
    }

    return h;
}

int StringDoubleHash::Find(const String & string) const
{
    unsigned int key = getKey(string);
    unsigned int h   = Iterate(key, string);

    if (strings[h] == NULL)
        return -1;

    return h;
}

void StringDoubleHash::Delete(unsigned int index)
{
    if (index >= size || strings[index] == NULL)
        return;

    delete strings[index];
    strings[index] = NULL;
    count--;

    if (count * 8 < size && size > 32)
        Shrink();
    else
    {
        // rehash the next strings until we find empty slot
        index = (index + 1) & mask;

        while (strings[index] != NULL)
        {
            if ((keys[index] & mask) != index)
            {
                unsigned int h = Iterate(keys[index], *strings[index]);

                if (h != (unsigned int) index)
                {
                    keys[h] = keys[index];
                    strings[h] = strings[index];
                    doubles[h] = doubles[index];

                    strings[index] = NULL;
                }
            }

            index = (index + 1) & mask;
        }
    }
}

void StringHash::Print()
{
    Print(stdout);
}

void StringHash::Print(const char * filename)
{
    FILE * output = fopen(filename, "wt");
    if (output == NULL)
        return;
    Print(output);
    fclose(output);
}

void StringHash::Print(FILE * output)
{
    for (unsigned int i = 0; i < size; i++)
        if (SlotInUse(i))
            strings[i]->WriteLine(output);
}

String StringHash::StringList(char separator)
{
    String list;

    for (unsigned int i = 0; i < size; i++)
        if (SlotInUse(i))
            list += *strings[i] + separator;

    list.SetLength(list.Length() - 1);

    return list;
}

int StringIntHash::GetCount(const String & key) const
{
    int index = Find(key);
    return index == -1 ?  0 : integers[index];
}

int StringIntHash::IncrementCount(const String & key)
{
    int index = Find(key);

    if (index != -1)
        return ++(integers[index]);

    SetInteger(key, 1);
    return 1;
}

int StringIntHash::IncrementCount(const String & key, int amount)
{
    int index = Find(key);

    if (index != -1)
        return (integers[index] += amount);

    SetInteger(key, amount);
    return amount;
}

int StringIntHash::DecrementCount(const String & key)
{
    int index = Find(key);

    if (index != -1)
        return --(integers[index]);

    SetInteger(key, -1);
    return -1;
}

void StringDoubleHash::Clear()
{
    for (unsigned int i = 0; i < size; i++)
        if (strings[i] != NULL)
        {
            delete strings[i];
            strings[i] = NULL;
        }

    count = 0;

    if (size > 256)
        SetSize(256);
}

StringHash & StringHash::operator = (const StringHash & rhs)
{
    Clear();

    for (int i = 0; i < rhs.Capacity(); i++)
        if (rhs.SlotInUse(i))
            Add(*(rhs.strings[i]), rhs.objects[i]);

    return *this;
}

StringIntHash & StringIntHash::operator = (const StringIntHash & rhs)
{
    Clear();

    for (int i = 0; i < rhs.Capacity(); i++)
        if (rhs.SlotInUse(i))
            Add(*(rhs.strings[i]), rhs.integers[i]);

    return *this;
}

bool StringIntHash::operator == (const StringIntHash & rhs) const
{
    if (Capacity() != rhs.Capacity()) return false;
    if (Entries() != rhs.Entries()) return false;
    for (int i = 0; i < rhs.Capacity(); i++)
    {
        if(rhs.SlotInUse(i) != SlotInUse(i))
        {
            return(false);
        }
        if (rhs.SlotInUse(i))
        {
            if(*(strings[i]) != *(rhs.strings[i]))
            {
                return(false);
            }
            if(rhs.integers[i] != integers[i])
            {
                return(false);
            }
        }
    }
    return(true);
}

StringDoubleHash & StringDoubleHash::operator = (const StringDoubleHash & rhs)
{
    Clear();

    for (int i = 0; i < rhs.Capacity(); i++)
        if (rhs.SlotInUse(i))
            Add(*(rhs.strings[i]), rhs.doubles[i]);

    return *this;
}

void StringHash::Swap(StringHash & s)
{
    String ** tstrings = s.strings;
    s.strings = strings;
    strings = tstrings;

    void ** tobjects = s.objects;
    s.objects = objects;
    objects = tobjects;

    unsigned int * tkeys = s.keys;
    s.keys = keys;
    keys = tkeys;

    unsigned int temp = s.count;
    s.count = count;
    count = temp;

    temp = s.size;
    s.size = size;
    size = temp;

    temp = s.mask;
    s.mask = mask;
    mask = temp;
}

