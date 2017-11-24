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

#include "StringMap.h"

int StringMap::alloc = 8;

StringMap::StringMap(int startsize)
{
    count = 0;
    size = (startsize + alloc) / alloc * alloc;
    strings = new ::String * [size];
    objects = new void * [size];
};

StringMap::~StringMap()
{
    for (int i = 0; i < count; i++)
        delete strings[i];
    delete [] strings;
    delete [] objects;
}

void StringMap::Grow(int newsize)
{
    if (newsize >= size)
    {
        if ((newsize >> 1) >= size)
            size = (newsize + alloc) / alloc * alloc;
        else
        {
            size = alloc;
            while (size <= newsize)
                size *= 2;
        }

        size = (newsize + alloc) / alloc * alloc;

        ::String ** newStrings = new ::String * [size];
        void     ** newObjects = new void * [size];

        for (int i = 0; i < count; i++)
        {
            newStrings[i] = strings[i];
            newObjects[i] = objects[i];
        }

        delete [] strings;
        delete [] objects;

        strings = newStrings;
        objects = newObjects;
    }
}

int StringMap::Add(const ::String & key, void * object)
{
    if (count == 0)
    {
        Grow(1);
        strings[0] = new ::String(key);
        objects[0] = object;
        return count++;
    }

    int left = 0;
    int right = count - 1;

    while (right > left)
    {
        int probe = (left + right) / 2;
        int test  = key.SlowCompare(*(strings[probe]));

        if (test == 0)
        {
            objects[probe] = object;
            return probe;
        }

        if (test < 0)
            right = probe - 1;
        else
            left  = probe + 1;
    }

    int insertAt = left;
    int test = key.SlowCompare(*(strings[insertAt]));

    if (test == 0)
    {
        objects[insertAt] = object;
        return insertAt;
    }

    if (test > 0) insertAt++;

    Grow(count + 1);

    if (insertAt < count)
    {
        for (int i = count; i > insertAt; i--)
        {
            strings[i] = strings[i - 1];
            objects[i] = objects[i - 1];
        }
    }

    strings[insertAt] = new ::String(key);
    objects[insertAt] = object;
    count++;

    return insertAt;
}

int StringMap::Find(const ::String & s,  void *(*create_object)())
{
    if (!count)
        return create_object == NULL ? -1 : Add(s, create_object());

    int left = 0;
    int right = count - 1;

    while (right > left)
    {
        int probe = (left + right) / 2;
        int test  = s.SlowCompare(*(strings[probe]));

        if (test == 0)
            return probe;

        if (test < 0)
            right = probe - 1;
        else
            left  = probe + 1;
    }

    int position = left;
    int test = s.SlowCompare(*(strings[left]));

    if (test == 0)
        return position;

    if (create_object == NULL)
        return -1;

    if (test > 0)
        position++;

    Grow(count + 1);

    if (position < count)
    {
        for (int i = count; i > position; i--)
        {
            strings[i] = strings[i - 1];
            objects[i] = objects[i - 1];
        }
    }

    strings[position] = new ::String(s);
    objects[position] = create_object();
    count++;

    return position;
}

int StringMap::Find(const ::String & s) const
{
    if (!count) return -1;

    int left = 0;
    int right = count - 1;

    while (right > left)
    {
        int probe = (left + right) / 2;
        int test  = s.SlowCompare(*(strings[probe]));

        if (test == 0)
            return probe;

        if (test < 0)
            right = probe - 1;
        else
            left  = probe + 1;
    }

    int position = left;
    int test = s.SlowCompare(*(strings[left]));

    if (test == 0)
        return position;

    return -1;
}

int StringMap::FindStem(const ::String & stem) const
{
    if (!count) return -1;

    int left = 0;
    int right = count - 1;

    while (right > left)
    {
        int probe = (left + right) / 2;
        int test  = strings[probe]->SlowCompareToStem(stem);

        if (test == 0)
        {
            if ((left  < probe && strings[probe-1]->SlowCompareToStem(stem) == 0) ||
                    (right > probe && strings[probe+1]->SlowCompareToStem(stem) == 0))
                return -2;

            return probe;
        }

        if (test > 0)
            right = probe - 1;
        else
            left  = probe + 1;
    }

    if (strings[left]->SlowCompareToStem(stem) == 0)
        return left;

    return -1;
}

int StringMap::FindFirstStem(const ::String & stem) const
{
    if (!count) return -1;

    int left = 0;
    int right = count - 1;

    while (right > left)
    {
        int probe = (left + right) / 2;
        int test  = strings[probe]->SlowCompareToStem(stem);

        if (test == 0)
        {
            while (left < probe && strings[probe-1]->SlowCompareToStem(stem) == 0)
                probe--;

            return probe;
        }

        if (test > 0)
            right = probe - 1;
        else
            left  = probe + 1;
    }

    if (strings[left]->SlowCompareToStem(stem) == 0)
        return left;

    return -1;
}

void * StringMap::CreateMap()
{
    return (void *) new StringMap();
}

void StringMap::Clear()
{
    for (int i = 0; i < count; i++)
        delete strings[i];
    count = 0;
}

void StringMap::Delete(int index)
{
    count--;

    delete strings[index];

    for (int i = index; i < count; i++)
    {
        strings[i] = strings[i+1];
        objects[i] = objects[i+1];
    }
}

// StringIntMap class
//

int StringIntMap::alloc = 8;

StringIntMap::StringIntMap(int startsize)
{
    count = 0;
    size = (startsize + alloc) / alloc * alloc;
    strings = new ::String * [size];
    integers = new int[size];
};

StringIntMap::~StringIntMap()
{
    for (int i = 0; i < count; i++)
        delete strings[i];
    delete [] strings;
    delete [] integers;
}

void StringIntMap::Grow(int newsize)
{
    if (newsize >= size)
    {
        if ((newsize >> 1) >= size)
            size = (newsize + alloc) / alloc * alloc;
        else
        {
            size = alloc;
            while (size <= newsize)
                size *= 2;
        }

        ::String ** newStrings = new ::String * [size];
        int       * newIntegers = new int [size];

        for (int i = 0; i < count; i++)
        {
            newStrings[i] = strings[i];
            newIntegers[i] = integers[i];
        }

        delete [] strings;
        delete [] integers;

        strings = newStrings;
        integers = newIntegers;
    }
}

int StringIntMap::Add(const ::String & key, int integer)
{
    if (count == 0)
    {
        Grow(1);
        strings[0] = new ::String(key);
        integers[0] = integer;
        return count++;
    }

    int left = 0;
    int right = count - 1;

    while (right > left)
    {
        int probe = (left + right) / 2;
        int test  = key.SlowCompare(*(strings[probe]));

        if (test == 0)
        {
            integers[probe] = integer;
            return probe;
        }

        if (test < 0)
            right = probe - 1;
        else
            left  = probe + 1;
    }

    int insertAt = left;
    int test = key.SlowCompare(*(strings[insertAt]));

    if (test == 0)
    {
        integers[insertAt] = integer;
        return insertAt;
    }

    if (test > 0) insertAt++;

    Grow(count + 1);

    if (insertAt < count)
    {
        for (int i = count; i > insertAt; i--)
        {
            strings[i] = strings[i - 1];
            integers[i] = integers[i - 1];
        }
    }

    strings[insertAt] = new ::String(key);
    integers[insertAt] = integer;
    count++;

    return insertAt;
}

int StringIntMap::Find(const ::String & s, int defaultValue)
{
    if (!count)
        return Add(s, defaultValue);

    int left = 0;
    int right = count - 1;

    while (right > left)
    {
        int probe = (left + right) / 2;
        int test  = s.SlowCompare(*(strings[probe]));

        if (test == 0)
            return probe;

        if (test < 0)
            right = probe - 1;
        else
            left  = probe + 1;
    }

    int position = left;
    int test = s.SlowCompare(*(strings[left]));

    if (test == 0)
        return position;

    if (test > 0)
        position++;

    Grow(count + 1);

    if (position < count)
    {
        for (int i = count; i > position; i--)
        {
            strings[i] = strings[i - 1];
            integers[i] = integers[i - 1];
        }
    }

    strings[position] = new ::String(s);
    integers[position] = defaultValue;
    count++;

    return position;
}

int StringIntMap::Find(const ::String & s) const
{
    if (!count) return -1;

    int left = 0;
    int right = count - 1;

    while (right > left)
    {
        int probe = (left + right) / 2;
        int test  = s.SlowCompare(*(strings[probe]));

        if (test == 0)
            return probe;

        if (test < 0)
            right = probe - 1;
        else
            left  = probe + 1;
    }

    int position = left;
    int test = s.SlowCompare(*(strings[left]));

    if (test == 0)
        return position;

    return -1;
}

int StringIntMap::FindStem(const ::String & stem) const
{
    if (!count) return -1;

    int left = 0;
    int right = count - 1;

    while (right > left)
    {
        int probe = (left + right) / 2;
        int test  = strings[probe]->SlowCompareToStem(stem);

        if (test == 0)
        {
            if ((left  < probe && strings[probe-1]->SlowCompareToStem(stem) == 0) ||
                    (right > probe && strings[probe+1]->SlowCompareToStem(stem) == 0))
                return -2;

            return probe;
        }

        if (test > 0)
            right = probe - 1;
        else
            left  = probe + 1;
    }

    if (strings[left]->SlowCompareToStem(stem) == 0)
        return left;

    return -1;
}

void StringIntMap::Clear()
{
    for (int i = 0; i < count; i++)
        delete strings[i];
    count = 0;
}

int StringIntMap::GetCount(const ::String & key) const
{
    int index = Find(key);
    return index == -1 ?  0 : integers[index];
}

int StringIntMap::IncrementCount(const ::String & key)
{
    int index = Find(key);

    if (index != -1)
        return ++(integers[index]);

    SetInteger(key, 1);
    return 1;
}

int StringIntMap::DecrementCount(const ::String & key)
{
    int index = Find(key);

    if (index != -1)
        return --(integers[index]);

    SetInteger(key, -1);
    return -1;
}

void StringIntMap::Delete(int index)
{
    count--;

    delete strings[index];

    for (int i = index; i < count; i++)
    {
        strings[i] = strings[i+1];
        integers[i] = integers[i+1];
    }
}



