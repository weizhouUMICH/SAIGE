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

#ifndef __STRINGMAP_H__
#define __STRINGMAP_H__

#include "StringBasics.h"

class StringMap
{
protected:
    ::String ** strings;
    void     ** objects;
    int         count, size;

public:
    static int alloc;

    StringMap(int startsize = 0);
    virtual ~StringMap();

    void Grow(int newsize);
    void Clear();
    int  Length() const
    {
        return count;
    }

    void * Object(int i) const
    {
        return objects[i];
    }
    void * Object(const ::String & key) const
    {
        int index = Find(key);
        return (index >= 0) ? objects[index] : NULL;
    }
    void * Object(const ::String & key, void *(*create_object)())
    {
        return objects[Find(key, create_object)];
    }

    void SetObject(int i, void * object)
    {
        objects[i] = object;
    }
    void SetObject(const ::String & key, void * object)
    {
        Add(key, object);
    }

    int Add(const ::String & s, void * object = NULL);
    int Find(const ::String & s, void *(*create_object)() = NULL);
    int Find(const ::String & s) const;
    int FindStem(const ::String & stem) const;
    int FindFirstStem(const ::String & stem) const;

    StringMap & operator = (const StringMap & rhs);

    const ::String & operator [](int i) const
    {
        return *(strings[i]);
    }
    ::String & operator [](int i)
    {
        return *(strings[i]);
    }
    ::String & String(int i)
    {
        return *(strings[i]);
    }

    static void * CreateMap();

    void Delete(int index);
};

class StringIntMap
{
protected:
    ::String ** strings;
    int       * integers;
    int         count, size;

public:
    static int alloc;

    StringIntMap(int startsize = 0);
    virtual ~StringIntMap();

    void Grow(int newsize);
    void Clear();
    int  Length() const
    {
        return count;
    }

    int Integer(int i) const
    {
        return integers[i];
    }
    int Integer(const ::String & key) const
    {
        int index = Find(key);
        return (index >= 0) ? (int) integers[index] : -1;
    }

    void SetInteger(int i, int value)
    {
        integers[i] = value;
    }
    void SetInteger(const ::String & key, int value)
    {
        Add(key, value);
    }

    int Add(const ::String & s, int i);
    int Find(const ::String & s, int defaultValue);
    int Find(const ::String & s) const;
    int FindStem(const ::String & stem) const;

    StringIntMap & operator = (const StringIntMap & rhs);

    const ::String & operator [](int i) const
    {
        return *(strings[i]);
    }
    ::String & operator [](int i)
    {
        return *(strings[i]);
    }
    ::String & String(int i)
    {
        return *(strings[i]);
    }

    static void * CreateMap();

    int IncrementCount(const ::String & key);
    int DecrementCount(const ::String & key);
    int GetCount(const ::String & key) const;
    int GetCount(int index) const
    {
        return integers[index];
    }

    void Delete(int index);
};

#endif

