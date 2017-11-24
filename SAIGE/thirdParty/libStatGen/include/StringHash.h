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

#ifndef __STRINGHASH_H__
#define __STRINGHASH_H__

#include "StringBasics.h"
#include "Constant.h"
#include "Hash.h"


class StringHashBase
{
public:
    inline void setCaseSensitive(bool caseSensitive) {myCaseSensitive = caseSensitive;}
    StringHashBase()
        : myCaseSensitive(false)
    {}

    virtual ~StringHashBase() {}

    // Make pure virtual
    virtual void SetSize(int newsize) = 0;

protected:
    inline bool stringsEqual(const String& string1, const String& string2) const
    {
        if(myCaseSensitive)
        {
            // Case sensitive is faster.
            return(string1.FastCompare(string2) == 0);
        }
        // Case insensitive - slow compare - convert to same case.
        return(string1.SlowCompare(string2) == 0);
    }

    inline unsigned int getKey(const String& string) const
    {
        if(myCaseSensitive)
        {
            return(hash(string.uchar(), string.Length(), 0));
        }
        // Case insensitive.
        return(hash_no_case(string.uchar(), string.Length(), 0));
    }

    bool myCaseSensitive;
};


class StringHash : public StringHashBase
{
protected:
    String        ** strings;
    void          ** objects;
    unsigned int      * keys;
    unsigned int count, size;
    unsigned int        mask;

public:
    StringHash(int startsize = 32);
    virtual ~StringHash();

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
    void * Object(const String & key) const
    {
        int index = Find(key);

        return index >= 0 ? objects[index] : NULL;
    }
    void * Object(const String & key, void *(*create_object)())
    {
        int index = Find(key, create_object);

        return objects[index];
    }

    void SetObject(int i, void * object)
    {
        objects[i] = object;
    }
    void SetObject(const String & key, void * object)
    {
        Add(key, object);
    }

    int Add(const String & s, void * object = NULL);
    int Find(const String & s, void *(*create_object)() = NULL);
    int Find(const String & s) const;

    StringHash & operator = (const StringHash & rhs);

    const String & operator [](int i) const
    {
        return *(strings[i]);
    }
    String & operator [](int i)
    {
        return *(strings[i]);
    }
//      String & String(int i) { return *(strings[i]); }

    static void * CreateHash();

    void Delete(unsigned int index);
    void Delete(const String & key)
    {
        Delete(Find(key));
    }

    bool SlotInUse(int index) const
    {
        return strings[index] != NULL;
    }

    void Print();
    void Print(FILE * file);
    void Print(const char * filename);

    String StringList(char separator = ',');

    // Initialize hash with the contents of a file
    void ReadLinesFromFile(FILE * file);
    void ReadLinesFromFile(const char * filename);

    void ReadLinesFromFile(IFILE & file);

    void Swap(StringHash & s);

private:

    unsigned int Iterate(unsigned int key, const String & string) const
    {
        unsigned int h = key & mask;

        while (strings[h] != NULL &&
               (keys[h] != key ||
                (!stringsEqual(*(strings[h]), string))))
            h = (h + 1) & mask;

        return h;
    }

    void Insert(unsigned int where, unsigned int key, const String & string)
    {
        strings[where] = new String;
        *(strings[where]) = string;
        keys[where] = key;

        count++;
    }
};

class StringIntHash : public StringHashBase
{
protected:
    String        ** strings;
    int           * integers;
    unsigned int      * keys;
    unsigned int count, size;
    unsigned int        mask;

public:
    StringIntHash(int startsize = 32);
    virtual ~StringIntHash();

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

    int Integer(int i) const
    {
        return integers[i];
    }
    int Integer(const String & key) const
    {
        int index = Find(key);

        return index >= 0 ? integers[index] : -1;
    }

    void SetInteger(int i, int value)
    {
        integers[i] = value;
    }
    void SetInteger(const String & key, int value)
    {
        Add(key, value);
    }

    int IncrementCount(const String & key);
    int IncrementCount(const String & key, int amount);
    int DecrementCount(const String & key);
    int GetCount(const String & key) const;
    int GetCount(int index) const
    {
        return integers[index];
    }

    int Add(const String & s, int integer);
    int Find(const String & s, int defaultValue);
    int Find(const String & s) const;

    StringIntHash & operator = (const StringIntHash & rhs);
    bool operator == (const StringIntHash & rhs) const;

    const String & operator [](int i) const
    {
        return *(strings[i]);
    }
    String & operator [](int i)
    {
        return *(strings[i]);
    }
//      String & String(int i) { return *(strings[i]); }

    void Delete(unsigned int index);
    void Delete(const String & key)
    {
        Delete(Find(key));
    }

    bool SlotInUse(int index) const
    {
        return strings[index] != NULL;
    }

private:

    unsigned int Iterate(unsigned int key, const String & string) const
    {
        unsigned int h = key & mask;

        while (strings[h] != NULL &&
               (keys[h] != key ||
                (!stringsEqual(*(strings[h]), string))))
            h = (h + 1) & mask;

        return h;
    }

    void Insert(unsigned int where, unsigned int key, const String & string)
    {
        strings[where] = new String;
        *(strings[where]) = string;
        keys[where] = key;

        count++;
    }
};

class StringDoubleHash : public StringHashBase
{
protected:
    String        ** strings;
    double         * doubles;
    unsigned int      * keys;
    unsigned int count, size;
    unsigned int        mask;

public:
    StringDoubleHash(int startsize = 32);
    virtual ~StringDoubleHash();

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

    double Double(int i) const
    {
        return doubles[i];
    }
    double Double(const String & key) const
    {
        int index = Find(key);

        return index >= 0 ? doubles[index] : _NAN_;
    }

    void SetDouble(int i, double value)
    {
        doubles[i] = value;
    }
    void SetDouble(const String & key, double value)
    {
        Add(key, value);
    }

    int Add(const String & s, double value);
    int Find(const String & s, double defaultValue);
    int Find(const String & s) const;

    StringDoubleHash & operator = (const StringDoubleHash & rhs);

    const String & operator [](int i) const
    {
        return *(strings[i]);
    }
    String & operator [](int i)
    {
        return *(strings[i]);
    }
//      String & String(int i) { return *(strings[i]); }

    void Delete(unsigned int index);
    void Delete(const String & key)
    {
        Delete(Find(key));
    }

    bool SlotInUse(int index) const
    {
        return strings[index] != NULL;
    }

private:

    unsigned int Iterate(unsigned int key, const String & string) const
    {
        unsigned int h = key & mask;

        while (strings[h] != NULL &&
               (keys[h] != key ||
                (!stringsEqual(*(strings[h]), string))))
            h = (h + 1) & mask;

        return h;
    }

    void Insert(unsigned int where, unsigned int key, const String & string)
    {
        strings[where] = new String;
        *(strings[where]) = string;
        keys[where] = key;

        count++;
    }
};


#endif
