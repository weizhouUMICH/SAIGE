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

#ifndef __STRING_ARRAY_H__
#define __STRING_ARRAY_H__

#include "StringBasics.h"

class StringArray
{
protected:
    String ** strings;
    int size, count;

public:
    static int alloc;
    static bool lazyMemoryManagement;

    StringArray(int startsize = 0);
    StringArray(StringArray & original);
    virtual ~StringArray();

    // Each line in a file is parsed into a separate array element
    //

    void Read(FILE * f);
    void Write(FILE * f);
    void WriteLine(FILE * f);
    void Read(const char * filename);
    void Write(const char * filename);
    void WriteLine(const char * filename);

    void Read(IFILE & f);

    // Write all strings to the screen
    void Print();
    void PrintLine();

    // Write all strings to a file
    void Print(FILE * f);
    void PrintLine(FILE * f);

    void Grow(int newsize);
    void Clear();

    int Length() const
    {
        return count;
    }
    int Dimension(int newcount);
    int CharLength();

    String & operator [](int i)
    {
        return *(strings[i]);
    }
    const String & operator [](int i) const
    {
        return *(strings[i]);
    }

    // These functions divide a string into tokens and append these to the
    // array. Return value is the new array length
    //

    int AddColumns(const String & s, char ch = '\t');
    int AddColumns(const String & s, char ch, int maxColumns);
    int AddTokens(const String & s, char ch);
    int AddTokens(const String & s, const String & separators = " \t\r\n");

    int ReplaceColumns(const String & s, char ch = '\t')
    {
        Clear();
        return AddColumns(s, ch);
    }
    int ReplaceTokens(const String & s, const String & separators = " \t\r\n")
    {
        Clear();
        return AddTokens(s, separators);
    }

    // These functions add, insert or remove a single array element
    //

    int  Add(const String & s);
    void InsertAt(int position, const String & s);
    void Delete(int position);

    // These functions manipulate a string as a stack
    //

    String & Last() const;
    int      Push(const String & s)
    {
        return Add(s);
    }
    String   Pop();

    // Linear search (N/2 comparisons on average) for a single element
    // If searching is required, StringMaps are a better option
    //

    int Find(const String & s) const;
    int FastFind(const String & s) const;
    int SlowFind(const String & s) const;

    // Alphetically orders strings
    //
    void Sort();

    // Trims strings to remove whitespace
    void Trim();

    StringArray & operator = (const StringArray & rhs);

    bool operator == (const StringArray & rhs) const;
    bool operator != (const StringArray & rhs) const
    {
        return !(*this == rhs);
    }

    void Swap(StringArray & s);

private:
    static int ComparisonForSort(const void * a, const void * b);
};

#endif

