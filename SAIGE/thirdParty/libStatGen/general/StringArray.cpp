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

#include "StringArray.h"
#include "InputFile.h"
#include "Sort.h"
#include "Error.h"

#include <string.h>

int   StringArray::alloc = 32;
bool  StringArray::lazyMemoryManagement = false;

StringArray::StringArray(int startsize)
{
    count = startsize;
    size = (startsize + alloc) / alloc * alloc;
    strings = new String * [size];
    for (int i = 0; i < count; i++)
        strings[i] = new String;
    for (int i = count; i < size; i++)
        strings[i] = NULL;
};

StringArray::StringArray(StringArray & rhs)
{
    count = rhs.count;
    size = (rhs.count + alloc) / alloc * alloc;
    strings = new String * [size];

    for (int i = 0; i < count; i++)
        strings[i] = new String(rhs[i]);;
    for (int i = count; i < size; i++)
        strings[i] = NULL;
}

StringArray::~StringArray()
{
    for (int i = 0; i < size; i++)
        if (strings[i] != NULL)
            delete strings[i];
        else
            break;

    delete [] strings;
}

int StringArray::CharLength()
{
    int charlen = 0;
    for (int i = 0; i < count; i++)
        charlen += strings[i]->Length();
    return charlen;
}

void StringArray::Read(const char * filename)
{
    IFILE f = ifopen(filename, "rb");
    if (f == NULL) return;
    Read(f);
    ifclose(f);
}

void StringArray::Write(const char * filename)
{
    FILE * f = fopen(filename, "wt");
    if (f == NULL) return;
    Write(f);
    fclose(f);
}

void StringArray::WriteLine(const char * filename)
{
    FILE * f = fopen(filename, "wt");
    if (f == NULL) return;
    WriteLine(f);
    fclose(f);
}

void StringArray::Read(FILE * f)
{
    while (!feof(f))
    {
        Grow(count + 1);
        if (strings[count] == NULL)
            strings[count] = new String;
        strings[count]->ReadLine(f);
        count++;
    }
}

void StringArray::Write(FILE * f)
{
    for (int i = 0; i < count; i++)
        strings[i]->WriteLine(f);
}

void StringArray::WriteLine(FILE * f)
{
    for (int i = 0; i < count; i++)
        fprintf(f, "%s%c", (const char *)(*strings[i]), i == count-1 ? '\n' : '\t');
}

void StringArray::Read(IFILE & f)
{
    while (!ifeof(f))
    {
        Grow(count + 1);
        if (strings[count] == NULL)
            strings[count] = new String;
        strings[count]->ReadLine(f);
        if (ifeof(f) && strings[count]->Length()==0)
        {
            return;
        }
        count++;
    }
}

void StringArray::Grow(int newsize)
{
    if (newsize >= size)
    {
        int oldsize = size;

        if ((newsize >> 1) >= size)
            size = (newsize + alloc) / alloc * alloc;
        else
        {
            size = alloc;
            while (size <= newsize)
                size *= 2;
        }
        String ** tmp = new String * [size];
        for (int i = 0; i < oldsize; i++)
            tmp[i] = strings[i];
        for (int i = oldsize; i < size; i++)
            tmp[i] = NULL;
        delete [] strings;
        strings = tmp;
    }
}

void StringArray::Clear()
{
    if (!lazyMemoryManagement)
    {
        for (int i = 0; i < size; i++)
        {
            if (strings[i] != NULL)
            {
                delete strings[i];
                strings[i] = NULL;
            }
            else
            {
                break;
            }
        }
    }
    count = 0;
}

int StringArray::AddColumns(const String & s, char ch)
{
    if (s.Length() > 0)
        for (int pos = 0; pos <= s.Length(); pos++)
        {
            int oldpos = pos;
            pos = s.FindChar(ch, pos);
            if (pos == -1) pos = s.Length();
            Grow(count + 1);

            if (strings[count] == NULL)
            {
                strings[count] = new String(pos - oldpos);
            }
            strings[count]->SetLength(pos - oldpos);
            memcpy((char *) *strings[count++], ((const char *) s) + oldpos, pos - oldpos);
        }

    return count;
}

int StringArray::AddColumns(const String & s, char ch, int maxColumns)
{
    maxColumns += count;

    if (s.Length() > 0)
        for (int pos = 0; pos <= s.Length() && maxColumns != count; pos++)
        {
            int oldpos = pos;
            pos = s.FindChar(ch, pos);
            if (pos == -1) pos = s.Length();
            Grow(count + 1);

            if (strings[count] == NULL)
                strings[count] = new String(pos - oldpos);
            strings[count]->SetLength(pos - oldpos);
            memcpy((char *) *strings[count++], ((const char *) s) + oldpos, pos - oldpos);
        };

    return count;
}

int StringArray::AddTokens(const String & s, char ch)
{
    for (int pos = 0; pos < s.Length(); pos++)
    {
        while (pos < s.Length() && s[pos] == ch) pos++;
        int oldpos = pos;

        while (pos < s.Length() && s[pos] != ch) pos++;

        if (oldpos < s.Length())
        {
            Grow(count + 1);
            if (strings[count] == NULL)
            {
                strings[count] = new String(pos - oldpos);
            }
            strings[count]->SetLength(pos - oldpos);
            memcpy((char *) *strings[count++], (const char *) s + oldpos, pos - oldpos);
        }
    }

    return count;
}

int StringArray::AddTokens(const String & s, const String & separators)
{
    for (int pos = 0; pos < s.Length(); pos++)
    {
        while (pos < s.Length() && separators.FindChar(s[pos]) != -1) pos++;
        int oldpos = pos;

        while (pos < s.Length() && separators.FindChar(s[pos]) == -1) pos++;

        if (oldpos < s.Length())
        {
            Grow(count + 1);
            if (strings[count] == NULL)
                strings[count] = new String(pos - oldpos);
            strings[count]->SetLength(pos - oldpos);
            memcpy((char *) *strings[count++], ((const char *) s) + oldpos, pos - oldpos);
        }
    }

    return count;
}

int StringArray::Dimension(int newcount)
{
    if (newcount > count)
    {
        Grow(newcount);
        for (int i = count; i < newcount; i++)
        {
            if (strings[i] == NULL)
                strings[i] = new String;
            else
                strings[i]->Clear();
        }
        count = newcount;
    }
    else if (newcount < count)
    {
        if (!lazyMemoryManagement)
        {
            for (int i = newcount; i < size; i++)
            {
                if (strings[i] != NULL)
                {
                    delete strings[i];
                    strings[i] = NULL;
                }
                else
                {
                    break;
                }
            }
        }
        count = newcount;
    }

    return count;
}

int StringArray::Find(const String & s) const
{
    for (int i = 0; i < count; i++)
        if (*(strings[i]) == s)
            return i;
    return -1;
}

int StringArray::FastFind(const String & s) const
{
    for (int i = 0; i < count; i++)
        if (strings[i]->FastCompare(s) == 0)
            return i;
    return -1;
}

int StringArray::SlowFind(const String & s) const
{
    for (int i = 0; i < count; i++)
        if (strings[i]->SlowCompare(s) == 0)
            return i;
    return -1;
}

int StringArray::Add(const String & s)
{
    Grow(count + 1);
    if (strings[count] == NULL)
    {
        strings[count] = new String(s);
    }
    else
    {
        *strings[count] = s;
    }
    return ++count;
}

void StringArray::InsertAt(int position, const String & s)
{
    Grow(count + 1);

    String * newString = strings[count];
    if (newString == NULL)
        newString = new String(s);
    else
        *newString = s;

    for (int i = count; i > position; i--)
        strings[i] = strings[i - 1];
    strings[position] = newString;
    count++;
}

String & StringArray::Last() const
{
    if (!count) error("StringArray: Null String Access");
    return *(strings[count - 1]);
}

void StringArray::Delete(int index)
{
    String * oldString = strings[index];

    count--;
    for (; index < count; index++)
        strings[index] = strings[index + 1];
    strings[count] = oldString;
}

StringArray & StringArray::operator = (const StringArray & rhs)
{
    Dimension(rhs.count);
    for (int i = 0; i < rhs.count; i++)
        *strings[i] = *rhs.strings[i];
    return *this;
}

bool StringArray::operator == (const StringArray & rhs) const
{
    if (count != rhs.count) return false;
    for (int i = 0; i < rhs.count; i++)
        if (*strings[i] != *rhs.strings[i])
            return false;
    return true;
}

void StringArray::Sort()
{
    QuickSort(strings, count, sizeof(String *), ComparisonForSort);
}

int StringArray::ComparisonForSort(const void * a, const void * b)
{
    String * string1 = *(String **) a;
    String * string2 = *(String **) b;

    return Compare(*string1, *string2);
}

String StringArray::Pop()
{
    String result = *(strings[count - 1]);

    Dimension(count - 1);

    return result;
}

void StringArray::Trim()
{
    for (int i = 0; i < count; i++)
        strings[i]->Trim();
}

void StringArray::Print()
{
    Print(stdout);
}

void StringArray::Print(FILE * output)
{
    for (int i = 0; i < count; i++)
        fprintf(output, "%s\n", (const char *)(*strings[i]));
}

void StringArray::PrintLine()
{
    PrintLine(stdout);
}

void StringArray::PrintLine(FILE * output)
{
    for (int i = 0; i < count; i++)
        fprintf(output, "%s%c", (const char *)(*strings[i]), i == count - 1 ? '\n' : '\t');
}

void StringArray::Swap(StringArray & s)
{
    String ** temp = s.strings;
    s.strings = strings;
    strings = temp;

    int swap = s.size;
    s.size = size;
    size = swap;

    swap = s.count;
    s.count = count;
    count = swap;
}

