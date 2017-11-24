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

#ifndef __BASICSTRING_H__
#define __BASICSTRING_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <vector>
#include <iostream>

using std::vector;

#define READBUF     128
#define READBUFSTR  "128"

#ifdef __PREFIX_STRING__
#define String      BasicString
#endif

#include "InputFile.h"

class String
{
private:
    void NewString(int startsize);

protected:
    char * buffer;
    int  len, size;

public:
    static int  alloc;
    static bool caseSensitive;

    explicit String(int startsize = 0)
    {
        NewString(startsize);
    }
    String(const char * s);
    String(const String & s);
    explicit String(char ch, int count = 1);

    ~String()
    {
        if(buffer != NULL)
            delete [] buffer;
    }

    String & Clear()
    {
        len = buffer[0] = 0;
        return *this;
    }

    String & Copy(const String & s);
    String & Copy(const String & s, int start, int count);
    String & Copy(const char * s);

    bool     IsEmpty() const
    {
        return len == 0;
    }
    String & ToUpper();
    String & ToLower();
    String   AsUpper();
    String   AsLower();
    String   Capitalize();
    String & Reverse();

    String & LeftClip(int clipAmount);
    String & RightClip(int clipAmount);

    String & operator = (char ch);
    String   operator + (char ch) const;
    String & operator += (char ch);

    String & operator = (const String & rhs);
    String   operator + (const String & rhs) const;
    String & operator += (const String & rhs);

    String & operator = (const char * rhs);
    String   operator + (const char * rhs) const;
    String & operator += (const char * rhs);

    String & operator = (int rhs);
    String   operator + (int rhs) const;
    String & operator += (int rhs);

    String & operator = (double rhs);
    String   operator + (double rhs) const;
    String & operator += (double rhs);

    void appendFullFloat(float rhs);

    String & operator = (unsigned int rhs);
    String   operator + (unsigned int rhs) const;
    String & operator += (unsigned int rhs);
    String   operator *(unsigned int rhs) const;
    String & operator *= (unsigned int rhs);

    int Compare(const String & rhs) const;
    int FastCompare(const String & rhs) const;
    int SlowCompare(const String & rhs) const;

    int Compare(const char * rhs) const;
    int FastCompare(const char * rhs) const;
    int SlowCompare(const char * rhs) const;

    int CompareToStem(const String & stem) const;
    int FastCompareToStem(const String & stem) const;
    int SlowCompareToStem(const String & stem) const;

    int CompareToStem(const char * stem) const;
    int FastCompareToStem(const char * stem) const;
    int SlowCompareToStem(const char * stem) const;

    int MatchesBeginningOf(const String & stem) const;
    int FastMatchesBeginningOf(const String & stem) const;
    int SlowMatchesBeginningOf(const String & stem) const;

    int MatchesBeginningOf(const char * stem) const;
    int FastMatchesBeginningOf(const char * stem) const;
    int SlowMatchesBeginningOf(const char * stem) const;

    int operator == (const String & rhs) const
    {
        return Compare(rhs) == 0;
    }
    int operator != (const String & rhs) const
    {
        return Compare(rhs) != 0;
    }
    int operator < (const String & rhs) const
    {
        return Compare(rhs)  < 0;
    }
    int operator > (const String & rhs) const
    {
        return Compare(rhs)  > 0;
    }
    int operator >= (const String & rhs) const
    {
        return Compare(rhs) >= 0;
    }
    int operator <= (const String & rhs) const
    {
        return Compare(rhs) <= 0;
    }

    int operator == (const char * rhs) const
    {
        return Compare(rhs) == 0;
    }
    int operator != (const char * rhs) const
    {
        return Compare(rhs) != 0;
    }
    int operator < (const char * rhs) const
    {
        return Compare(rhs)  < 0;
    }
    int operator > (const char * rhs) const
    {
        return Compare(rhs)  > 0;
    }
    int operator <= (const char * rhs) const
    {
        return Compare(rhs) <= 0;
    }
    int operator >= (const char * rhs) const
    {
        return Compare(rhs) >= 0;
    }

    operator const char *() const
    {
        return buffer;
    }
    const char *c_str() const
    {
        return (const char *) buffer;
    }
    operator char *()
    {
        return buffer;
    }

    operator int () const
    {
        return atoi(buffer);
    }
    operator double() const
    {
        return atof(buffer);
    }
    operator long double() const;

    char operator [](int i)  const
    {
        return buffer[i];
    }
    char & operator [](int i)
    {
        return buffer[i];
    }

    char & Last()
    {
        return buffer[len - 1];
    }
    char & First()
    {
        return buffer[0];
    }

    void Grow(int newSize);
    void Swap(String & s);

    char * LockBuffer(int size = -1);
    String & UnlockBuffer();

    String & Read();
    // Return the status.  A negative number indicates an error/EOF.
    int ReadLine();
    void     WriteLine();
    void     Write();

    String & Read(FILE * f);
    // Return the status.  A negative number indicates an error/EOF.
    int ReadLine(FILE * f);
    void     WriteLine(FILE * f);
    void     Write(FILE * f);

    String & Read(IFILE & f);

    // Read a line using getc
    // Return the status.  A negative number indicates an error/EOF.
    int ReadLine(IFILE & f);

    String Left(int count) const;
    String Right(int count) const;
    String Mid(int start, int end) const;
    String SubStr(int start, int count) const;
    String SubStr(int start) const;

    int FindChar(char ch, int start = 0) const;
    int FastFindChar(char ch, int start = 0) const;
    int SlowFindChar(char ch, int start = 0) const;

    int FindLastChar(char ch) const;
    int FastFindLastChar(char ch) const;
    int SlowFindLastChar(char ch) const;

    // Since there is no longer implicit conversion
    // from char to String, declare this method that
    // takes a character rather than a String reference.
    int Find(char ch, int start = 0) const
    {
        return(FindChar(ch, start));
    }
    int Find(const String & str, int start = 0) const;
    int FastFind(const String & str, int start = 0) const;
    int SlowFind(const String & str, int start = 0) const;

    String & Filter(const String & s);
    String & Filter(const char * s);

    String & ExcludeCharacters(const String & s);
    String & ExcludeCharacters(const char * s);

    int Length() const
    {
        return len;
    }
    int BufferSize() const
    {
        return size;
    }

    int SetLength(int newlen);
    int Dimension(int newlen)
    {
        return SetLength(newlen);
    }

    String & Add(const String & s)
    {
        return *this += s;
    }
    String & Add(char ch)
    {
        return *this += ch;
    }

    String   RightToLeft();
    String & Invert();
    String & Invert(const String & s);

    String & Trim();
    String & Trim(char character);
    vector<String> *Split(char splitChar);

    long   AsInteger() const;
    bool   AsInteger(long& intValue) const;
    bool   AsInteger(int& intValue) const;
    double AsDouble() const
    {
        return (double) *this;
    }
    long double AsLongDouble() const
    {
        return (long double) *this;
    }

    int    printf(const char * format, ...);
    int    vprintf(const char * format, va_list arglist);

    int    catprintf(const char * format, ...);
    int    vcatprintf(const char * format, va_list arglist);

    // Replacement vsnprintf and snprint functions for
    // problematic architectures...

    static int  my_snprintf(char * buffer, int bufsize, const char * format, ...);
    static int  my_vsnprintf(char * buffer, int bufsize, const char * format, va_list args);
    static void my_vsnprintf_close_file();
    static void check_vsnprintf();

    // Check string contents
    bool   IsNumber();

    // Explicit conversions
    const unsigned char * uchar() const
    {
        return (unsigned char *) buffer;
    }
    const signed char * schar() const
    {
        return (signed char *) buffer;
    }

    static FILE * my_vsnprintf_file;

    // Utility functions
    void Fill(char ch, int length = -1);

private:

    static int vsnprintfChecked;
};

inline int Compare(const String & s1, const String & s2)
{
    return s1.Compare(s2);
}

inline int Compare(const String & s1, const char * s2)
{
    return s1.Compare(s2);
}

inline int Compare(const char * s1, const String & s2)
{
    return -s2.Compare(s1);
}

inline int FastCompare(const String & s1, const String & s2)
{
    return s1.FastCompare(s2);
}

inline int FastCompare(const String & s1, const char * s2)
{
    return s1.FastCompare(s2);
}

inline int FastCompare(const char * s1, const String & s2)
{
    return -s2.FastCompare(s1);
}

inline int SlowCompare(const String & s1, const String & s2)
{
    return s1.SlowCompare(s2);
}

inline int SlowCompare(const String & s1, const char * s2)
{
    return s1.SlowCompare(s2);
}

inline int SlowCompare(const char * s1, const String & s2)
{
    return -s2.SlowCompare(s1);
}

String operator + (char lhs, const String & rhs);
String operator + (const char * lhs, const String & rhs);
String operator + (int lhs, const String & rhs);
String operator + (unsigned int lhs, const String & rhs);

std::ostream& operator << (std::ostream& os, const String& s);

/// Write to a file using streaming.
/// \param stream file to write to - IFILE is a pointer to an InputFile object
/// \param str string containing what should be written to the file.
inline InputFile& operator << (InputFile& stream, const String& str)
{
    unsigned int numExpected = str.Length();
    unsigned int numWritten = 
        stream.ifwrite(str.c_str(), numExpected);
    if(numExpected != numWritten)
    {
        std::cerr << "Failed to stream to IFILE, expected " 
                  << numExpected << " but only wrote "
                  << numWritten << std::endl;
    }
    return(stream);
}


#endif



