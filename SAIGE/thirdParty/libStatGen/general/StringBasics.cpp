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

#include "StringBasics.h"
#include "Error.h"
#include "Constant.h"
#include "MathConstant.h"

#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <sstream>

#define SWP(A,B) {int tmp=a; a=b; b=tmp;}

#ifdef  _MSC_VER
#ifndef snprintf
#define vsnprintf _vsnprintf
#define snprintf  _snprintf
#endif
#endif

// If natural ordering is defined, comparisons will
// order strings including numbers correctly
// (eg, ... "8", "9", "10" ...) rather than using
// ASCII ordering (... "10", "8", "9", ...)
#define NATURAL_ORDERING         1

int  String::alloc = 8;
bool String::caseSensitive = true;

void String::NewString(int startsize)
{
    len = 0;
    size = (startsize + alloc) / alloc * alloc;
    buffer = new char [size];
    buffer[0] = 0;
}

String::String(const char * s)
{
    int clen = s == NULL ? 0 : strlen(s);
    NewString(clen);
    if (clen)
    {
        len  = clen;
        memcpy(buffer, s, len + 1);
    }
}

String::String(char ch, int count)
{
    NewString(count);
    memset(buffer, ch, count);
    buffer[count] = 0;
    len = count;
}

String::String(const String & s)
{
    len = s.len;
    size = (s.len + alloc) / alloc * alloc;;
    buffer = new char [size];
    memcpy(buffer, s.buffer, len + 1);
}

void String::Grow(int newSize)
{
    if (newSize >= size)
    {
        if ((newSize >> 1) >= size)
            size = (newSize + alloc) / alloc * alloc;
        else
        {
            size = alloc;
            while (size <= newSize)
                size *= 2;
        }

        char * tmp = new char [size];
        // // len + 1 due to terminating NUL which is not counted in len
        // memcpy(tmp, buffer, len + 1);
        memcpy(tmp, buffer, len);
        tmp[len] = '\0';
        delete [] buffer;
        buffer = tmp;
    }
}

void String::Swap(String & s)
{
    char * temp = s.buffer;
    s.buffer = buffer;
    buffer = temp;

    int swap = s.size;
    s.size = size;
    size = swap;

    swap = s.len;
    s.len = len;
    len = swap;
}

String & String::Copy(const String & s)
{
    Grow(s.len);
    len = s.len;
    memcpy(buffer, s.buffer, len + 1);
    return *this;
}

String & String::Copy(const String & s, int start, int n)
{
    if (s.len <= start) return Clear();
    if (s.len < start + n) n = s.len - start;
    Grow(n);
    memcpy(buffer, s.buffer + start, n);
    buffer[len = n] = 0;
    return *this;
}

String & String::Copy(const char * s)
{
    if (s == NULL)
    {
        len = 0;
        buffer[0] = 0;
    }
    else
    {
        int clen = strlen(s);
        Grow(clen);
        len = clen;
        memcpy(buffer, s, len + 1);
    }
    return *this;
}

String & String::ToUpper()
{
    for (int i = 0; i < len; i++)
        buffer[i] = (char) toupper(buffer[i]);
    return *this;
}

String & String::ToLower()
{
    for (int i = 0; i < len; i++)
        buffer[i] = (char) tolower(buffer[i]);
    return *this;
}

String String::AsUpper()
{
    String temp;
    temp = *this;
    return temp.ToUpper();
}

String String::AsLower()
{
    String temp;
    temp = *this;
    return temp.ToLower();
}

String String::Capitalize()
{
    String temp;
    temp = *this;
    temp.buffer[0] = (char) toupper(temp.buffer[0]);
    return temp;
}

String & String::operator = (const String & rhs)
{
    Copy(rhs);
    return *this;
}

String & String::operator = (const char * rhs)
{
    Copy(rhs);
    return * this;
}

String & String::operator += (const String & rhs)
{
    Grow(len + rhs.len);
    memcpy(buffer + len, rhs.buffer, rhs.len + 1);
    len += rhs.len;
    return *this;
}

String & String::operator += (const char * rhs)
{
    if (rhs != NULL)
    {
        int clen = strlen(rhs);
        Grow(len + clen);
        memcpy(buffer + len, rhs, clen + 1);
        len += clen;
    }
    return *this;
}

String String::operator + (const String & rhs) const
{
    String result(len + rhs.len);
    memcpy(result.buffer, buffer, len);
    memcpy(result.buffer + len, rhs.buffer, rhs.len + 1);
    result.len = len + rhs.len;
    return result;
}

String String::operator + (const char * rhs) const
{
    if (rhs != NULL)
    {
        int clen = strlen(rhs);
        String result(len + clen);
        memcpy(result.buffer, buffer, len);
        memcpy(result.buffer + len, rhs, clen + 1);
        result.len = len + clen;
        return result;
    }
    return *this;
}

String & String::operator = (char ch)
{
    if (ch)
    {
        Grow(1);
        buffer[0] = ch;
        buffer[1] = 0;
        len = 1;
    }
    else
        len  = buffer[0] = 0;
    return *this;
}

String & String::operator += (char ch)
{
    if (ch)
    {
        Grow(len + 1);
        buffer[len] = ch;
        buffer[++len] = 0;
    }
    return *this;
}

String String::operator + (char ch) const
{
    String result(*this);
    result += ch;
    return result;
}

String & String::operator = (int rhs)
{
    Clear();

    if (rhs < 0)
    {
        Add('-');
        *this += (unsigned int) -rhs;
    }
    else
        *this = (unsigned int) rhs;
    return *this;
}

String & String::operator = (unsigned int rhs)
{
    Clear();

    unsigned long long base = 10;
    int   digits = 1;

    while (rhs >= base)
    {
        base *= 10;
        digits++;
    }

    Grow(digits);

    while (base /= 10)
    {
        char ch = char(rhs / base);
        rhs = rhs - ch * base;
        buffer[len++] = char(ch + '0');
    }
    buffer[len] = 0;
    return *this;
};

String String::operator + (int rhs) const
{
    String result(*this);
    result += rhs;
    return result;
};

String String::operator + (unsigned int rhs) const
{
    String result(*this);
    result += rhs;
    return result;
};

String & String::operator += (int rhs)
{
    String temp;
    temp = rhs;
    return *this += temp;
}

String & String::operator += (unsigned int rhs)
{
    String temp;
    temp = rhs;
    return *this += temp;
}

String & String::operator *= (unsigned int rhs)
{
    if (rhs == 0)
        Clear();
    else
    {
        String original(*this);

        Grow(len * rhs);

        for (unsigned int i = 1; i < rhs; i++)
            *this += original;
    }
    return *this;
}

String & String::operator = (double rhs)
{
    LockBuffer(32);
    sprintf(buffer, "%.3f", rhs);
    UnlockBuffer();
    return *this;
}

String String::operator + (double rhs) const
{
    String result(*this);
    result += rhs;
    return result;
}

String & String::operator += (double rhs)
{
    String temp;
    temp = rhs;
    return *this += temp;
}


void String::appendFullFloat(float rhs)
{
    std::ostringstream os;
    os << rhs;
    *this += os.str().c_str();
}

char * String::LockBuffer(int min)
{
    if (min > 0) Grow(min);
    return buffer;
}

String & String::UnlockBuffer()
{
    for (len = 0; len < size; len++)
        if (buffer[len] == 0)
            return *this;
    error("BasicString - direct access overflowed buffer");
    return *this;
}

int String::Compare(const String & s) const
{
    if (caseSensitive)
        return String::FastCompare(s);
    else
        return String::SlowCompare(s);
}

int String::Compare(const char * s) const
{
    return caseSensitive ? FastCompare(s) : SlowCompare(s);
}

int String::FastCompare(const String & s) const
{
    for (int i = 0; i <= len; i++)
        if (buffer[i] - s.buffer[i])
        {
#ifdef NATURAL_ORDERING
            int d = i;
            while (isdigit(buffer[d]) && isdigit(s.buffer[d]))
                d++;
            if (isdigit(buffer[d]))
                return 1;
            if (isdigit(s.buffer[d]))
                return -1;
#endif
            return buffer[i] - s.buffer[i];
        }
    return 0;
}

int String::FastCompare(const char * s) const
{
    if (s == NULL)
        return -len;

    for (int i = 0; i <= len; i++)
        if (buffer[i] - s[i])
        {
#ifdef NATURAL_ORDERING
            int d = i;
            while (isdigit(buffer[d]) && isdigit(s[d]))
                d++;
            if (isdigit(buffer[d]))
                return 1;
            if (isdigit(s[d]))
                return -1;
#endif
            return buffer[i] - s[i];
        }
    return 0;
}

int String::SlowCompare(const String & s) const
{
    for (int i = 0; i <= len; i++)
        if (toupper(buffer[i]) - toupper(s.buffer[i]))
        {
#ifdef NATURAL_ORDERING
            int d = i;
            while (isdigit(buffer[d]) && isdigit(s[d]))
                d++;
            if (isdigit(buffer[d]))
                return 1;
            if (isdigit(s.buffer[d]))
                return -1;
#endif
            return toupper(buffer[i]) - toupper(s.buffer[i]);
        }
    return 0;
}

int String::SlowCompare(const char * s) const
{
    if (s == NULL)
        return -len;

    for (int i = 0; i <= len; i++)
        if (toupper(buffer[i]) - toupper(s[i]))
        {
#ifdef NATURAL_ORDERING
            int d = i;
            while (isdigit(buffer[d]) && isdigit(s[d]))
                d++;
            if (isdigit(buffer[d]))
                return 1;
            if (isdigit(s[d]))
                return -1;
#endif
            return toupper(buffer[i]) - toupper(s[i]);
        }
    return 0;
}

int String::ReadLine(FILE * f)
{
    len = 0;
    buffer[len] = 0;

    if (f == NULL) return -1;

    int  clen = 0;
    char check[2] = {0, 0};

    int     step = 128;
    String  format("%128[^\n\r]%1[\n\r]");

    int returnValue = 1;

    int io = 0;

    while (check[0] != '\n' && check[0] != '\r')
    {
        if (clen)
        {
            step *= 2;
            format.printf("%%%d%s", step, "[^\n\r]%1[\n\r]");
        }
        clen += step;

        io = fscanf(f, format, LockBuffer(clen) + len, check);
        UnlockBuffer();
        // Avoid getting stuck on zero length lines (system specific!)
        if (io == 0 && check[0] != '\n' && check[0] != '\r')
            io = fscanf(f, "%1[\n\r]", check);
        if (io == 0 || io == EOF)
        {
            // Set return value to indicate error/EOF
            returnValue = -1;
            break;
        }
    }

    if (check[0] == '\n') io = fscanf(f, "%*1[\r]");
    if (check[0] == '\r') io = fscanf(f, "%*1[\n]");

    return returnValue;
}


String & String::Read(FILE * f)
{
    len = 0;
    buffer[len] = 0;

    if (f == NULL) return *this;

    int  clen = 0;
    char check[2] = {'G', 0};

    while (strchr(WHITESPACE, check[0]) == NULL)
    {
        clen += READBUF;
        int io = fscanf(f, " %" READBUFSTR "[^" WHITESPACE "]"
                        "%1[" WHITESPACE "]", LockBuffer(clen) + len, check);
        if (io == 0 || io == EOF) break;
        UnlockBuffer();
    }

    return *this;
}

String & String::Read()
{
    return Read(stdin);
}

String & String::Read(IFILE & f)
{
    len = 0;
    buffer[len] = 0;

    if (f == NULL) return *this;

    bool leading = true;

    while (true)
    {
        int ch = ifgetc(f);

        if (ch == -1) break;

        if (strchr(WHITESPACE, ch) != NULL)
        {
            if (leading)
            {
                continue;
            }
            else
            {
                break;
            }
        }

        if (len + 1 == size)
            Grow(len + 1);

        buffer[len++] = (char) ch;
        buffer[len] = 0;

        leading = false;
    }

    return *this;
}

int String::ReadLine()
{
    static int last = 0;
    int ch;

    len = 0;
    buffer[len] = 0;

    while (true)
    {
        ch = getchar();

        if (ch == EOF)
        {
            break;
        }

        if (ch == 10)
        {
            if (last == 13)
            {
                last = 0;
                continue;
            }
            else
            {
                last = 10;
                break;
            }
        }

        if (ch == 13)
        {
            if (last == 10)
            {
                last = 0;
                continue;
            }
            else
            {
                last = 13;
                break;
            }
        }

        if (len + 1 == size)
        {
            Grow(len + 1);
        }

        last = ch;
        buffer[len++] = (char) last;
        buffer[len] = 0;
    }

    if ((ch == EOF) && (len == 0))
    {
        // Indicate error/EOF if nothing was read.
        return -1;
    }

    // Return success.
    return 1;
}


// Read line using getc.

#if defined(_WIN32)
int String::ReadLine(IFILE & f)
{
    static int last = 0;
    int ch;

    len = 0;
    buffer[len] = 0;

    while (true)
    {
        ch = f->ifgetc();

        if (ch == EOF)
        {
            break;
        }

        if (ch == 10)
        {
            if (last == 13)
            {
                last = 0;
                continue;
            }
            else
            {
                last = 10;
                break;
            }
        }

        if (ch == 13)
        {
            if (last == 10)
            {
                last = 0;
                continue;
            }
            else
            {
                last = 13;
                break;
            }
        }

        if (len + 1 == size)
        {
            Grow(len + 1);
        }

        last = ch;
        buffer[len++] = (char) last;
        buffer[len] = 0;
    }

    if ((ch == EOF) && (len == 0))
    {
        // Indicate error/EOF if nothing was read.
        return -1;
    }
    return 1;
}
#else
int String::ReadLine(IFILE & f)
{
    int ch;
    char *ptr = buffer;
    char *endBuffer = buffer + size;
    len = 0;

    while ( ((ch = f->ifgetc()) != EOF) && (ch != '\n'))
    {
      if (ptr >= endBuffer - 1)
        {
            // resize: 1 byte for the next character, 1 byte
            //  for the NUL at the end.
            Grow(len + 2);
            endBuffer = buffer + size;
            ptr = buffer + len;
        }

        *ptr++ = ch;
        len++;
    }

    // NB: assumes that buffer is always allocated.
    buffer[len] = 0;

    if ((ch == EOF) && (len == 0))
    {
        // Indicate error/EOF if nothing was read.
        return -1;
    }
    return 1;
}
#endif

void String::Write(FILE * f)
{
    fprintf(f, "%s", buffer);
}

void String::Write()
{
    Write(stdout);
}

void String::WriteLine()
{
    WriteLine(stdout);
}

void String::WriteLine(FILE * f)
{
    if (f == NULL) return;
    fprintf(f, "%s\n", buffer);
}

std::ostream& operator << (std::ostream& os, const String& s)
{
    return os << s.c_str();
}

String String::Left(int n) const
{
    if (n < 0) n = 0;
    if (len < n) n = len;
    String result(n);
    memcpy(result.buffer, buffer, n);
    result.buffer[result.len = n] = 0;
    return result;
}

String String::Right(int n) const
{
    if (n < 0) n = 0;
    if (len < n) n = len;
    String result(n);
    memcpy(result.buffer, buffer + len - n, n);
    result.buffer[result.len = n] = 0;
    return result;
}

String String::SubStr(int start, int n) const
{
    if (start < 0)
    {
        n += start;
        start = 0;
    };
    n = min(len - start, n);
    n = max(n, 0);
    String result(n);
    if (start > len) return result;
    memcpy(result.buffer, buffer + start, n);
    result.buffer[result.len = n] = 0;
    return result;
}

String String::SubStr(int start) const
{
    return SubStr(start, len - start);
}

String String::Mid(int start, int end) const
{
    return SubStr(start, end - start + 1);
}

int String::FindChar(char ch, int start) const
{
    return caseSensitive ? FastFindChar(ch, start) : SlowFindChar(ch, start);
}

int String::FastFindChar(char ch, int start) const
{
    for (; start < len; start++)
        if (buffer[start] == ch)
            return start;
    return -1;
}

int String::SlowFindChar(char ch, int start) const
{
    ch = (char) toupper(ch);
    for (; start < len; start++)
        if (toupper(buffer[start]) == ch)
            return start;
    return -1;
}

int String::FindLastChar(char ch) const
{
    return caseSensitive ? FastFindLastChar(ch) : SlowFindLastChar(ch);
}

int String::FastFindLastChar(char ch) const
{
    for (int start = len-1; start >= 0; start--)
        if (buffer[start] == ch)
            return start;
    return -1;
}

int String::SlowFindLastChar(char ch) const
{
    ch = (char) toupper(ch);
    for (int start = len-1 ; start >= 0; start--)
        if (toupper(buffer[start]) == ch)
            return start;
    return -1;
}

int String::Find(const String & pattern, int start) const
{
    return caseSensitive ? FastFind(pattern, start) : SlowFind(pattern, start);
}

// TODO -- We should have a better string search algorithm

int String::FastFind(const String & pattern, int start) const
{
    for (int i ; start <= len - pattern.Length(); start++)
        if (buffer[start] == pattern[0])
        {
            for (i = 1; i < pattern.Length(); i++)
                if (pattern[i] != buffer[start + i])
                    break;
            if (i == pattern.Length()) return start;
        }
    return -1;
}

int String::SlowFind(const String & pattern, int start) const
{
    int firstchar = toupper(pattern[0]);

    for (int i ; start <= len - pattern.Length(); start++)
        if (toupper(buffer[start]) == firstchar)
        {
            for (i = 1; i < pattern.Length(); i++)
                if (toupper(pattern[i]) != toupper(buffer[start + i]))
                    break;
            if (i == pattern.Length()) return start;
        }
    return -1;
}

int String::SetLength(int newlen)
{
    if (newlen > len)
    {
        Grow(newlen);
        memset(buffer + len, ' ', newlen - len);
    }
    buffer[newlen] = 0;
    return len = newlen;
}

String & String::Filter(const String & s)
{
    int to = 0;
    for (int from = 0; from < len; from++)
        if (s.FindChar(buffer[from]) != -1)
            buffer[to++] = buffer[from];
    buffer[len = to] = 0;
    return *this;
}

String & String::Filter(const char * s)
{
    String filter(s);
    return Filter(filter);
}

String & String::ExcludeCharacters(const String & s)
{
    int to = 0;
    for (int from = 0; from < len; from++)
        if (s.FindChar(buffer[from]) == -1)
            buffer[to++] = buffer[from];
    buffer[len = to] = 0;
    return *this;
}

String & String::ExcludeCharacters(const char * s)
{
    String excluded(s);
    return ExcludeCharacters(excluded);
}

String operator + (const char * lhs, const String & rhs)
{
    String result(lhs);
    result += rhs;
    return result;
}

String operator + (char lhs, const String & rhs)
{
    String result(lhs);
    result += rhs;
    return result;
}

String operator + (int lhs, const String & rhs)
{
    String result;
    result = lhs;
    result += rhs;
    return result;
}

String operator + (unsigned int lhs, const String & rhs)
{
    String result;
    result = lhs;
    result += rhs;
    return result;
}

long String::AsInteger() const
{
    long returnValue = 0;
    if(!AsInteger(returnValue))
    {
        // This is not an integer, but nothing to do but return a value.
    }
    return(returnValue);
}


// Check that the string is an integer when converting it.
// If the entire string is an integer, return true, if not, return false.
bool String::AsInteger(long& intValue) const
{
    long integer = 0;
    int  base = 10;
    int  pos = 0;
    int  sign = 1;
    bool isInt = true;

    // If this is no value for this integer, return false.
    if (pos == len)
    {
        return(false);
    }

    if (buffer[pos] == '-')
    {
        sign = -1, pos++;
    }

    if ((len > pos + 2) && (buffer[pos] == '0') &&
            ((buffer[pos+1] == 'x') || (buffer[pos+1] == 'X')))
    {
        base = 16, pos += 2;
    }

    // If this is no value for this integer, return false.
    if (pos == len)
    {
        return(false);
    }

    for (;  pos < len; pos++)
    {
        char digit = (char) toupper(buffer[pos]);

        if (digit >= '0' && digit <= '9')
        {
            integer = integer * base + digit - '0';
        }
        else if (digit >= 'A' && digit <= 'F' && base == 16)
        {
            integer = integer * base + digit - 'A' + 10;
        }
        else
        {
            isInt = false;
            break;
        }
    }

    intValue = sign*integer;

    return(isInt);
}


// Check that the string is an integer when converting it.
// If the entire string is an integer, return true, if not, return false.
bool String::AsInteger(int& intValue) const
{
    int integer = 0;
    int  base = 10;
    int  pos = 0;
    int  sign = 1;
    bool isInt = true;

    // If this is no value for this integer, return false.
    if (pos == len)
    {
        return(false);
    }

    if (buffer[pos] == '-')
    {
        sign = -1, pos++;
    }

    if ((len > pos + 2) && (buffer[pos] == '0') &&
            ((buffer[pos+1] == 'x') || (buffer[pos+1] == 'X')))
    {
        base = 16, pos += 2;
    }

    // If this is no value for this integer, return false.
    if (pos == len)
    {
        return(false);
    }

    for (;  pos < len; pos++)
    {
        char digit = (char) toupper(buffer[pos]);

        if (digit >= '0' && digit <= '9')
        {
            integer = integer * base + digit - '0';
        }
        else if (digit >= 'A' && digit <= 'F' && base == 16)
        {
            integer = integer * base + digit - 'A' + 10;
        }
        else
        {
            isInt = false;
            break;
        }
    }

    intValue = sign*integer;

    return(isInt);
}


String & String::Invert()
{
    for (int i = 0, j = len - 1; i < j; i++, j--)
    {
        char tmp = buffer[i];
        buffer[i] = buffer[j];
        buffer[j] = tmp;
    }
    return *this;
}

String String::RightToLeft()
{
    String result(*this);
    result.Invert();
    return result;
}

String & String::Invert(const String & s)
{
    Copy(s);
    return Invert();
}

int String::CompareToStem(const String & stem) const
{
    if (caseSensitive)
        return String::FastCompareToStem(stem);
    else
        return String::SlowCompareToStem(stem);
}

int String::FastCompareToStem(const String & stem) const
{
    for (int i = 0; i < stem.len; i++)
        if (buffer[i] - stem.buffer[i])
        {
#ifdef NATURAL_ORDERING
            int d = i;
            while (isdigit(buffer[d]) && isdigit(stem.buffer[d]) && d < stem.len)
                d++;
            if (isdigit(buffer[d]) && d < stem.len)
                return 1;
            if (isdigit(stem.buffer[d]))
                return -1;
#endif
            return buffer[i] - stem.buffer[i];
        }
    return 0;
}

int String::SlowCompareToStem(const String & stem) const
{
    for (int i = 0; i < stem.len; i++)
        if (toupper(buffer[i]) - toupper(stem.buffer[i]))
        {
#ifdef NATURAL_ORDERING
            int d = i;
            while (isdigit(buffer[d]) && isdigit(stem.buffer[d]) && d < stem.len)
                d++;
            if (isdigit(buffer[d]) && d < stem.len)
                return 1;
            if (isdigit(stem.buffer[d]))
                return -1;
#endif
            return toupper(buffer[i]) - toupper(stem.buffer[i]);
        }
    return 0;
}

int String::CompareToStem(const char * stem) const
{
    if (caseSensitive)
        return String::FastCompareToStem(stem);
    else
        return String::SlowCompareToStem(stem);
}

int String::FastCompareToStem(const char * stem) const
{
    for (int i = 0; stem[i] != 0; i++)
        if (buffer[i] - stem[i])
            return buffer[i] - stem[i];
    return 0;
}

int String::SlowCompareToStem(const char * stem) const
{
    for (int i = 0; stem[i] != 0; i++)
        if (toupper(buffer[i]) - toupper(stem[i]))
            return toupper(buffer[i]) - toupper(stem[i]);
    return 0;
}

int String::MatchesBeginningOf(const String & stem) const
{
    if (caseSensitive)
        return String::FastMatchesBeginningOf(stem);
    else
        return String::SlowMatchesBeginningOf(stem);
}

int String::FastMatchesBeginningOf(const String & stem) const
{
    for (int i = 0; i < len; i++)
        if (buffer[i] - stem.buffer[i])
            return buffer[i] - stem.buffer[i];
    return 0;
}

int String::SlowMatchesBeginningOf(const String & stem) const
{
    for (int i = 0; i < len; i++)
        if (toupper(buffer[i]) - toupper(stem.buffer[i]))
            return toupper(buffer[i]) - toupper(stem.buffer[i]);
    return 0;
}

int String::MatchesBeginningOf(const char * stem) const
{
    if (caseSensitive)
        return String::FastMatchesBeginningOf(stem);
    else
        return String::SlowMatchesBeginningOf(stem);
}

int String::FastMatchesBeginningOf(const char * stem) const
{
    for (int i = 0; i < len; i++)
        if (buffer[i] - stem[i])
            return buffer[i] - stem[i];
    return 0;
}

int String::SlowMatchesBeginningOf(const char * stem) const
{
    for (int i = 0; i < len; i++)
        if (toupper(buffer[i]) - toupper(stem[i]))
            return toupper(buffer[i]) - toupper(stem[i]);
    return 0;
}

String & String::Trim(char character)
{
    int first = 0;
    while (buffer[first] && buffer[first] == character)
        first++;

    int last = len - 1;
    while (last >= 0 && buffer[last] == character)
        last--;

    int out = 0;
    while (first <= last)
        buffer[out++] = buffer[first++];

    buffer[len = out] = 0;

    return *this;
}

String & String::Trim()
{
    int first = 0;
    while (buffer[first] && isspace(buffer[first]))
        first++;

    int last = len - 1;
    while (last >= 0 && isspace(buffer[last]))
        last--;

    int out = 0;
    while (first <= last)
        buffer[out++] = buffer[first++];

    buffer[len = out] = 0;

    return *this;
}

vector<String> *String::Split(char splitChar)
{
    vector<String> *result = new vector<String>;
    String word;

    for (int i = 0; i<Length(); i++)
    {
        if ((*this)[i]==splitChar)
        {
            result->push_back(word);
            word.Clear();
        }
        else
            word.Add((*this)[i]);
    }
    if (word.Length()>0) result->push_back(word);
    return result;
}


#define VSNPRINTF_NOT_CHECKED    0
#define VSNPRINTF_IS_OK          1
#define VSNPRINTF_NOT_OK         2

int String::vsnprintfChecked = 0;

int String::printf(const char * format, ...)
{
    va_list  ap;
    va_start(ap, format);

    vprintf(format, ap);

    va_end(ap);
    return len;
}

int String::catprintf(const char * format, ...)
{
    va_list  ap;
    va_start(ap, format);

    vcatprintf(format, ap);

    va_end(ap);
    return len;
}

int String::vprintf(const char * format, va_list ap)
{
    check_vsnprintf();

    while (true)
    {
        int bytes_needed;
#ifdef va_copy
        va_list arguments;
        va_copy(arguments, ap);
#else
        va_list & arguments = ap;
#endif

        if (vsnprintfChecked == VSNPRINTF_IS_OK)
            bytes_needed = vsnprintf(buffer, size, format, arguments);
        else
            bytes_needed = my_vsnprintf(buffer, size, format, arguments);

#ifdef va_copy
        va_end(arguments);
#endif

        if (bytes_needed >= size)
            Grow(bytes_needed);
        else if (bytes_needed == -1)
            Grow(size * 2);
        else
        {
            return len = bytes_needed;
        }
    }
}

void String::check_vsnprintf()
{
    if (vsnprintfChecked == VSNPRINTF_NOT_CHECKED)
    {
        char temp[100];

        memset(temp, 0, 100);
        int check = snprintf(temp, 5, "%5s", "VSNPRINTF");

        if (temp[6] != 0 || temp[7] != 0 || (check != 9 && check != -1))
            /*
            error("This program requires a working version of vsnprintf\n"
                  "However, vsnprintf in the current library seems buggy\n\n"
                  "Recompiling this program with the -D__REPLACE_SNPRINTF__ flag\n"
                  "may solve this problem.\n\n");
            */
            vsnprintfChecked = VSNPRINTF_NOT_OK;
        else
            vsnprintfChecked = VSNPRINTF_IS_OK;
    }
}

int String::vcatprintf(const char * format, va_list ap)
{
    check_vsnprintf();

    if (len == size)
        Grow(size * 2);

    while (true)
    {
        int bytes_needed;
#ifdef va_copy
        va_list arguments;
        va_copy(arguments, ap);
#else
        va_list & arguments = ap;
#endif

        if (vsnprintfChecked == VSNPRINTF_IS_OK)
            bytes_needed = len + vsnprintf(buffer + len, size - len, format, arguments);
        else
            bytes_needed = len + my_vsnprintf(buffer + len, size - len, format, arguments);

#ifdef va_copy
        va_end(arguments);
#endif

        if (bytes_needed >= size)
            Grow(bytes_needed);
        else if (bytes_needed < len)
            Grow(size * 2);
        else
        {
            return len = bytes_needed;
        }
    }
}

FILE * String::my_vsnprintf_file = NULL;

int String::my_vsnprintf(char * buffer, int bufsize, const char * format, va_list args)
{
    if (my_vsnprintf_file == NULL)
    {
        my_vsnprintf_file = tmpfile();
        atexit(my_vsnprintf_close_file);
    }

    rewind(my_vsnprintf_file);

    int len = vfprintf(my_vsnprintf_file, format, args);

    rewind(my_vsnprintf_file);

    if (len < bufsize)
        buffer[bufsize = len] = 0;
    int numRead = fread(buffer, 1, bufsize, my_vsnprintf_file);
    if(numRead != bufsize)
    {
        std::cerr
            << "Warning, StringBasics failed reading stream in my_vsnprintf\n";
    }
    return len;
}

int String::my_snprintf(char * buffer, int bufsize, const char * format, ...)
{
    va_list  ap;
    va_start(ap, format);

    int bytes = my_vsnprintf(buffer, bufsize, format, ap);

    va_end(ap);

    return bytes;
}

void String::my_vsnprintf_close_file()
{
    fclose(my_vsnprintf_file);
}

bool String::IsNumber()
{
    int  pos = 0;
    bool digits = false;

    // Skip leading sign
    if (buffer[pos] == '-' || buffer[pos] == '+')
        pos++;

    // Check integer portion
    while (buffer[pos] >= '0' && buffer[pos] <= '9')
        pos++, digits = true;

    // Skip decimal point
    if (buffer[pos] == '.')
    {
        pos++;

        // Check fractional portion
        while (buffer[pos] >= '0' && buffer[pos] <= '9')
            pos++, digits = true;
    }

    if (!digits) return false;

    // Check exponent
    if (buffer[pos] == 'E' || buffer[pos] == 'e')
    {
        pos++;

        // Skip leading sign
        if (buffer[pos] == '-' || buffer[pos] == '+')
            pos++;

        digits = false;

        // Check exponent digits
        while (buffer[pos] >= '0' && buffer[pos] <= '9')
            pos++, digits = true;
    }

    return (pos == len) && digits;
}

void String::Fill(char ch, int length)
{
    if (length >= 0)
        SetLength(length);

    for (int i = 0; i < len; i++)
        buffer[i] = ch;
}

String & String::Reverse()
{
    for (int i = 0, j = len - 1; i < j; i++, j--)
    {
        int tmp = buffer[i];
        buffer[i] = buffer[j];
        buffer[j] = tmp;
    }

    return *this;
}

// String::LeftClip() trims the string so only characters after clipPoint remain

String & String::LeftClip(int clipAmount)
{
    if (clipAmount == 0)
        return *this;

    if (clipAmount > Length())
    {
        len = 0;
        return *this;
    }

    // Use memory move, because the two blocks can overlap
    memmove(buffer, buffer + clipAmount, len - clipAmount);
    buffer[len -= clipAmount] = 0;

    return *this;
}

String & String::RightClip(int clipAmount)
{
    if (clipAmount == 0) return *this;

    if (clipAmount > Length())
    {
        len = 0;
        return *this;
    }

    len -= clipAmount;
    buffer[len] = 0;

    return *this;
}

// Implementation of long double convertors is in flux across different platforms

#ifdef __GNUC__
String::operator long double() const
{
    return strtold(buffer, NULL);
}
#else
#ifdef __BORLANDC__
String::operator long double() const
{
    return _strtold(buffer, NULL);
}
#else
String::operator long double() const
{
    return atof(buffer);
}
#endif
#endif


