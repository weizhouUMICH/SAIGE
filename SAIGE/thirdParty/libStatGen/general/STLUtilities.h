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

#ifndef _STLUTILITIES_H
#define _STLUTILITIES_H
#include <assert.h>
#include <iomanip>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string.h>
#include <vector>


///
/// This file is inspired by the poor quality of string support in
/// STL for what should be trivial capabiltiies, for example setting
/// or appending the ASCII representation of a floating point
/// or integer number to a string.
///
/// This file uses variadic templates to implement a type safe
/// version (subset) of C-library printf.
///
/// Therefore, -std=c++0x is a required option on g++
///

namespace STLUtilities
{
///
/// use std streams API to do float conversion to string,
/// then append it.
///
inline std::string &append(std::string &s, float f)
{
    std::ostringstream buffer;
    buffer << f;
    s += buffer.str();
    return s;
}

///
/// use std streams API to do double conversion to string,
/// then append it.
///
inline std::string &append(std::string &s, double f)
{
    std::ostringstream buffer;
    buffer << f;
    s += buffer.str();
    return s;
}

///
/// The rest we can handle readily ourselves.
/// Unlike std::string operator +, this operator treats c as
/// a character and appends the ASCII character c.
///
inline std::string &append(std::string &s, char c)
{
    s += c;
    return s;
}

///
/// Similar to signed char case, but this time for unsigned.
///
inline std::string &append(std::string &s, unsigned char c)
{
    s += c;
    return s;
}

///
/// Now append a full C style NUL terminated string to
/// the std::string.
///
inline std::string &append(std::string &s, const char *rhs)
{
    s += rhs;
    return s;
}

///
/// Prevent the generic template from picking up std::string
///
inline std::string &append(std::string &s, std::string &rhs)
{
    s += rhs;
    return s;
}

///
/// iterate over the provided vector, appending all elements with
/// an optional separator
///
template<typename T> std::string &append(std::string &s, std::vector<T> v, std::string separator="")
{
    for (typename T::iterator i=v.begin(); i!=v.end(); i++)
    {
        if (i!=v.begin()) s += separator;
        s << *i;
    }
    return s;
}

///
/// This template handles the rest of the cases for
/// integral types.  Not user friendly if you pass
/// in a type T that is for example a std::vector.
///
template<typename T> std::string &append(std::string &s, T i)
{
    char digits[20];
    char *p = digits;
    bool negative = false;

    if (i<0)
    {
        negative = true;
        i = -i;
    }

    do
    {
        *p++ = '0' + i % 10;
        i = i/10;
    }
    while (i);

    if (negative) s += '-';

    do
    {
        s += *--p;
    }
    while (p > digits);

    return s;
}

inline std::string &operator <<(std::string &s, char c)
{
    return append(s, c);
}

inline std::string &operator <<(std::string &s, unsigned char c)
{
    return append(s, c);
}

inline std::string &operator <<(std::string &s, uint64_t i)
{
    return append(s, i);
}

inline std::string &operator <<(std::string &s, int64_t i)
{
    return append(s, i);
}

template<typename T> std::string &operator <<(std::string &s, T val)
{
    return append(s, val);
}

template<typename S> std::string &append(std::string &s, std::vector<std::string> v, S delimeter, bool itemize = false)
{
    bool showDelimeter = false;
    for (std::vector<std::string>::iterator i=v.begin(); i!=v.end(); i++)
    {
        if (showDelimeter) s << delimeter;
        else showDelimeter = true;
        if (itemize) s << (i - v.begin()) << ": ";
        s << *i;
    }
    return s;
}

template<typename T, typename S> std::string &append(std::string &s, std::vector<T> v, S delimeter, bool itemize = false)
{
    bool showDelimeter = false;
    for (typename std::vector<T>::iterator i=v.begin(); i!=v.end(); i++)
    {
        if (showDelimeter) s << delimeter;
        else showDelimeter = true;
        if (itemize) s << (i - v.begin()) << ": ";
        s << *i;
    }
    return s;
}

//
// Split the string input into words delimited by the character
// delimiter.  For a given number of input delimiters, result.size()
// will not change, regardless of the data in between the delimiters.
//
// Refactor this to pre-allocate the word that we place data into,
// then we have minimal data copy.
//
int Tokenize(std::vector<std::string> &result, const char *input, char delimiter);




//
// Variadic templates necessary for reasonable printf implementation
// are only supported as an experimental feature that in theory is
// subject to changes in the future draft standard for C++.
// 
// Only defined when the g++ option -std=c++0x is used.
//
//
#if defined(__GXX_EXPERIMENTAL_CXX0X__)
//
// problems in compatibility with stdio printf/fprintf:
//   - variable field width (%*d) not implemented
//   - L 0 fills to the left of the number through precision width
//     (ie printf("%20.6d",12) yields 000012 in a 20 width field)
//
// What is different from C-style printf:
//   type safe
//   no truncation of type data
//   no vairable width fields
//

inline void fprintf(std::ostream &stream, const char* s)
{
    while (*s)
    {
        if (*s == '%' && *++s != '%')
            throw std::runtime_error("invalid format string: missing arguments");
        stream << *s++;
    }
}

template<typename T, typename... Args>
void fprintf(std::ostream &stream, const char* s, const T& value, const Args&... args)
{
    while (*s)
    {
        if (*s == '%' && *++s != '%')
        {
            bool leftJustify = false;
            bool zeroPad = false;
            int fieldWidth = 0;
            int precision = 3;
            char fillChar = ' ';
            if (*s && *s == '-')
            {
                leftJustify = true;
                s++;
            }

            if (*s && *s == '0')
            {
                fillChar = '0';
                zeroPad = true;
                s++;
            }

            while (*s && isdigit(*s))
            {
                fieldWidth *= 10;
                fieldWidth += (*s - '0');
                s++;
            }

            if (*s && *s == '.')
            {
                precision = 0;
                s++;
                while (*s && isdigit(*s))
                {
                    precision *= 10;
                    precision += (*s - '0');
                    s++;
                }
                s++;
            }

            while (*s)
            {
                switch (*s)
                {
                    case 's':
                        s++;
                        stream << std::setw(fieldWidth) << (leftJustify ? std::left : std::right) << value;
                        break;
                    case 'p':
                    case 'x':
                    case 'X':
                        s++;
                        stream << std::setw(fieldWidth) << std::setfill(fillChar) << (leftJustify ? std::left : std::right) << std::hex << value;
                        break;
                    case 'l':
                    case 'L':
                        s++;
                        continue;
                    case 'f':
                    case 'd':
                    case 'h':
                    case 'j':
                    case 't':
                    case 'z':
                        s++;
                        stream << std::setw(fieldWidth) << std::setfill(fillChar) << (leftJustify ? std::left : std::right) << std::dec << value;
                        break;
                    default:
                        throw std::runtime_error("Unrecognized printf conversion character");
                        break;
                }
                break;
            }
            fprintf(stream, s, args...);
            return;
        }
        stream << *s++;
    }
    throw std::runtime_error("extra arguments provided to printf");
}

template<typename T, typename... Args>
void printf(const char* s, const T& value, const Args&... args)
{
    fprintf(std::cout, s, value, args...);
}

template<typename... Args>
void sprintf(std::string &buffer, const char *fmt, const Args&... args)
{
    std::ostringstream stream;

    fprintf((std::ostream &) stream, fmt, args...);

    // stream.str() returns a const std::string &, so we
    // can't do a std::swap()
    buffer = stream.str();
}
#endif


} // end namespace STLUtilities

#endif
