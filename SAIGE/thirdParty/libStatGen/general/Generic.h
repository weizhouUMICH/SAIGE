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

#if !defined(_GENERIC_H)
#define _GENERIC_H

#include <stdint.h>

#include <list>
#include <iostream>
#include <ostream>
#include <utility>
#include <vector>
#include <string>

template <typename T>
inline T abs(T x)
{
    return (x < 0) ? -x : x;
}

//
// this is safe for signed/unsigned:
//
template <typename T>
inline T absDiff(T x, T y)
{
    return (x < y) ? (y - x) : (x - y);
}

//
//
template <typename T>
inline T in(T x, T y, T z)
{
    return (x >= y && x < z);
}

//
// These overloaded operators and functions are largely
// for diagnostic debug printing.  The underlying problem
// is that gdb is unable to decipher any STL use, let alone
// complex STL use.  printf debugging is a poor second choice,
// but these functions at least make it practical to do rapidly.
//

//
// Write a std::pair to a stream
//
template <typename A, typename B>
std::ostream &operator << (std::ostream &stream, std::pair<A, B> p)
{
    stream << "(" << p.first << ", " << p.second << ")";
    return stream;
}

//
// generic vector print -- in normal use, you should
// be able to simply do foostream << somevector, and get
// sane results, provided that the vector elements themselves
// can be written to the stream.
//
// Example code is in Generic.cpp
//
template <typename T>
std::ostream &operator << (std::ostream &stream, std::vector<T> const &v)
{

    typename std::vector<T>::const_iterator i;
    for (i = v.begin(); i != v.end(); i++)
    {
        stream << (i - v.begin()) << ": " << *i << std::endl;
    }
    return stream;
}

//
// same overload as above, except for std::list
//
template <typename T>
std::ostream &operator << (std::ostream &stream, std::list<T> const &l)
{

    typename std::list<T>::const_iterator i;
    int j = 0;
    for (i = l.begin(); i != l.end(); i++, j++)
    {
        stream << j << ": " << *i << std::endl;
    }
    return stream;
}

template <typename TITLE, typename ITEM, typename EXPECT, typename GOT>
void check(int &returnCode, TITLE title, ITEM item, EXPECT expect, GOT got)
{
    if (expect!=got)
    {
        std::cout << "Test " << title << ": expect " << item << " = '" << expect << "', but got '" << got << "'." << std::endl;
        returnCode += 1;
    }
}

//
// specialization of template below:
// load a set of lines from a file into a vector of strings.
//
inline std::istream &operator >> (std::istream &stream, std::vector<std::string> &vec)
{
    std::string val;
    while (true)
    {
        if (!stream.good()) break;
        getline(stream, val);
        stream >> val;
        vec.push_back(val);
    }
    return stream;
}


//
// read values from a stream, appending to the provided
// vec.  stops when the stream is consumed.
//
template<typename T>
std::istream &operator >> (std::istream &stream, std::vector<T> &vec)
{
    T val;
    while (true)
    {
        if (!stream.good()) break;
        stream >> val;
        vec.push_back(val);
    }
    return stream;
}


#if 0
//
// generic vector of iterators print
//
template <typename T>
std::ostream &operator << (
    std::ostream &stream,
    std::vector<
    std::pair< std::vector<typename T>::iterator , std::vector< typename T>::iterator >
    > v
)
{

    typename IteratorType i;
    typename std::vector<T>::iterator i;
    for (i = v.begin(); i != v.end(); i++)
    {
        stream << *i << std::endl;
    }
    return stream;
}

#endif


//
// These are packed set/get functions for dealing with
// packed 1, 2 and 4 bit unsigned values inside of arbitrary
// arrays of data (char */std::vector<char> whatever).
//
template<typename T>
inline uint32_t PackedAccess_1Bit(T byteSequence, uint32_t bitIndex)
{
    return (((byteSequence)[bitIndex>>3] >> (bitIndex&0x7)) & 0x1);
}

template<typename T>
inline void PackedAssign_1Bit(T byteSequence, uint32_t bitIndex, uint32_t value)
{
    (byteSequence)[bitIndex>>3] =
        ((byteSequence)[bitIndex>>3]
        & ~(1<<(bitIndex&0x07)))
        | ((value&0x01)<<(bitIndex&0x7));
}

inline size_t Packed1BitElementCount2Bytes(uint32_t i)
{
    return (size_t)(i+7)/8;
}

template<typename T>
inline uint32_t PackedAccess_2Bit(T byteSequence, uint32_t index)
{
    return (((byteSequence)[index>>2] >> ((index&0x3)<<1)) & 0x3);
}

template<typename T>
inline void PackedAssign_2Bit(T byteSequence, uint32_t index, uint32_t value)
{
    (byteSequence)[index>>2] =
        ((byteSequence)[index>>2]
        & ~(3<<((index&0x03)<<1)))
        | ((value&0x03)<<((index&0x3)<<1));
}

inline size_t Packed2BitElementCount2Bytes(uint32_t i)
{
    return (size_t)(i+3)/4;
}

template<typename T>
inline uint32_t PackedAccess_4Bit(T byteSequence, uint32_t index)
{
    return (((byteSequence)[index>>1] >> ((index&0x1)<<2)) & 0xf);
}

template<typename T>
inline void PackedAssign_4Bit(T byteSequence, uint32_t index, uint32_t value)
{
    (byteSequence)[index>>1] =
        ((byteSequence)[index>>1]
        & ~(7<<((index&0x01)<<2)))
        | ((value&0x0f)<<((index&0x1)<<2));
}

inline size_t Packed4BitElementCount2Bytes(uint32_t i)
{
    return (size_t)(i+1)/2;
}

#endif
