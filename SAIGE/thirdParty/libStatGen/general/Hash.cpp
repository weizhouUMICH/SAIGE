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

#include "Hash.h"

#include <ctype.h>

// ********************************************************
//
// This code is based on the original by Robert Jenkins.
//
// http://burtleburtle.net/bob/hash/doobs.html
//
// ********************************************************

#define MIX_INTEGERS(a,b,c) \
    { \
    a -= b; a -= c; a ^= (c>>13); \
    b -= c; b -= a; b ^= (a<<8);  \
    c -= a; c -= b; c ^= (b>>13); \
    a -= b; a -= c; a ^= (c>>12); \
    b -= c; b -= a; b ^= (a<<16); \
    c -= a; c -= b; c ^= (b>>5);  \
    a -= b; a -= c; a ^= (c>>3);  \
    b -= c; b -= a; b ^= (a<<10); \
    c -= a; c -= b; c ^= (b>>15); \
    }

#define ui   (unsigned int)

unsigned int hash(const unsigned char * key, unsigned int length, unsigned int initval)
{
    unsigned int a = 0x9e3779b9;
    unsigned int b = 0x9e3779b9;
    unsigned int c = initval;
    unsigned int len = length;

    /*---------------------------------------- handle most of the key */
    while (len >= 12)
    {
        a += (key[0] +(ui(key[1])<<8) +(ui(key[2])<<16) +(ui(key[3])<<24));
        b += (key[4] +(ui(key[5])<<8) +(ui(key[6])<<16) +(ui(key[7])<<24));
        c += (key[8] +(ui(key[9])<<8) +(ui(key[10])<<16)+(ui(key[11])<<24));
        MIX_INTEGERS(a,b,c);
        key += 12;
        len -= 12;
    }

    /*------------------------------------- handle the last 11 bytes */
    c += length;
    switch (len)             /* all the case statements fall through */
    {
        case 11:
            c+=(ui(key[10])<<24);
        case 10:
            c+=(ui(key[9])<<16);
        case 9 :
            c+=(ui(key[8])<<8);
            /* the first byte of c is reserved for the length */

        case 8 :
            b+=(ui(key[7])<<24);
        case 7 :
            b+=(ui(key[6])<<16);
        case 6 :
            b+=(ui(key[5])<<8);
        case 5 :
            b+=key[4];

        case 4 :
            a+=(ui(key[3])<<24);
        case 3 :
            a+=(ui(key[2])<<16);
        case 2 :
            a+=(ui(key[1])<<8);
        case 1 :
            a+=key[0];
            /* case 0: nothing left to add */
    }
    MIX_INTEGERS(a,b,c);

    /*-------------------------------------------- report the result */
    return c;
}

unsigned int hash_no_case(const unsigned char * key, unsigned int length, unsigned int initval)
{
    unsigned int a = 0x9e3779b9;
    unsigned int b = 0x9e3779b9;
    unsigned int c = initval;
    unsigned int len = length;

    /*---------------------------------------- handle most of the key */
    while (len >= 12)
    {
        a += (toupper(key[0]) +(ui(toupper(key[1]))<<8) +(ui(toupper(key[2]))<<16) +(ui(toupper(key[3]))<<24));
        b += (toupper(key[4]) +(ui(toupper(key[5]))<<8) +(ui(toupper(key[6]))<<16) +(ui(toupper(key[7]))<<24));
        c += (toupper(key[8]) +(ui(toupper(key[9]))<<8) +(ui(toupper(key[10]))<<16)+(ui(toupper(key[11]))<<24));
        MIX_INTEGERS(a,b,c);
        key += 12;
        len -= 12;
    }

    /*------------------------------------- handle the last 11 bytes */
    c += length;
    switch (len)             /* all the case statements fall through */
    {
        case 11:
            c+=(ui(toupper(key[10]))<<24);
        case 10:
            c+=(ui(toupper(key[9]))<<16);
        case 9 :
            c+=(ui(toupper(key[8]))<<8);
            /* the first byte of c is reserved for the length */

        case 8 :
            b+=(ui(toupper(key[7]))<<24);
        case 7 :
            b+=(ui(toupper(key[6]))<<16);
        case 6 :
            b+=(ui(toupper(key[5]))<<8);
        case 5 :
            b+=toupper(key[4]);

        case 4 :
            a+=(ui(toupper(key[3]))<<24);
        case 3 :
            a+=(ui(toupper(key[2]))<<16);
        case 2 :
            a+=(ui(toupper(key[1]))<<8);
        case 1 :
            a+=toupper(key[0]);
            /* case 0: nothing left to add */
    }
    MIX_INTEGERS(a,b,c);

    /*-------------------------------------------- report the result */
    return c;
}
