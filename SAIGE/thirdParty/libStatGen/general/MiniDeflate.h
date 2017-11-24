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

#ifndef __MINIDEFLATE_H__
#define __MINIDEFLATE_H__

#include <stdio.h>

// MiniDeflate reads and writes files in a simple Deflate like format
// A quick overview of this format follows, at the bottom of this file
//

// Performance tuning constants
//

// Hash table size is HASH_SIZE (a prime)
#define HASH_SIZE    4093
// Hash table depth is HASH_DEPTH (a power of 2)
#define HASH_DEPTH   8
// Matches that are not at least OKAY_MATCH chars are added to hash table
#define OKAY_MATCH   32
// Buffer size for FILE I/O
#define BUFFER_SIZE  (32 * 1024)

class MiniDeflate
{
public:
    MiniDeflate();
    ~MiniDeflate();

    void Deflate(FILE * output, void * input, size_t bytes);
    void Inflate(FILE * input, void * ouput, size_t bytes);

private:
    unsigned char *  buffer;
    unsigned char *  hash_keys;
    unsigned char ** hash_values;

    // Inline functions used during file compression
    inline void EvaluateMatch(unsigned char * in, int len, int hash,
                              unsigned char * & best_pos, int & best_match);
    inline void QuoteLiterals(unsigned char * & in, int literal,
                              unsigned char * & out, int & buffer_len,
                              FILE * output);
    inline void OutputLiterals(unsigned char * & in, int literal,
                               unsigned char * & out, int & buffer_len,
                               FILE * output);
    inline void CiteLiteral(unsigned char * & out, int literal,
                            unsigned char * & in, int & buffer_len,
                            FILE * input);
};

// Format specification for deflate files
//
// A compressed file is a sequence of bytes {0 .. N}.
// Each byte is a sequence of bits [0 .. 7] with 0 as the Most Significant Bit.
//
// The following tokens are recognized:
//
// Literal quotes -- refer to unique strings
//
//   BYTE0    BYTE1     BYTE2       Description
//     0       HI        LO         Quote of 31 bytes of more
//                                  Followed by (HI << 8 + LO + 31) quoted chars
//    0:4|LEN                       Quote of up to 1-15 bytes
//                                  Followed by LEN quoted chars
//
// String matches -- refer to previous strings in the input stream
//
//   BYTE0    BYTE1     BYTE2     BYTE3   BYTE4     Description
//  1:4|OFF   OFF1     OFF2:2|0    HI      LO       Long match of > 66 bytes
//                                                  Offset of OFF|OFF1|OFF2 + 1
//                                                  Length of HI|LO + 66
//  1:4|OFF   OFF1     OFF2:2|LEN                   Distant match of < 66 bytes
//                                                  Offset of OFF|OFF1|OFF2 + 1
//                                                  Length of LEN + 2
//  LEN|OFF   OFF1                                  Nearby short match
//                                                  Offset OFF|OFF1 + 1
//                                                  Length LEN
//

// NOTE: When partitioning bytes, I use the notation X:n|Y so that
// X takes the n MSB bits of byte and Y takes the remaining bits.


#endif


