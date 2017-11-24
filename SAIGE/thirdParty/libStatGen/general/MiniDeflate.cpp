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

#include "MiniDeflate.h"

// Convenient constants and macros
//
#define EMPTY_KEY    123
#define uchar        unsigned char

#ifndef min
#define min(a,b)   (((a)<(b))?(a):(b))
#endif

MiniDeflate::MiniDeflate()
{
    buffer = new uchar [BUFFER_SIZE + 5];
    hash_keys = new uchar [HASH_SIZE];
    hash_values = new uchar * [HASH_SIZE * HASH_DEPTH];
}

MiniDeflate::~MiniDeflate()
{
    delete [] buffer;
    delete [] hash_keys;
    delete [] hash_values;
}

void MiniDeflate::EvaluateMatch(unsigned char * in, int len, int hash,
                                unsigned char * & best_pos, int & best_match)
{
    int max = min(len, 0xFFFF + 66);

    for (int i = HASH_DEPTH; i > 0; i--)
        // Check each possible match (up to HASH_DEPTH)
    {
        uchar * pos = hash_values[hash * HASH_DEPTH + ((hash_keys[hash] + i) % HASH_DEPTH)];

        if (pos == NULL || in - pos >= 0x4001) break;

        int match = 0;

        while (match < max && pos[match] == in[match])
            match++;

        if (match > best_match)
        {
            best_match = match;
            best_pos = pos;
        }
    }

    // If string seems pretty unique, add to hash table
    if (best_match < OKAY_MATCH)
    {
        int delta = hash_keys[hash] = (uchar)((hash_keys[hash] + 1) & 7);
        hash_values[hash * 8 + delta] = in;
    }
}

void MiniDeflate::QuoteLiterals(unsigned char * & in, int literal,
                                unsigned char * & out, int & buffer_len,
                                FILE * output)
{
    if (buffer_len < 0)
    {
        fwrite(buffer, out - buffer, 1, output);
        buffer_len = BUFFER_SIZE;
        out = buffer;
    }

    while (buffer_len < literal)
    {
        literal -= buffer_len;
        while (buffer_len--)
        {
            *out = *in;
            in++;
            out++;
        }
        fwrite(buffer, BUFFER_SIZE, 1, output);
        buffer_len = BUFFER_SIZE;
        out = buffer;
    }

    while (literal--)
    {
        *out = *in;
        in++;
        out++;
        buffer_len--;
    }
}

void MiniDeflate::OutputLiterals(unsigned char * & in, int literal,
                                 unsigned char * & out, int & buffer_len,
                                 FILE * output)
{
    while (literal > 0)
        if (literal < 16)
        {
            *out = (char) literal;
            out++;
            buffer_len--;
            QuoteLiterals(in, literal, out, buffer_len, output);
            break;
        }
        else if (literal < 31)
        {
            *out = 15;
            out++;
            buffer_len--;
            QuoteLiterals(in, 15, out, buffer_len, output);
            *out = (uchar)(literal - 15);
            out++;
            buffer_len--;
            QuoteLiterals(in, literal - 15, out, buffer_len, output);
            break;
        }
        else
        {
            int length = min(literal, 0xFFFF + 31);
            literal   -= length;
            length    -= 31;

            *out = 0;
            out++;
            *out = (uchar)(length >> 8);
            out++;
            *out = (uchar)(length & 0xFF);
            out++;
            buffer_len -= 3;

            QuoteLiterals(in, length + 31, out, buffer_len, output);
        }
}


void MiniDeflate::Deflate(FILE * output, void * void_input, size_t len)
{
    uchar    * in  = (uchar *) void_input;
    uchar    * out = (uchar *) buffer;
    int buffer_len = BUFFER_SIZE;

    for (int i = 0; i < HASH_SIZE; i++) hash_keys[i] = EMPTY_KEY;

    uchar * in2 = in;

    while (len > 2)
    {
        // Hash the current input value
        int hash = ((in[0] << 16) | (in[1] << 8) | in[2]) % HASH_SIZE;

        if (hash_keys[hash] != EMPTY_KEY)
            // Possible matches in hash table
        {
            int best_match = 0;
            uchar * best_pos = NULL;

            EvaluateMatch(in, len, hash, best_pos, best_match);

            // If there are no decent matches
            if (best_match < 3)
            {
                in++;
                len--;
                continue;
            }

            // Try look ahead if match isn't great
            while (best_match < OKAY_MATCH && len > 3)
            {
                // Peek to see if we could get a better match
                int next_hash = ((in[1] << 16) | (in[2] << 8) | in[3]) % HASH_SIZE;

                if (hash_keys[next_hash] == EMPTY_KEY) break;

                int next_match = 0;
                uchar * next_pos = NULL;

                EvaluateMatch(in + 1, len - 1, next_hash, next_pos, next_match);

                // Didn't find a better match
                if (next_match <= best_match + 1) break;

                // Found a better match, so try again
                in++;
                len--;
                best_match = next_match;
                best_pos = next_pos;
            }

            int best_offset = in - best_pos - 1;

            // This is where we output stuff
            // Check if we have some literals to output first
            OutputLiterals(in2, in - in2, out, buffer_len, output);

            in2 = in += best_match;
            len -= best_match;

            if (best_match < 17 && best_offset < 0x1000)
            {
                *out = (uchar)(((best_match - 1) << 4) | (best_offset >> 8));
                out++;
                *out = (uchar)(best_offset & 0xFF);
                out++;
                buffer_len -= 2;
            }
            else if (best_match < 66)
            {
                *out = (uchar)(16 | (best_offset >> 10));
                out++;
                *out = (uchar)((best_offset >> 2) & 0xFF);
                out++;
                *out = (uchar)((best_offset << 6) | (best_match - 2));
                out++;
                buffer_len -= 3;
            }
            else
            {
                *out = (uchar)(16 | (best_offset >> 10));
                out++;
                *out = (uchar)((best_offset >> 2) & 0xFF);
                out++;
                *out = (uchar)(best_offset << 6);
                out++;
                best_match -= 66;
                *out = (uchar)(best_match >> 8);
                out++;
                *out = (uchar)(best_match & 0xFF);
                out++;
                buffer_len -= 5;
            }

            if (buffer_len <= 0)
            {
                fwrite(buffer, out - buffer, 1, output);
                buffer_len = BUFFER_SIZE;
                out = buffer;
            }
        }
        // Never seen this sequence before
        else
        {
            hash_keys[hash] = 0;
            for (int i = 1; i < HASH_DEPTH; i++) hash_values[hash * 8 + i] = NULL;
            hash_values[hash * 8] = in;
            in++;
            len--;
        }
    }

    // Check if we have some trailing literals to output
    in += len;
    OutputLiterals(in2, in - in2, out, buffer_len, output);

    // Flush output
    if (out != buffer) fwrite(buffer, out - buffer, 1, output);
}

void MiniDeflate::CiteLiteral(unsigned char * & out, int literal,
                              unsigned char * & in, int & buffer_len,
                              FILE * input)
{
    while (buffer_len < literal)
    {
        literal -= buffer_len;
        while (buffer_len--)
        {
            *out = *in;
            in++;
            out++;
        }
        buffer_len = fread(buffer + 5, 1, BUFFER_SIZE, input);
        in = buffer + 5;
    }

    while (literal--)
    {
        *out = *in;
        in++;
        out++;
        buffer_len--;
    }
}


void MiniDeflate::Inflate(FILE * input, void * void_output, size_t len)
{
    uchar    * out = (uchar *) void_output;
    uchar    * in  = (uchar *) buffer + 5;
    int buffer_len = BUFFER_SIZE;

    buffer_len = fread(buffer + 5, 1, BUFFER_SIZE, input);

    while (len)
    {
        int match_len = *in >> 4;

        // Matching a literal
        if (match_len == 0)
        {
            match_len = *in & 0x0F;
            in++, buffer_len--;

            // If match_len == 0 then string is longer than 30 characters
            // Strings of 16 - 30 characters are encoded as two short strings
            if (match_len == 0)
            {
                match_len = (in[0] << 8) + in[1] + 31;
                in += 2;
                buffer_len -= 2;
            }

            CiteLiteral(out, match_len, in, buffer_len, input);
            len -= match_len;
        }
        // Long match, 14 bit offset
        else if (match_len == 1)
        {
            int offset = (((in[0] & 0x0F) << 10) | (in[1] << 2) | (in[2] >> 6)) + 1;
            match_len  = (in[2] & 0x3F) + 2;
            in += 3;
            buffer_len -= 3;

            if (match_len == 2)
            {
                match_len = ((in[0] << 8) | in[1]) + 66;
                in  += 2;
                buffer_len -= 2;
            }

            uchar * match_pos = out - offset;
            len -= match_len;
            while (match_len--)
            {
                *out = *match_pos;
                out++, match_pos++;
            }
        }
        // Typical short match
        else
        {
            int offset = (((in[0] & 0x0F) << 8) | in[1]) + 1;
            in += 2;
            buffer_len -= 2;

            uchar * match_pos = out - offset;
            len -= ++match_len;
            while (match_len--)
            {
                *out = *match_pos;
                out++, match_pos++;
            }
        }

        if (buffer_len < 5)
        {
            uchar * in2 = (uchar *) buffer + 5 - buffer_len;
            while (in2 != buffer + 5)
            {
                *in2 = *in;
                in2++;
                in++;
            }

            in = buffer + 5 - buffer_len;
            buffer_len += fread(buffer + 5, 1, BUFFER_SIZE, input);
        }
    }

    if (buffer_len) fseek(input, -buffer_len, SEEK_CUR);
}



