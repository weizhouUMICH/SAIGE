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

#ifndef _GENOME_SEQUENCE_HELPERS_H
#define _GENOME_SEQUENCE_HELPERS_H

#if !defined(MD5_DIGEST_LENGTH)
#define MD5_DIGEST_LENGTH 16
#endif

#include "MemoryMapArray.h"

#include <stdint.h>


//
// ChromosomeInfo represents the per chromosome information
// necessary to write out SAM/BAM records.  In addition, it
// contains a single internal index used to point to the vector
// offset where the chromosome bases start.
//
// This is mildly non-optimal for larger collections of chromosomes
// or contigs - one use case described having millions of contigs,
// in which case, this structure alone would take a gigabyte or more.
//
struct ChromosomeInfo
{
    static const int  MAX_GENOME_INFO_STRING=128;

    void constructorClear()
    {
        memset(this,0, sizeof(*this));
    }
    void setChromosomeName(const char *n)
    {
        strncpy(name, n, sizeof(name)-1);
        name[sizeof(name)-1] = '\0';
    }
    genomeIndex_t   start;                              // internal offset to combined genome vector
    genomeIndex_t   size;                               // SAM SQ:LN value
    char            md5[2*MD5_DIGEST_LENGTH + 1];       // 32 chars plus NUL, SAM SQ:M5 value
    char            name[MAX_GENOME_INFO_STRING];       // SAM SQ:SN value
    char            assemblyID[MAX_GENOME_INFO_STRING]; // SAM SQ:AS value
    char            uri[MAX_GENOME_INFO_STRING];        // SAM SQ:UR value
    char            species[MAX_GENOME_INFO_STRING];    // SAM SQ:SP value

    // handy setting methods:
    void setAssemblyID(const char *newID)
    {
        strncpy(assemblyID, newID, sizeof(assemblyID)-1);
        name[sizeof(name)-1] = '\0';
    }
    void setSpecies(const char *newSpecies)
    {
        strncpy(species, newSpecies, sizeof(species)-1);
        species[sizeof(species)-1] = '\0';
    }
    void setURI(const char *newURI)
    {
        strncpy(uri, newURI, sizeof(uri)-1);
        uri[sizeof(uri)-1] = '\0';
    }
};

class genomeSequenceMmapHeader : public MemoryMapArrayHeader
{
    friend class GenomeSequence;
    friend std::ostream &operator << (std::ostream &, genomeSequenceMmapHeader &);
private:
    uint32_t    _chromosomeCount;
    bool        _colorSpace;

    ChromosomeInfo  _chromosomes[0];

public:
    //
    // getHeaderSize is special in that it must not access any
    // member variables, since it is called when the header has
    // not been created yet.
    //
    static size_t  getHeaderSize(int chromosomeCount)
    {
        return sizeof(genomeSequenceMmapHeader) + sizeof(ChromosomeInfo[1]) * chromosomeCount;
    }
    //
    // below methods return TRUE if it failed, false otherwise (primarily
    // a length check).
    //
};

std::ostream &operator << (std::ostream &stream, genomeSequenceMmapHeader &h);

//
// define the genomeSequence array type:
//
// NB the access/set routines use the encoded base values in the range
// 0-15, not the corresponding base pair character.
//
inline uint8_t genomeSequenceAccess(void *base, genomeIndex_t index)
{
    if ((index&1)==0)
    {
        return ((uint8_t *) base)[index>>1] & 0xf;
    }
    else
    {
        return (((uint8_t *) base)[index>>1] & 0xf0) >> 4;
    }
};
inline void genomeSequenceSet(void *base, genomeIndex_t index, uint8_t v)
{
    if ((index&1)==0)
    {
        ((uint8_t *) base)[index>>1] = (((uint8_t *) base)[index>>1] & 0xf0) | v;
    }
    else
    {
        ((uint8_t *) base)[index>>1] = (((uint8_t *) base)[index>>1] & 0x0f) | v<<4;
    }
}

inline size_t mmapGenomeSequenceElementCount2Bytes(genomeIndex_t i)
{
    return sizeof(uint8_t) * i / 2;
}

class PackedRead
{
    void set(int index, int val)
    {
        packedBases[index>>1] =
            (packedBases[index>>1]             // original value
             & ~(7<<((index&0x01)<<2)))         // logical AND off the original value
            | ((val&0x0f)<<((index&0x1)<<2));  // logical OR in the new value
    }
public:
    std::vector<uint8_t> packedBases;
    uint32_t    length;
    int size()
    {
        return length;
    }
    void clear()
    {
        packedBases.clear();
        length=0;
    }
    void set(const char *rhs, int padWithNCount = 0);
    uint8_t operator [](int index)
    {
        return (packedBases[index>>1] >> ((index&0x1)<<2)) & 0xf;
    }
};

#endif
