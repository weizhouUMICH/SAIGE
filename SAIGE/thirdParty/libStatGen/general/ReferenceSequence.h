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

#ifndef __REFERENCESEQUENCE_H
#define __REFERENCESEQUENCE_H

#include "BaseAsciiMap.h"
#include "Generic.h"
#include "PackedVector.h"

#include <algorithm>
#include <iostream>

//
// The namespace Sequence is for templated algorithms that are by and large
// independent of the underlying storage mechanism for the bases.
//
// They are written in such a way as to assume that the array operator []
// will return an ASCII representation of a base space nucleotide.
//
// In theory, this set of templates will work with a variety of combinations
// of means for representing bases - String, std::string, and others.
//
// The containers are expected to allow for an overidden [] operator,
// and provide a size() method to return the number of bases in the
// container.
//
// In containers where sequence data is placed, they must in addition
// have a clear() method as well as a push_back() method as done
// in std::string containers.
//
namespace Sequence {

//
// wordMatch(Sequnece, Index, Word) - compare a word to a sequence of bases
//  Sequence is a generic container with a large set of bases
//  Index is the starting point to start the comparison
//  Word is the small sequence of bases to match
//
template<typename sequenceType, typename sequenceIndexType, typename wordType>
bool wordMatch(sequenceType &sequence,  sequenceIndexType index, wordType &word) 
{
    if( (index + word.size()) >= sequence.size() ) return false;

    for(size_t i = 0; i < word.size(); i++) {
        if( sequence[index + i] != word[i]) return false;
    }
    return true;
}

//
// printNearbyWords(output, sequence, index, word, deviation) searches
// for 'deviation' bases on either side of the index into sequence
// and prints all occurrences where word appears.
//
template<typename sequenceType, typename sequenceIndexType, typename wordType>
void printNearbyWords(std::ostream &output, sequenceType &sequence,  sequenceIndexType index, wordType &word, int deviation)
{
    for (sequenceIndexType i = index - deviation; i < index + deviation; i++)
    {
        if (wordMatch(sequence, i, word))
        {
            output << "word '"
            << word
            << "' found "
            << i - index
            << " away from position "
            << index
            << "."
            << std::endl;
        }
    }
}

//
// getString(sequence, index, baseCount, word) - populate word with the 'baseCount'
// bases that occur at the 'index' starting position in sequence.
//
template<typename sequenceType, typename sequenceIndexType, typename wordType>
void getString(sequenceType &sequence, sequenceIndexType index, int baseCount, wordType &word)
{
    word.clear();
    for (sequenceIndexType i=0; i < (sequenceIndexType) baseCount; i++)
    {
        word.push_back(sequence[index + i]);
    }
}

//
// getHighLightedString() is a debugging aid for printing "highlighted"
// subsets of bases, where the highlighting is done via turning the
// base into a lower case ASCII equivalent.
//
template<typename sequenceType, typename sequenceIndexType, typename wordType>
void getHighLightedString(
        sequenceType &sequence,
        sequenceIndexType index,
        int baseCount,
        wordType &word,
        sequenceIndexType highLightStart,
        sequenceIndexType highLightEnd)
{
    word.clear();
    for (sequenceIndexType i=0; i < (sequenceIndexType) baseCount; i++)
    {
        char base = sequence[index+i];

        if(in(index+i, highLightStart, highLightEnd))
            base = tolower(base);

        word.push_back(base);
    }
}

//
// printBaseContext() outputs a base at location 'index' along with 'baseCount'
// bases on either side of that base (default 30).
//
template<typename sequenceType, typename sequenceIndexType>
void printBaseContext(std::ostream &output, sequenceType &sequence, sequenceIndexType index, int baseCount = 30)
{
    output << "index: " << index << std::endl;
    for (sequenceIndexType i=index-baseCount; i<=index+baseCount; i++)
        output << sequence[i];
    output << std::endl;
    for (sequenceIndexType i=index-baseCount; i<index; i++)
        std::cout << " ";
    output << "^";
    output << std::endl;
}

template<typename sequenceType, typename sequenceIndexType>
void getMismatchHatString(sequenceType &sequence, sequenceIndexType location, std::string &result, std::string &read)
{
    result = "";
    for (uint32_t i=0; i < read.size(); i++)
    {
        if (read[i] == sequence[location+i])
            result.push_back(' ');
        else
            result.push_back('^');
    }
}

template<typename sequenceType, typename sequenceIndexType>
void getMismatchString(sequenceType &sequence, sequenceIndexType location, std::string &result, std::string &read)
{
    result = "";
    for (uint32_t i=0; i < read.size(); i++)
    {
        if (read[i] == sequence[location+i])
            result.push_back(toupper(read[i]));
        else
            result.push_back(tolower(read[i]));
    }
}

/// Return the mismatch count, disregarding CIGAR strings
///
/// \param read is the sequence we're counting mismatches in
/// \param location is where in the genmoe we start comparing
/// \param exclude is a wildcard character (e.g. '.' or 'N')
///
/// \return number of bases that don't match the reference, except those that match exclude

template<typename sequenceType, typename sequenceIndexType, typename readType>
int getMismatchCount(sequenceType &sequence, sequenceIndexType location, readType &read, char exclude='\0')
{
    int mismatchCount = 0;
    for (uint32_t i=0; i<read.size(); i++)
        if (read[i]!=exclude) mismatchCount += read[i]!=sequence[location + i];
    return mismatchCount;
}

/// brute force sumQ - no sanity checking
///
/// \param read shotgun sequencer read string
/// \param qualities phred quality string of same length
/// \param location the alignment location to check sumQ
template<typename sequenceType, typename sequenceIndexType, typename readType, typename qualityType>
int getSumQ(sequenceType &sequence, sequenceIndexType location, readType &read, qualityType &qualities)
{
    int sumQ = 0;
    for (uint32_t i=0; i<read.size(); i++)
        sumQ += (read[i]!=sequence[location + i] ? (qualities[i]-33) : 0);
    return sumQ;
};



template<typename sequenceType, typename sequenceIndexType, typename readType, typename qualityType>
sequenceIndexType simpleLocalAligner(
        sequenceType &sequence,
        sequenceIndexType index,
        readType &read,
        qualityType &quality,
        int windowSize)
{
    int bestScore = 1000000; // either mismatch count or sumQ
    sequenceIndexType bestMatchLocation = -1;
    for (int i=-windowSize; i<windowSize; i++)
    {
        int newScore;

        // 'in' is a template, and types of all three args have to match
        if( not in((size_t) index + i, (size_t) 0, sequence.size())) continue;

        if (quality.size() == 0)
        {
            newScore = getMismatchCount(sequence, index + i, read);
        }
        else
        {
            newScore = getSumQ(sequence, index + i, read, quality);
        }
        if (newScore < bestScore)
        {
            bestScore = newScore;
            bestMatchLocation = index + i;
        }
    }
    return bestMatchLocation;
}


}

typedef uint32_t PackedVectorIndex_t;

//
// this is actually a base space implementation, since
// the load/set methods do indirection from the symbol
// value through a symbol value to symbol lookup table.
//
class PackedSequenceData : public PackedVector4Bit_t
{

    public:
    inline char operator[](PackedVectorIndex_t baseIndex)
    {
        return BaseAsciiMap::int2base[ (static_cast<PackedVector4Bit_t>(*this))[baseIndex]];
    }
    inline void set(PackedVectorIndex_t baseIndex, char value)
    {
        this->PackedVector4Bit_t::set(baseIndex, BaseAsciiMap::base2int[(uint32_t) value]);
    }
    inline void push_back(char value)
    {
        this->PackedVector4Bit_t::push_back(BaseAsciiMap::base2int[(uint32_t) value]);
    }
};

std::ostream &operator << (std::ostream &stream, PackedSequenceData &v)
{
    for(size_t i=0; i<v.size(); i++) {
        stream << i << ": " << v[i] << std::endl;
    }
    return stream;
}


//
// Load a fasta format file from filename into the buffer
// provided by the caller.
// While parsing the fasta file, record each chromosome name,
// its start location, and its size.
//
// NB: the caller must implement the logic to determine how
// large the sequence data is.  There is no correct way to do
// this, because we can't reliably estimate here how much sequence
// data is contained in a compressed file.
//
// To safely pre-allocate space in sequenceData, use the reserve() method
// before calling this function.
//
template<typename SequenceDataType, typename ChromosomeNameType>
bool loadFastaFile(const char *filename,
        std::vector<SequenceDataType> &sequenceData,
        std::vector<ChromosomeNameType> &chromosomeNames)
{
    InputFile inputStream(filename, "r", InputFile::DEFAULT);

    if(!inputStream.isOpen()) {
        std::cerr << "Failed to open file " << filename << "\n";
        return true;
    }

    int whichChromosome = -1;

    //
    // chromosomeNames is cheap to clear, so do it here.
    //
    // NB: I explicitly choose not to clear the sequence data
    // container, this allows the caller to pre-allocate based
    // on their knowledge of the size of the expected genome.
    //

    chromosomeNames.clear();

    char ch;
    while((ch = inputStream.ifgetc()) != EOF) {
        switch (ch)
        {
            case '\n':
            case '\r':
                break;
            case '>':
            {
                std::string chromosomeName = "";
                //
                // pull out the chromosome new name
                //
                while (!isspace((ch = inputStream.ifgetc())) && ch != EOF)
                {
                    chromosomeName += ch;  // slow, but who cares
                }
                //
                // eat the rest of the line
                //
                do {
                    ch = inputStream.ifgetc();
                } while(ch != EOF && ch != '\n' && ch != '\r');
                //
                // save the Chromosome name and index into our
                // header so we can use them later.
                //
                chromosomeNames.push_back(chromosomeName);

                whichChromosome++;

                sequenceData.resize(whichChromosome+1);

                break;
            }
            default:
                // we get here for sequence data.
                //
                // save the base value
                // Note: invalid characters come here as well, but we
                // let ::set deal with mapping them.

                sequenceData[whichChromosome].push_back(toupper(ch));
#if 0
                if (isColorSpace())
                {
//
// anything outside these values represents an invalid base
// base codes: 0-> A,    1-> C,     2-> G,      3-> T
// colorspace: 0-> blue, 1-> green, 2-> oragne, 3->red
//
                    const char fromBase2CS[] =
                    {
                        /* 0000 */ 0,   // A->A
                        /* 0001 */ 1,   // A->C
                        /* 0010 */ 2,   // A->G
                        /* 0011 */ 3,   // A->T
                        /* 0100 */ 1,   // C->A
                        /* 0101 */ 0,   // C->C
                        /* 0110 */ 3,   // C->G
                        /* 0111 */ 2,   // C->T
                        /* 1000 */ 2,   // G->A
                        /* 1001 */ 3,   // G->C
                        /* 1010 */ 0,   // G->G
                        /* 1011 */ 1,   // G->T
                        /* 1100 */ 3,   // T->A
                        /* 1101 */ 2,   // T->C
                        /* 1110 */ 1,   // T->G
                        /* 1111 */ 0,   // T->T
                    };
                    //
                    // we are writing color space values on transitions,
                    // so we don't write a colorspace value when we
                    // get the first base value.
                    //
                    // On second and subsequent bases, write based on
                    // the index table above
                    //
                    char thisBase = base2int[(int)(fasta[fastaIndex])];
                    if (lastBase>=0)
                    {
                        char color;
                        if (lastBase>3 || thisBase>3) color=4;
                        else color = fromBase2CS[(int)(lastBase<<2 | thisBase)];
                        // re-use the int to base, because ::set expects a base char (ATCG), not
                        // a color code (0123).  It should only matter on final output.
                        set(header->elementCount++, int2base[(int) color]);
                    }
                    lastBase = thisBase;
                }
                else
                {
                    set(header->elementCount++, toupper(fasta[fastaIndex]));
                }
#endif
                break;
        }
    }
    return false;
}

#endif
