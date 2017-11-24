/*
 *  Copyright (C) 2010-2012  Regents of the University of Michigan
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

#include "assert.h"
#include "ctype.h"
#include "stdio.h"
#include "Error.h"


#include "Generic.h"
#include "GenomeSequence.h"

#include <algorithm>
#include <istream>
#include <fstream>
#include <sstream>
#include <stdexcept>

#if defined(_WIN32)
#include <io.h>
#ifndef R_OK
#define R_OK 4
#endif
#endif

// not general use:
#include "CSG_MD5.h"

//
// given a read in a string, pack it into the vector of
// bytes coded as two bases per byte.
//
// The goal is to allow us to more rapidly compare against
// the genome, which is itself packed 2 bases per byte.
//
// Unfortunately, the match position may be odd or even,
// so as a result, we also need to be able to prepad
// with 'N' bases to get the byte alignment the same, so
// padWithNCount may be a positive number indicating how
// many N bases to prepend.
//
void PackedRead::set(const char *rhs, int padWithNCount)
{
    clear();

    // pad this packed read with 'N' bases two at a time
    while (padWithNCount>1)
    {
        packedBases.push_back(
            BaseAsciiMap::base2int[(int) 'N'] << 4 |
            BaseAsciiMap::base2int[(int) 'N']
        );
        padWithNCount -= 2;
        length+=2;
    }

    // when we have only one base, pack one 'N' base with
    // the first base in rhs if there is one.
    if (padWithNCount)
    {
        // NB: *rhs could be NUL, which is ok here - just keep
        // the length straight.
        packedBases.push_back(
            BaseAsciiMap::base2int[(int) *rhs] << 4 |
            BaseAsciiMap::base2int[(int) 'N']
        );
        // two cases - have characters in rhs or we don't:
        if (*rhs)
        {
            length+=2;  // pad byte plus one byte from rhs
            rhs++;
        }
        else
        {
            length++;
        }
        padWithNCount--;    // should now be zero, so superfluous.
        assert(padWithNCount==0);
    }

    // pad pairs of bases from rhs, two at a time:
    while (*rhs && *(rhs+1))
    {
        packedBases.push_back(
            BaseAsciiMap::base2int[(int) *(rhs+1)] << 4 |
            BaseAsciiMap::base2int[(int) *(rhs+0)]
        );
        rhs+=2;
        length+=2;
    }

    // if there is an odd base left at the end, put it
    // in a byte all its own (low 4 bits == 0):
    if (*rhs)
    {
        packedBases.push_back(
            BaseAsciiMap::base2int[(int) *(rhs+0)]
        );
        length++;
    }
    return;
}

std::string GenomeSequence::IntegerToSeq(unsigned int n, unsigned int wordsize) const
{
    std::string sequence("");
    for (unsigned int i = 0; i < wordsize; i ++)
        sequence += "N";

    unsigned int clearHigherBits = ~(3U << (wordsize<<1));  // XXX refactor - this appears several places

    if (n > clearHigherBits)
        error("%d needs to be a non-negative integer < clearHigherBits\n", n);

    for (unsigned int i = 0; i < wordsize; i++)
    {
        sequence[wordsize-1-i] = BaseAsciiMap::int2base[n & 3];
        n >>= 2;
    }
    return sequence;
}



GenomeSequence::GenomeSequence()
{
    constructorClear();
}

void GenomeSequence::constructorClear()
{
    _debugFlag = 0;
    _progressStream = NULL;
    _colorSpace = false;
    _createOverwrite = false;
}

void GenomeSequence::setup(const char *referenceFilename)
{
    setReferenceName(referenceFilename);

    if (_progressStream) *_progressStream << "open and prefetch reference genome " << referenceFilename << ": " << std::flush;

    if (open(false))
    {
        std::cerr << "Failed to open reference genome " << referenceFilename << std::endl;
        std::cerr << errorStr << std::endl;
        exit(1);
    }

    prefetch();
    if (_progressStream) *_progressStream << "done." << std::endl << std::flush;
}

GenomeSequence::~GenomeSequence()
{
    // free up resources:
    _umfaFile.close();
}

//
// mapped open.
//
// if the file exists, map in into memory, and fill in a few useful
// fields.
//

bool GenomeSequence::open(bool isColorSpace, int flags)
{
    bool rc;

    if (isColorSpace)
    {
        _umfaFilename = _baseFilename + "-cs.umfa";
    }
    else
    {
        _umfaFilename = _baseFilename + "-bs.umfa";
    }

    if(access(_umfaFilename.c_str(), R_OK) != 0)
    {
        // umfa file doesn't exist, so try to create it.
        if(create(isColorSpace))
        {
            // Couldon't access or create the umfa.
            std::cerr << "GenomeSequence::open: failed to open file "
                      << _umfaFilename
                      << " also failed creating it."
                      << std::endl;
            return true;
        }
    }

    rc = genomeSequenceArray::open(_umfaFilename.c_str(), flags);
    if (rc)
    {
        std::cerr << "GenomeSequence::open: failed to open file "
                  << _umfaFilename
                  << std::endl;
        return true;
    }

    _colorSpace = header->_colorSpace;

    return false;
}

void GenomeSequence::sanityCheck(MemoryMap &fasta) const
{
    unsigned int i;

    unsigned int genomeIndex = 0;
    for (i=0; i<fasta.length(); i++)
    {
        switch (fasta[i])
        {
            case '>':
                while (fasta[i]!='\n' && fasta[i]!='\r') i++;
                break;
            case '\n':
            case '\r':
                break;
            default:
                assert(BaseAsciiMap::base2int[(int)(*this)[genomeIndex]] == 
                       BaseAsciiMap::base2int[(int) fasta[i]]);
                genomeIndex++;
                break;
        }
    }
}

#define HAS_SUFFIX(str, suffix) ((strlen(suffix) < str.size()) && (str.substr(str.size() - strlen(suffix)) == suffix))
//
// referenceFilename is  either a fasta or a UM fasta (.fa or .umfa)
// filename.  In both cases, the suffix gets removed and the
// base name is kept for later use depending on context.
// @return always return false
//
bool GenomeSequence::setReferenceName(std::string referenceFilename)
{

    if (HAS_SUFFIX(referenceFilename, ".fa"))
    {
        _referenceFilename = referenceFilename;
        _baseFilename = _referenceFilename.substr(0, referenceFilename.size() - 3);
    }
    else if (HAS_SUFFIX(referenceFilename, ".umfa"))
    {
        _baseFilename = referenceFilename.substr(0, referenceFilename.size() - 5);
    }
    else if (HAS_SUFFIX(referenceFilename, "-cs.umfa"))
    {
        _baseFilename = referenceFilename.substr(0, referenceFilename.size() - 8);
    }
    else if (HAS_SUFFIX(referenceFilename, "-bs.umfa"))
    {
        _baseFilename = referenceFilename.substr(0, referenceFilename.size() - 8);
    }
    else
    {
        _baseFilename = referenceFilename;
    }
    _fastaFilename = _baseFilename + ".fa";

    if (HAS_SUFFIX(referenceFilename, ".fasta"))
    {
        _referenceFilename = referenceFilename;
        _baseFilename = _referenceFilename.substr(0, referenceFilename.size() - 6);
        _fastaFilename = _baseFilename + ".fasta";        
    }
    
    return false;
}

//
// this works in lockstep with ::create to populate
// the per chromosome header fields size and md5
// checksum.
//
// It relies on header->elementCount being set to
// the data length loaded so far ... not the ultimate
// reference length.
//
bool GenomeSequence::setChromosomeMD5andLength(uint32_t whichChromosome)
{
    if (whichChromosome>=header->_chromosomeCount) return true;

    ChromosomeInfo *c = &header->_chromosomes[whichChromosome];
    c->size = header->elementCount - c->start;

    MD5_CTX md5Context;
    uint8_t md5Signature[MD5_DIGEST_LENGTH];

    //
    // it's easier to ensure we do this right if we just do it
    // in one big chunk:
    //
    char *md5Buffer = (char *) malloc(c->size);

    MD5Init(&md5Context);

    for (genomeIndex_t i = 0; i < c->size; i ++)
    {
        md5Buffer[i] = (*this)[c->start + i];
    }
    MD5Update(&md5Context, (unsigned char *) md5Buffer, c->size);
    MD5Final((unsigned char *) &md5Signature, &md5Context);
    free(md5Buffer);
    for (int i=0; i<MD5_DIGEST_LENGTH; i++)
    {
        // cheesy but works:
        sprintf(c->md5+2*i, "%02x", md5Signature[i]);
    }
    // redundant, strictly speaking due to sprintf NUL terminating
    // it's output strings, but put it here anyway.
    c->md5[2*MD5_DIGEST_LENGTH] = '\0';

    return false;
}

//
// Given a buffer with a fasta format contents, count
// the number of chromosomes in it and return that value.
//
static bool getFastaStats(const char *fastaData, size_t fastaDataSize, uint32_t &chromosomeCount, uint64_t &baseCount)
{
    chromosomeCount = 0;
    baseCount = 0;
    bool atLineStart = true;

    //
    // loop over the fasta file, essentially matching for the
    // pattern '^>.*$' and counting them.
    //
    for (size_t fastaIndex = 0; fastaIndex < fastaDataSize; fastaIndex++)
    {
        switch (fastaData[fastaIndex])
        {
            case '\n':
            case '\r':
                atLineStart = true;
                break;
            case '>':
            {
                if (!atLineStart) break;
                chromosomeCount++;
                //
                // eat the rest of the line
                //
                while (fastaIndex < fastaDataSize && fastaData[fastaIndex]!='\n' && fastaData[fastaIndex]!='\r')
                {
                    fastaIndex++;
                }
                break;
            }
            default:
                baseCount++;
                atLineStart = false;
                break;
        }

    }
    return false;
}

class PackedSequenceData : public std::vector<uint8_t>
{
    std::vector<uint8_t> m_packedBases;

    size_t    m_baseCount;

    void set(size_t index, uint8_t value) {
        m_packedBases[index>>1] =
            (m_packedBases[index>>1]            // original value
            & ~(7<<((index&0x01)<<2)))          // logical AND off the original value
            | ((value&0x0f)<<((index&0x1)<<2)); // logical OR in the new value
    }

public:

    void reserve(size_t baseCount) {m_packedBases.reserve(baseCount/2);}
    size_t size() {return m_baseCount;}
    void clear() {m_packedBases.clear(); m_baseCount = 0;}
    uint8_t operator [](size_t index)
    {
        return (m_packedBases[index>>1] >> ((index&0x1)<<2)) & 0xf;
    }
    void push_back(uint8_t base);
};

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
bool loadFastaFile(const char *filename,
        std::vector<PackedSequenceData> &sequenceData,
        std::vector<std::string> &chromosomeNames)
{
    InputFile inputStream(filename, "r", InputFile::DEFAULT);

    if(!inputStream.isOpen()) {
        std::cerr << "Failed to open file " << filename << "\n";
        return true;
    }

    int whichChromosome = -1;
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

                break;
            }
            default:
                // we get here for sequence data.
                //
                // save the base value
                // Note: invalid characters come here as well, but we
                // let ::set deal with mapping them.
                break;
        }
    }
    return false;
}

//
// recreate the umfa file from a reference fasta format file
//
// The general format of a FASTA file is best described
// on wikipedia at http://en.wikipedia.org/wiki/FASTA_format
//
// The format parsed here is a simpler subset, and is
// described here http://www.ncbi.nlm.nih.gov/blast/fasta.shtml
//
bool GenomeSequence::create(bool isColor)
{
    setColorSpace(isColor);

    if (_baseFilename=="")
    {
        std::cerr << "Base reference filename is empty." << std::endl;
        return true;
    }

    if (isColorSpace())
    {
        _umfaFilename = _baseFilename + "-cs.umfa";
    }
    else
    {
        _umfaFilename = _baseFilename + "-bs.umfa";
    }

    if (!_createOverwrite && access(_umfaFilename.c_str(), R_OK) == 0)
    {
        std::cerr << "Output file '" << _umfaFilename << "' exists or is not writable - please remove." << std::endl;
        return true;
    }

    MemoryMap fastaFile;

    if (fastaFile.open(_fastaFilename.c_str()))
    {
        std::cerr << "failed to open input fasta file '" << _fastaFilename << "'." << std::endl;
        return true;
    }

    std::cerr << "Creating FASTA "
              << (isColorSpace() ? "color space " : "")
              << "binary cache file '"
              << _umfaFilename
              << "'."
              << std::endl;

    std::cerr << std::flush;

    //
    // simple ptr to fasta data -- just treat the memory map
    // as an array of fastaDataSize characters...
    //
    const char *fasta = (const char *) fastaFile.data;
    size_t fastaDataSize = fastaFile.length();

    uint32_t    chromosomeCount = 0;
    uint64_t    baseCount = 0;
    getFastaStats(fasta, fastaDataSize, chromosomeCount, baseCount);

    if (genomeSequenceArray::create(_umfaFilename.c_str(), baseCount, chromosomeCount))
    {
        std::cerr << "failed to create '"
                  << _umfaFilename
                  << "'."
                  << std::endl;
        perror("");
        return true;
    }
    header->elementCount = 0;
    header->_colorSpace = isColorSpace();
    header->setApplication(_application.c_str());
    header->_chromosomeCount = chromosomeCount;
    //
    // clear out the variable length chromosome info array
    //
    for (uint32_t i=0; i<header->_chromosomeCount; i++) header->_chromosomes[i].constructorClear();

    std::string chromosomeName;

    //
    // for converting the reference to colorspace, the first base is always 5 (in base space it is 'N')
    signed char lastBase = BaseAsciiMap::base2int[(int) 'N'];
    bool terminateLoad = false;
    int percent = -1, newPercent;
    uint32_t whichChromosome = 0;
    for (uint64_t fastaIndex = 0; fastaIndex < fastaDataSize; fastaIndex++)
    {
        if (_progressStream)
        {
            newPercent = (int) (1.0 * fastaIndex / fastaDataSize) * 100;
            if (newPercent>percent)
            {
                *_progressStream << "\r" << newPercent << "% ";
                *_progressStream << std::flush;
                percent = newPercent;
            }
        }
        switch (fasta[fastaIndex])
        {
            case '\n':
            case '\r':
                break;
            case '>':
            {
                chromosomeName = "";
                fastaIndex++;       // skip the > char
                //
                // pull out the chromosome new name
                //
                while (!isspace(fasta[fastaIndex]))
                {
                    chromosomeName += fasta[fastaIndex++];  // slow, but who cares
                }
                //
                // eat the rest of the line
                //
                while (fasta[fastaIndex]!='\n' && fasta[fastaIndex]!='\r')
                {
                    fastaIndex++;
                }
                //
                // save the Chromosome name and index into our
                // header so we can use them later.
                //
                ChromosomeInfo *c = &header->_chromosomes[whichChromosome];
                c->setChromosomeName(chromosomeName.c_str());
                c->start = header->elementCount;
                // c->size gets computed at the next '>' line or at the EOF

                if (whichChromosome>0)
                {
                    //
                    // compute md5 checksum for the chromosome that we just
                    // loaded (if there was one) - note that on the last
                    // chromosome, we have to duplicate this code after
                    // the end of this loop
                    //
                    setChromosomeMD5andLength(whichChromosome - 1);
                }
                whichChromosome++;
                if (whichChromosome > header->_chromosomeCount)
                {
                    std::cerr << "BUG: Exceeded computed chromosome count ("
                              << header->_chromosomeCount
                              << ") - genome is now truncated at chromosome "
                              << header->_chromosomes[header->_chromosomeCount-1].name
                              << " (index "
                              << header->_chromosomeCount
                              << ")."
                              << std::endl;
                    terminateLoad = true;
                }
                break;
            }
            default:
                // save the base pair value
                // Note: invalid characters come here as well, but we
                // let ::set deal with mapping them.
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
                    char thisBase = 
                        BaseAsciiMap::base2int[(int)(fasta[fastaIndex])];
                    if (lastBase>=0)
                    {
                        char color;
                        if (lastBase>3 || thisBase>3) color=4;
                        else color = fromBase2CS[(int)(lastBase<<2 | thisBase)];
                        // re-use the int to base, because ::set expects a base char (ATCG), not
                        // a color code (0123).  It should only matter on final output.
                        set(header->elementCount++, 
                            BaseAsciiMap::int2base[(int) color]);
                    }
                    lastBase = thisBase;
                }
                else
                {
                    set(header->elementCount++, toupper(fasta[fastaIndex]));
                }
                break;
        }

        //
        // slightly awkward exit handling when we exceed the fixed
        // number of chromosomes
        //
        if (terminateLoad) break;
    }

    //
    // also slightly awkward code to handle the last dangling chromosome...
    // all we should need to do is compute the md5 checksum
    //
    if (whichChromosome==0)
    {
        fastaFile.close();
        throw std::runtime_error("No chromosomes found - aborting!");
    }
    else
    {
        setChromosomeMD5andLength(whichChromosome-1);
    }

    fastaFile.close();

    if (_progressStream) *_progressStream << "\r";

    std::cerr << "FASTA binary cache file '"
              << _umfaFilename
              << "' created."
              << std::endl;

    //
    // leave the umfastaFile open in case caller wants to use it
    //
    return false;
}

int GenomeSequence::getChromosomeCount() const
{
    return header->_chromosomeCount;
}

//return chromosome index: 0, 1, ... 24;
int GenomeSequence::getChromosome(genomeIndex_t position) const
{
    if (position == INVALID_GENOME_INDEX) return INVALID_CHROMOSOME_INDEX;

    if (header->_chromosomeCount == 0)
        return INVALID_CHROMOSOME_INDEX;

    int start = 0;
    int stop = header->_chromosomeCount - 1;

    // eliminate case where position is in the last chromosome, since the loop
    // below falls off the end of the list if it in the last one.

    if (position > header->_chromosomes[stop].start)
        return (stop);

    while (start <= stop)
    {
        int middle = (start + stop) / 2;

        if (position >= header->_chromosomes[middle].start && position < header->_chromosomes[middle + 1].start)
            return middle;

        if (position == header->_chromosomes[middle + 1].start)
            return (middle + 1);

        if (position > header->_chromosomes[middle + 1].start)
            start = middle + 1;

        if (position < header->_chromosomes[middle].start)
            stop = middle - 1;
    }

    return -1;
}

//
// Given a chromosome name and 1-based chromosome index, return the
// genome index (0 based) into sequence for it.
//
// NB: the header->chromosomes array contains zero based genome positions
//
genomeIndex_t GenomeSequence::getGenomePosition(
    const char *chromosomeName,
    unsigned int chromosomeIndex) const
{
    genomeIndex_t i = getGenomePosition(chromosomeName);
    if (i == INVALID_GENOME_INDEX) return INVALID_GENOME_INDEX;
    return i + chromosomeIndex - 1;
}

genomeIndex_t GenomeSequence::getGenomePosition(
    int chromosome,
    unsigned int chromosomeIndex) const
{
    if (chromosome<0 || chromosome >= (int) header->_chromosomeCount) return INVALID_GENOME_INDEX;

    genomeIndex_t i = header->_chromosomes[chromosome].start;
    if (i == INVALID_GENOME_INDEX) return INVALID_GENOME_INDEX;
    return i + chromosomeIndex - 1;
}

//
// return the genome index (0 based) of the start of the named
// chromosome.  If none is found, INVALID_GENOME_INDEX is returned.
//
// XXX may need to speed this up - and smarten it up with some
// modest chromosome name parsing.... e.g. '%d/X/Y' or 'chr%d/chrX/chrY' or
// other schemes.
//
genomeIndex_t GenomeSequence::getGenomePosition(const char *chromosomeName) const
{
    int chromosome = getChromosome(chromosomeName);
    if (chromosome==INVALID_CHROMOSOME_INDEX) return INVALID_GENOME_INDEX;
    return header->_chromosomes[chromosome].start;
}

int GenomeSequence::getChromosome(const char *chromosomeName) const
{
    unsigned int i;
    for (i=0; i<header->_chromosomeCount; i++)
    {
        if (strcmp(header->_chromosomes[i].name, chromosomeName)==0)
        {
            return i;
        }
    }
    return INVALID_CHROMOSOME_INDEX;
}

//
// Given a read, reverse the string and swap the base
// pairs for the reverse strand equivalents.
//
void GenomeSequence::getReverseRead(std::string &read)
{
    std::string newRead;
    if (read.size()) for (int32_t i=(int) read.size() - 1; i>=0; i--)
        {
            newRead.push_back(BasePair(read[i]));
        }
    read = newRead;
}

void GenomeSequence::getReverseRead(String& read)
{
    int i = 0;
    int j = read.Length()-1;
    char temp;
    while (i < j)
    {
        temp = read[j];
        read[j] = read[i];
        read[i] = temp;
    }
}

#define ABS(x) ( (x) > 0 ? (x) : -(x) )
int GenomeSequence::debugPrintReadValidation(
    std::string &read,
    std::string &quality,
    char   direction,
    genomeIndex_t   readLocation,
    int sumQuality,
    int mismatchCount,
    bool recurse
) 
{
    int validateSumQ = 0;
    int validateMismatchCount = 0;
    int rc = 0;
    std::string genomeData;

    for (uint32_t i=0; i<read.size(); i++)
    {
        if (tolower(read[i]) != tolower((*this)[readLocation + i]))
        {
            validateSumQ += quality[i] - '!';
            // XXX no longer valid:
            if (direction=='F' ? i<24 : (i >= (read.size() - 24))) validateMismatchCount++;
            genomeData.push_back(tolower((*this)[readLocation + i]));
        }
        else
        {
            genomeData.push_back(toupper((*this)[readLocation + i]));
        }
    }
    assert(validateSumQ>=0);
    if (validateSumQ != sumQuality && validateMismatchCount == mismatchCount)
    {
        printf("SUMQ: Original Genome: %s  test read: %s : actual sumQ = %d, test sumQ = %d\n",
               genomeData.c_str(),
               read.c_str(),
               validateSumQ,
               sumQuality
              );
        rc++;
    }
    else if (validateSumQ == sumQuality && validateMismatchCount != mismatchCount)
    {
        printf("MISM: Original Genome: %s  test read: %s : actual mismatch %d test mismatches %d\n",
               genomeData.c_str(),
               read.c_str(),
               validateMismatchCount,
               mismatchCount
              );
        rc++;
    }
    else if (validateSumQ != sumQuality && validateMismatchCount != mismatchCount)
    {
        printf("BOTH: Original Genome: %s  test read: %s : actual sumQ = %d, test sumQ = %d, actual mismatch %d test mismatches %d\n",
               genomeData.c_str(),
               read.c_str(),
               validateSumQ,
               sumQuality,
               validateMismatchCount,
               mismatchCount
              );
        rc++;
    }

    if (recurse && ABS(validateMismatchCount - mismatchCount) > (int) read.size()/2)
    {
        printf("large mismatch difference, trying reverse strand: ");
        std::string reverseRead = read;
        std::string reverseQuality = quality;
        getReverseRead(reverseRead);
        reverse(reverseQuality.begin(), reverseQuality.end());
        rc = debugPrintReadValidation(reverseRead, reverseQuality, readLocation, sumQuality, mismatchCount, false);
    }
    return rc;
}
#undef ABS


bool GenomeSequence::wordMatch(unsigned int index, std::string &word) const
{
    for (uint32_t i = 0; i<word.size(); i++)
    {
        if ((*this)[index + i] != word[i]) return false;
    }
    return true;
}

bool GenomeSequence::printNearbyWords(unsigned int index, unsigned int deviation, std::string &word) const
{
    for (unsigned int i = index - deviation; i < index + deviation; i++)
    {
        if (wordMatch(i, word))
        {
            std::cerr << "word '"
                      << word
                      << "' found "
                      << i - index
                      << " away from position "
                      << index
                      << "."
                      << std::endl;
        }
    }
    return false;
}

void GenomeSequence::dumpSequenceSAMDictionary(std::ostream &file) const
{
    for (unsigned int i=0; i<header->_chromosomeCount; i++)
    {
        file
        <<  "@SQ"
        << "\tSN:" << header->_chromosomes[i].name  // name
        << "\tLN:" << header->_chromosomes[i].size  // number of bases
        << "\tAS:" << header->_chromosomes[i].assemblyID // e.g. NCBI36.3
        << "\tM5:" << header->_chromosomes[i].md5
        << "\tUR:" << header->_chromosomes[i].uri
        << "\tSP:" << header->_chromosomes[i].species // e.g. Homo_sapiens
        << std::endl;
    }
}

void GenomeSequence::dumpHeaderTSV(std::ostream &file) const
{
    file << "# Reference: " << _baseFilename << std::endl;
    file << "# SN: sample name - must be unique" << std::endl;
    file << "# AS: assembly name" << std::endl;
    file << "# SP: species" << std::endl;
    file << "# LN: chromosome/contig length" << std::endl;
    file << "# M5: chromosome/contig MD5 checksum" << std::endl;
    file << "# LN and M5 are only printed for informational purposes." << std::endl;
    file << "# Karma will only set those values when creating the index." << std::endl;
    file << "SN" << "\t" << "AS" << "\t" << "SP" << "\t" << "UR" << "\t" << "LN" << "\t" << "M5" << std::endl;
    for (unsigned int i=0; i<header->_chromosomeCount; i++)
    {
        file
        << header->_chromosomes[i].name  // name
        << "\t" << header->_chromosomes[i].assemblyID // e.g. NCBI36.3
        << "\t" << header->_chromosomes[i].uri
        << "\t" << header->_chromosomes[i].species // e.g. Homo_sapiens
        << "\t" << header->_chromosomes[i].size  // number of bases
        << "\t" << header->_chromosomes[i].md5
        << std::endl;
    }
}


void GenomeSequence::getString(std::string &str, int chromosome, uint32_t index, int baseCount) const
{
    //
    // calculate the genome index for the lazy caller...
    //
    genomeIndex_t genomeIndex = header->_chromosomes[chromosome].start + index - 1;

    getString(str, genomeIndex, baseCount);
}

void GenomeSequence::getString(String &str, int chromosome, uint32_t index, int baseCount) const
{
    std::string string;
    this-> getString(string, chromosome, index, baseCount);
    str = string.c_str();
}

void GenomeSequence::getString(std::string &str, genomeIndex_t index, int baseCount) const
{
    str.clear();
    if (baseCount > 0)
    {
        for (int i=0; i<baseCount; i++)
        {
            str.push_back((*this)[index + i]);
        }
    }
    else
    {
        // if caller passed negative basecount, give them
        // the read for the 3' end
        for (int i=0; i< -baseCount; i++)
        {
            str.push_back(BaseAsciiMap::base2complement[(int)(*this)[index + i]]);
        }
    }
}

void GenomeSequence::getString(String &str, genomeIndex_t index, int baseCount) const
{
    std::string string;
    getString(string, index, baseCount);
    str = string.c_str();
}

void GenomeSequence::getHighLightedString(std::string &str, genomeIndex_t index, int baseCount, genomeIndex_t highLightStart, genomeIndex_t highLightEnd) const
{
    str.clear();
    if (baseCount > 0)
    {
        for (int i=0; i<baseCount; i++)
        {
            char base = (*this)[index + i];
            if (in(index+i, highLightStart, highLightEnd))
                base = tolower(base);
            str.push_back(base);
        }
    }
    else
    {
        // if caller passed negative basecount, give them
        // the read for the 3' end
        for (int i=0; i< -baseCount; i++)
        {
            char base = BaseAsciiMap::base2complement[(int)(*this)[index + i]];
            if (in(index+i, highLightStart, highLightEnd))
                base = tolower(base);
            str.push_back(base);
        }
    }
}

void GenomeSequence::print30(genomeIndex_t index) const
{
    std::cout << "index: " << index << "\n";
    for (genomeIndex_t i=index-30; i<index+30; i++)
        std::cout << (*this)[i];
    std::cout << "\n";
    for (genomeIndex_t i=index-30; i<index; i++)
        std::cout << " ";
    std::cout << "^";
    std::cout << std::endl;
}

void GenomeSequence::getMismatchHatString(std::string &result, const std::string &read, genomeIndex_t location) const
{
    result = "";
    for (uint32_t i=0; i < read.size(); i++)
    {
        if (read[i] == (*this)[location+i])
            result.push_back(' ');
        else
            result.push_back('^');
    }
}

void GenomeSequence::getMismatchString(std::string &result, const std::string &read, genomeIndex_t location) const
{
    result = "";
    for (uint32_t i=0; i < read.size(); i++)
    {
        if (read[i] == (*this)[location+i])
            result.push_back(toupper(read[i]));
        else
            result.push_back(tolower(read[i]));
    }
}

genomeIndex_t GenomeSequence::simpleLocalAligner(std::string &read, std::string &quality, genomeIndex_t index, int windowSize) const
{
    int bestScore = 1000000; // either mismatch count or sumQ
    genomeIndex_t bestMatchLocation = INVALID_GENOME_INDEX;
    for (int i=-windowSize; i<windowSize; i++)
    {
        int newScore;

        if (i<0 && ((uint32_t) -i) > index) continue;
        if (index + i + read.size() >= getNumberBases()) continue;

        if (quality=="")
        {
            newScore = this->getMismatchCount(read, index + i);
        }
        else
        {
            newScore = this->getSumQ(read, quality, index + i);
        }
        if (newScore < bestScore)
        {
            bestScore = newScore;
            bestMatchLocation = index + i;
        }
    }
    return bestMatchLocation;
}

std::ostream &operator << (std::ostream &stream, genomeSequenceMmapHeader &h)
{
    stream << (MemoryMapArrayHeader &) h;
    stream << "chromosomeCount: " << h._chromosomeCount << std::endl;
    stream << "isColorSpace: " << h._colorSpace << std::endl;
    stream << "chromosomeCount: " << h._chromosomeCount << std::endl;
    uint64_t totalSize = 0;
    for (uint32_t i=0; i < h._chromosomeCount; i++)
    {
        totalSize += h._chromosomes[i].size;
        stream << "Chromosome Index " << i << " name: " << h._chromosomes[i].name << std::endl;
        stream << "Chromosome Index " << i << " whole genome start: " << h._chromosomes[i].start << std::endl;
        stream << "Chromosome Index " << i << " whole genome size: " << h._chromosomes[i].size << std::endl;
        stream << "Chromosome Index " << i << " md5 checksum: " << h._chromosomes[i].md5 << std::endl;
        stream << "Chromosome Index " << i << " assemblyID: " << h._chromosomes[i].assemblyID << std::endl;
        stream << "Chromosome Index " << i << " species: " << h._chromosomes[i].species << std::endl;
        stream << "Chromosome Index " << i << " URI: " << h._chromosomes[i].uri << std::endl;
    }
    stream << "Total Genome Size: " << totalSize << " bases."<< std::endl;
    if (totalSize != h.elementCount)
    {
        stream << "Total Genome Size: does not match elementCount!\n";
    }

    stream << std::endl;
    return stream;
}

void GenomeSequence::getChromosomeAndIndex(std::string &s, genomeIndex_t i) const
{
    int whichChromosome = 0;

    whichChromosome = getChromosome(i);

    if (whichChromosome == INVALID_CHROMOSOME_INDEX)
    {
        s = "invalid genome index";     // TODO include the index in error
    }
    else
    {
        std::ostringstream buf;
        genomeIndex_t   chromosomeIndex = i - getChromosomeStart(whichChromosome) + 1;
        buf << header->_chromosomes[whichChromosome].name  << ":"  << chromosomeIndex;
#if 0
        buf << " (GenomeIndex " << i << ")";
#endif
        s = buf.str();
    }
    return;
}


void GenomeSequence::getChromosomeAndIndex(String& s, genomeIndex_t i) const
{
    std::string ss;
    getChromosomeAndIndex(ss, i);
    s = ss.c_str();
    return;
}


//
// This is intended to be a helper routine to get dbSNP files
// loaded.  In some cases, we will load into an mmap() file (ie
// when we are creating it), in others, we will simply be loading
// an existing dbSNP file into RAM (when the binary file does not
// exist or when we are running with useMemoryMapFlag == false.
//
// Assume that dbSNP exists, is writable, and is the right size.
//
// Using the dbSNPFilename given, mark each dbSNP position
// with a bool true.
//
// Return value:
//   True: if populateDBSNP() succeed
//   False: if not succeed
bool GenomeSequence::populateDBSNP(
    mmapArrayBool_t &dbSNP,
    IFILE inputFile) const
{
    assert(dbSNP.getElementCount() == getNumberBases());

    if(inputFile == NULL)
    {
        // FAIL, file not opened.
        return(false);
    }

    std::string chromosomeName;
    std::string position;
    genomeIndex_t chromosomePosition1;  // 1-based
    uint64_t    ignoredLineCount = 0;

    // Read til the end of the file.
    char* postPosPtr = NULL;
    while(!inputFile->ifeof())
    {
        chromosomeName.clear();
        position.clear();
        // Read the chromosome
        if(inputFile->readTilTab(chromosomeName) <= 0)
        {
            // hit either eof or end of line, check if
            // it is a header.
            if(chromosomeName.size()>0 && chromosomeName[0]=='#')
            {
                // header, so just continue.
                continue;
            }
            // Not the header, so this line is poorly formatted.
            ++ignoredLineCount;
            // Continue to the next line.
            continue;
        }

        // Check if it is a header line.
        if(chromosomeName.size()>0 && chromosomeName[0]=='#')
        {
            // did not hit eof or end of line, 
            // so discard the rest of the line.
            inputFile->discardLine();
            continue;
        }

        // Not a header, so read the position.
        if(inputFile->readTilTab(position) > 0)
        {
            // Additional data on the line, so discard it.
            inputFile->discardLine();
        }

        // Convert the position to a string.
        chromosomePosition1 = strtoul(position.c_str(), &postPosPtr, 0);


        if(postPosPtr == position.c_str())
        {
            ++ignoredLineCount;
            continue;
        }

        // 1-based genome index.
        genomeIndex_t genomeIndex = 
            getGenomePosition(chromosomeName.c_str(), chromosomePosition1);

        // if the genome index is invalid, ignore it
        if((genomeIndex == INVALID_GENOME_INDEX) || 
           (genomeIndex > getNumberBases()))
        {
            ignoredLineCount++;
            continue;
        }

        dbSNP.set(genomeIndex, true);
    }

    if (ignoredLineCount > 0)
    {
        std::cerr << "GenomeSequence::populateDBSNP: ignored " << ignoredLineCount << " SNP positions due to invalid format of line." << std::endl;
        return false;
    }
    return true;
}

bool GenomeSequence::loadDBSNP(
    mmapArrayBool_t &dbSNP,
    const char *inputFileName) const
{
    //
    // the goal in this section of code is to allow the user
    // to either specify a valid binary version of the SNP file,
    // or the original text file that it gets created from.
    //
    // To do this, we basically open, sniff the error message,
    // and if it claims it is not a binary version of the file,
    // we go ahead and treat it as the text file and use the
    // GenomeSequence::populateDBSNP method to load it.
    //
    // Further checking is really needed to ensure users don't
    // mix a dbSNP file for a different reference, since it is really
    // easy to do.
    //
    if (strlen(inputFileName)!=0)
    {
        std::cerr << "Load dbSNP file '" << inputFileName << "': " << std::flush;

        if (dbSNP.open(inputFileName, O_RDONLY))
        {
            //
            // failed to open, possibly due to bad magic.
            //
            // this is really awful ... need to have a return
            // code that is smart enough to avoid this ugliness:
            //
            if (dbSNP.getErrorString().find("wrong type of file")==std::string::npos)
            {
                std::cerr << "Error: " << dbSNP.getErrorString() << std::endl;
                exit(1);
            }
            //
            // we have a file, assume we can load it as a text file
            //
            IFILE inputFile = ifopen(inputFileName, "r");
            if(inputFile == NULL)
            {
                std::cerr << "Error: failed to open " << inputFileName << std::endl;
                exit(1);
            }

            std::cerr << "(as text file) ";

            // anonymously (RAM resident only) create:
            dbSNP.create(getNumberBases());

            // now load it into RAM
            populateDBSNP(dbSNP, inputFile);
            ifclose(inputFile);

        }
        else
        {
            std::cerr << "(as binary mapped file) ";
        }

        std::cerr << "DONE!" << std::endl;
        return false;
    }
    else
    {
        return true;
    }
}


#if defined(TEST)

void simplestExample(void)
{
    GenomeSequence  reference;
    genomeIndex_t   index;

    // a particular reference is set by:
    // reference.setFastaName("/usr/cluster/share/karma/human_g1k_v37_12CS.fa")
    //
    // In the above example, the suffix .fa is stripped and replaced with .umfa,
    // which contains the actual file being opened.
    //
    if (reference.open())
    {
        perror("GenomeSequence::open");
        exit(1);
    }


    index = 1000000000; // 10^9
    //
    // Write the base at the given index.  Here, index is 0 based,
    // and is across the whole genome, as all chromosomes are sequentially
    // concatenated, so the allowed range is
    //
    // 0.. (reference.getChromosomeStart(last) + reference.getChromosomeSize(last))
    //
    // (where int last = reference.getChromosomeCount() - 1;)
    //
    std::cout << "base[" << index << "] = " << reference[index] << std::endl;

    //
    // Example for finding chromosome and one based chromosome position given
    // and absolute position on the genome in 'index':
    //
    int chr = reference.getChromosome(index);
    genomeIndex_t   chrIndex = index - reference.getChromosomeStart(chr) + 1;   // 1-based

    std::cout << "genome index " << index << " corresponds to chromosome " << chr << " position " << chrIndex << std::endl;

    //
    // Example for finding an absolute genome index position when the
    // chromosome name and one based position are known:
    //
    const char *chromosomeName = "5";
    chr = reference.getChromosome(chromosomeName);     // 0-based
    chrIndex = 100000;                      // 1-based

    index = reference.getChromosomeStart(chr) + chrIndex - 1;

    std::cout << "Chromosome '" << chromosomeName << "' position " << chrIndex << " corresponds to genome index position " << index << std::endl;

    reference.close();
}

void testGenomeSequence(void)
{
    GenomeSequence reference;

#if 0
    std::string referenceName = "someotherreference";
    if (reference.setFastaName(referenceName))
    {
        std::cerr << "failed to open reference file "
                  << referenceName
                  << std::endl;
        exit(1);
    }
#endif

    std::cerr << "open and prefetch the reference genome: ";

    // open it
    if (reference.open())
    {
        exit(1);
    }
    std::cerr << "done!" << std::endl;

    //
    // For the human genome, genomeIndex ranges from 0 to 3.2x10^9
    //
    genomeIndex_t   genomeIndex;    // 0 based
    unsigned int chromosomeIndex;   // 1 based
    unsigned int chromosome;        // 0..23 or so
    std::string chromosomeName;

    //
    // Here we'll start with a chromosome name, then obtain the genome
    // index, and use it to find the base we want:
    //
    chromosomeName = "2";
    chromosomeIndex = 1234567;
    // this call is slow (string search for chromsomeName):
    genomeIndex = reference.getGenomePosition(chromosomeName.c_str(), chromosomeIndex);
    assert(genomeIndex!=INVALID_GENOME_INDEX);
    std::cout << "Chromosome " << chromosomeName << ", index ";
    std::cout << chromosomeIndex << " contains base " << reference[genomeIndex];
    std::cout << " at genome index position " << genomeIndex << std::endl;

    //
    // now reverse it - given a genomeIndex from above, find the chromosome
    // name and index:
    //

    // slow (binary search on genomeIndex):
    chromosome = reference.getChromosome(genomeIndex);
    unsigned int newChromosomeIndex;
    // not slow:
    newChromosomeIndex = genomeIndex - reference.getChromosomeStart(chromosome) + 1;

    assert(chromosomeIndex == newChromosomeIndex);

    // more testing... at least test and use PackedRead:
    //
    PackedRead pr;

    pr.set("ATCGATCG", 0);
    assert(pr.size()==8);
    assert(pr[0]==BaseAsciiMap::base2int[(int) 'A']);
    assert(pr[1]==BaseAsciiMap::base2int[(int) 'T']);
    assert(pr[2]==BaseAsciiMap::base2int[(int) 'C']);
    assert(pr[3]==BaseAsciiMap::base2int[(int) 'G']);
    pr.set("ATCGATCG", 1);
    assert(pr.size()==9);
    pr.set("", 0);
    assert(pr.size()==0);
    pr.set("", 1);
    assert(pr.size()==1);
    pr.set("", 2);
    assert(pr.size()==2);
    pr.set("", 3);
    assert(pr.size()==3);
    assert(pr[0]==BaseAsciiMap::base2int[(int) 'N']);
    assert(pr[1]==BaseAsciiMap::base2int[(int) 'N']);
    assert(pr[2]==BaseAsciiMap::base2int[(int) 'N']);
    pr.set("C", 1);
    assert(pr.size()==2);
    assert(pr[0]==BaseAsciiMap::base2int[(int) 'N']);
    assert(pr[1]==BaseAsciiMap::base2int[(int) 'C']);

}

//
// After I build libcsg, I compile and run this test code using:
//
// g++ -DTEST -o try GenomeSequence.cpp -L. -lcsg -lm -lz -lssl
// you also may need -fno-rtti
// ./try
//
int main(int argc, const char **argv)
{
#if 1
    simplestExample();
#else
    testGenomeSequence();
#endif
    exit(0);
}

#endif
