#ifndef _CHROMOSOME_H_
#define _CHROMOSOME_H_

#include "GenomeSequence.h"

class Chromosome{
public:
    explicit Chromosome(GenomeSequence* gs, unsigned int chrosomeIndex);
    explicit Chromosome(GenomeSequence* gs, const char* chromosomeName);
    explicit Chromosome(const char* genomseSequenceFileName, unsigned int chromosomeIndex, bool isColorSpace);
    explicit Chromosome(const std::string& genomseSequenceFileName, unsigned int chromosomeIndex, bool isColorSpace);
    genomeIndex_t Length() const
    {
        return chromosomeSize;
    }
    // 0-based index 
    inline char operator[](genomeIndex_t index) const
    {
        index += offset;
        return (*gs)[index];
    }
    const char* Name() const {
        return gs->getChromosomeName(this->chromosomeIndex);
    }
private:
    GenomeSequence* gs;
    int chromosomeIndex;
    genomeIndex_t offset;           // chromosome index 0 corresponds (*gs)[offset]
    genomeIndex_t chromosomeSize;   // return the length of the chromosome
};

#endif /* _CHROMOSOME_H_ */
