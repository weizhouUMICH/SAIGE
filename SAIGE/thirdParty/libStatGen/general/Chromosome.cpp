#include <cassert>
#include "Chromosome.h"

Chromosome::Chromosome(GenomeSequence* gs, unsigned int chromosomeIndex)
{
    assert(gs);
    assert(chromosomeIndex < (unsigned int)gs->getChromosomeCount());

    this->gs = gs;
    this->chromosomeIndex = chromosomeIndex;
    this->offset = gs->getChromosomeStart((int)chromosomeIndex);
    this->chromosomeSize = gs->getChromosomeSize((int)chromosomeIndex);
}

Chromosome::Chromosome(GenomeSequence* gs, const char* chromosomeName)
{
    assert(gs);
    this->gs = gs;

    this->chromosomeIndex = gs->getChromosome(chromosomeName);
    assert(chromosomeIndex != INVALID_CHROMOSOME_INDEX);

    this->offset = gs->getChromosomeStart((int)chromosomeIndex);
    this->chromosomeSize = gs->getChromosomeSize((int)chromosomeIndex);
}

Chromosome::Chromosome(const char* genomseSequenceFileName, unsigned int chromosomeIndex, bool isColorSpace) 
{
    std::string s(genomseSequenceFileName);
    if (this->gs) delete gs;
    gs = new GenomeSequence;
    assert(gs);
    gs->setReferenceName(s);
    assert(!gs->open(isColorSpace));
    this->chromosomeIndex = chromosomeIndex;
    this->offset = gs->getChromosomeStart((int)chromosomeIndex);
    this->chromosomeSize = gs->getChromosomeSize((int)chromosomeIndex);
}

Chromosome::Chromosome(const std::string& genomseSequenceFileName, unsigned int chromosomeIndex, bool isColorSpace) 
{
    if (this->gs) delete gs;
    gs = new GenomeSequence;
    assert(gs);
    gs->setReferenceName(genomseSequenceFileName);
    assert(!gs->open(isColorSpace));
    this->chromosomeIndex = chromosomeIndex;
    this->offset = gs->getChromosomeStart((int)chromosomeIndex);
    this->chromosomeSize = gs->getChromosomeSize((int)chromosomeIndex);
}
