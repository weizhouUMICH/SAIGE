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

#include <gtest/gtest.h>

#include "GenomeSequence.h"

const char* RM_BS_REFERENCE = "rm -f ./phiX-bs.umfa";
const char* RM_CS_REFERENCE = "rm -f ./phiX-cs.umfa";
const char* REFERENCE_NAME = "./phiX.fa";

TEST(GenomeSequenceTest, staticLookupTest)
{
    GenomeSequence s;
    // quick sanity check...
    EXPECT_EQ(GenomeSequence::int2base[GenomeSequence::base2int[(int) 'A']], 'A');
    EXPECT_EQ(GenomeSequence::int2base[GenomeSequence::base2int[(int) 'a']], 'A');
    EXPECT_EQ(GenomeSequence::int2base[GenomeSequence::base2int[(int) 'T']], 'T');
    EXPECT_EQ(GenomeSequence::int2base[GenomeSequence::base2int[(int) 't']], 'T');
    EXPECT_EQ(GenomeSequence::int2base[GenomeSequence::base2int[(int) 'C']], 'C');
    EXPECT_EQ(GenomeSequence::int2base[GenomeSequence::base2int[(int) 'c']], 'C');
    EXPECT_EQ(GenomeSequence::int2base[GenomeSequence::base2int[(int) 'G']], 'G');
    EXPECT_EQ(GenomeSequence::int2base[GenomeSequence::base2int[(int) 'g']], 'G');
    EXPECT_EQ(GenomeSequence::int2base[GenomeSequence::base2int[(int) 'N']], 'N');
    EXPECT_EQ(GenomeSequence::int2base[GenomeSequence::base2int[(int) 'n']], 'N');
    EXPECT_EQ(GenomeSequence::int2base[GenomeSequence::base2int[(int) 'M']], 'M');
    EXPECT_EQ(GenomeSequence::int2base[GenomeSequence::base2int[(int) 'm']], 'M');

    EXPECT_EQ(GenomeSequence::base2int[(int) 'N'], 4);
    EXPECT_EQ(GenomeSequence::base2int[(int) 'n'], 4);
    EXPECT_EQ(GenomeSequence::base2int[(int) 'A'], 0);
    EXPECT_EQ(GenomeSequence::base2int[(int) 'a'], 0);
    EXPECT_EQ(GenomeSequence::base2int[(int) 'T'], 3);
    EXPECT_EQ(GenomeSequence::base2int[(int) 't'], 3);
    EXPECT_EQ(GenomeSequence::base2int[(int) 'C'], 1);
    EXPECT_EQ(GenomeSequence::base2int[(int) 'c'], 1);
    EXPECT_EQ(GenomeSequence::base2int[(int) 'G'], 2);
    EXPECT_EQ(GenomeSequence::base2int[(int) 'g'], 2);
}


TEST(GenomeSequenceTest, testBaseSpaceReference)
{
    GenomeSequence s;
    int exitCode = system(RM_BS_REFERENCE);
    EXPECT_EQ(exitCode, 0);

    s.setReferenceName(REFERENCE_NAME);
    bool rc = s.create(false);
    EXPECT_EQ(rc, false);
    EXPECT_EQ(s[0], 'G');
    EXPECT_EQ(s[1], 'A');
    EXPECT_EQ(s[2], 'G');
    EXPECT_EQ(s[s.getNumberBases()-3], 'G');
    EXPECT_EQ(s[s.getNumberBases()-2], 'C');
    EXPECT_EQ(s[s.getNumberBases()-1], 'A');
    EXPECT_EQ(s[s.getNumberBases()], 'N');  // check bounds checker

    s.close();
}

TEST(GenomeSequenceTest, testColorSpaceReference)
{
    GenomeSequence s;
    int exitCode = system(RM_CS_REFERENCE);
    EXPECT_EQ(exitCode, 0);

    s.setReferenceName(REFERENCE_NAME);
    bool rc = s.create(true);

    // NB: I did not calculate these expected values, I just
    // read them from the converted genome and set them here.
    // So in theory, they should be checked by hand to ensure
    // that they are correct.
    EXPECT_EQ(rc, false);
    EXPECT_EQ(s[0], 'N');   // in color space, first symbol is unknown
    EXPECT_EQ(s[1], '2');
    EXPECT_EQ(s[2], '2');
    EXPECT_EQ(s[s.getNumberBases()-3], '1');
    EXPECT_EQ(s[s.getNumberBases()-2], '3');
    EXPECT_EQ(s[s.getNumberBases()-1], '1');
    EXPECT_EQ(s[s.getNumberBases()], 'N');  // check bounds checker

    s.close();
}

#if 0
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
    if (reference.open()) {
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
    if (reference.setFastaName(referenceName)) {
        std::cerr << "failed to open reference file "
                  << referenceName
                  << std::endl;
        exit(1);
    }
#endif

    std::cerr << "open and prefetch the reference genome: ";

    // open it
    if (reference.open()) {
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
    assert(pr[0]==GenomeSequence::base2int[(int) 'A']);
    assert(pr[1]==GenomeSequence::base2int[(int) 'T']);
    assert(pr[2]==GenomeSequence::base2int[(int) 'C']);
    assert(pr[3]==GenomeSequence::base2int[(int) 'G']);
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
    assert(pr[0]==GenomeSequence::base2int[(int) 'N']);
    assert(pr[1]==GenomeSequence::base2int[(int) 'N']);
    assert(pr[2]==GenomeSequence::base2int[(int) 'N']);
    pr.set("C", 1);
    assert(pr.size()==2);
    assert(pr[0]==GenomeSequence::base2int[(int) 'N']);
    assert(pr[1]==GenomeSequence::base2int[(int) 'C']);

}

#endif
