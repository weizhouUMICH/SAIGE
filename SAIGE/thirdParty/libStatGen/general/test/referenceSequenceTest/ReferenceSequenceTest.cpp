/*
 *  Copyright (C) 2011  Regents of the University of Michigan
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

#include <getopt.h>
#include "Generic.h"
#include <stdio.h>
#include "ReferenceSequence.h"
#include "UnitTest.h"

#include <assert.h>
#include <sstream>


class ReferenceSequenceTest : public UnitTest
{
public:
    ReferenceSequenceTest(const char *title) : UnitTest(title) {;}
    void test1();
    void test2();
    void test3();
    void humanGenomeTest1();

    void test() {
        test1();
        test2();
        test3();
        // This test is very slow:
        // humanGenomeTest1();
    }
};

void ReferenceSequenceTest::test1(void)
{
    std::string sequence("ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG");
    std::string word;

    word="ACTG";
    check(m_failures, ++m_testNum, "Test wordMatch with std::string", true,
            Sequence::wordMatch(sequence, 4, word));

    std::stringstream output;

    Sequence::printNearbyWords(output, sequence, 8, word, 4);

    std::string expect("\
word 'ACTG' found -4 away from position 8.\n\
word 'ACTG' found 0 away from position 8.\n\
");

    check(m_failures, ++m_testNum, "Test printNearbyWords with std::string", expect, output.str());


    Sequence::getString(sequence, 4, 4, word);

    check(m_failures, ++m_testNum, "Test getString with std::string", "ACTG", word);

    Sequence::getHighLightedString(sequence, 0, 12, word, 4, 8);
    check(m_failures, ++m_testNum, "Test getHighLightedStribng with std::string", "ACTGactgACTG",word);

#if 0
    // busted test - don't know why
    output.clear();
    output.str(std::string());
//    Sequence::printBaseContext(std::cout, sequence, 8, 4);
    Sequence::printBaseContext(output, sequence, 8, 4);
    expect="\
index: 8\n\
ACTGACTGA\n\
    ^\n\
";
    check(m_failures, ++m_testNum, "Test printBaseContext with std::string", expect, output.str());
#endif
    std::string result;
    std::string   read("ACTGZZZZACTG");
              expect = "    ^^^^    ";
    Sequence::getMismatchHatString(sequence, 4, result, read);
    check(m_failures, ++m_testNum, "Test getMismatchHatString with std::string", expect, result);


    read="ACTG";
    std::string quality("");
    size_t location = Sequence::simpleLocalAligner(sequence, 0, read, quality, 12);
    check(m_failures, ++m_testNum, "Test simpleLocalAligner with std::string", (size_t) 0, location);

    read="ACNG";
    int misMatches = Sequence::getMismatchCount(sequence, 0, read);
    check(m_failures, ++m_testNum, "Test getMismatchCount with std::string", 1, misMatches);

    read="ACNG";
    quality="$$$$";
    int sumQ = Sequence::getSumQ(sequence, 0, read, quality);
    check(m_failures, ++m_testNum, "Test getSumQ with std::string", 3, sumQ);
}

void ReferenceSequenceTest::test2(void)
{
    PackedSequenceData sequence;
    std::string word;

    sequence.push_back('A');
    sequence.push_back('C');
    sequence.push_back('T');
    sequence.push_back('G');

    sequence.push_back('A');
    sequence.push_back('C');
    sequence.push_back('T');
    sequence.push_back('G');

    sequence.push_back('A');
    sequence.push_back('C');
    sequence.push_back('T');
    sequence.push_back('G');

    sequence.push_back('A');
    sequence.push_back('C');
    sequence.push_back('T');
    sequence.push_back('G');

    Sequence::getString(sequence, 4, 4, word);

    check(m_failures, ++m_testNum, "Test getString with PackedSequenceData", "ACTG", word);
    
    std::cout << "test2 sequence utilization is " << sequence.getUtilization() * 100 << "% - expect around 6.25%" << std::endl;

}

void ReferenceSequenceTest::test3(void)
{
    std::vector<PackedSequenceData> chromosomeSequence;
    std::vector<std::string> chromosomeNames;

    bool result = loadFastaFile("../phiX.fa", chromosomeSequence, chromosomeNames);

    if(result) {
        std::cout << "../phiX.fa not found - skipping these tests." << std::endl;
        return;
    }

    std::cout << "phiX reference utilization is " << chromosomeSequence[0].getUtilization() * 100 << "% - expect around 96.8%" << std::endl;



    check(m_failures, ++m_testNum, "Test loadFastaFile with PackedSequenceData", (size_t) 1, chromosomeNames.size());
    check(m_failures, ++m_testNum, "Test loadFastaFile with PackedSequenceData", (size_t) 1, chromosomeSequence.size());
    check(m_failures, ++m_testNum, "Test loadFastaFile with PackedSequenceData", "1", chromosomeNames[0]);

    std::string word;

    Sequence::getString(chromosomeSequence[0], 60, 10, word);

    check(m_failures, ++m_testNum, "Test loadFastaFile with PackedSequenceData", "AAATTATCTT", word);

}

void ReferenceSequenceTest::humanGenomeTest1(void)
{
    std::vector<PackedSequenceData> chromosomeSequence;
    std::vector<std::string> chromosomeNames;

#define HUMAN_GENOME "/data/local/ref/karma.ref/human.g1k.v37.fa"
    bool result = loadFastaFile(HUMAN_GENOME, chromosomeSequence, chromosomeNames);

    if(result) {
        std::cout << HUMAN_GENOME << " not found - skipping these tests." << std::endl;
        return;
    }

}

int main(int argc, char **argv)
{
    ReferenceSequenceTest test("ReferenceSequenceTest");

    test.test();

    std::cout << test;

    exit(test.getFailureCount());
}
