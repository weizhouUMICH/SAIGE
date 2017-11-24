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
#include "GreedyTupleAligner.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <gtest/gtest.h>

// design a mock class for GenomeSequence
class MockGenomeSequence{
public:
    MockGenomeSequence(std::string sequence) :
        sequence(sequence) {};
    char operator[] (int index) const {
        if (index < 0 || index >= (int)sequence.size()) {
            std::cerr << "exceeds boundary! at " << __FILE__ << ":" <<__LINE__ << std::endl;
            return 'N';
        }
        return sequence[index];
    }
    char& operator[] (int index) {
        if (index < 0 || index >= (int)sequence.size()) {
            std::cerr << "exceeds boundary! at " << __FILE__ << ":" <<__LINE__ << std::endl;
            return sequence[0]; 
        }
        return sequence[index];
    }

    int getNumberBases() const{
        return sequence.size();
    }
private:
    std::string sequence;
};

void printRefQueryCigar(const std::vector<char>& prettyPrintReference, 
                        const std::vector<char>& prettyPrintQuery,
                        CigarRoller& cs,
                        std::ostream& out){
    out << " ref  = ";
    for(std::vector<char>::const_iterator i=prettyPrintReference.begin(); i<prettyPrintReference.end(); i++) out << *i;
    out << std::endl;
    out <<"query = ";
    for(std::vector<char>::const_iterator i=prettyPrintQuery.begin(); i<prettyPrintQuery.end(); i++) out << *i;
    out << std::endl;
    out << "cigar = ";
    out << cs.getString();
    out << " ";
}

void runAlign(const char *query, const char *ref, CigarRoller& cs, int& matchPosition) {
    Weight wt;
    MockGenomeSequence reference (ref);
    GreedyTupleAligner<const char *, MockGenomeSequence, int> ga(wt);
    ga.Align(query, strlen(query), reference, 0, strlen(ref), cs, matchPosition);
}

void computePrettyOutput (std::vector<char>& prettyPrintReference,
                          const char* ref, 
                          std::vector<char>& prettyPrintQuery,
                          const char* query, 
                          const int matchPosition, 
                          const CigarRoller& expectedCigar) 
{
    const char* pRef = ref;
    const char* pQuery = query;
    for (int index = 0; index < matchPosition; index++){
        prettyPrintReference.push_back(*pRef++);
        prettyPrintQuery.push_back(' ');
    }
    for (int i = 0; i< expectedCigar.size(); i++) {
        switch( expectedCigar[i].operation) {
        case CigarRoller::mismatch:
        case CigarRoller::match:
            for (unsigned int j = 0; j < expectedCigar[i].count; j++){
                prettyPrintReference.push_back(*pRef++);
                    prettyPrintQuery.push_back(*pQuery++);
            }
            break;
            case CigarRoller::del:
                for (unsigned int j = 0; j < expectedCigar[i].count; j++){
                    prettyPrintReference.push_back(*pRef++);
                    prettyPrintQuery.push_back(' ');
                }
                break;
            case CigarRoller::insert:
                for (unsigned int j = 0; j < expectedCigar[i].count; j++){
                    prettyPrintReference.push_back(' ');
                    prettyPrintQuery.push_back(*pQuery++);
                }
                break;
            default:
                break;
            }
        } /// end for
        while (*pRef !='\0')
            prettyPrintReference.push_back(*pRef++);
        while (*pQuery !='\0')
            prettyPrintReference.push_back(*pQuery++);
}

bool verifyAlign(const char *query, const char *ref, const char *expectedCigarString, std::ostream& out = std::cout) {
    out.seekp(std::ios_base::beg);

    CigarRoller cs;
    int matchPosition;
    CigarRoller expectedCigar(expectedCigarString);

    runAlign(query, ref, cs, matchPosition);


    if (matchPosition < 0) {
        fprintf(stderr, "No match in %s, %d \n", __FILE__, __LINE__);
        return false;
    } 

    std::vector<char> prettyPrintReference,prettyPrintQuery;
    computePrettyOutput( prettyPrintReference, ref, 
                         prettyPrintQuery, query, 
                         matchPosition, 
                         expectedCigar); 

    const char *str = cs.getString();

    if((unsigned int)expectedCigar.getExpectedQueryBaseCount() != strlen(query)) {
        printRefQueryCigar(prettyPrintReference, prettyPrintQuery, cs, out);
        out << std::endl;
        out << "Expected Cigar string length " << expectedCigar.getExpectedQueryBaseCount() << " does not match the length of the query " << strlen(query) << ".  Please fix test case." << std::endl;
        return false;
    }
    if((unsigned int)cs.getExpectedQueryBaseCount() != strlen(query)) {
        printRefQueryCigar(prettyPrintReference, prettyPrintQuery, cs, out);
        out << std::endl;
        out << "Query Length of " << strlen(query) << " does not match computed cigar string length of " << cs.getExpectedQueryBaseCount() << std::endl;
        return false;
    }
    if (strcmp(expectedCigarString, str) == 0) {
        // printf("[Correct Answer = %s] \n", expectedCigarString) ;
        return true;
    } else {
        printRefQueryCigar(prettyPrintReference, prettyPrintQuery, cs, out);
        out << "[Correct Answer = " << expectedCigarString << "] --------------------- Wrong!" << std::endl;
        return false;
    }
    return true;
}

TEST(GreedyTupleAlignerTest, AlignToShortReference) {
    std::stringstream ss(std::stringstream::out);
#if 0
    //exact align
    EXPECT_TRUE( verifyAlign("12345", "123456789", "5M") ) << ss.str().c_str();
    EXPECT_TRUE( verifyAlign("23456", "123456789", "5M") ) << ss.str().c_str();
    //mismatch
    EXPECT_TRUE( verifyAlign("123B567", "123456789", "7M") ) << ss.str().c_str();
    EXPECT_TRUE( verifyAlign("234D678", "123456789", "7M") ) << ss.str().c_str();
    // del
    EXPECT_TRUE( verifyAlign("123467890","1234567890", "4M1D5M") ) << ss.str().c_str();
    EXPECT_TRUE( verifyAlign("123467890","B1234567890B", "4M1D5M") ) << ss.str().c_str();
    // ins
    EXPECT_TRUE( verifyAlign("12345067890","1234567890", "5M1I5M") ) << ss.str().c_str();
    EXPECT_TRUE( verifyAlign("12345067890","BBBB1234567890BBBB", "5M1I5M") ) << ss.str().c_str();
    // soft clip
    EXPECT_TRUE( verifyAlign("1234", "1235", "3M1S") ) << ss.str().c_str();
    // The following will treat as two mismatches
    // EXPECT_TRUE( verifyAlign("123456700", "123456789", "7M2S", ss) ) << ss.str().c_str();
#endif    
    EXPECT_TRUE( verifyAlign("1023456700", "123456789","1I7M2S") ) << ss.str().c_str();
}

TEST(GreedyTupleTestAligner, AlignToLongReference) {
    std::stringstream ss(std::stringstream::out);

    EXPECT_TRUE( verifyAlign("TTAGAATGCTATTGTGTTTGGAGATTTGAGGAAAGTGGGCGTGAAGACTTAGTGTTCATTTCCTCAACCTCTCTCTGTGTGAACATACGTCATCGGTCAGAAATTGGG","CCGAGATTGTGCCATTGCACTCCTGCCTGGGTAACAGAGTCAGACCCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAGATTAGGTTTTATAGATGGAAAATTCACAGCTCTCTCCAGATCAGAAATCTCCAAGAGTAAATTAGTGTCTTAAAGGGGTTGTAATAACTTTCCTATGTGACTAAGTGCATTATTAATCAATTTTTCTATGATCAAGTACTCCTTTACATACCTGCTAATACAATTTTTGATATGAAATCAGTCCTAGAGGGAATCAATGTAAGATACAGACTTGATGAGTGCTTGCAGTTTTTTATTGACAATCTGAAGAATGACTTGACTCTAAATTGCAGCTCAAGGCTTAGAATGCTATTGTGTTTGGAGATTTGAGGAAAGTGGGCGTGAAGACTTAGTGTTCATTTCCTCAACCTCTCTCTGTGTGAACATACAGGAATCAAATCTGTCTAGCCTCTCTTTTTGGCAAGGTTAAGAACAATTCCACTTCATCCTAATCCCAATGATTCCTGCCGACCCTCTTCCAAAAACTATTTAAAGACATGTTCTTCAAAGTTATATTTGTCTTTCCTTCAGGGAGAAAAAGAATACCAATCACTTATAATATGGAAACTAGCAGAAATGGGTCACATAAGTCATCTGTCAGAAATTGGGAAAATAGAGTAGGTCAGTCTTTCCAGTCATGGTACTTTTACCTTCAATCA", "88M200D20M") ) << ss.str().c_str();

    // This reads cigar string is the wrong length.  By that I mean that the number
    // of matches listed in the cigar should be the same as the length of
    // the original read.
    //
    EXPECT_TRUE( verifyAlign("GTGAAACTCCATCTCAAAAATAAGTAAATAAATAAATACATACATAGGCACAGTGCAGTTGTTAGTCAGAATTAGGTCACACTGGATTAGGGTGAGTACTTAATGCAACAGGTCTGGGG","GTGCCAGAGTTTAATTAATAGGATAAGGTTATGAGTCAGACTGTGTACCCCAAAAAAGATATGTTGAACTCCTAAGCCCCTGAACCACAGAATGGGATCCTATTCAGAAATAGGCACAGTGTCCGGGCACCATGGCTCACACTGGTAATCCCAGCACTCTGGGAGGCTGAGGTGGGTGCATCACCTGAGGTCAGGAGTTTGAGACCAGCCTGGCCAACATGGTGAAACCCCATCTCTACTAAAAATACAAACAGAACAGTTAGCCAGGTGTGGTGGTGGGCACCTGTAATCCCAGCTACTTGGGAGGCTGAGACAGGAGAATGGCTTGAACCCAGGAGGTGGAGGTTGCAGTGAGCCGAGATCGTGCCATTGCACTTCAGCCTGGGCCACAAGAGTGAAACTCCATCTCAAAAATAAGTAAATAAATAAATACATACGTAGGCACAGTGCAGTTGTTGTTAGTTAGAATTAGGTCACACTGGATTAGGGTGAGTCCTTAATCCAACAGGTCTGGTGTCCTTACAAATAGACAAATACACAGAAGGAACATGGCCACATGGAGATACAGACACACCAAAACATCATATTGAGATGTGGGCAAAGATTGGAGAGACACTTCTCCAAGTCAAGGAACATCTGGGACTACCCAGAAACTGTAAGAGGCAGAGAAAGGTCCTTCCCTGTAGGCTTTAGAGGAACATGGCCCTGCCAACATCTTGATCTTGGATTTCCAGCCTCCAGCATGTGAGACAAGTTTCTGGGTTTTTTTGGAGACAGAGTCTCACTCTTGTCACCCAGGCTGGAGTGCAGTGGCATGAACTTGGCTCACTGCAACCTCCTCCCAGGATCAAGGTATTGTCCTGCCTCAGCCTCCCGAGTAGCTGGGATGACAGGGGCCCGCCACCACGCCAGCTCATTTTTGTATTTTTTACTAGAGAAGGGGTTTCACCATGTTGGCCAGGCTGGTCTTGAACTCCTGACCTCAAGTGATCCACCCGCCTTGGCCTCCCAAAGTGCTAGGATTACAGGTGTGAGCCACTGCGCCTGGCAAGTTTCTGTTGCCTTAAGCCACTCTTTCTGTGGTAATTTGTTATCATGGCCCTAAGAAATGACTAGAGAGAGAAAGCAAATCCCTTTGTTTCTGCATTTACTGAAACAGATGAATAGATTTCTAGCTCCCTTGGGGTCTGAACTTTTAAAAGAGAGATTTCTTATACATATGATAATCATGATATTGT", "63M3D56M") ) << ss.str().c_str();

    EXPECT_TRUE( verifyAlign("ATATTGTTTTTTTCAATGCATATCAAAACAATGTTTACAATATACTACAGCCTAAGTGTGCAATAGCATTATGTGTAGAAATGCACATACCATAATTAGTTTTTTTTTTTGAAAAAACT","GTTCCAAAAGATTATATTTGTTAGGTTAGAGAATTTTAACTTATTTATATAATGGAGATTTTCTAATACTGAGAATACCTTAATTCTTATTGTAAGCCTACTTAACAGTGACAAAATGTTATTATAACGTGGTATTGAAATTAATATGATAGTATTTTATATGGATATTTGCATATGCAATTGACATATATTGTATATACAATATATAACTGTGTATTATATATTATATTTATATAATGTTATATTGTATATGAATATATTTGAATTATATGTATATACATATATATAGGCATTCATCAGAAATATTGCAGGTTTGGTTTCAGACGACTATAATAAAGTGAATATTGCAATAAAGCGAATCACAAGAAATTATTGTTTTTTTCAATGCATATCAAAACAATGTTTACAATATACTACAGCCTAAGTGTGCAATAGCATTATGTGTAGAAATGCACATACCATAATTAGTTTTTTTTTTGAAAAAACTGTTAATGATTATCTGAGCCTTCAGTGAGTTGTAATCTTTTCATGGTGGAGGATCATACCTCTACGTTGATGTCAGCTGACTGATCAGGGTAGTAGTTGCTGAAGGCTTGGGTGGCTGTGGCAATTTCTTAAAATAGGATAACAATGGCATTTACCACATTAATTGACTCCTTCTTTCACAAAAGATTTCTCTGTCTCATGCAATGCTGTTTGACAGCATTTTCCCCACAGTAGAATTTCTTTAAAAATTGG", "109M249D10M") ) << ss.str().c_str();

    EXPECT_TRUE( verifyAlign("CCAGACTATCTCAAGCAATCAACAGATTTAATGTAAGGAGTGTCAAAATCTGAATGATGCTTTTTGCAGAAATAGAAAATCCCTTTCTAATATTTTTATATTTTTGAC","TTATCGAGGCTGGCGGATTTTGTGAGGCCAGGAGTTCAAGACCAGCCTGGCCAACATGGCAAAACTCTGTCTCTACTAAAAATACAAAAGTTAGCTGGGCATAGTGGCACATGTCTATAGTTCTAGCTATGTGGGAGGCTGAGACACGAGAATCGCTTGAACCCAGGAGGTGGAGGTTGTGGTGAGCCGAGCTCACGCCATTGCTCTCCAGCCTGGGCAACAGAGCAAGACTGTCTCAAAAACAAAAACAAAAAACACAAAAACTACAAGACTTTTATGAAATAACTTAAGGAAGATATAAATAAATGGAAAGATATCCCATGTTCCTGACTTGGAAGACTTAATTTTGTTAAGATGTCCATACTATCTCAAGCAATCAACAGATTTAATGTAAGGAGTGTCAAAATCTGAATGATGCTTTTTGCAGAAATAGAAAATCCCTTTCTAATATTTTTATGTAATCTCAAGGGACCCCAAATAGCCAAAAGAATCCTGAAAAAGTAGAATAAAGCTGGAGGACTCATGATTCCTGATTTGAAAACTTACTACCAGATACAATAATCAAAACAGTTCCGTGCTTGTCATAAAGACAAACATATAGACCAATGGAACAGAATAGAGATTACAGGGACAAATCCTCATATATATGGTCAAATGATTTTTGACCAGTGCCAAGATCATTCATGGGTGAAAAGACAATCTTTTCAATAAAAGAT", "99M200D9M") ) << ss.str().c_str();

    // this is by no mean works, since our local aligner generates the output
    // 54M200D34M47D8M7D2M1I5M4S , since local aligner prefer gapped alignment.
    EXPECT_FALSE( verifyAlign("ATGAGGTCAGGAGATGGAGACCATCCTGGCTAACATGGTGAAACCCCATCTCTAAAAAAAGTGTAACAGAGGTGCATACTCAAAACTACAAAAGTCTCGTGAAAGGAA","CAAGAAAAAGAAATAAAATACATTTTAGTAGGAAAGGAAGAAGTTAAATTGTCTCCATTTGGTGACAACATGAGCTTATATGCAGAAAACCTAAAGACTCTACCAAAAAAACTGCTGGAACTGATAAATGAATTCGGTGAAGTCCTAGGGTATAAAATCAATGTACAAAATAAGTGGTGTTTCTATATTCTAATAAATTATTCAAAAGGGAAATTAAGAAATCAATCCCATTTTCAATAGCAACAACAAAAAAAATGACAATGCCAAAGTATAAATTTAACCAAGAAGCTACAAGAGTCTGGGCGCAGTGGCTTATGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGTGGATCATGAGGTCAGGAGATGGAGACCATCCTGGCTAACATGGTGAAACCCCATCTCTACTAAAAAAAAATAATAATAATAATAATAATAATAATAATTAGCCGGGCGTGGTGGTAGGCATCTGTAGTCCCAGCTACTCGGGAGGCTGAGGTAGGAGAATGGCATGAACCTGGGAGGTGGAGCTTGCAGTGAGCAGAGATCACACCCCTGCACTCCAGCCTGGGTGACAGGGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAGTGTAACAGAGTTGCATACTGAAAACTATAAAATTCTGATGAAAGAAAATGAAGCAACAAATAATTAATATAATAAAAAAGTCCATACTATCCAAAAT", "54M200D54M") ) << ss.str().c_str();

}

#if 0
    //
    // In this test case, there is no indel that I am aware of, however, the
    // aligner produces quite a long cigar string with many insertions
    // and deletions.  For now, I filter this case out later, but it would
    // be nice if it would limit itself to one or a small number of indels.
    //
    EXPECT_TRUE( verifyAlign("GAATCAATACGCTCGGGATGCAGCGCCTAGCCGTTGGTTTGAGAATGGTTCTCTAGAGTTATCTTCACCCTCTACCTTGTGTGGCACTATTTCTTCTATGACCTTGAC","TGGCTCAAGACCTGACCTTGTGCACGTCTTGGATGCCAGTTCTATTCCCCTCACAGGCCATATGAATCCTGTCCTTTCTGCCTCAAAATGCCCATCCAGAGCCTCTACATTGATTAGCTTTTCCCTCCCTTCCAGAAAAGTTCAAAGGCTACCTCCTCCTTGAAGCCTTCACAAATACCTTAATCTAACTGTTTATAACCCTCTGCCATCTTAGCACTGTGGAAAATACATAAACTTGGGGTTAGAACATCATCAGTTTTAATGTAAGCACCATCCTTTCTAGCTGTACAGCCCTCCTGAGCCTTAGTTCTTACATCTTGAAGATGGAACCAGCTCAACAAGCATAGGGATGTAGCAAGAATCAAGACACTGTAGATGCAGCACCTGGCCGGTGGTAAGAGCTTGGTTATCACAAGTTATCTTCACCCTCTACCTTGTGTGGCACTATTTCTTCTATGACCTTGACTGCTCTCTGCTCTGATCTGGAAGTTCGCTGGGAAAAGGTGTCCCCTTTTTTATTACCTACCGGGAGAACCATGAGTGATGCTCTACTTGTAGTATATATACCCTGAGATGATTATTCTTAAAGACTAGTTCTCATGACTTGAGAGTTTGCTCTGTGTTAGGTACCATTCTAACACTGGATGTTGACTATGTATGTTATTTAATACTTCCATCAACCCCATAATGTAGGGAGAATCATTATGCCCATTTTA", "108M") ) << ss.str().c_str();
#endif


TEST(GreedyTupleTestAligner, EquivalentAlignment) {
    std::stringstream ss(std::stringstream::out);
    // The test cases here are different to what we expected, 
    // but hte alignment are both valid. See the example:
    // ref:  AAAT        TGGGG    (4M5D5M)
    // ref:  AAA        TTGGGG    (3M5D6M)
    // read: AAAT CCCC  TTGGGG

    // We obtain 89M200D19M, the correct is 87M200D21M,
    // the 2 matched bases shifted to the right
    // if you check the output, you will see the 2 bases can be put before '200D' or after '200D'
    // and they are essentially the same
    EXPECT_FALSE( verifyAlign("TTTTCTTTTCAAAAATTTAAAAGTGACATACAAAATTATATGTGTATGTACAACAAAAGCTTAACTATAACACCTTGTTACATACTTTGGAATTGAAAGGCAGGAATG","CAGCACCCTAATTCACTATGCCCTAAGCTTCAAGGGCTTCAGAGTAAGCTCTCAGTGGAGTCTGATTGGAATCCCTCTTCGCCAGCTTGTGAGGTATGGGGCTAGGTTCCACAATATTCCCTTTGAGGGAGTAGATCTTCCAGCCTTCTGGGGCATGCTCTGAAAGTCCTCTTTGCAGAAGTAGCTCTTTAAAATCATATTCTCTTTCCAATTTGACCTCTTTTTTTATCCTTGTTCTGTCCATGCTGTCCAAAGCATCTTGGACTAAGTTTTGACTTTTTTTTTAAGTGCTGCATTTCCATTTGACATTTTACCTTTGTAAATTTCTATTTTTTTACCTTTGTGACTTATTAAAATATTTTCTTTTCAAAAATTTAAAAGTGACATACAAAATTATATGTGTATGTACAACAAAAGCTTAACTATAACACCTTGTTACATACTTTGCTATCCAGGCCACTGATCCTTTCTTACATAGTAAGTCAGCTATAGTTCATTAGCTTACAGTTTTTAGATACAAGTCTTAATCCATCCCTTCTCCTTTTGTATTCTTTACTTTCTGCAATATTTAAGACTTTTTGCGTTCTGACTAAAAGAAACCACCTGAAATTGGCATATGCAACTGTTCATGAATGAGAACTCGCATGGAATTGAAAGGCAGGAATGCAGCTTGACCTTAGAATGGATTTGATCCAGGAACTAGAAGGTGGGTAGGA", "87M200D21M") ) << ss.str().c_str();

    // for the same reason as before 
    EXPECT_FALSE( verifyAlign("AACAGTTGAGAGGTACTAAAATTGAGTTTTCTTGAAAAATATATTTAATCTAAAGTACTGAAAATTTGGGGGAAAATGCTTAAGGTCATATTCCTTTTTTGAAAAGAT","TCATCTTTCTCCCATACTGGCTGTTTCCTGCCCTCAAACACTGGACTCCAAGTTCTTCAGCTTGTGGACTCTTGGACCTACAACCAGTGGTCTGCCAGGGCCCTTTGGGCCTTCGGCCACAGACTGATGGCTACACTGTCGGCTCCCCTACTTTTGAGGTTTTGTGTCTTGGACTGGCTTTCTTGCTCCTCAGCTTGCAGACAGCCTACTGTAGGACTTCACTTTGTGACTATTTGAGTCAATACTCCTTAATAAACACCCTTTCATATATACATATATCCTATTAGTCCTGTTCCTCTAGAGAACCCTAATACAGTGTTGTACATTGAAATAAATATAATTATTCTGGTTTTGGTTGAACAGTTGAGAGGTACTAAAATTGAGTTTTCTTGAAAAATATATTTAATCTAAAGTACTGAAAATTTGGGGGAAAATGCTTCTGTAAATCCTAAGTTATTATTTCTTCAACTATATTCTGTAGTTAATTTCTCCAGCAATTCTTAATTTCAGCACAAATTAGCCACTGTTTGAATTAGGAATACTGAATCGTCTCCATTGCAGTGCAGTTAATAAGTCATTTCTTGATGAAGTAGTCCATGTAGGACTTGAAATCTTGTCTTTTTCATGATACATTATCATAAGGTCATATTCCTTTTTTGAAAAGATTGATGATACTATTCTGAAAGACACTAGTAGAGTTAGGCTTGGTTTTATGA", "80M200D28M") ) << ss.str().c_str();

}

