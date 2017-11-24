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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "SmithWaterman.h"

// put TEST below here, so that makedepend will see the .h, so that we
// can get a clean dependency for SmithWaterman.o, so that we can at least
// compile the header when we change it.

#if defined(TEST)

#include <getopt.h>
#include "Generic.h"

// g++ -g -o testSW -DTEST SmithWaterman.cpp
//
// Smith-Waterman - test code uses a 256x256 array of int16
//
int swat(
    bool        showAllCases,
    const char *A,
    const char *qualities,
    const char *B,
    int direction,
    const char *expectedCigarString,
    int expectedSumQ
)
{
    int allowedInsertDelete = 1024;
    int errors = 0;

    //            read length 256
    //                 reference length 1024
    SmithWaterman<256, 1024, uint16_t, const char *, const char *, const char *, uint32_t, uint32_t > sw(&A, &qualities, &B, strlen(A), strlen(B), allowedInsertDelete, direction);

    //
    // now we align the read:
    //
    sw.populateH();
    sw.populateAlignment();

    int sumQ = sw.getSumQ();

    CigarRoller cigar;
    cigar.clear();
    sw.rollCigar(cigar);

    const char *cigarStr = cigar.getString();

    //
    // now we pretty print the results
    //

    bool badCigar = false, badQuality = false;

    if (strcmp(cigarStr, expectedCigarString)!=0)
    {
        badCigar = true;
        errors ++;
    }

    if (sumQ != expectedSumQ)
    {
        badQuality = true;
        errors ++;
    }



    if (showAllCases || errors>0)
    {
        cout << "=============" << endl;
        cout << "        Read: " << A << endl;
        cout << "   Reference: " << B << endl;
        cout << "   Direction: " << direction << endl;
        cout << "Max Cell: " << sw.maxCostValue << " located at " << sw.maxCostPosition << endl;
        cout << "M: " << sw.m << "  N: " << sw.n << endl;
        cout << "Cigar String: " << cigarStr ;
        if (badCigar)
            cout << " (EXPECTED: " << expectedCigarString << ")";
        cout << endl;
        cout << "        sumQ:" << sumQ;
        if (badQuality)
            cout << " (EXPECTED: " << expectedSumQ << ")";
        cout << endl;

        if (strlen(B) < 100 || showAllCases)
            sw.printH(false);

        for (vector<pair<int,int> >::iterator i = sw.alignment.begin(); i != sw.alignment.end(); i++) cout << *i << endl;

        cout << "=============" << endl << endl;
    }

    delete cigarStr;

    return errors;
}


// test with Sequence 1 = ACACACTA
//           Sequence 2 = AGCACACA
int main(int argc, const char **argv)
{
    int errors = 0;

    bool showAllCasesFlag = false;
    int opt;

    while ((opt = getopt(argc, (char **) argv, "v")) != -1)
    {
        switch (opt)
        {
            case 'v':
                showAllCasesFlag = true;
                break;
            default:
                cerr << "usage: testSW [-v]" << std::endl;
                exit(1);
        }
    }


    // CIGAR explanation - for backward SW runs, the corresponding
    // CIGAR string is generated from the back of the string to the
    // front.  Recall that the soft clipping is only done at the
    // "end" of the string, taking direction into account.

    // forwards - simple
    errors += swat(showAllCasesFlag, "1234", "\"#$-", "1235", 1, "3M1S", 0);

    // backwards - simple
    errors += swat(showAllCasesFlag, "1234", "\"#$-", "1235", -1, "4M", 12);

    // backwards - soft left clip
    errors += swat(showAllCasesFlag, "1234", "\"#$-", "0234", -1, "1S3M", 0);

    // delete in read (arg 1) - forward
    errors += swat(showAllCasesFlag, "123467890", "\"#$%^&*()-", "1234567890", +1, "4M1D5M", 50);

    // insert in read (arg 1) - forward
    errors += swat(showAllCasesFlag, "1234556789", "\"#$%^&*()-", "1234567890", +1, "5M1I4M", 50);

    // delete in read (arg 1) - backward
    errors += swat(showAllCasesFlag, "X123467890", "#\"#$%^&*()-", "1234567890", -1, "1S4M1D5M", 50);

    // insert in read (arg 1) - backward
    errors += swat(showAllCasesFlag, "1234556789", "\"#$%^&*()-", "0123456789", -1, "4M1I5M", 50);

    // insert in read (arg 1) - backward
    errors += swat(showAllCasesFlag, "X1223456789", "00000000000", "00123456789", -1, "1S1M1I8M", 50);

    // insert in read (arg 1) - backward
    errors += swat(showAllCasesFlag, "XY1223456789", "000000000000", "000123456789", -1, "2S1M1I8M", 50);

    // forward - soft right clip of 2 bases - sumQ should be 0
    errors += swat(showAllCasesFlag, "123456700", "\"#$%^&*()-", "123456789", +1, "7M2S", 0);

    // insert in read (arg 1) - forward w/2 mismatches at end
    errors += swat(showAllCasesFlag, "1023456700", "\"#$%^&*()-", "123456789", +1, "1M1I6M2S", 50);

    // insert in read (arg 1) - forward w/2 mismatches at end
    errors += swat(showAllCasesFlag, "CTCCACCTCCCGGTT", "111111111111111", "TCCACCTCCCAGGTT", -1, "1S10M1D4M", 50);

    //
    errors += swat(showAllCasesFlag, "1234", "0000", "12345", +1, "4M", 0);

    //
    errors += swat(showAllCasesFlag, "1234X", "00000", "12345", +1, "4M1S", 0);

    //
    errors += swat(showAllCasesFlag, "4321", "0000", "7654321", -1, "4M", 0);

    //
    errors += swat(showAllCasesFlag, "X4321", "00000", "7654321", -1, "1S4M", 0);

    //
    errors += swat(showAllCasesFlag, "X432A10", "0000000", "76543210", -1, "1S3M1I2M", 50);

    //
    errors += swat(showAllCasesFlag, "1345", "0000", "12345", -1, "1M1D3M", 50);

    errors += swat(showAllCasesFlag, "45689", "00000", "1234567890", -1, "3M1D2M", 50);

//  errors += swat(showAllCasesFlag, "AATAATTTTTTATATACAGATCGCTGTAGAGTGTAGTTATAGTATGATTCCAACTTTTATTTCTTTCATGACTAATTATATGTATACATGTGCCATGTTGGTGTGCTG", "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", "TCCAATGTAGGGCTGTTATAAACAGTGTTGATACATATGTTTTTGTATAAGTCTTTGTTGAATACATGCTTTCATTTTTGTAGGGTATATGTCCAGGAATTAAATTTTTGCATTATTGGGGAAGTTCAAACGTAGATCAGTAGATGTTCCCAAATGATTTTCAGGATATGTATCCATGTAAATTCCTACCAGCAATGCAGGAGAATTCCAATTGCCCATGTTCTAATCAGAATATTGTTATATCCTAAGACTAATTTTAAATATTCTGATGGGTGTAGAGTGGAGGCATAGTATGATTTCAACTTGTATTTCTTTCATGACTAATTATCTTCTATGTTAATTGTTATTTTGTATGTTTATTGCAAAGTGCCTATCCAGAATTTTTGTCTATAATTTTGTTGTGCTGTCTCTTGCTTTATGAATTTTATAGGATTCTTAATATTATAATTGAGTTATCTTTCTTTTTTATTATTATTATTATACTTTAAGTTTTAGGGTATATGTGCACAACGTGCAAGTTTGTCACATATGTATACATGTGCCATGTTGGTGTGCTGCACCCATTAACTCATCATTTAGCATTAGGTATATCTCCTAATGCTATCCCTTCCTCCTCCCCCCACCCCACAACAGTCCCCGGTGTGTGATGTTCCCCTGCCTTTGTCCTCTTTCTTATACTTGCATGAGCAATCTCCTCAAACTGATACTTGCCTTTTTTGTCCTTGGTGTGGTTTGGCTCTGTGTTCCCACCCAAATCTTCATAATACCCATGTGCCAAGGGTGGGACTGGGTGGAGGTAATTGGGTCATGGGGATGGTTTCCCTCATACTATTATGATAGTGAGTGTTTTCACGAGACCTGATGGTTTTATAACTGTGTGGCATTTCCCTTGCTTCCACTCACTCCATCCTGCCACCCTGTGAAGAAGGTGCCTGCTTCTCCTTTGGTTACTGCTATGATTGTAAGTTTCCTGAGGCCTCCCCAGCAACGCAAAACTGTGAATCAATTAAACCTTTTTCCTTTATAAATTACTAAGTCTTGGGTATTTCTTCATAGTGTTGTGAGCATAGACTAAAACAGTAAGTTGTTACCAGGAGTGGGGTACTGCTGTAAGATAACTGAGAATGTGAAAGTGACTTAGGAACTAGGTAATGAGCAGAGGTTGGAACAGTTTAAAAGGCTCAGAAGAAGACAGAAAGATGTGGGAAAGTTTGGA", -1, "77M200D31M", 50);

    errors += swat(showAllCasesFlag, "TTAGAATGCTATTGTGTTTGGAGATTTGAGGAAAGTGGGCGTGAAGACTTAGTGTTCATTTCCTCAACCTCTCTCTGTGTGAACATACGTCATCGGTCAGAAATTGGG", "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", "CCGAGATTGTGCCATTGCACTCCTGCCTGGGTAACAGAGTCAGACCCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAGATTAGGTTTTATAGATGGAAAATTCACAGCTCTCTCCAGATCAGAAATCTCCAAGAGTAAATTAGTGTCTTAAAGGGGTTGTAATAACTTTCCTATGTGACTAAGTGCATTATTAATCAATTTTTCTATGATCAAGTACTCCTTTACATACCTGCTAATACAATTTTTGATATGAAATCAGTCCTAGAGGGAATCAATGTAAGATACAGACTTGATGAGTGCTTGCAGTTTTTTATTGACAATCTGAAGAATGACTTGACTCTAAATTGCAGCTCAAGGCTTAGAATGCTATTGTGTTTGGAGATTTGAGGAAAGTGGGCGTGAAGACTTAGTGTTCATTTCCTCAACCTCTCTCTGTGTGAACATACAGGAATCAAATCTGTCTAGCCTCTCTTTTTGGCAAGGTTAAGAACAATTCCACTTCATCCTAATCCCAATGATTCCTGCCGACCCTCTTCCAAAAACTATTTAAAGACATGTTCTTCAAAGTTATATTTGTCTTTCCTTCAGGGAGAAAAAGAATACCAATCACTTATAATATGGAAACTAGCAGAAATGGGTCACATAAGTCATCTGTCAGAAATTGGGAAAATAGAGTAGGTCAGTCTTTCCAGTCATGGTACTTTTACCTTCAATCA", -1, "88M200D20M", 50);

    // prefix TTAGAATGCTATTGTGTTTGGAGATTTGAGGAAAGTGGGCGTGAAGACTTAGTGTTCATTTCCTCAACCTCTCTCTGTGTGAACATAC
    // suffix GTCATCTGTCAGAAATTGGGA
    cout << endl << "Total Errors found: " << errors << endl;
}
#endif
