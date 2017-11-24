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
#include <string>
#include "STLUtilities.h"

//
#define _STLUTILITIES_BENCHMARK_
// This can turn on the benchmark of STLUtilities class and String class
//
#ifdef _STLUTILITIES_BENCHMARK_
#include "Performance.h"
#include "Random.h"
#include "StringBasics.h"
#endif /* _STLUTILITIES_BENCHMARK_ */


#include <gtest/gtest.h>

TEST(STLUtilitiesTest, tSTLUtilitiesTest)
{
#if 0
    std::string test;
    std::string::iterator result;

    test = "445566";
    result = trimSequence(test, '5', true);
    EXPECT_EQ(result - test.begin() , 2);

    test = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, 'A', true);
    EXPECT_TRUE(result == test.begin());
#endif


    using namespace STLUtilities; // for overloaded std::string << operator

    std::string toot;

    toot << "double: " << 0.123456 << " LL: " << -5LL << " UL: " << 999UL << " char: " << '!';

    EXPECT_TRUE(toot == "double: 0.123456 LL: -5 UL: 999 char: !");


    // same result as above using different methods
    toot.clear();

    append(toot, "double: ");

    append(toot, 0.123456);

    append(toot, " LL: ");

    append(toot, -5LL);

    append(toot, " UL: ");

    append(toot, 999UL);

    append(toot, " char: ");

    append(toot, (char) (' ' + 1));

    EXPECT_TRUE(toot == "double: 0.123456 LL: -5 UL: 999 char: !");


    toot.clear();

    std::vector<int> v;
    v.push_back(1);
    v.push_back(2);
    v.push_back(3);
    v.push_back(4);
    v.push_back(5);

    toot = "array: ";
    append(toot, v, "\t", true);

    EXPECT_TRUE(toot == "array: 0: 1\t1: 2\t2: 3\t3: 4\t4: 5");

    std::vector<std::string> tokens;

    Tokenize(tokens, "ab\tcd\tefg\thi\tjk", '\t');
    EXPECT_EQ(tokens.size(), 5U);
    EXPECT_TRUE(tokens[0] == "ab");
    EXPECT_TRUE(tokens[1] == "cd");
    EXPECT_TRUE(tokens[2] == "efg");
    EXPECT_TRUE(tokens[3] == "hi");
    EXPECT_TRUE(tokens[4] == "jk");


    Tokenize(tokens, "ab\tcd\tefg\thi\tjk\t", '\t');
    EXPECT_EQ(tokens.size(), 6U);
    EXPECT_TRUE(tokens[5] == "");

    // a single tab splits two empty fields, so should see two tokens here:
    Tokenize(tokens, "\t", '\t');
    EXPECT_EQ(tokens.size(), 2U);
    EXPECT_TRUE(tokens[0] == "");
    EXPECT_TRUE(tokens[1] == "");


    Tokenize(tokens, "bahbah", '\t');
    EXPECT_EQ(tokens.size(), 1U);
    EXPECT_TRUE(tokens[0] == "bahbah");

    //
    // no data on the line is the same as a single empty field.
    // the reason is we don't want to have a file with a single
    // column of data, but two separate values for .size().  Better
    // to let the caller simply say 'if tokens[0]==""'
    //
    Tokenize(tokens, "", '\t');
    EXPECT_EQ(tokens.size(), 1U);
    EXPECT_TRUE(tokens[0] == "");

#if 0
    toot = "";
    append(toot, tokens, '\t');
    std::cout << toot << std::endl;

    exit(0);
#endif
}

//
// Variadic templates necessary for reasonable printf implementation
// are only supported as an experimental feature that in theory is
// subject to changes in the future draft standard for C++.
// 
// Only defined when the g++ option -std=c++0x is used.
//
#if defined(__GXX_EXPERIMENTAL_CXX0X__)
TEST(STLUtilitiesTestPrintf, tSTLUtilitiesTestPrintf)
{
    using namespace STLUtilities; // for overloaded std::string << operator

    std::string result;

    sprintf(result, "Hello, world!");
    EXPECT_TRUE(result=="Hello, world!");

    sprintf(result, "n = %20.5lXXX", 123ULL);
    EXPECT_TRUE(result=="n =                   7bXX");

    sprintf(result, "hello, world! %20sand boo", "well then");
    EXPECT_TRUE(result=="hello, world!            well thenand boo");

    sprintf(result, "addr = %08xXXX", 1234);
    EXPECT_TRUE(result=="addr = 000004d2XXX");

    sprintf(result, "Hello, world!! Imagine: %d!", 2345.1234);
    EXPECT_TRUE(result=="Hello, world!! Imagine: 2345.12!");

}
#endif

#ifdef _STLUTILITIES_BENCHMARK_

//
// Compare StringBasics.h String with std::string and STLUtilities append methods
//
// NB: these are mostly inline with the exception of String::operator +(char c), which
// is a function call, so as currently implemented, String is at a disadvantage.
//
// However, all of these potentially suffer from limitations in g++, as I've noticed
// that it at times is unable to inline more than a few levels of nested functions
// deep even if all are trivially short and inlined.
//

TEST(STLUtilitiesTestPrintf, Benchmark1) {
    using namespace STLUtilities; // for overloaded std::string << operator

    std::string result;
    Random random;
    int range = 'z' - 'A';
    unsigned int round = 1e6;
    Timing timing;
    random.Reset(0);
    timing.start();
    for (unsigned int i = 0; i < round; i++)
    {
        result << (char)('A' + ( random.NextInt() % range));
    }
    timing.end();
    std::cout << "STLUtilities " << round << " times takes " << timing.interval() << " second " << std::endl;

    String s;
    random.Reset(0);
    timing.start();
    for (unsigned int i = 0; i < round; i++)
    {
        s += (char)('A' + ( random.NextInt() % range));
    }
    timing.end();
    std::cout << "String " << round << " times takes " << timing.interval() << " second " << std::endl;
    EXPECT_EQ(result, s.c_str());

    std::string st;
    random.Reset(0);
    timing.start();
    for (unsigned int i = 0; i < round; i++)
    {
        st += (char)('A' + ( random.NextInt() % range));
    }
    timing.end();
    std::cout << "std::string " << round << " times takes " << timing.interval() << " second " << std::endl;
    EXPECT_EQ(result, st);

}

TEST(STLUtilitiesTestPrintf, Benchmark2) {
    using namespace STLUtilities; // for overloaded std::string << operator

    std::string result;
    Random random;
    unsigned int round = 1e6;
    Timing timing;
    timing.start();
    for (unsigned int i = 0; i < round; i++)
    {
        result = "";
        for(int j=0; j<15; j++) result << (char) 'A';
    }
    timing.end();
    std::cout << "STLUtilities " << round << " times takes " << timing.interval() << " second " << std::endl;

    String s;
    timing.start();
    for (unsigned int i = 0; i < round; i++)
    {
        s = "";
        for(int j=0; j<15; j++) s += 'A';
    }
    timing.end();
    std::cout << "String " << round << " times takes " << timing.interval() << " second " << std::endl;
    EXPECT_EQ(result, s.c_str());

    std::string st;
    timing.start();
    for (unsigned int i = 0; i < round; i++)
    {
        st = "";
        for(int j=0; j<15; j++) st += 'A';
    }
    timing.end();
    std::cout << "std::string " << round << " times takes " << timing.interval() << " second " << std::endl;
    EXPECT_EQ(result, st);
}

#endif /* _STLUTILITIES_BENCHMARK_ */
