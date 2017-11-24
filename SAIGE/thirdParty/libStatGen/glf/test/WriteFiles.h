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
#ifndef __WRITE_FILES_H__
#define __WRITE_FILES_H__

#include <stdio.h>
#include "GlfFile.h"


void testWrite();
void testHeaderWrite();
void testWriteCopiedHeader();

class TestWrite
{
public:
    void testWrite();
private:
    void writeHeader(GlfFile& glfOut, int headerNum);
    void writeRefSection1(GlfFile& glfOut);
    void writeSec1Record1(GlfFile& glfOut);
    void writeSec1Record2(GlfFile& glfOut);
    void writeEndMarker(GlfFile& glfOut);
    void writeRefSection2(GlfFile& glfOut);
    void writeSec2Record1(GlfFile& glfOut);
    
    void readHeader(GlfFile& glfIn, int headerNum);
    void readRefSection1(GlfFile& glfIn);
    void readSec1Record1(GlfFile& glfIn);
    void readSec1Record2(GlfFile& glfIn);
    void readEndMarker(GlfFile& glfIn);
    void readRefSection2(GlfFile& glfIn);
    void readSec2Record1(GlfFile& glfIn);
    
    void checkEOF(GlfFile& glfIn);
    
    // 1st file header values:
    static const std::string HEADER_TEXT1;
    
    // SEC1 values:
    static const std::string SEC1_REFNAME;
    static const uint32_t SEC1_REFLEN = 200;
    
    // SEC1REC1 values:
    static const uint8_t SEC1REC1_RECTYPE = 1;
    static const uint8_t SEC1REC1_REFBASE = 4;
    static const uint32_t SEC1REC1_OFFSET = 99;
    static const uint32_t SEC1REC1_MINLK = 55;
    static const uint32_t SEC1REC1_READDEPTH = 31;
    static const uint8_t SEC1REC1_RMSMAPQ = 25;
    
    // SEC1REC2 values:
    static const uint8_t SEC1REC2_RECTYPE = 2;
    static const uint8_t SEC1REC2_REFBASE = 1;
    static const uint32_t SEC1REC2_OFFSET = 6;
    static const uint32_t SEC1REC2_MINLK = 44;
    static const uint32_t SEC1REC2_READDEPTH = 66;
    static const uint8_t SEC1REC2_RMSMAPQ = 32;
    static const uint8_t SEC1REC2_LKHOM1 = 98;
    static const uint8_t SEC1REC2_LKHOM2 = 86;
    static const uint8_t SEC1REC2_LKHET = 73;
    static const int16_t SEC1REC2_INDELLEN1 = 2;
    static const int16_t SEC1REC2_INDELLEN2 = -3;
    static const std::string SEC1REC2_INDELSEQ1;
    static const std::string SEC1REC2_INDELSEQ2;
    
    // SEC2 values
    static const std::string SEC2_REFNAME;
    static const uint32_t SEC2_REFLEN = 102;
    
    // SEC2REC1 values:
    static const uint8_t SEC2REC1_RECTYPE = 1;
    static const uint8_t SEC2REC1_REFBASE = 2;
    static const uint32_t SEC2REC1_OFFSET = 50;
    static const uint32_t SEC2REC1_MINLK = 55;
    static const uint32_t SEC2REC1_READDEPTH = 31;
    static const uint8_t SEC2REC1_RMSMAPQ = 25;

    // 2nd file header.
    static const std::string HEADER_TEXT2;
    // 3rd file header.
    static const std::string HEADER_TEXT3;
};
#endif
