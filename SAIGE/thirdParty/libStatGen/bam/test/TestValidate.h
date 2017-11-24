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

#include "SamFile.h"

void validateRead1(SamRecord& samRecord);
void validateRead2(SamRecord& samRecord);
void validateRead3(SamRecord& samRecord);
void validateRead4(SamRecord& samRecord);
void validateRead5(SamRecord& samRecord);
void validateRead6(SamRecord& samRecord);
void validateRead7(SamRecord& samRecord);
void validateRead8(SamRecord& samRecord);
void validateRead9(SamRecord& samRecord);
void validateRead10(SamRecord& samRecord);


void validateHeader(SamFileHeader& samHeader);
void validateHeaderFields(SamFileHeader& samHeader);
void validateHeaderString(SamFileHeader& samHeader);

class TestValidate
{
public:
    static const int READ1_POS = 1010;
    static const int READ1_ALIGN_END = 1016;
    static const int READ1_UNCLIP_START = 1010;
    static const int READ1_UNCLIP_END = 1016;
    static const int READ1_ALIGN_LEN = 7;
    static const std::string READ1_CIGAR;
    static const std::string READ1_SEQ;
    static const std::string READ1_QUAL;

    static const int READ2_POS = 1011;

    static const int READ6_POS = 1750;
    static const int READ6_ALIGN_END = 1754;
    static const int READ6_UNCLIP_START = 1745;
    static const int READ6_UNCLIP_END = 1754;
    static const int READ6_ALIGN_LEN = 5;
    static const std::string READ6_CIGAR;
    static const std::string READ6_SEQ;
    static const std::string READ6_QUAL;

    static const int READ7_POS = 1750;
    static const int READ7_ALIGN_END = 1754;
    static const int READ7_UNCLIP_START = 1747;
    static const int READ7_UNCLIP_END = 1758;
    static const int READ7_ALIGN_LEN = 5;
    static const std::string READ7_CIGAR;
    static const std::string READ7_SEQ;
    static const std::string READ7_QUAL;

};
