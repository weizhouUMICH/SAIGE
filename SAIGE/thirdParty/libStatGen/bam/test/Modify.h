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

void testModify();

class modify
{
public:
    void testModify(const char* filename);
   
private:
    void modifyPosition();

    void modifyCigar();

    void modifyFlag();

    // Open and read the first record.
    void openAndRead1Rec();

    void modifyTags();

    // Variables.
    std::string myFilename;

    // Rather than passing around all these variables, just store them in the class.
    SamFile samIn;
    SamFileHeader samHeader;
    SamRecord samRecord;
};
