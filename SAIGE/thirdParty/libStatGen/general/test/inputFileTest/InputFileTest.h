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
#include "InputFile.h"


class IFILE_Test : public InputFile
{
public:
    void test();

    static const int TEST_FILE_SIZE;
    static const int BGZF_TEST_FILE_SIZE;
    static const std::string TEST_FILE_CONTENTS;

private:
    void testAll(const char* extension);
    void test_readFromFile(const char* extension);
    void test_readTilChar(const char* extension);

    // Tested together because they are used to test each other.
    void test_ifeof_ifrewind(const char* extension);

    // Tested together to verify they can be successfully be called after the
    // other has been called.
    void test_ifread_ifgetc(const char* extension);

    void test_ifclose(const char* extension);

    void test_ifseek(const char* extension);

    void test_noExistRead(const char *extension);

    void openFile(const char* extension);
    void openLargeFile(const char* extension);
    void openNoExistFile(const char* extension);

    // Buffer used for reading into.
    static const int MAX_TEST_BUFFER_SIZE = 100;
    char myTestBuffer[MAX_TEST_BUFFER_SIZE];

};
