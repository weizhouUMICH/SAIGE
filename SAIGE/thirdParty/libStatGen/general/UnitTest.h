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

#ifndef __UNITTEST_H
#define __UNITTEST_H

#include <string>
#include <iostream>

class UnitTest
{
protected:
    std::string m_title;
    int m_failures;
    int m_testNum;

public:
    UnitTest(const char *title) : m_title(title), m_failures(0), m_testNum(0) {;};
    void test();
    int getPassCount() {return m_testNum - m_failures;}
    int getFailureCount() {return m_failures;}
    const std::string getTitle() const {return m_title;}
};

std::ostream &operator << (std::ostream &stream, UnitTest &test)
{
    stream << test.getTitle() << " PASS: " << test.getPassCount() <<
        "  FAIL: " << test.getFailureCount() << std::endl;
    return stream;
}

#endif
