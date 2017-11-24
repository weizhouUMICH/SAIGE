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

#include <iostream>

#include "String.h"
#include <vector>

#ifdef OBSOLETE

std::vector<csg::string> *csg::string::split(char splitChar)
{
    std::vector<csg::string> *result = new std::vector<csg::string>;
    csg::string word;

    for (size_t i = 0; i<size(); i++)
    {
        if ((*this)[i]==splitChar)
        {
            result->push_back(word);
            word.clear();
        }
        else
            word.push_back((*this)[i]);
    }
    if (word.size()>0) result->push_back(word);
    return result;
}


#if defined(TEST)

int main(int argc, const char **argv)
{
    csg::string string("abcdef:abcdefghijk");

    std::vector<csg::string>    *result = string.split(':');

    for (int i=0; i<result->size(); i++)
    {
        std::cout << i << "\t" << (*result)[i] << std::endl;
    }
    delete result;  // suck

}
#endif

#endif
