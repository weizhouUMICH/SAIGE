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

#include "Generic.h"

#if defined(TEST)

#include <iostream>
#include <list>
#include <utility>
#include <vector>

//
// g++ -g -o testGeneric -DTEST Generic.cpp
//
int main(int argc, const char **argv)
{
    std::vector<int> a;
    std::vector< std::pair<int, int> > b;
    std::pair<int, int> c;

    std::vector<int>::iterator i;

    a.push_back(0);
    a.push_back(1);
    a.push_back(2);
    a.push_back(3);

    std::cout << a;

    c.first = 10;
    c.second = 20;
    b.push_back(c);
    b.push_back(c);
    b.push_back(c);

    std::cout << b;

    i = a.begin();

    std::list<std::pair<int, int> > l;

    l.push_back(c);

    std::cout << l;

//  std::cout << "iterator i: " << i << std::endl;

//  std::cout << argv;
}

#endif
