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

#include "InplaceMerge.h"

#if defined(TEST)
#include "iostream"
#include "Generic.h"

int main(int argc, const char **argv)
{
    int a[] = {1,2,3,4,5};
    int b[] = {2,4,6,7,10};
    int c[] = {3,5,8,10,11};

    std::vector<int> z(15);

    std::copy(a, a+5, z.begin());
    std::copy(b, b+5, z.begin() + 5);
    std::copy(c, c+5, z.begin() + 10);

    std::vector<int> indeces, counts;

    indeces.push_back(0);
    indeces.push_back(5);
    indeces.push_back(10);

    counts.push_back(5);
    counts.push_back(5);
    counts.push_back(5);

    inplace_merge(indeces, counts, 0, 3, z);

    std::cout << z;
}

#endif
