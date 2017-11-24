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
#include "MemoryMapArray.h"

void MemoryMapArrayHeader::debugPrint(FILE *f)
{
    time_t local = creationDate;
    fprintf(f, "typeCookie = %08x\n", typeCookie);
    fprintf(f, "typeVersion = %08x\n", typeVersion);
    fprintf(f, "contentCookie = %08x\n", contentCookie);
    fprintf(f, "contentVersion = %08x\n", contentVersion);
    fprintf(f, "Created on %s", asctime(localtime(&local)));
    fprintf(f, "Created by user %s on host %s for application '%s'.\n",
            creationUser,
            creationHost,
            application);
}

std::ostream &operator << (std::ostream &stream, MemoryMapArrayHeader &h)
{
    time_t local = h.creationDate;
    stream << "typeCookie = " << h.typeCookie << "\n";
    stream << "typeVersion = " << h.typeVersion << "\n";
    stream << "contentCookie = " << h.contentCookie << "\n";
    stream << "contentVersion = " << h.contentVersion << "\n";
    stream << "headerSize = " << h.headerSize << "\n";
    stream << "elementCount = " << h.elementCount << "\n";

    stream << "Created on " << asctime(localtime(&local)) << "\n";
    stream << "Created by user " << h.creationUser << " on host " << h.creationHost << " for application '" << h.application << "'.\n";
    return stream;
}

#if defined(TEST)
#include <assert.h>
#include <stdlib.h>

void test32()
{
    mmapArrayUint32_t   test;

    unlink("twinkypie");
    assert(test.create("twinkypie", 11)==0);
    test.set(0,0);
    test.set(1,1);
    test.set(2,2);
    test.set(3,3);
    test.set(4,4);
    test.set(5,5);
    test.set(6,6);
    test.set(7,7);
    test.set(8,8);
    test.set(9,9);
    test.set(10,10);
    assert(test[0]==0);
    assert(test[10]==10);
    test.close();
    assert(test.open("twinkypie")==0);
    assert(test[0]==0);
    assert(test[10]==10);
    test.close();
    unlink("twinkypie");
}

void testbool()
{
    mmapArrayBool_t   test;

    unlink("twinkypie");
    assert(test.create("twinkypie", 11)==0);
    test.set(0,0);
    test.set(1,1);
    test.set(2,0);
    test.set(3,1);
    test.set(4,0);
    test.set(5,1);
    test.set(6,0);
    test.set(7,1);
    test.set(8,0);
    test.set(9,0);
    test.set(10,1);
    assert(test[0]==0);
    assert(test[1]==1);
    assert(test[10]==1);
    test.close();
    assert(test.open("twinkypie")==0);
    assert(test[0]==0);
    assert(test[10]==1);
    test.close();
    unlink("twinkypie");
}

void test2bit()
{
    mmapArray2Bit_t   test;

    unlink("twinkypie");
    assert(test.create("twinkypie", 11)==0);
    test.set(0,0);
    test.set(1,1);
    test.set(2,2);
    test.set(3,3);
    test.set(4,3);
    test.set(5,2);
    test.set(6,1);
    test.set(7,0);
    test.set(8,2);
    test.set(9,1);
    test.set(10,3);
    test.setApplication("testing 2 bit values!");
    assert(test[0]==0);
    assert(test[1]==1);
    assert(test[2]==2);
    assert(test[3]==3);
    assert(test[4]==3);
    assert(test[5]==2);
    assert(test[6]==1);
    assert(test[7]==0);
    assert(test[8]==2);
    assert(test[9]==1);
    assert(test[10]==3);
    test.close();
    assert(test.open("twinkypie")==0);
    test.debugPrint(stdout);
    test.close();
    unlink("twinkypie");
}

void test4bit()
{
    mmapArray4Bit_t   test;

    unlink("twinkypie");
    assert(test.create("twinkypie", 11)==0);
    test.set(0,0);
    test.set(1,1);
    test.set(2,2);
    test.set(3,3);
    test.set(4,4);
    test.set(5,5);
    test.set(6,6);
    test.set(7,7);
    test.set(8,8);
    test.set(9,9);
    test.set(10,10);
    test.setApplication("testing 4 bit values!");
    assert(test[0]==0);
    assert(test[1]==1);
    assert(test[7]==7);
    assert(test[10]==10);
    test.close();
    assert(test.open("twinkypie")==0);
    assert(test[0]==0);
    assert(test[1]==1);
    assert(test[7]==7);
    assert(test[10]==10);
    test.debugPrint(stdout);
    test.close();
    unlink("twinkypie");
}

int main(int argc, char **argv)
{

    test32();
    testbool();
    test2bit();
    test4bit();
    exit(0);
}

#endif
