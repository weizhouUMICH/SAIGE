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
#include "TestPosList.h"
#include <assert.h>
#include <stdexcept>


void testPosList()
{
    TestPosList posListTest;
    posListTest.testPosList();
}

TestPosList::TestPosList()
{
}


TestPosList::~TestPosList()
{
}


void TestPosList::testPosList()
{
    assert(myPosList.size() == 24);
    
    for(int i = 0; i < 24; i++)
    {
        assert(myPosList.at(i).size() == 100);
    }

    bool caught = false;
    try
    {
        myPosList.at(24);
    }
    catch(std::out_of_range& oor)
    {
        caught = true;
    }

    assert(caught == true);

    //////////////////////////////
    // Test accessing
    for(int i = 0; i < 24; i++)
    {
        for(int j = 0; j < 100; j++)
        {
            assert(!hasPosition(i, j));
        }
    }

    //////////////////////////////
    // Test setting all
    for(int i = 0; i < 24; i++)
    {
        for(int j = 0; j < 100; j++)
        {
            addPosition(i, j);
        }
    }
    for(int i = 0; i < 24; i++)
    {
        for(int j = 0; j < 100; j++)
        {
            assert(hasPosition(i, j));
        }
    }




    //////////////////////////////
    // Test accessing out of range
    assert(!hasPosition(-1, 0));
    assert(!hasPosition(0, -1));
    assert(!hasPosition(100, 0));
    assert(!hasPosition(0, 1000));

    //////////////////////////////
    // Test adding more to ref 4,
    // but skipping positions.
    for(int j = 300; j < 350; j++)
    {
        addPosition(4, j);
    }
    for(int j = 0; j < 100; j++)
    {
        assert(hasPosition(4, j));
    }
    for(int j = 100; j < 300; j++)
    {
        assert(!hasPosition(4, j));
    }
    for(int j = 300; j < 350; j++)
    {
        assert(hasPosition(4, j));
    }

    // Test adding a new reference, 30,
    // position 16.
    addPosition(30, 16);

    // Check the size now.
    assert(myPosList.size() == 31);
    
    for(int i = 0; i < 24; i++)
    {
        if(i != 4)
        {
            assert(myPosList.at(i).size() == 100);
        }
        else
        {
            assert(myPosList.at(i).size() == 350);
        }
    }

    for(int i = 24; i < 31; i++)
    {
        assert(myPosList.at(i).size() == 350);
    }

    //////////////////////////////
    // Test accessing
    for(int i = 24; i < 30; i++)
    {
        for(int j = 0; j < 350; j++)
        {
            assert(!hasPosition(i, j));
        }
    }
    for(int j = 0; j < 350; j++)
    {
        if(j != 16)
        {
            assert(!hasPosition(30, j));
        }
        else
        {
            assert(hasPosition(30, 16));
        }
    }
}
