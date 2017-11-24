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
#include "TestPileup.h"

void testPileup()
{
    TestPileup pileupTest;
    pileupTest.testPileupPosition();
}

void TestPileupElement::analyze()
{
    assert(strcmp(getChromosome(), "") == 0);
    assert(getRefPosition() == 14000);
}


void TestPileup::testPileupPosition()
{
    assert(pileupPosition(14000) == 0);
    assert(pileupHead == 14000);
    assert(pileupStart == 14000);
    assert(pileupTail == 14000);

    bool caught = false;
    try
    {
        pileupPosition(13999);
    }
    catch(std::exception& e)
    {
        caught = true;
        assert(strcmp(e.what(), "Overflow on the pileup buffer: specifiedPosition = 13999, pileup buffer start position: 14000, pileup buffer end position: 15024") == 0);
    }
    assert(caught);

    caught = false;
    try
    {
        pileupPosition(15025);
    }
    catch(std::exception& e)
    {
        caught = true;
        assert(strcmp(e.what(), "Overflow on the pileup buffer: specifiedPosition = 15025, pileup buffer start position: 14000, pileup buffer end position: 15024") == 0);
    }
    assert(caught);
}
