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
#include "NonOverlapRegions.h"
#include "NonOverlapRegionsTest.h"
#include <assert.h>
#include <iostream>

int main(int argc, char ** argv)
{

   NonOverlapRegionsTest myTest;

   myTest.test();
}

void NonOverlapRegionsTest::test()
{
    testPos();
    testChrom();
}


void NonOverlapRegionsTest::testChrom()
{
    NonOverlapRegions reg;

    // Assert that the regions are empty.
    assert(reg.myRegions.size() == 0);
    // Verify no regions.
    for(int i = 0; i < 30; i++)
    {
        assert(reg.inRegion("a", i) == false);
        assert(reg.inRegion("3", i) == false);
    }
    // The chromosomes checked for were added.
    assert(reg.myRegions.size() == 2);
    assert(reg.myRegions["a"].myRegions.size() == 0);
    assert(reg.myRegions["a"].myRegionIter == 
           reg.myRegions["a"].myRegions.end());
    assert(reg.myRegions["3"].myRegions.size() == 0);
    assert(reg.myRegions["3"].myRegionIter == 
           reg.myRegions["3"].myRegions.end());

    // Add a region.
    reg.add("3", 13, 15);
    // Verify regions.
    assert(reg.myRegions.size() == 2);
    for(int i = 0; i < 30; i++)
    {
        assert(reg.inRegion("a", i) == false);
        if((i >= 13) && (i < 15))
        {
            assert(reg.inRegion("3", i) == true);
        }
        else
        {
            assert(reg.inRegion("3", i) == false);
        }
    }
    
    // Add a region.
    reg.add("a", 1, 5);
    // Verify regions.
    assert(reg.myRegions.size() == 2);
    for(int i = 0; i < 30; i++)
    {
        if((i >= 1) && (i < 5))
        {
            assert(reg.inRegion("a", i) == true);
        }
        else
        {
            assert(reg.inRegion("a", i) == false);
        }
        if((i >= 13) && (i < 15))
        {
            assert(reg.inRegion("3", i) == true);
        }
        else
        {
            assert(reg.inRegion("3", i) == false);
        }
    }
    

}

void NonOverlapRegionsTest::testPos()
{
    NonOverlapRegionPos reg;
    std::list< std::pair<int32_t, int32_t> >::iterator iter;

    // Assert that the regions are empty.
    assert(reg.myRegions.empty());
    assert(reg.myRegionIter == reg.myRegions.end());
    assert(reg.myTmpIter == reg.myRegions.end());
    // Verify regions.
    for(int i = 0; i < 30; i++)
    {
        assert(reg.inRegion(i) == false);
    }


    // Add a region
    reg.add(13, 15);
    // Verify regions.
    assert(reg.myRegions.size() == 1);
    assert(reg.myRegionIter->first == 13);
    assert(reg.myRegionIter->second == 15);
    for(int i = 0; i < 30; i++)
    {
        if((i >= 13) && (i < 15))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }

    // Insert before this.
    reg.add(4,6);
    assert(reg.myRegions.size() == 2);
    assert(reg.myRegionIter->first == 4);
    assert(reg.myRegionIter->second == 6);
    iter = reg.myRegions.begin();
    assert(iter->first == 4);
    assert(iter->second == 6);
    ++iter;
    assert(iter->first == 13);
    assert(iter->second == 15);
    ++iter;
    assert(iter == reg.myRegions.end());
    for(int i = 0; i < 30; i++)
    {
        if(((i >= 4) && (i < 6)) || ((i >= 13) && (i < 15)))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }

    // Insert at the end.
    reg.add(22,26);
    assert(reg.myRegions.size() == 3);
    assert(reg.myRegionIter->first == 22);
    assert(reg.myRegionIter->second == 26);
    iter = reg.myRegions.begin();
    assert(iter->first == 4);
    assert(iter->second == 6);
    ++iter;
    assert(iter->first == 13);
    assert(iter->second == 15);
    ++iter;
    assert(iter->first == 22);
    assert(iter->second == 26);
    ++iter;
    assert(iter == reg.myRegions.end());
    for(int i = 0; i < 30; i++)
    {
        if(((i >= 4) && (i < 6)) || ((i >= 13) && (i < 15)) || 
           ((i >= 22) && (i < 26)))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }

    // Insert in the middle.
    reg.add(8,9);
    assert(reg.myRegions.size() == 4);
    assert(reg.myRegionIter->first == 8);
    assert(reg.myRegionIter->second == 9);
    iter = reg.myRegions.begin();
    assert(iter->first == 4);
    assert(iter->second == 6);
    ++iter;
    assert(iter->first == 8);
    assert(iter->second == 9);
    ++iter;
    assert(iter->first == 13);
    assert(iter->second == 15);
    ++iter;
    assert(iter->first == 22);
    assert(iter->second == 26);
    ++iter;
    assert(iter == reg.myRegions.end());
    for(int i = 0; i < 30; i++)
    {
        if(((i >= 4) && (i < 6)) || ((i >= 8) && (i < 9)) ||
           ((i >= 13) && (i < 15)) || ((i >= 22) && (i < 26)))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }

    // Insert start does not overlap, but the end does.
    reg.add(20,24);
    assert(reg.myRegions.size() == 4);
    assert(reg.myRegionIter->first == 20);
    assert(reg.myRegionIter->second == 26);
    iter = reg.myRegions.begin();
    assert(iter->first == 4);
    assert(iter->second == 6);
    ++iter;
    assert(iter->first == 8);
    assert(iter->second == 9);
    ++iter;
    assert(iter->first == 13);
    assert(iter->second == 15);
    ++iter;
    assert(iter->first == 20);
    assert(iter->second == 26);
    ++iter;
    assert(iter == reg.myRegions.end());
    for(int i = 0; i < 30; i++)
    {
        if(((i >= 4) && (i < 6)) || ((i >= 8) && (i < 9)) ||
           ((i >= 13) && (i < 15)) || ((i >= 20) && (i < 26)))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }

    // Add another region
    reg.add(18,19);
    assert(reg.myRegions.size() == 5);
    assert(reg.myRegionIter->first == 18);
    assert(reg.myRegionIter->second == 19);
    iter = reg.myRegions.begin();
    assert(iter->first == 4);
    assert(iter->second == 6);
    ++iter;
    assert(iter->first == 8);
    assert(iter->second == 9);
    ++iter;
    assert(iter->first == 13);
    assert(iter->second == 15);
    ++iter;
    assert(iter->first == 18);
    assert(iter->second == 19);
    ++iter;
    assert(iter->first == 20);
    assert(iter->second == 26);
    ++iter;
    assert(iter == reg.myRegions.end());
    for(int i = 0; i < 30; i++)
    {
        if(((i >= 4) && (i < 6)) || ((i >= 8) && (i < 9)) ||
           ((i >= 13) && (i < 15)) || ((i >= 18) && (i < 19)) ||
           ((i >= 20) && (i < 26)))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }

    // Start is not in, but overlap two others (ending not at the end).
    reg.add(12,19);
    assert(reg.myRegions.size() == 4);
    assert(reg.myRegionIter->first == 12);
    assert(reg.myRegionIter->second == 19);
    iter = reg.myRegions.begin();
    assert(iter->first == 4);
    assert(iter->second == 6);
    ++iter;
    assert(iter->first == 8);
    assert(iter->second == 9);
    ++iter;
    assert(iter->first == 12);
    assert(iter->second == 19);
    ++iter;
    assert(iter->first == 20);
    assert(iter->second == 26);
    ++iter;
    assert(iter == reg.myRegions.end());
    for(int i = 0; i < 30; i++)
    {
        if(((i >= 4) && (i < 6)) || ((i >= 8) && (i < 9)) ||
           ((i >= 12) && (i < 19)) || ((i >= 20) && (i < 26)))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }

    // Completely in region to left.
    reg.add(5,6);
    assert(reg.myRegions.size() == 4);
    assert(reg.myRegionIter->first == 4);
    assert(reg.myRegionIter->second == 6);
    iter = reg.myRegions.begin();
    assert(iter->first == 4);
    assert(iter->second == 6);
    ++iter;
    assert(iter->first == 8);
    assert(iter->second == 9);
    ++iter;
    assert(iter->first == 12);
    assert(iter->second == 19);
    ++iter;
    assert(iter->first == 20);
    assert(iter->second == 26);
    ++iter;
    assert(iter == reg.myRegions.end());
    for(int i = 0; i < 30; i++)
    {
        if(((i >= 4) && (i < 6)) || ((i >= 8) && (i < 9)) ||
           ((i >= 12) && (i < 19)) || ((i >= 20) && (i < 26)))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }

    // Completely in region to right.
    reg.add(22,24);
    assert(reg.myRegions.size() == 4);
    assert(reg.myRegionIter->first == 20);
    assert(reg.myRegionIter->second == 26);
    iter = reg.myRegions.begin();
    assert(iter->first == 4);
    assert(iter->second == 6);
    ++iter;
    assert(iter->first == 8);
    assert(iter->second == 9);
    ++iter;
    assert(iter->first == 12);
    assert(iter->second == 19);
    ++iter;
    assert(iter->first == 20);
    assert(iter->second == 26);
    ++iter;
    assert(iter == reg.myRegions.end());
    for(int i = 0; i < 30; i++)
    {
        if(((i >= 4) && (i < 6)) || ((i >= 8) && (i < 9)) ||
           ((i >= 12) && (i < 19)) || ((i >= 20) && (i < 26)))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }

    // Add region to right.
    reg.add(28,29);
    assert(reg.myRegions.size() == 5);
    assert(reg.myRegionIter->first == 28);
    assert(reg.myRegionIter->second == 29);
    iter = reg.myRegions.begin();
    assert(iter->first == 4);
    assert(iter->second == 6);
    ++iter;
    assert(iter->first == 8);
    assert(iter->second == 9);
    ++iter;
    assert(iter->first == 12);
    assert(iter->second == 19);
    ++iter;
    assert(iter->first == 20);
    assert(iter->second == 26);
    ++iter;
    assert(iter->first == 28);
    assert(iter->second == 29);
    ++iter;
    assert(iter == reg.myRegions.end());
    for(int i = 0; i < 30; i++)
    {
        if(((i >= 4) && (i < 6)) || ((i >= 8) && (i < 9)) ||
           ((i >= 12) && (i < 19)) || ((i >= 20) && (i < 26)) ||
           ((i >= 28) && (i < 29)))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }

    // Add region to left, start is in the region, and end extends past.
    reg.add(8,10);
    assert(reg.myRegions.size() == 5);
    assert(reg.myRegionIter->first == 8);
    assert(reg.myRegionIter->second == 10);
    iter = reg.myRegions.begin();
    assert(iter->first == 4);
    assert(iter->second == 6);
    ++iter;
    assert(iter->first == 8);
    assert(iter->second == 10);
    ++iter;
    assert(iter->first == 12);
    assert(iter->second == 19);
    ++iter;
    assert(iter->first == 20);
    assert(iter->second == 26);
    ++iter;
    assert(iter->first == 28);
    assert(iter->second == 29);
    ++iter;
    assert(iter == reg.myRegions.end());
    for(int i = 0; i < 30; i++)
    {
        if(((i >= 4) && (i < 6)) || ((i >= 8) && (i < 10)) ||
           ((i >= 12) && (i < 19)) || ((i >= 20) && (i < 26)) ||
           ((i >= 28) && (i < 29)))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }

    // Add region  start is in the region, and end extends past and overlaps
    // the next region.
    reg.add(5,9);
    assert(reg.myRegions.size() == 4);
    assert(reg.myRegionIter->first == 4);
    assert(reg.myRegionIter->second == 10);
    iter = reg.myRegions.begin();
    assert(iter->first == 4);
    assert(iter->second == 10);
    ++iter;
    assert(iter->first == 12);
    assert(iter->second == 19);
    ++iter;
    assert(iter->first == 20);
    assert(iter->second == 26);
    ++iter;
    assert(iter->first == 28);
    assert(iter->second == 29);
    ++iter;
    assert(iter == reg.myRegions.end());
    for(int i = 0; i < 30; i++)
    {
        if(((i >= 4) && (i < 10)) ||
           ((i >= 12) && (i < 19)) || ((i >= 20) && (i < 26)) ||
           ((i >= 28) && (i < 29)))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }

    // Add region  start is in the region, and end extends past and overlaps
    // the next region.
    reg.add(10,11);
    assert(reg.myRegions.size() == 5);
    assert(reg.myRegionIter->first == 10);
    assert(reg.myRegionIter->second == 11);
    iter = reg.myRegions.begin();
    assert(iter->first == 4);
    assert(iter->second == 10);
    ++iter;
    assert(iter->first == 10);
    assert(iter->second == 11);
    ++iter;
    assert(iter->first == 12);
    assert(iter->second == 19);
    ++iter;
    assert(iter->first == 20);
    assert(iter->second == 26);
    ++iter;
    assert(iter->first == 28);
    assert(iter->second == 29);
    ++iter;
    assert(iter == reg.myRegions.end());
    for(int i = 0; i < 30; i++)
    {
        if(((i >= 4) && (i < 10)) || ((i >= 10) && (i < 11)) ||
           ((i >= 12) && (i < 19)) || ((i >= 20) && (i < 26)) ||
           ((i >= 28) && (i < 29)))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }

    // Add region start is in the region, and end extends past and overlaps
    // the next 2 regions.
    reg.add(10,24);
    assert(reg.myRegions.size() == 3);
    assert(reg.myRegionIter->first == 10);
    assert(reg.myRegionIter->second == 26);
    iter = reg.myRegions.begin();
    assert(iter->first == 4);
    assert(iter->second == 10);
    ++iter;
    assert(iter->first == 10);
    assert(iter->second == 26);
    ++iter;
    assert(iter->first == 28);
    assert(iter->second == 29);
    ++iter;
    assert(iter == reg.myRegions.end());
    for(int i = 0; i < 30; i++)
    {
        if(((i >= 4) && (i < 10)) || ((i >= 10) && (i < 26)) ||
           ((i >= 28) && (i < 29)))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }

    // Add region start outside of a region and ends at the end.
    reg.add(2,30);
    assert(reg.myRegions.size() == 1);
    assert(reg.myRegionIter->first == 2);
    assert(reg.myRegionIter->second == 30);
    iter = reg.myRegions.begin();
    assert(iter->first == 2);
    assert(iter->second == 30);
    ++iter;
    assert(iter == reg.myRegions.end());
    for(int i = 0; i < 50; i++)
    {
        if(((i >= 2) && (i < 30)))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }

    // Add invalid region (start = end)
    reg.add(40,40);
    assert(reg.myRegions.size() == 1);
    iter = reg.myRegions.begin();
    assert(iter->first == 2);
    assert(iter->second == 30);
    ++iter;
    assert(iter == reg.myRegions.end());
    for(int i = 0; i < 50; i++)
    {
        if(((i >= 2) && (i < 30)))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }

    // Add invalid region (start < end)
    reg.add(40, 38);
    assert(reg.myRegions.size() == 1);
    iter = reg.myRegions.begin();
    assert(iter->first == 2);
    assert(iter->second == 30);
    ++iter;
    assert(iter == reg.myRegions.end());
    for(int i = 0; i < 50; i++)
    {
        if(((i >= 2) && (i < 30)))
        {
            assert(reg.inRegion(i) == true);
        }
        else
        {
            assert(reg.inRegion(i) == false);
        }
    }
}
