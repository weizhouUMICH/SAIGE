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

//////////////////////////////////////////////////////////////////////////

#include "NonOverlapRegions.h"
#include <iostream>

NonOverlapRegions::NonOverlapRegions()
    : myRegions()
{
}


NonOverlapRegions::~NonOverlapRegions()
{
    myRegions.clear();
}


void NonOverlapRegions::add(const char* chrom, int32_t start, int32_t end)
{
    // Add the region.
    myRegions[chrom].add(start, end);
}


bool NonOverlapRegions::inRegion(const char* chrom, int32_t pos)
{
    // Return whether or not the position was found within a region.
    // Note, this will create a NonOverlapRegion for this chrom if it
    // did not already exist, but it won't have any regions.
    return(myRegions[chrom].inRegion(pos));
}


NonOverlapRegionPos::NonOverlapRegionPos()
    : myRegions()
{
    myRegionIter = myRegions.begin();
    myTmpIter = myRegions.begin();
}


NonOverlapRegionPos::NonOverlapRegionPos(const NonOverlapRegionPos& reg)
    : myRegions()
{
    myRegionIter = myRegions.begin();
    myTmpIter = myRegions.begin();
}


NonOverlapRegionPos::~NonOverlapRegionPos()
{
    myRegionIter = myRegions.begin();
    myTmpIter = myRegions.begin();
    myRegions.clear();
}


void NonOverlapRegionPos::add(int32_t start, int32_t end)
{
    // Check to see if the start/end are valid in relation.
    if(start >= end)
    {
        std::cerr << "NonOverlapRegionPos::add: Invalid Range, "
                  << "start must be < end, but " << start << " >= " << end 
                  << std::endl;
        return;
    }

    bool added = false;
    // Locate the correct position in the region list for this start/end.
    if(inRegion(start))
    {
        // Check if the region end needs to be updated.
        if(end > myRegionIter->second)
        {
            myRegionIter->second = end;
        }
        added = true;
    }
    else
    {
        // Check to see if we are at the end.
        if(myRegionIter != myRegions.end())
        {
            // Not at the end.
            // Check to see if the region overlaps the current region.
            if(end >= myRegionIter->first)
            {
                // Overlaps, so update the start.
                myRegionIter->first = start;
                // Check if the end needs to be updated.
                if(myRegionIter->second < end)
                {
                    myRegionIter->second = end;
                }
                added = true;
            }
        }
    }

    // If we already added the record, check to see if the end of the
    // new region overlaps any additional regions (know that myRegionIter is
    // not at the end.
    if(added)
    {
        // Check to see if any other regions were overlapped by this record.
        myTmpIter = myRegionIter;
        ++myTmpIter;
        while(myTmpIter != myRegions.end())
        {
            // If the region starts before the end of this one, consume it.
            if(myTmpIter->first <= end)
            {
                if(myTmpIter->second > myRegionIter->second)
                {
                    // Update this region with the new end.
                    myRegionIter->second = myTmpIter->second;
                }
                
                myTmpIter = myRegions.erase(myTmpIter);
            }
            else
            {
                // This region is not overlapped by the new region, so stop.
                break;
            }
        }
    }
    else
    {
        // Add the region.
        myRegionIter = myRegions.insert(myRegionIter, 
                                         std::make_pair(start, end));
    }
}


bool NonOverlapRegionPos::inRegion(int32_t pos)
{
    // Return whether or not the position was found within a region.
    // If it is found within the region, myRegionIter will point to the region
    // otherwise myRegionIter will point to the region after the position 
    // or to the end if the position is after the last region.

    // Determine if it needs to search to the left
    //   a) it is at the end
    //   b) the region starts after the position.
    if(myRegionIter == myRegions.end())
    {
        // If the iterator is at the end, search to the left.
        return(findLeft(pos));
    }
    else if(pos < myRegionIter->first)
    {
        // Not at the end, so search left if the position is less
        // than this region's start.
        return(findLeft(pos));
    }
    else
    {
        return(findRight(pos));
    }
}


bool NonOverlapRegionPos::findRight(int32_t pos)
{
    // Keep looping until the end or until the position is found.
    while(myRegionIter != myRegions.end())
    {
        // Check to see if we have passed the position.
        if(pos < myRegionIter->first)
        {
            // stop here, position comes before this region,
            // so myRegionIter is pointing to just after it,
            // but was not found in the region.
            return(false);
        }
        else if(pos < myRegionIter->second)
        {
            // the position is in the region, so return true.
            return(true);
        }
        else
        {
            // The position is after this region, so increment.
            ++myRegionIter;
        }
    }
    // exited because we are at the end of the regions and the position was
    // not found.
    return(false);
}


bool NonOverlapRegionPos::findLeft(int32_t pos)
{
    if(myRegionIter == myRegions.end())
    {
        if(myRegionIter == myRegions.begin())
        {
            // There is nothing in this list, so just return.
            return(false);
        }
        // Move 1 lower than the end.
        --myRegionIter;
    }

    while(myRegionIter->first > pos)
    {
        // This region is before our position, so move to the previous region
        // unless this is the first region in the list.
        if(myRegionIter == myRegions.begin())
        {
            // Did not find the position and the iter is at the element
            // just after the position.
            return(false);
        }
        // Not yet to the beginning of the list, so decrement.
        --myRegionIter;
    }

    // At this point, the regions iter points to a region that starts
    // before the position.
    // Determine if the position is in the region by checking if it is 
    // less than the end of the region.
    if(pos < myRegionIter->second)
    {
        // in the region.
        return(true);
    }

    // This region ends before this position.  The iterator needs to point
    // to the region after the position, so increment it.
    ++myRegionIter;
    return(false);
}
