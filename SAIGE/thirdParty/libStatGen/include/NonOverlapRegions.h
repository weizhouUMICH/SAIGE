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

#ifndef __NONOVERLAP_REGIONS_H__
#define __NONOVERLAP_REGIONS_H__

#include <map>
#include <string>
#include <list>
#include <stdint.h>

/// This class contains a list of non-overlapping regions, just positions, not
/// including chromosomes (see NonOverlapRegions for chromosomes and positions).
///  When regions are added that overlap, it merges them.  After adding regions,
/// you can check to see if a position is found in one of the regions.  It is
/// designed to work fastest if you make calls in sequential order.
class NonOverlapRegionPos
{
public:
    friend class NonOverlapRegionsTest;
    NonOverlapRegionPos();
    /// Copy constructor, does not copy, but initializes with an empty
    /// region list.
    NonOverlapRegionPos(const NonOverlapRegionPos& reg);

    ~NonOverlapRegionPos();

    /// End position is not included in the region.
    /// If this region overlaps another region(s), they will be merged into
    /// one region.
    void add(int32_t start, int32_t end);

    /// Return whether or not the position was found within a region.
    /// If it is found within the region, myRegionIter will point to the region
    /// otherwise myRegionIter will point to the region after the position
    /// or to the end if the position is after the last region.
    bool inRegion(int32_t pos);

private:
    // True if pos found in the region pointed to by myRegionIter or to
    // the right of myRegionIter.  If the position is found in a region,
    // myRegionIter will point to the region containing the position.
    // If the position is not found in a region, myRegionIter will point
    // to the region after the position, or to the end if the position is
    // after the last region.
    bool findRight(int32_t pos);

    // True if pos found in the region pointed to by myRegionIter or to
    // the left of myRegionIter.  If the position is found in a region,
    // myRegionIter will point to the region containing the position.
    // If the position is not found in a region, myRegionIter will point
    // to the region after the position, or to the end if the position is
    // after the last region.
    bool findLeft(int32_t pos);


    std::list< std::pair<int32_t, int32_t> > myRegions;
    std::list< std::pair<int32_t, int32_t> >::iterator myRegionIter;
    std::list< std::pair<int32_t, int32_t> >::iterator myTmpIter;
};


/// This class contains a list of non-overlapping regions.  When regions are
/// added that overlap, it merges them.  After adding regions, you can check
/// to see if a position is found in one of the regions.  It is designed to
/// work fastest if you make calls in sequential order.
class NonOverlapRegions
{
public:
    friend class NonOverlapRegionsTest;

    NonOverlapRegions();
    ~NonOverlapRegions();

    /// End position is not included in the region.
    /// If this region overlaps another region(s), they will be merged into 
    /// one region.
    void add(const char* chrom, int32_t start, int32_t end);
 
    /// Return whether or not the position was found within a region.
    /// If it is found within the region, myRegionIter will point to the region
    /// otherwise myRegionIter will point to the region after the position 
    /// or to the end if the position is after the last region.
    bool inRegion(const char* chrom, int32_t pos);

private:
    // Copy Constructor - unimplimented.
    NonOverlapRegions(const NonOverlapRegions& reg);

    std::map<std::string, NonOverlapRegionPos> myRegions;
};

#endif
