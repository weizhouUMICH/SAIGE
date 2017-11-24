/*
 *  Copyright (C) 2011  Regents of the University of Michigan,
 *                      Hyun Min Kang, Matthew Flickenger, Matthew Snyder,
 *                      and Goncalo Abecasis
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

 #include <string>

#include "VcfRecordFilter.h"
#include "VcfHelper.h"


bool VcfRecordFilter::ourParseFilter = false;

bool VcfRecordFilter::read(IFILE filePtr)
{
    static const std::string fieldStopCharsNoParse = "\n\t";
    static const int tabPos = 1;

    static std::string fieldStopChars = fieldStopCharsNoParse;
    fieldStopChars += FILTER_DELIM;

    // The start of the first character in stopChars that means there is more
    // filter info in the format field, so continue reading the format field.
    static const int contPos = fieldStopCharsNoParse.length();

    // Clear out any previously set values.
    reset();
    
    if(ifeof(filePtr))
    {
        // End of file, just return false.
        return(false);
    }
    

    // Check how much  the filter should be parsed.
    int stopPos = 0;
    if(!ourParseFilter)
    {
        // Do not need to parse the filter, so just read until the tab.
        stopPos = filePtr->readTilChar(fieldStopCharsNoParse, myFilterString);
    }
    else
    {
        // Parse the filter as we go.
        stopPos = contPos;
        std::string* nextFilter;
        while(stopPos >= contPos)
        {
            nextFilter = &(myFilterVector.getNextEmpty());
            stopPos = filePtr->readTilChar(fieldStopChars, *nextFilter);
        }
    }
    return(stopPos == tabPos);
}


void VcfRecordFilter::reset()
{
    myFilterString.clear();
    myFilterVector.reset();
}


bool VcfRecordFilter::passedAllFilters()
{
    static std::string pass("PASS");

    return(pass.compare(getString()) == 0);
}


const std::string& VcfRecordFilter::getString()
{
    // Check if the filter string needs to be set.
    if(myFilterString.size() == 0)
    {
        // Filter string is not yet set, so set it.
        for(int i = 0; i < myFilterVector.size(); i++)
        {
            if(i != 0)
            {
                myFilterString += ';';
            }
            myFilterString += myFilterVector.get(i);
        }
    }
    return(myFilterString);
}

int VcfRecordFilter::getNumFilters()
{
    // Check if the filter has been parsed.
    if(myFilterVector.size() == 0)
    {
        // Filter is not parsed, so parse the filter.
        VcfHelper::parseString(myFilterString, FILTER_DELIM, myFilterVector);
    }

    return(myFilterVector.size());
}


const std::string& VcfRecordFilter::getString(int index)
{
    // Check if the filter has been parsed yet.
    if(myFilterVector.size() == 0)
    {
        // Filter has not yet been parsed, so parse it.
        VcfHelper::parseString(myFilterString, FILTER_DELIM, myFilterVector);
    }
    return(myFilterVector.get(index));
}


void VcfRecordFilter::setFilter(const char* filter)
{
    reset();

    myFilterString = filter;
}


void VcfRecordFilter::addFilter(const char* filter)
{
    // If both are empty, add to both.
    if((myFilterString.size() == 0) && 
       (myFilterVector.size() == 0))
    {
        myFilterString += filter;
        myFilterVector.getNextEmpty() = filter;
    }
    else
    {
        // Either filter (or both) have contents, so append
        // as appropriate.
        if(myFilterString.size() != 0)
        {
            // String is set, so add the filter to the string.
            myFilterString += ';';
            myFilterString += filter;
        }
        
        if(myFilterVector.size() != 0)
        {
            // Vector is set, so add the filter to the vector.
            myFilterVector.getNextEmpty() = filter;
        }
    }
}
