/*
 *  Copyright (C) 2013  Regents of the University of Michigan
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

#include "VcfRecordDiscardRules.h"

void VcfRecordDiscardRules::reset()
{
    myExcludeIDs.clear();
    myIncludeIDs.clear();
    myNumDiscarded = 0;
}


bool VcfRecordDiscardRules::setExcludeIDs(const char* filename)
{
    return(setIDs(myExcludeIDs, filename));
}


bool VcfRecordDiscardRules::setIncludeIDs(const char* filename)
{
    return(setIDs(myIncludeIDs, filename));
}


bool VcfRecordDiscardRules::discardForID(std::string& myID)
{
    if(!myExcludeIDs.empty())
    {
        if(myExcludeIDs.find(myID) != myExcludeIDs.end())
        {
            // The ID is in the exclude list,
            // so return true, discard the record.
            // increment the discard counter.
            ++myNumDiscarded;
            return(true);
        }
    }
    else if(!myIncludeIDs.empty())
    {
        if(myIncludeIDs.find(myID) == myIncludeIDs.end())
        {
            // The ID is not in the include list,
            // so return false, discard the record.
            // increment the discard counter.
            ++myNumDiscarded;
            return(true);
        }
    }
    return(false);
}


bool VcfRecordDiscardRules::setIDs(IDList& idlist, const char* filename)
{
    // Open the file nad read in all the exclude ids.
    IFILE idFile = ifopen(filename, "r");
    if(idFile == NULL)
    {
        return(false);
    }
    std::string line;
    while(idFile->readLine(line) == 0)
    {
        idlist.insert(line);
        line.clear();
    }
    if(!line.empty())
    {
        idlist.insert(line);
        line.clear();
    }
    return(true);
}
