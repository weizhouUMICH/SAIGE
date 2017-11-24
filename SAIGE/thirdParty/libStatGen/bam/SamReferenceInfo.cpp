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

#include "SamReferenceInfo.h"

SamReferenceInfo::SamReferenceInfo()
    : myReferenceContigs(),
      myReferenceHash(),
      myReferenceLengths()
{
    clear();
}


SamReferenceInfo::~SamReferenceInfo()
{
    clear();
}

// Add reference sequence name and reference sequence length.
void SamReferenceInfo::add(const char* referenceSequenceName, 
                           int32_t referenceSequenceLength)
{
        myReferenceHash.Add(referenceSequenceName, 
                            myReferenceContigs.Length());
        myReferenceContigs.Push(referenceSequenceName);
        myReferenceLengths.Push(referenceSequenceLength);
}


int SamReferenceInfo::getReferenceID(const String & referenceName, 
                                     bool addID)
{
    if (referenceName == "*")
        return -1;
    
    int id = myReferenceHash.Find(referenceName);

    if (id >= 0)
        return myReferenceHash.Integer(id);
    
    if(!addID)
    {
        // Don't add the id, so return NO_REF_ID
        return(NO_REF_ID);
    }

    id = myReferenceContigs.Length();
    myReferenceContigs.Push(referenceName);
    myReferenceLengths.Push(0);
    myReferenceHash.Add(referenceName, id);

    return id;
}


int SamReferenceInfo::getReferenceID(const char* referenceName, 
                                     bool addID)
{
    String referenceNameString = referenceName;

    return(getReferenceID(referenceNameString, addID));
}


const String & SamReferenceInfo::getReferenceLabel(int id) const
{
    static String noname("*");

    if ((id < 0) || (id >= myReferenceContigs.Length()))
    {
        return noname;
    }

    return myReferenceContigs[id];
}


int32_t SamReferenceInfo::getNumEntries() const
{
    // The number of entries is the size of referenceLengths.
    return(myReferenceLengths.Length());
}


const char* SamReferenceInfo::getReferenceName(int index) const
{
    if((index >= 0) && (index < getNumEntries()))
    {
        return(myReferenceContigs[index].c_str());
    }
    
    // Out of range, return blank
    return("");
}
   

int32_t SamReferenceInfo::getReferenceLength(int index) const
{
    if((index >= 0) && (index < getNumEntries()))
    {
        return(myReferenceLengths[index]);
    }
    
    // Out of bounds, return 0
    return(0);
}

void SamReferenceInfo::clear()
{
    myReferenceContigs.Clear();
    myReferenceHash.Clear();
    myReferenceLengths.Clear();
}


SamReferenceInfo& SamReferenceInfo::operator = (const SamReferenceInfo &newInfo)
{
    clear();
    // Copy Reference contigs, hash, lengths.
    myReferenceContigs = newInfo.myReferenceContigs;
    myReferenceHash = newInfo.myReferenceHash;
    myReferenceLengths = newInfo.myReferenceLengths;
    return(*this);
}


bool SamReferenceInfo::operator== (const SamReferenceInfo& rhs) const
{
    // Hash may be different, but if Contigs are the same, the hashes will
    // contain the same basic info (maybe just at different indices.
    return((myReferenceContigs == rhs.myReferenceContigs) &&
           (myReferenceLengths == rhs.myReferenceLengths));
}
