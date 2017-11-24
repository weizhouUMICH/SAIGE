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

#include "SamHeaderRecord.h"

// Constructor
SamHeaderRecord::SamHeaderRecord()
    : myTagHash(),
      myTags(),
      myNumActiveTags(0)
{
}


// Destructor
SamHeaderRecord::~SamHeaderRecord()
{
    reset();
}


// Set the fields from the passed in line.
// Return true if successfully set.
bool SamHeaderRecord::setFields(const StringArray& tokens)
{
    bool status = true;
   
    // Loop through the tags for this type.
    // The tags start in column 1 since column 0 contains the type.
    for(int columnIndex = 1; columnIndex < tokens.Length(); columnIndex++)
    {
        // Validate that the tag is at least 3 characters. Two for the token,
        // one for the ':'.
        if((tokens[columnIndex].Length() < 3) || 
           (tokens[columnIndex][2] != ':'))
        {
            // Continue to the next tag, this one is too small/invalid.
            status = false;
            std::cerr << "ERROR: Poorly formatted tag in header: " 
                      << tokens[columnIndex] << std::endl;
            continue;
        }
      
        // Get the tag from the token.
        char tag[3];
        tag[0] = tokens[columnIndex][0];
        tag[1] = tokens[columnIndex][1];
        tag[2] = 0;

        // The tag value is the rest of the substring.
        String tagValue = (tokens[columnIndex]).SubStr(3);

        // Set the tag.      
        status &= setTag(tag, tagValue.c_str());
    }

    status &= isValid();

    return(status);
}


// Check to see if the record is valid.
bool SamHeaderRecord::isValid()
{
    bool status = true;
    // Check that the required tags are set. If they aren't, return false.
    for(unsigned int reqIndex = 0; reqIndex < myRequiredTags.size(); reqIndex++)
    {
        // Check to see if the required tag at this index exists and has
        // a value.
        int index = myTagHash.Integer(myRequiredTags[reqIndex].c_str());
        if((index < 0) || !(myTags[index]->hasValue()))
        {
            // Did not find the tag, stet status to false.
            std::cerr << "ERROR: Missing required tag: " 
                      << myRequiredTags[reqIndex] << "." << std::endl;
            status = false;
        }
    }
    return(status);
}


// Return the value associated with the specified tag.
const char* SamHeaderRecord::getTagValue(const char* tag) const
{
    // Look up the tag in myTags.
    int index = myTagHash.Integer(tag);
    if(index < 0)
    {
        // The tag was not found in the hash, so return "".
        return("");
    }

    // The tag was found in the hash, so return the tag value found at the 
    // index associated with the tag.
    return(myTags[index]->getValue());
}


// Set the value of the specified tag to the specified value.
// Set value to NULL in order to delete the tag.
// Returns whether or not it was successful.
bool SamHeaderRecord::setTag(const char* tag, const char* value)
{
    // Lookup the tag in the hash.
    int vectorIndex = myTagHash.Integer(tag);
    if(vectorIndex < 0)
    {
        // The tag was not found in the hash, so create a new one.
        SamHeaderTag* tagPtr = new SamHeaderTag(tag, value);
      
        if(tagPtr == NULL)
        {
            // Failed to allocate the tag, return false.
            std::cerr << "Failed to allocate space (new) for a SamHeaderTag.\n";
            return(false);
        }

        // Add the new tag to the back of the tag values.
        vectorIndex = myTags.size();
        myTags.push_back(tagPtr);

        // If the value is not null, increment the number of active tags.
        if(value[0] != 0)
        {
            ++myNumActiveTags;
        }

        // Add the tag to the hash.
        int hashIndex = myTagHash.Add(tag, vectorIndex);
      
        if((myTagHash.Integer(hashIndex) != vectorIndex) ||
           (myTagHash[hashIndex] != tag))
        {
            // Failed to add the tag, so return false.
            std::cerr << "Failed to add tag, " << tag
                      << ", to the hash." << std::endl;
            return(false);
        }
        return(true);
    }
    else if((unsigned int)vectorIndex < myTags.size())
    {
        // Found the tag in the hash.  So, update the tag if it
        // is not the key.
        if(myKeyTag != tag)
        {
            // Not the key, so update the tag.
            // If the new value is null and the old one is not, decrement the
            // number of active tags.
            if((value[0] == 0) && ((myTags[vectorIndex]->getValue())[0] != 0))
            {
                // Tag was deleted since the new value is blank but the old
                // value was not.
                --myNumActiveTags;
            }
            else if((value[0] != 0) && 
                    ((myTags[vectorIndex]->getValue())[0] == 0))
            {
                // Tag was added since the old value was blank and the new value
                // is not.
                ++myNumActiveTags;
            }

            // Just modifying a tag, so this does not affect the number 
            // of active tags.
            return(myTags[vectorIndex]->setValue(value));
        }
        else if(strcmp(value, myTags[vectorIndex]->getValue()) == 0)
        {
            // The new key value is the same as the previous value, so
            // it is not a change, return true.
            return(true);
        }
        else
        {
            // Can't modify the key tag's value since that will
            // screw up the hash.
            std::cerr << "Can't modify the key tag, " << tag << " from "
                      << myTags[vectorIndex]->getValue() << " to " 
                      << value << std::endl;
            return(false);
        }
    }

    // Got an invalid index from the hash.  This is not supposed to happen.
    // so return false.
    std::cerr << "Invalid tag index found: " << vectorIndex 
              << ", but max index is " << myTags.size() << " for tag: " 
              << tag << std::endl;
    return(false);
}


// Reset this header record to an empty state.
void SamHeaderRecord::reset()
{
    // Delete the tag hash.
    myTagHash.Clear();

    // Loop through deleting all the tags in the vector.
    for(unsigned int vectorIndex = 0; 
        vectorIndex < myTags.size();
        vectorIndex++)
    {
        delete myTags[vectorIndex];
        myTags[vectorIndex] = NULL;
    }
    // Clear the tag vector.
    myTags.clear();

    myNumActiveTags = 0;
}


// Appends the string representation of this header record
// to the passed in string.
bool SamHeaderRecord::appendString(std::string& header)
{
    // Track whether or not the header type has been written.
    // Only write the header type if at least one of the tags has
    // an associated value.
    bool writtenHeader = false;

    if(isActiveHeaderRecord() && isValid())
    {
        // Loop through all the entries in the tag vector.
        for(unsigned int vectorIndex = 0; 
            vectorIndex < myTags.size(); 
            vectorIndex++)
        {
            if(!writtenHeader && (myTags[vectorIndex]->hasValue()))
            {
                // The tag has a value and the header type has not yet been written,
                // so write it.
                header += "@";
                header += myTypeString;
                writtenHeader = true;
            }
            myTags[vectorIndex]->getTagString(header);
        }

        // If a header has been written, add a new line character.
        if(writtenHeader)
        {
            header += "\n";
            return(true);
        }
    }

    // Nothing was written, return false.
    return(false);
}


// Add the key tag with the specified value.
bool SamHeaderRecord::addKey(const char* value)
{
    if(myKeyTag.size() == 0)
    {
        return(false);
    }
    return(setTag(myKeyTag.data(), value));
}


// Return the value associated with the specified tag.
const char* SamHeaderRecord::getKeyValue() const
{
    // Look up the tag in myTags.
    int index = myTagHash.Integer(myKeyTag.c_str());
    if(index < 0)
    {
        // The tag was not found in the hash, so return "".
        return("");
    }

    // The tag was found in the hash, so return the tag value found at the 
    // index associated with the tag.
    return(myTags[index]->getValue());
}


// This header is active if there is at least one tag set.
bool SamHeaderRecord::isActiveHeaderRecord()
{
    return(myNumActiveTags != 0);
}


// Return the type of this header record.
const char* SamHeaderRecord::getTypeString()
{
    return(myTypeString.c_str());
}


// Return the type of this header record.
SamHeaderRecord::SamHeaderRecordType SamHeaderRecord::getType()
{
    return(myType);
}


void SamHeaderRecord::addRequiredTag(const char* requiredTag)
{
    myRequiredTags.push_back(requiredTag);
}


void SamHeaderRecord::internalCopy(SamHeaderRecord& newRec) const
{
    newRec.myTagHash = myTagHash;

    newRec.myTags.clear();

    // Loop through copying the tags.
    for(unsigned int vectorIndex = 0; 
        vectorIndex < myTags.size();
        vectorIndex++)
    {
        if(myTags[vectorIndex] != NULL)
        {
            newRec.myTags.push_back(new SamHeaderTag(*(myTags[vectorIndex])));
        }
    }
    newRec.myRequiredTags = myRequiredTags;
    newRec.myNumActiveTags = myNumActiveTags;
}
