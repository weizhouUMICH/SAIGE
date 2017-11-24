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

#include "SamHeaderTag.h"


SamHeaderTag::SamHeaderTag(const char* tag, const char* value)
{
    setTag(tag, value);
}


SamHeaderTag::SamHeaderTag(const SamHeaderTag& oldTag)
{
    setTag(oldTag.myTag.c_str(), oldTag.myValue.c_str());
}


SamHeaderTag::~SamHeaderTag()
{
}


// Add this tag to the passed in tag string.
// NOTE: does not clear tagString.
bool SamHeaderTag::getTagString(std::string& tagString)
{
    if(myValue.length() != 0)
    {
        // There is a value associated with this tag, so add it to the string.
        tagString += "\t";
        tagString += myTag;
        tagString += ":";
        tagString += myValue;
        return(true);
    }
    // This tag has no associated value, return false.
    return(false);
}


// Set this tag to the passed in tag and value.
bool SamHeaderTag::setTag(const char* tag, const char* value)
{
    myTag = tag;
    myValue = value;
    return(true);
}


// Set the value associated with this tag to the passed in value.
bool SamHeaderTag::setValue(const char* value)
{
    myValue = value;
    return(true);
}


// Return the tag for this tag.
const char* SamHeaderTag::getTag()
{
    return(myTag.c_str());
}


// Return the value associated with this tag.
const char* SamHeaderTag::getValue()
{
    return(myValue.c_str());
}


// Return true if there is a non-blank value associated with this tag.
bool SamHeaderTag::hasValue()
{
    return(myValue.size() != 0);
}
